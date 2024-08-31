MODULE ConvPar_GF2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- This convective parameterization is build to attempt                 !
!  a smooth transition to cloud resolving scales as proposed            !
!  by Arakawa et al (2011, ACP). The scheme is  described               !
!  in the paper Grell and Freitas (ACP, 2014).                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- Implemented in GEOS5 GCM by Saulo Freitas (July 2016)                !
!- Use the following references for this implementation:                !
!- Freitas et al (2018, JAMES/AGU, https://doi.org/10.1029/2017MS001251)!
!- Freitas et al (2020, GMD/EGU,   https://doi.org/10.5194/gmd-2020-38) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
USE module_gate
USE MAPL
USE ConvPar_GF_SharedParams
USE GEOSmoist_Process_Library, only : sigma, SH_MD_DP, ICE_FRACTION, make_DropletNumber, make_IceNumber

 IMPLICIT NONE
 PRIVATE
 PUBLIC  gf2020_interface &
        ,use_memory,convection_tracer &
        ,downdraft                    &
        ,use_rebcb, vert_discr, satur_calc, clev_grid, apply_sub_mp        &
        ,sgs_w_timescale, lightning_diag, tau_ocea_cp,  tau_land_cp        &
        ,autoconv, bc_meth,overshoot,use_wetbulb                           &
        ,c1,c0_deep,qrc_crit_lnd,qrc_crit_ocn                              &
        ,lambau_deep,lambau_shdn,c0_mid                                    &
        ,cum_min_edt_land  ,cum_min_edt_ocean                              &
        ,cum_max_edt_land  ,cum_max_edt_ocean, cum_hei_down_land           &
        ,cum_hei_down_ocean,cum_hei_updf_land, cum_hei_updf_ocean          &
        ,use_momentum_transp,cum_entr_rate,min_entr_rate                   &
        ,zero_diff , nmp, lsmp, cnmp,moist_trigger,frac_modis,max_tq_tend  &
        ,cum_fadj_massflx, cum_use_excess, cum_ave_layer, adv_trigger      &
        ,evap_fix,output_sound,use_cloud_dissipation      &
        ,use_smooth_tend,GF_convpar_init,beta_sh,c0_shal                   &
        ,use_linear_subcl_mf,cap_maxs,entrversion

 LOGICAL :: wrtgrads = .false.
 INTEGER :: nrec = 0, ntimes = 0
 REAL    :: int_time = 0.
 !-
 !- number of microphysics schemes in the host model
 INTEGER ,PARAMETER  :: nmp = 2, lsmp = 1, cnmp = 2

 INTEGER :: USE_MEMORY        =-1 != -1/0/1/2 .../10    !-

 INTEGER :: CONVECTION_TRACER = 0 != 0/1:  turn ON/OFF the "convection" tracer

 INTEGER :: CLEV_GRID         = 1 != 0/1/2: interpolation method to define environ state at the cloud levels (at face layer), default = 0
                                  != clev_grid = 0 default method
                                  != clev_grid = 1 interpolation method based on Tiedtke (1989)
                                  != clev_grid = 2 for GATE soundings only

 INTEGER :: USE_REBCB         = 1 != 0/1: turn ON/OFF rainfall evap below cloud base, default = 0

 INTEGER :: VERT_DISCR        = 1 != 0/1: 1=new vert discretization, default = 0

 INTEGER :: SATUR_CALC        = 1 != 0/1: 1=new saturation specific humidity calculation, default = 0

 INTEGER :: SGS_W_TIMESCALE   = 0 != 0/1: vertical velocity for tau_ecmwf, default = 0

 INTEGER :: LIGHTNING_DIAG    = 0 != 0/1: do LIGHTNING_DIAGgnostics based on Lopez (2016, MWR)

 INTEGER :: APPLY_SUB_MP      = 0 != 0/1: subsidence transport applied the to grid-scale/anvil ice/liq mix
                                  !=      ratio and cloud fraction

 INTEGER :: USE_WETBULB       = 0 != 0/1

 REAL    :: OVERSHOOT         = 0.!= 0, 1

 REAL    :: MIN_ENTR_RATE     = 1.0/40000.0 ! minimum allowed entrainment rate [m-1]

 INTEGER :: AUTOCONV          = 1     != 1, 3 or 4 autoconversion formulation: (1) Kessler, (3) Kessler with temp dependence, (4) Sundvisqt

 INTEGER :: USE_MOMENTUM_TRANSP = 1   != 0/1:  turn ON/OFF conv transp of momentum
 REAL    ::  LAMBAU_DEEP        = 0.0 != default= 2.0 lambda parameter for deep/congestus convection momentum transp
 REAL    ::  LAMBAU_SHDN        = 2.0 != default= 2.0 lambda parameter for shallow/downdraft convection momentum transp

 INTEGER :: DOWNDRAFT           = 1   != 0/1:  turn ON/OFF downdrafts, default = 1

 REAL    :: MAX_TQ_TEND         = 100.   != max T,Q tendency allowed (100 K/day)

 INTEGER :: ZERO_DIFF           = 0      != to get the closest solution of the stable version Dec 2019 for single-moment

 INTEGER :: USE_SMOOTH_TEND     = 0      != 0 => OFF, > 0 produces smoother tendencies (e.g.: for 1=> makes average between k-1,k,k+1)
 !---                                              deep, shallow, congestus
 REAL,   DIMENSION(maxiens) :: CUM_HEI_DOWN_LAND =(/0.30,  0.20,  0.20/)!= [0.2,0.8] height of the max Z Downdraft , default = 0.50
 REAL,   DIMENSION(maxiens) :: CUM_HEI_DOWN_OCEAN=(/0.30,  0.20,  0.20/)!= [0.2,0.8] height of the max Z Downdraft , default = 0.50

 REAL,   DIMENSION(maxiens) :: CUM_HEI_UPDF_LAND =(/0.35,  0.10,  0.10/)!= [0.2,0.8] height of the max Z Updraft   , default = 0.35
 REAL,   DIMENSION(maxiens) :: CUM_HEI_UPDF_OCEAN=(/0.35,  0.10,  0.10/)!= [0.2,0.8] height of the max Z Updraft   , default = 0.35

 REAL,   DIMENSION(maxiens) :: CUM_MIN_EDT_LAND  =(/0.10,  0.00,  0.10/)!= minimum evap fraction allowed over the land  ,default= 0.1
 REAL,   DIMENSION(maxiens) :: CUM_MIN_EDT_OCEAN =(/0.10,  0.00,  0.10/)!= minimum evap fraction allowed over the ocean ,default= 0.1

 REAL,   DIMENSION(maxiens) :: CUM_MAX_EDT_LAND  =(/0.90,  0.00,  0.90/)!= maximum evap fraction allowed over the land  ,default= 0.9
 REAL,   DIMENSION(maxiens) :: CUM_MAX_EDT_OCEAN =(/0.90,  0.00,  0.90/)!= maximum evap fraction allowed over the ocean ,default= 0.9

 REAL,   DIMENSION(maxiens) :: CUM_FADJ_MASSFLX  =(/1.00,  1.00,  1.00/)!= multiplicative factor for tunning the mass flux at cloud base
                                                                        != default = 1.0
 INTEGER,DIMENSION(maxiens) :: CUM_USE_EXCESS    =(/1,     1,     1   /)!= use T,Q excess sub-grid scale variability

 INTEGER :: MOIST_TRIGGER  = 0 != relative humidity effects on the cap_max trigger function
 INTEGER :: FRAC_MODIS     = 0 != use fraction liq/ice content derived from MODIS/CALIPO sensors
 INTEGER :: ADV_TRIGGER    = 0 !=  1=> Kain (2004),  2=> moisture adv trigger (Ma & Tan, 2009, Atmos Res)
 INTEGER :: EVAP_FIX       = 1 != fix total evap > column rainfall

 INTEGER :: OUTPUT_SOUND   = 0

 REAL    :: tau_ocea_cp    = 6.*3600.
 REAL    :: tau_land_cp    = 6.*3600.

 REAL    :: use_cloud_dissipation = 0.
 REAL    :: use_gustiness         = 0.
 REAL    :: use_random_num        = 0.
 REAL    :: dcape_threshold       = 0.
 REAL    :: beta_sh               = 2.2
 INTEGER :: use_linear_subcl_mf   = 1
 REAL    :: cap_maxs              = 50.

 !------------------- internal variables  -------------------

 !-- turn ON/OFF deep/shallow/mid plumes
 INTEGER, PARAMETER :: ON=1, OFF=0

 REAL    ::  HEI_DOWN_LAND     != [0.2,0.8] height of the max Z Downdraft , default = 0.50
 REAL    ::  HEI_DOWN_OCEAN    != [0.2,0.8] height of the max Z Downdraft , default = 0.50
 REAL    ::  HEI_UPDF_LAND     != [0.2,0.8] height of the max Z Updraft   , default = 0.35
 REAL    ::  HEI_UPDF_OCEAN    != [0.2,0.8] height of the max Z Updraft   , default = 0.35
 REAL    ::  MIN_EDT_LAND      != default= 0.1 - minimum evap fraction allowed over the land
 REAL    ::  MIN_EDT_OCEAN     != default= 0.1 - minimum evap fraction allowed over the ocean
 REAL    ::  MAX_EDT_LAND      != default= 0.9 - maximum evap fraction allowed over the land
 REAL    ::  MAX_EDT_OCEAN     != default= 0.9 - maximum evap fraction allowed over the ocean
 REAL    ::  FADJ_MASSFLX      != default= 1.0 - multiplicative factor for the mass flux at cloud base
 INTEGER ::  USE_EXCESS        != default= 1   - use T,Q excess sub-grid scale variability

 !-- General internal controls for the diverse options in GF
 INTEGER            :: ENTRVERSION    = 1       != entr formulations

 LOGICAL, PARAMETER :: COUPL_MPHYSICS = .TRUE.  != coupling with cloud microphysics (do not change  to false)

 LOGICAL, PARAMETER :: MELT_GLAC      = .TRUE.  != turn ON/OFF ice phase/melting

 LOGICAL, PARAMETER :: FEED_3DMODEL   = .TRUE.  != set "false" to not feedback the AGCM with the
                                                != heating/drying/transport conv tendencies
 LOGICAL            :: USE_C1D        = .FALSE. != turn ON/OFF the 'c1d' detrainment approach, don't change this.

 LOGICAL            :: FIRST_GUESS_W  = .FALSE. != use it to calculate a 1st guess of the updraft vert velocity

 INTEGER, PARAMETER :: LIQ_ICE_NUMBER_CONC = 1 !

 !-rainfall evaporation(1) orig (2) mix orig+new (3) new
 INTEGER, PARAMETER :: aeroevap = 1

 INTEGER, PARAMETER ::               &
                       maxens  = 1,  & ! 1  ensemble one on cap_max
                       maxens2 = 1,  & ! 1  ensemble two on precip efficiency
                       maxens3 = 16, & !16 ensemble three done in cup_forcing_ens16 for G3d
                       ensdim  = maxens*maxens2*maxens3,&
                       ens4    = 1
 !
 !- proportionality constant to estimate pressure
 !- gradient of updraft (Zhang and Wu, 2003, JAS)
!REAL, PARAMETER ::    pgcon=-0.55
 REAL, PARAMETER ::    pgcon= 0.0
 !- numerical constraints
 REAL, PARAMETER ::       &
  xmbmaxshal  =  0.05,    &  ! kg/m2/s
  mintracer   =  tiny(1.),&  ! kg/kg - tiny(x)
  smallerQV   =  1.e-16  ,&  ! kg/kg
  PI = 3.1415926536

 INTEGER :: whoami_all, JCOL

CONTAINS
!---------------------------------------------------------------------------------------------------
  SUBROUTINE GF2020_INTERFACE(  mxp,myp,mzp,LONS,LATS,DT_MOIST                    &
                               ,PLE, PLO, ZLE, ZLO, PK, MASS, KH          &
                               ,T1, TH1, Q1, U1,V1,W1,BYNCY,QLCN,QICN,QLLS,QILS   &
                               ,CNPCPRATE                                         &
                               ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQCDT ,ENTLAM     &
                               ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC, CLCN,CLLS     &
                               ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                               ,RADSW   ,RADLW ,DQDT_BL  ,DTDT_BL                 &
                               ,FRLAND  ,AREA  ,T2M ,Q2M     &
                               ,TA      ,QA    ,SH    ,EVAP  ,PHIS                &
                               ,KPBLIN  ,CNVFRC,SRFTYPE                           &
                               ,STOCHASTIC_SIG, SIGMA_DEEP, SIGMA_MID             &
                               ,DQDT_GF,DTDT_GF,DUDT_GF,DVDT_GF                   &
                               ,CNV_TOPP_DP,CNV_TOPP_MD,CNV_TOPP_SH               &
                               ,MUPDP,MUPSH,MUPMD,MDNDP                           &
                               ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD,WQT_DC           &
                               ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
                               ,DTDTDYN,DQVDTDYN                                  &
                               ,REVSU, entr3d, entr_dp, entr_md, entr_sh, PRFIL   &
                               ,TPWI,TPWI_star,LIGHTN_DENS                        &
                               ,CNV_TR)

    IMPLICIT NONE

    INTEGER                          ,INTENT(IN)   :: mxp,myp,mzp

    REAL                             ,INTENT(IN)   :: DT_moist

    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(IN)   :: PLE,ZLE

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(IN)   :: ZLO, PLO, PK, MASS, KH,       &
                                                      T1,TH1,Q1,U1,V1,W1,BYNCY,QLCN,QICN,QLLS,QILS, &
                                                      CLLS,CLCN

    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(IN)   :: PLE_DYN_IN

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(IN)   ::  QV_DYN_IN, U_DYN_IN, V_DYN_IN, T_DYN_IN, &
                                                      DTDTDYN, DQVDTDYN

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(IN)   :: RADSW, RADLW, DQDT_BL, DTDT_BL

    REAL   ,DIMENSION(mxp,myp)       ,INTENT(IN)   :: FRLAND, AREA, &
                                                      T2M, Q2M, TA, QA, SH, EVAP, PHIS,  &
                                                      LONS, LATS, &
                                                      STOCHASTIC_SIG

    REAL   ,DIMENSION(mxp,myp)       ,INTENT(IN)   :: KPBLIN, CNVFRC, SRFTYPE


    REAL   ,DIMENSION(mxp,myp)       ,INTENT(IN)   :: TPWI, TPWI_star

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(INOUT):: CNV_TR

    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(OUT)  :: CNV_MFC, WQT_DC

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(OUT)  :: CNV_MF0 , CNV_PRC3 , CNV_MFD , CNV_DQCDT,  &
                                                      CNV_UPDF, CNV_CVW  , CNV_QC  , ENTLAM


    REAL   ,DIMENSION(mxp,myp)       ,INTENT(OUT)  :: CNPCPRATE, LIGHTN_DENS

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(OUT)  :: REVSU, entr3d, entr_dp, entr_md, entr_sh

    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(OUT)  :: PRFIL

    !- Tendencies
    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(OUT)  :: DQDT_GF,DTDT_GF,DUDT_GF,DVDT_GF

    !-for debug/diagnostoc purposes
    REAL   ,DIMENSION(mxp,myp)       ,INTENT(OUT)  :: SIGMA_DEEP, SIGMA_MID
    REAL   ,DIMENSION(mxp,myp)       ,INTENT(OUT)  :: CNV_TOPP_DP,CNV_TOPP_MD,CNV_TOPP_SH
    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(OUT)  :: MUPDP,MDNDP,MUPSH,MUPMD
    REAL   ,DIMENSION(mxp,myp)       ,INTENT(OUT)  :: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
    REAL   ,DIMENSION(mxp,myp)       ,INTENT(OUT)  :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC

    INTEGER      :: ims,ime, jms,jme, kms,kme,    &
                    its,ite, jts,jte, kts,kte,    &
                    mynum

    REAL,  DIMENSION(mzp , mxp, myp ) :: up       &
                                        ,vp       &
                                        ,wp       &
                                        ,rvap     &
                                        ,temp     &
                                        ,press    &
                                        ,zm3d     &
                                        ,zt3d     &
                                        ,dm3d     &
                                        ,ec3d     &
                                        ,curr_rvap&
                                        ,buoy_exc &
                                        ,khloc



    REAL,  DIMENSION(mzp , mxp, myp ) ::          &
                                          gsf_t   & ! grid-scale forcing for temp
                                        , gsf_q   & ! advection forcing for rv
                                        ,advf_t   & ! advection forcing for temp
                                        ,sgsf_t   & ! sub-grid scale forcing for temp
                                        ,sgsf_q   & ! sub-grid scale forcing for rv
                                        ,SRC_T    & ! temp tendency      from convection
                                        ,SRC_Q    & ! rv tendency        from convection
                                        ,SRC_CI   & ! cloud/ice tendency from convection
                                        ,SRC_U    & ! U tendency         from convection
                                        ,SRC_V    & ! V tendency         from convection
                                        ,SRC_NI   & ! Ice     number tendency from convection
                                        ,SRC_NL   & ! Droplet number tendency from convection
                                        ,SRC_BUOY & ! buoyancy tendency from downdrafts
                                        ,REVSU_GF & ! evaporation_or_sublimation of_convective_precipitation kg kg-1 s-1
                                        ,PRFIL_GF   ! ice_or_liq convective_precipitation flux: kg m2 s-1 (deep only)


    REAL,  DIMENSION(nmp, mzp , mxp, myp ) ::     &
                                         mp_ice   &
                                        ,mp_liq   &
                                        ,mp_cf

    REAL,  DIMENSION(nmp, mzp , mxp, myp ) ::     &
                                         SUB_MPQI & ! subsidence transport applied to ice mix ratio
                                        ,SUB_MPQL & ! subsidence transport applied to cloud mix ratio
                                        ,SUB_MPCF   ! subsidence transport applied to cloud fraction


    REAL,  ALLOCATABLE, DIMENSION(:,:,:,:) :: SRC_CHEM ! tracer mixing ratio tendencies from the parameterized convection

    CHARACTER(LEN=100)          :: AER_CHEM_MECH

    REAL,    DIMENSION(mxp,myp) :: CONPRR

    REAL,    DIMENSION(mxp,myp) ::  aot500  ,temp2m  ,sfc_press &
                                   ,sflux_r ,sflux_t ,topt      &
                                   ,xland   ,dx2d    ,water_bud &
                                   ,col_sat
    INTEGER, DIMENSION(mxp,myp) :: kpbl,do_this_column
    INTEGER, DIMENSION(mzp)     :: flip
    INTEGER :: k,i,j,iens,ispc
    !- for convective transport
    INTEGER, DIMENSION(mxp,myp,maxiens) ::  &
               ierr4d                       &
              ,jmin4d                       &
              ,klcl4d                       &
              ,k224d                        &
              ,kbcon4d                      &
              ,ktop4d                       &
              ,kstabi4d                     &
              ,kstabm4d

    REAL,DIMENSION(mxp,myp,maxiens)     ::  &
               cprr4d                       &
              ,xmb4d                        &
              ,edt4d                        &
              ,pwav4d                       &
              ,sigma4d
    REAL,DIMENSION(mxp,mzp,myp,maxiens) ::  &
               entr5d                       &
              ,pcup5d                       &
              ,up_massentr5d                &
              ,up_massdetr5d                &
              ,dd_massentr5d                &
              ,dd_massdetr5d                &
              ,zup5d                        &
              ,zdn5d                        &
              ,prup5d                       &
              ,prdn5d                       &
              ,clwup5d                      &
              ,tup5d                        &
              ,conv_cld_fr5d

    !-----------local var in GEOS-5 data structure
    REAL,  DIMENSION(mxp, myp, mzp) :: DZ, AIR_DEN
    REAL,  DIMENSION(mxp, myp, mzp) :: ZLO_N,PLO_N,PK_N,MASS_N
    REAL,  DIMENSION(mxp,myp,0:mzp) :: ZLE_N
    INTEGER :: status,alloc_stat ,wantrank=-99999
    REAL    :: tem1,src_cnvtr,snk_cnvtr,tau_cp
    INTEGER, PARAMETER :: itest=1!3
    REAL :: RL, RI, disp_factor,x1,x2
    INTEGER :: mtp
    real :: dtqw,qsatup

   WQT_DC = 0.0
   CNV_MFC = 0.0
   CNV_MF0 = 0.0
   CNV_PRC3 = 0.0
   CNV_MFD = 0.0
   CNV_DQCDT = 0.0
   CNV_UPDF = 0.0
   CNV_CVW = 0.0
   CNV_QC = 0.0
   ENTLAM = 0.0
   CNPCPRATE = 0.0
   LIGHTN_DENS = 0.0
   REVSU = 0.0
   PRFIL = 0.0

   DQDT_GF = 0.0
   DTDT_GF = 0.0
   DUDT_GF = 0.0
   DVDT_GF = 0.0

   SIGMA_DEEP = 0.0
   SIGMA_MID = 0.0
   MUPDP = 0.0
   MDNDP = 0.0
   MUPSH = 0.0
   MUPMD = 0.0
   MFDP = 0.0
   MFSH = 0.0
   MFMD = 0.0
   ERRDP = 0.0
   ERRSH = 0.0
   ERRMD = 0.0
   AA0 = 0.0
   AA1 = 0.0
   AA2 = 0.0
   AA3 = 0.0
   AA1_BL = 0.0
   AA1_CIN = 0.0
   TAU_BL = 0.0
   TAU_EC = 0.0

    !-- some initialization
    do_this_column =0
    ierr4d         =0
    jmin4d         =0
    klcl4d         =0
    k224d          =0
    kbcon4d        =0
    ktop4d         =0
    kstabi4d       =0
    kstabm4d       =0
    xmb4d          =0.
    cprr4d         =0.
    edt4d          =0.
    pwav4d         =0.
    sigma4d        =0.
    pcup5d         =0.
    entr5d         =0.                      
    up_massentr5d  =0.
    up_massdetr5d  =0.
    dd_massentr5d  =0.
    dd_massdetr5d  =0.
    zup5d          =0.
    zdn5d          =0.
    prup5d         =0.
    prdn5d         =0.
    clwup5d        =0.
    tup5d          =0.
    conv_cld_fr5d  =0.
    SRC_NI         =0.
    SRC_NL         =0.
    SRC_T          =0.
    SRC_Q          =0.
    SRC_CI         =0.
    SRC_U          =0.
    SRC_V          =0.
    CNPCPRATE      =0.
    SUB_MPQI       =0.
    SUB_MPQL       =0.
    SUB_MPCF       =0.
    LIGHTN_DENS    =0.
    SRC_BUOY       =0.
    REVSU_GF       =0.
    PRFIL_GF       =0.
    !-
    !---temporary settings for debugging purposes
    !- special setting for SCM runs
    if(mxp==1 .and. myp==1 .and. maxval(T2m) < 1.e-6) return

    !- special setting for SCM runs
    if(mxp>1 .and. myp>1) wrtgrads = .false.
    !call mpi_comm_rank(MPI_COMM_WORLD,WHOAMI_ALL,status)
    !----
    if(wrtgrads) call alloc_grads_arr(1,mzp,1,jl)
    !--------------------------------------------------------

    !- time counter
    ntimes = ntimes + 1
    !if(ntimes==1 .and. WHOAMI_ALL == 0) print *,'==> Applying GF convection scheme '
    mynum = -999
    call set_index_loops( ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          mxp,myp,mzp                   )

    !- define the vector "flip" to invert the z-axis orientation
    call flipz(flip,mzp)
    !
    mtp = size(CNV_Tracers)
    IF(.not.allocated(SRC_CHEM)) THEN
       allocate(SRC_CHEM(mtp, mzp, mxp, myp),stat=alloc_stat) !- tendency from convection
       IF(alloc_stat==0) SRC_CHEM=0.0
    ENDIF
    !
    !- 2-d input data
    aot500  (:,:) = 0.1  ! #
    !- as moist is called before surface, at the 1st time step all arrays
    !- from surface are zero
    if(maxval(T2m) < 1.e-6) then
       temp2m   (:,:) = T1(:,:,mzp) ! Kelvin
    else
       temp2m   (:,:) = T2M(:,:) ! or TA(:,:) ! Kelvin
    endif
    !- moisture flux from sfc
    sflux_r  (:,:) = EVAP(:,:) ! kg m-2 s-1
    !- sensible heat flux (sh) comes in W m-2, below it is converted to K m s-1
    sflux_t  (:,:) = SH  (:,:) /(1004. * PLE(:,:,mzp)/(287.04*T1(:,:,mzp)*(1.+0.608*Q1(:,:,mzp)))) ! K m s-1
    !- topography height  (m)
    topt     (:,:) = PHIS(:,:)/MAPL_GRAV
    !- land/ocean fraction: land if < 1 ,ocean if = 1
    xland    (:,:) = 1.0-FRLAND(:,:)
    !
    !- grid length for the scale awareness
    dx2d(:,:) = sqrt(AREA(:,:)) ! meters
    !- special setting for SCM runs
    if(mxp==1 .and. myp==1) dx2d = 100000.

    !-pbl heigth index
    DO j=1,myp
     DO i=1,mxp
       if (KPBLIN(i,j) /= 0.0) then
          kpbl(i,j) = flip(NINT(KPBLIN(i,j)))
       else
          kpbl(i,j) = 1
       endif
    ENDDO
    ENDDO
    !
    !- 3-d input data
    !- any var with index "1" (and w and pk) are already updated with dynamics
    !  tendencies and everything else (from physics) that was called before moist
    !
    IF(trim(GF_ENV_SETTING)=='CURRENT') then
     !- 1st setting: enviromental state is the one already modified by dyn + physics
     DZ      = -( ZLE(:,:,1:mzp) - ZLE(:,:,0:mzp-1) )
     AIR_DEN = PLO/(287.04*T1*(1.+0.608*Q1))
     DO j=1,myp
      DO i=1,mxp
       DO k=1,mzp
        temp  (k,i,j) = T1    (i,j,flip(k))
        press (k,i,j) = PLO   (i,j,flip(k))!Pa
        rvap  (k,i,j) = Q1    (i,j,flip(k))!
        up    (k,i,j) = U1    (i,j,flip(k))! already @ A-grid (m/s)
        vp    (k,i,j) = V1    (i,j,flip(k))! already @ A-grid (m/s)
        wp    (k,i,j) = W1    (i,j,flip(k))! m/s
        zt3d  (k,i,j) = ZLO   (i,j,flip(k))! mid -layer level
        zm3d  (k,i,j) = ZLE   (i,j,flip(k))! edge-layer level
        dm3d  (k,i,j) = MASS  (i,j,flip(k))
        khloc (k,i,j) = KH    (i,j,flip(k))
        curr_rvap(k,i,j) = Q1 (i,j,flip(k))!current rvap

        if (ENTRVERSION==0) then
          ! eq 6 of https://doi.org/10.1029/2021JD034881
          ec3d(k,i,j) = 0.71*max(0.5,W1(i,j,flip(k)))**(-1.17) * max(0.1,BYNCY(i,j,flip(k)))**(-0.36)
        else
          ec3d(k,i,j) = 1.0
        endif
        entr3d(i,j,flip(k)) = ec3d(k,i,j) 

        mp_ice(lsmp,k,i,j) = QILS  (i,j,flip(k))
        mp_liq(lsmp,k,i,j) = QLLS  (i,j,flip(k))
        mp_cf (lsmp,k,i,j) = CLLS  (i,j,flip(k))
        mp_ice(cnmp,k,i,j) = QICN  (i,j,flip(k))
        mp_liq(cnmp,k,i,j) = QLCN  (i,j,flip(k))
        mp_cf (cnmp,k,i,j) = CLCN  (i,j,flip(k))

       ENDDO
      ENDDO
     ENDDO
     !- sfc pressure (Pa)
     sfc_press(:,:) = PLE(:,:,mzp)
     !- Grid and sub-grid scale forcings for convection
     DO j=1,myp
      DO i=1,mxp
       DO k=1,mzp
         gsf_t (k,i,j) = 0.
         gsf_q (k,i,j) = 0.
        sgsf_t (k,i,j) = 0.
        sgsf_q (k,i,j) = 0.
        advf_t (k,i,j) = 0.
       ENDDO
      ENDDO
     ENDDO
    ELSEIF(trim(GF_ENV_SETTING)=='DYNAMICS') then
     !-2nd setting: environmental state is that one before any tendency
     !- is applied (i.e, at begin of each time step).
     !- Get back the model state, heights and others variables at time N
     !- (or at the beggining of current time step)
     !- In physics, the state vars (T,U,V,PLE) are untouched and represent the
     !- model state after dynamics phase 1. But, "Q" is modified by physics, so
     !- depending on what was called before this subroutine, "Q" may be already
     !- changed from what it was just after dynamics phase 1. To solve this issue,
     !- "Q" just after dynamics is saved in the var named "QV_DYN_IN" in "GEOS_AgcmGridComp.F90".
     MASS_N = ( PLE_DYN_IN(:,:,1:mzp)-PLE_DYN_IN(:,:,0:mzp-1) )*(1./MAPL_GRAV)
     call get_vars(mzp,mxp,myp,QV_DYN_IN,T_DYN_IN,PLE_DYN_IN,ZLE_N,ZLO_N,PLO_N,PK_N)
     DZ      = -( ZLE_N(:,:,1:mzp) - ZLE_N(:,:,0:mzp-1) )
     AIR_DEN = PLO_N/(287.04*T_DYN_IN*(1.+0.608*QV_DYN_IN))
     DO j=1,myp
      DO i=1,mxp
       DO k=1,mzp
        temp       (k,i,j) = T_DYN_IN   (i,j,flip(k))! (K)
        press      (k,i,j) = PLO_N      (i,j,flip(k))! (Pa) @ mid-layer level
        rvap       (k,i,j) = QV_DYN_IN  (i,j,flip(k))! water vapor mix ratio
        up         (k,i,j) = U_DYN_IN   (i,j,flip(k))! already @ A-grid (m/s)
        vp         (k,i,j) = V_DYN_IN   (i,j,flip(k))! already @ A-grid (m/s)
        wp         (k,i,j) = W1         (i,j,flip(k))! (m/s)
        zt3d       (k,i,j) = ZLO_N      (i,j,flip(k))! mid -layer level (m)
        zm3d       (k,i,j) = ZLE_N      (i,j,flip(k))! edge-layer level (m)
        dm3d       (k,i,j) = MASS_N     (i,j,flip(k))
        khloc      (k,i,j) = KH         (i,j,flip(k))
        curr_rvap  (k,i,j) = Q1         (i,j,flip(k)) ! current rvap (dyn+phys)

        if (ENTRVERSION==0) then
          ! eq 6 of https://doi.org/10.1029/2021JD034881
          ec3d(k,i,j) = 0.71*max(0.5,W1(i,j,flip(k)))**(-1.17) * max(0.1,BYNCY(i,j,flip(k)))**(-0.36)
        else
          ec3d(k,i,j) = 1.0
        endif
        entr3d(i,j,flip(k)) = ec3d(k,i,j) 

        mp_ice(lsmp,k,i,j) = QILS       (i,j,flip(k))
        mp_liq(lsmp,k,i,j) = QLLS       (i,j,flip(k))
        mp_cf (lsmp,k,i,j) = CLLS       (i,j,flip(k))
        mp_ice(cnmp,k,i,j) = QICN       (i,j,flip(k))
        mp_liq(cnmp,k,i,j) = QLCN       (i,j,flip(k))
        mp_cf (cnmp,k,i,j) = CLCN       (i,j,flip(k))
       ENDDO
      ENDDO
     ENDDO
     !- sfc pressure (Pa)
     sfc_press(:,:) = PLE_DYN_IN(:,:,mzp)
     !- Grid and sub-grid scale forcings for convection
     DO j=1,myp
      DO i=1,mxp
       DO k=1,mzp
         gsf_t (k,i,j) = DTDTDYN (i,j,flip(k)) + RADSW(i,j,flip(k)) + RADLW(i,j,flip(k))
         gsf_q (k,i,j) = DQVDTDYN(i,j,flip(k))
        sgsf_t (k,i,j) = DTDT_BL (i,j,flip(k))
        sgsf_q (k,i,j) = DQDT_BL (i,j,flip(k))
        advf_t (k,i,j) = DTDTDYN (i,j,flip(k))
       ENDDO
      ENDDO
     ENDDO
     !
    ELSE
     stop 'unknown GF_ENV_SETTING at convpar_gf2020.F90'
    ENDIF

    IF(CONVECTION_TRACER==1) THEN
      DO j=1,myp
       DO i=1,mxp
        !--- saturation CWV
         col_sat(i,j) = TPWI(i,j)/(1.e-6+TPWI_star(i,j))
         col_sat(i,j) = min(1.,max(0.,col_sat(i,j)))
        !--- temporary to hold only CWV in mm
        ! col_sat(i,j) = TPWI(i,j)
        !
        DO k=1,mzp
           buoy_exc(k,i,j)=CNV_TR (i,j,flip(k))
        ENDDO
       ENDDO
      ENDDO
      !print*,"buoy_exc1=",maxval(buoy_exc),minval(buoy_exc)
    ELSE
      buoy_exc = 0.0
    ENDIF

    !- call the driver routine to apply the parameterization
    CALL GF2020_DRV(mxp,myp,mzp,mtp,nmp         &
                     ,ims,ime, jms,jme, kms,kme   &
                     ,its,ite, jts,jte, kts,kte   &
                     ,flip        &
                     ,mynum       &
                     ,dt_moist    &
                     ,dx2d        &
                     ,stochastic_sig &
                     ,zm3d        &
                     ,zt3d        &
                     ,dm3d        &
                     !--- sfc inputs
                     ,lons        &
                     ,lats        &
                     ,aot500      &
                     ,temp2m      &
                     ,sflux_r     &
                     ,sflux_t     &
                     ,topt        &
                     ,xland       &
                     ,sfc_press   &
                     ,kpbl        &
                     ,CNVFRC      &
                     ,SRFTYPE     &
                     !--- atmos state
                     ,col_sat     &
                     ,up          &
                     ,vp          &
                     ,wp          &
                     ,ec3d        &
                     ,temp        &
                     ,press       &
                     ,rvap        &
                     ,mp_ice      &
                     ,mp_liq      &
                     ,mp_cf       &
                     ,curr_rvap   &
                     !---- forcings---
                     ,buoy_exc    &
                     ,gsf_t       &
                     ,gsf_q       &
                     ,advf_t      &
                     ,sgsf_t      &
                     ,sgsf_q      &
                     !---- output ----
                     ,CONPRR      &
                     ,LIGHTN_DENS &
                     ,SRC_T       &
                     ,SRC_Q       &
                     ,SRC_CI      &
                     ,SRC_U       &
                     ,SRC_V       &
                     ,SUB_MPQI    &
                     ,SUB_MPQL    &
                     ,SUB_MPCF    &
                     ,SRC_BUOY    &
                     ,SRC_CHEM    &
                     ,REVSU_GF    &
                     ,PRFIL_GF    &
                     !
                     !
                     ,do_this_column&
                     ,ierr4d       &
                     ,jmin4d       &
                     ,klcl4d       &
                     ,k224d        &
                     ,kbcon4d      &
                     ,ktop4d       &
                     ,kstabi4d     &
                     ,kstabm4d     &
                     ,cprr4d       &
                     ,xmb4d        &
                     ,edt4d        &
                     ,pwav4d       &
                     ,sigma4d      &
                     ,pcup5d       &
                     ,entr5d       &
                     ,up_massentr5d&
                     ,up_massdetr5d&
                     ,dd_massentr5d&
                     ,dd_massdetr5d&
                     ,zup5d        &
                     ,zdn5d        &
                     ,prup5d       &
                     ,prdn5d       &
                     ,clwup5d      &
                     ,tup5d        &
                     ,conv_cld_fr5d&
                     !-- for debug/diagnostic
                     ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC)


 IF(FEED_3DMODEL)THEN
      !-- update GEOS-5 model state with the feedback from cumulus convection
      !- to include the tendencies from the convection,  update the vars th1,q1,v1 and u1
      DO j=1,myp
        DO i=1,mxp
          IF(do_this_column(i,j) == 0) CYCLE
          !- conv precip rate: mm/s = kg m-2 s-1
          CNPCPRATE (i,j) =  CONPRR (i,j)
          IF(ITEST==0) CNPCPRATE(i,j) =  0.

          !--- sublimation/evaporation tendencies (kg/kg/s)
          DO k=1,mzp
          !--- sublimation/evaporation tendencies (kg/kg/s)
               REVSU (i,j,k) = REVSU_GF(flip(k),i,j)
          !--- preciptation fluxes (kg/kg/s)
               PRFIL (i,j,k) = PRFIL_GF(flip(k),i,j)
          ENDDO

        ENDDO
      ENDDO
!-----
      IF(USE_MOMENTUM_TRANSP > 0) THEN
        DO j=1,myp
          DO i=1,mxp
            IF(do_this_column(i,j) == 0) CYCLE
          ENDDO
        ENDDO
      ENDIF


      IF(APPLY_SUB_MP == 1) THEN
        DO j=1,myp
          DO i=1,mxp
            IF(do_this_column(i,j) == 0) CYCLE
          ENDDO
        ENDDO
      ENDIF

      IF(USE_TRACER_TRANSP==1) THEN

        DO j=1,myp
          DO i=1,mxp
            IF(do_this_column(i,j) == 0) CYCLE
            DO k=1,mzp
              !- update tracer mass mixing ratios
              DO ispc=1,mtp

                 CNV_Tracers(ispc)%Q(i,j,k) = CNV_Tracers(ispc)%Q(i,j,k) + DT_moist * SRC_CHEM(ispc,flip(k),i,j)

                 !-- final check for negative tracer mass mixing ratio
                 CNV_Tracers(ispc)%Q(i,j,k) = max(CNV_Tracers(ispc)%Q(i,j,k), mintracer)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      entr_dp = MAPL_UNDEF
      entr_md = MAPL_UNDEF
      entr_sh = MAPL_UNDEF

      DO IENS=1, maxiens
        if(icumulus_gf(IENS) == ON) then
           DO j=1,myp
            DO i=1,mxp
             if(ierr4d(i,j,IENS) .ne. 0) cycle
              DO k=mzp,flip(ktop4d(i,j,IENS))-1,-1

                !- Export entrainment rates used by GF
                if (IENS==DEEP) entr_dp(i,j,k) = entr5d(i,flip(k),j,IENS)
                if (IENS==MID ) entr_md(i,j,k) = entr5d(i,flip(k),j,IENS)
                if (IENS==SHAL) entr_sh(i,j,k) = entr5d(i,flip(k),j,IENS)

                !- special treatment for CNV_DQCDT: 'convective_condensate_source',  UNITS     = 'kg m-2 s-1',
                !- SRC_CI contains contributions from deep, shallow,... . So, do not need to be accumulated over  CNV_DQCDT
                !- note the SRC_CI has different array structure (k,i,j) _not_ (i,j,k)
                CNV_DQCDT(i,j,k)=  SRC_CI(flip(k),i,j) * DZ(i,j,k) * AIR_DEN(i,j,k) !units: kg[w]/(kg[air] s) * m * kg[air]/m^3 = kg[w]/(m^2 s)

                !
                !-'detraining_mass_flux',UNITS     = 'kg m-2 s-1',
                CNV_MFD (i,j,k) = CNV_MFD (i,j,k) + up_massdetr5d(i,flip(k),j,IENS)


                !-'cloud_base_mass_flux',    units = 'kg m-2 s-1',
                CNV_MF0 (i,j,k) = CNV_MF0 (i,j,k) + ( zup5d(i,flip(k),j,IENS) )

                !-convective mass flux [kg m^{-2} s^{-1}] - with downdrafts
           !!!  CNV_MFC (i,j,k) = CNV_MFC (i,j,k) + ( zup5d(i,flip(k),j,IENS) + &
           !!!                                        zdn5d(i,flip(k),j,IENS)*edt4d(i,j,IENS) )

                !---only updraft
                CNV_MFC (i,j,k) = CNV_MFC (i,j,k) + ( zup5d(i,flip(k),j,IENS) )

                ! deep convective total water flux. assumes .033 fractional area.
                qsatup = MAPL_EQsat(tup5d(i,flip(k),j,IENS),press(flip(k),i,j),dtqw) + clwup5d(i,flip(k),j,IENS)/0.033
                WQT_DC (i,j,k) = WQT_DC (i,j,k) + zup5d(i,flip(k),j,IENS)*(qsatup - rvap(flip(k),i,j))

                if(zup5d(i,flip(k),j,IENS) > 1.0e-6) then
                   !-'entrainment parameter',  UNITS     ='m-1',
                   ENTLAM  (i,j,k) =  ENTLAM (i,j,k) + (up_massentr5d(i,flip(k),j,IENS)/(DZ(i,j,k)*zup5d(i,flip(k),j,IENS)))

                   !-'updraft_vertical_velocity',            UNITS     = 'hPa s-1',
                   CNV_CVW (i,j,k) = -0.2 ! hPa/s =>  4 m/s
                endif

                !-'grid_mean_convective_condensate', UNITS     ='kg kg-1'
                CNV_QC  (i,j,k) = CNV_QC (i,j,k) + clwup5d(i,flip(k),j,IENS)
                !
                !
                !~ !----------------------------------------------------------------------------------------------------
                !- not using progno-cloud to calculate the precip from the convective column
                !- if CNV_PRC3 will be send to progno-cloud, set CNPCPRATE = zero
                !-'convective_precipitation_from_GF',UNITS     = 'kg m-2 s-1',
                !- JAN/17/2017 : the units above are wrong. The correct are kg[precip water]/kg[air]
                CNV_PRC3(i,j,k) = CNV_PRC3 (i,j,k) +                 (prup5d(i,flip(k),j,IENS) + &
                                                     edt4d(i,j,IENS)* prdn5d(i,flip(k),j,IENS) ) &
                                                     * DT_moist/(DZ(i,j,k)*AIR_DEN(i,j,k))

                !-'updraft_areal_fraction',
                if(zup5d(i,flip(k),j,IENS) > 1.0e-6) CNV_UPDF(i,j,k) = 0.033

             ENDDO
            ENDDO
           ENDDO
        endif
      ENDDO
  ENDIF


  !
  !--- cold pool/"convection tracer"
  IF(CONVECTION_TRACER==1) THEN
           DO j=1,myp
            DO i=1,mxp

              tau_cp=FRLAND(i,j)*tau_land_cp + (1.0-FRLAND(i,j))*tau_ocea_cp

              DO k=1,mzp

                !- sink term (exp decay 1h)
                snk_cnvtr =  dt_moist * abs(CNV_TR (i,j,k))/tau_cp

                !- downdraft convective mass flux [kg m^{-2} s^{-1}]
                ! iens =?
                !src_cnvtr =  edt4d(i,j,iens)*zdn5d(i,flip(k),j,iens)

                !- downdraft detrainment mass flux [kg m^{-2} s^{-1}]
                ! iens =?
                !src_cnvtr =  edt4d(i,j,iens)*dd_massdetr5d(i,flip(k),j,iens)


                !- source term
                !- downdraft detrainment of buoyancy [ J/kg s^{-1}]
                !- negative sign => source for updraft lifting
                src_cnvtr = - dt_moist * min(0.,SRC_BUOY(flip(k),i,j))

                !- 'continuity' equation = ADV + SRC - SINK
                CNV_TR (i,j,k) = CNV_TR (i,j,k)  + src_cnvtr  - snk_cnvtr

              ENDDO
            ENDDO
           ENDDO
          !print*,"buoy_exc2=",maxval(SRC_BUOY),minval(SRC_BUOY)
   ENDIF

!
  CNV_TOPP_DP = MAPL_UNDEF
  CNV_TOPP_MD = MAPL_UNDEF
  CNV_TOPP_SH = MAPL_UNDEF
  IF(maxval(icumulus_gf)>0) then
      DO IENS=1, maxiens
        if(icumulus_gf(IENS) == ON .and. IENS== DEEP) then
         DO j=1,myp
           DO i=1,mxp
            if(ierr4d(i,j,DEEP) /= 0) cycle
            CNV_TOPP_DP(i,j)     = press(ktop4d(i,j,DEEP),i,j)
            MFDP      (i,j)      =   xmb4d(i,j,DEEP)
            SIGMA_DEEP(i,j)      = sigma4d(i,j,DEEP)
            MUPDP     (i,j,1:mzp)=zup5d(i,flip(1):flip(mzp):-1,j,DEEP)
            MDNDP     (i,j,1:mzp)=zdn5d(i,flip(1):flip(mzp):-1,j,DEEP) * edt4d(i,j,IENS)
           ENDDO
         ENDDO
        elseif(icumulus_gf(IENS) == ON .and. IENS== SHAL) then
         DO j=1,myp
           DO i=1,mxp
            if(ierr4d(i,j,SHAL) /= 0) cycle
            CNV_TOPP_SH(i,j)=press(ktop4d(i,j,SHAL),i,j)
            MFSH (i,j)      =xmb4d(i,j,SHAL)
            MUPSH(i,j,1:mzp)=zup5d(i,flip(1):flip(mzp):-1,j,SHAL)
           ENDDO
         ENDDO
        elseif(icumulus_gf(IENS) == ON .and. IENS== MID) then
         DO j=1,myp
           DO i=1,mxp
             if(ierr4d(i,j,MID) /= 0) cycle
             CNV_TOPP_MD(i,j)     = press(ktop4d(i,j,MID),i,j)
             MFMD      (i,j)      = xmb4d(i,j,MID)
             SIGMA_MID (i,j)      = sigma4d(i,j,MID )
             MUPMD     (i,j,1:mzp)=zup5d(i,flip(1):flip(mzp):-1,j,MID)
           ENDDO
         ENDDO
        endif
      ENDDO
      !- for output purposes, ierr=0 (convection is ON) will be changed to 1
      where(ierr4d==0)ierr4d=1
      where(ierr4d>1)ierr4d=0
      DO j=1,myp
       DO i=1,mxp
        !- Tendencies
         DTDT_GF(i,j,1:mzp)=SRC_T(flip(1):flip(mzp):-1,i,j)
         DQDT_GF(i,j,1:mzp)=SRC_Q(flip(1):flip(mzp):-1,i,j)
         DUDT_GF(i,j,1:mzp)=SRC_U(flip(1):flip(mzp):-1,i,j)
         DVDT_GF(i,j,1:mzp)=SRC_V(flip(1):flip(mzp):-1,i,j)
        !- Error codes
         ERRDP(i,j)=float(ierr4d(i,j,DEEP))
         ERRSH(i,j)=float(ierr4d(i,j,SHAL))
         ERRMD(i,j)=float(ierr4d(i,j,MID ))
       ENDDO
      ENDDO
  ENDIF
  !
  if( allocated(src_chem))  deallocate(src_chem,stat=alloc_stat) !tendency   from convection

  !- for debugging purposes only
  if(wrtgrads) call alloc_grads_arr(1,mzp,2,jl)

  END SUBROUTINE GF2020_INTERFACE
!---------------------------------------------------------------------------------------------------

  SUBROUTINE GF2020_DRV(mxp,myp,mzp,mtp,nmp         &
              ,ims,ime, jms,jme, kms,kme              &
              ,its,ite, jts,jte, kts,kte              &
              ,flip                                   &
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
              ,cnvfrc                &
              ,srftype               &

              ,col_sat               &
              ,u                     &
              ,v                     &
              ,w                     &
              ,entr_c                &
              ,temp                  &
              ,press                 &
              ,rvap                  &
              ,mp_ice                     &
              ,mp_liq                &
              ,mp_cf                 &
              ,curr_rvap             &

              !---- forcings---
              ,buoy_exc              &
              ,rthften               &! gsf_t
              ,rqvften               &! gsf_q
              ,rth_advten            &!advf_t
              ,rthblten              &!sgsf_t
              ,rqvblten              &!sgsf_q
              !---- output ----
              ,conprr                &
              ,lightn_dens           &
              ,rthcuten              &
              ,rqvcuten              &
              ,rqccuten              &
              ,rucuten               &
              ,rvcuten               &
              ,sub_mpqi              &
              ,sub_mpql              &
              ,sub_mpcf              &
              ,rbuoycuten            &
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
              ,entr5d                &
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
              ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC)

   IMPLICIT NONE
  !include "mpif.h"
  !------------------------------------------------------------------------
   INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          mynum,mzp,mxp,myp,mtp,nmp

   REAL,    INTENT(IN) :: DT

   INTEGER, INTENT(IN), dimension(mzp) :: flip

   REAL,    DIMENSION(kts:kte,its:ite,jts:jte), INTENT(IN)  ::       &
                                                          zm,        &
                                                          zt,        &
                                                          u,         &
                                                          v,         &
                                                          w,         &
                                                          entr_c,    &
                                                          rvap,      &
                                                          temp,      &
                                                          press,     &
                                                          dm,        &
                                                          curr_rvap, &
                                                          buoy_exc

   INTEGER, DIMENSION(its:ite,jts:jte), INTENT(IN) :: kpbl
   REAL,    DIMENSION(its:ite,jts:jte), INTENT(IN) :: cnvfrc, srftype

   !-- intent (in)
   REAL,    DIMENSION(its:ite,jts:jte) ::             topt ,aot500 ,temp2m ,sfc_press &
                                                     ,sflux_r ,sflux_t                &
                                                     ,xland,lons,lats,dx2d,col_sat    &
                                                     ,stochastic_sig

   REAL,    DIMENSION(kts:kte,its:ite,jts:jte), INTENT(IN) :: rthften    &
                                                             ,rqvften    &
                                                             ,rth_advten &
                                                             ,rthblten   &
                                                             ,rqvblten

   REAL,    DIMENSION(its:ite,jts:jte),         INTENT(OUT) ::   CONPRR,LIGHTN_DENS

   REAL,    DIMENSION(kts:kte,its:ite,jts:jte), INTENT(OUT) :: &
                                                    rthcuten   &
                                                   ,rqvcuten   &
                                                   ,rqccuten   &
                                                   ,rucuten    &
                                                   ,rvcuten    &
                                                   ,rbuoycuten &
                                                   ,revsu_gf   &
                                                   ,prfil_gf


   REAL,    DIMENSION(nmp,kts:kte,its:ite,jts:jte), INTENT(IN)  :: &
                                                    mp_ice     &
                                                   ,mp_liq     &
                                                   ,mp_cf

   REAL,    DIMENSION(nmp,kts:kte,its:ite,jts:jte), INTENT(OUT) :: &
                                                    sub_mpqi   &
                                                   ,sub_mpql   &
                                                   ,sub_mpcf

   !-***** rchemcuten uses the GF data structure (ispc,k,i,j) *********
   REAL,    DIMENSION(mtp,kts:kte,its:ite,jts:jte), INTENT(OUT)  :: rchemcuten

   INTEGER, DIMENSION(its:ite,jts:jte), INTENT(INOUT) :: do_this_column

   REAL,DIMENSION(mxp,mzp,myp,maxiens), INTENT(OUT)   :: entr5d
 
!- for convective transport and cloud/radiation (OUT)
   INTEGER,DIMENSION(mxp,myp,maxiens), INTENT(INOUT) ::  &
               ierr4d                    &
              ,jmin4d                    &
              ,klcl4d                    &
              ,k224d                     &
              ,kbcon4d                   &
              ,ktop4d                    &
              ,kstabi4d                  &
              ,kstabm4d

   REAL,DIMENSION(mxp,myp,maxiens), INTENT(INOUT) ::     &
               cprr4d                    &
              ,xmb4d                     &
              ,edt4d                     &
              ,pwav4d                    &
              ,sigma4d
   REAL,DIMENSION(mxp,mzp,myp,maxiens), INTENT(INOUT) :: &
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

!----------------------------------------------------------------------
! LOCAL VARS

! basic environmental input includes
! omega (omeg), windspeed (us,vs,ws), and a flag (aaeq) to turn off
! convection for this call only and at that particular gridpoint
!

    REAL,   DIMENSION (kts:kte,its:ite,jts:jte)  :: Tpert_h,Tpert_v

    REAL,   DIMENSION (its:ite,kts:kte) ::                                      &
                               zo,temp_old,qv_old,PO,US,VS,WS,rhoi,phil         &
                              ,temp_new_dp,qv_new_dp,temp_new_sh,qv_new_sh,z2d  &
                              ,tkeg,rcpg,dhdt,temp_new_md,qv_new_md             &
                              ,temp_new_BL,qv_new_BL,dm2d,temp_tendqv,qv_curr   &
                              ,buoy_exc2d,revsu_gf_2d,prfil_gf_2d               &
                              ,temp_new,qv_new,Tpert_2d,temp_new_adv,qv_new_adv


    REAL,   DIMENSION (its:ite,kts:kte,maxiens) ::                              &
                               outt,outq,outqc,outu,outv,outbuoy,outnliq,outnice

    REAL,   DIMENSION (mtp,its:ite,kts:kte)         :: se_chem
    REAL,   DIMENSION (mtp,its:ite,kts:kte,maxiens) :: out_chem

    REAL,   DIMENSION (nmp,its:ite,kts:kte)         :: mpqi,mpql,mpcf
    REAL,   DIMENSION (nmp,its:ite,kts:kte,maxiens) :: outmpqi,outmpql,outmpcf

    REAL,   DIMENSION (its:ite)   :: ter11, xlandi,pbl,zws,ccn,psur &
                                    ,ztexec,zqexec,h_sfc_flux,le_sfc_flux,tsur&
                                    ,xlons,xlats,fixout_qv,cum_ztexec,cum_zqexec


    REAL,   DIMENSION (its:ite,kts:kte,1:ens4)      ::  omeg

    REAL,   DIMENSION (kts:kte) :: min_tend,distance
    INTEGER,DIMENSION (its:ite) :: kpbli,last_ierr

    INTEGER :: i,j,k,kr,n,itf,jtf,ktf,ispc,zmax,status

    REAL :: dp,dq,exner, dtdt,PTEN,PQEN,PAPH,ZRHO,PAHFS,PQHFL,ZKHVFL,PGEOH
    REAL :: fixouts,dt_inv

    REAL,   DIMENSION(mxp,myp,-1:5) :: dummy_precip
    INTEGER :: imemory,irun,jlx,kk,kss,plume,ii_plume


!----------------------------------------------------------------------
    !-do not change this
    itf=ite
    ktf=kte-1
    jtf=jte
    int_time = int_time + dt
    WHOAMI_ALL=mynum
!----------------------------------------------------------------------
   IF(abs(C1) > 0.) USE_C1D = .TRUE.

   !-- for the moisture adv trigger
   IF(ADV_TRIGGER == 2) then
      call prepare_temp_pertubations(kts,kte,ktf,its,ite,itf,jts,jte,jtf,dt,xland,topt,zm &
                                 ,temp,rqvften,rthblten,rthften,Tpert_h,Tpert_v &
                                 ,AA1_CIN,AA1_BL)

      !-- for output only
      AA0(:,:) = Tpert_h(2, :,:); AA1(:,:) = Tpert_h(5,:,:)
      AA2(:,:) = Tpert_h(10,:,:); AA3(:,:) = Tpert_v(2,:,:)
   ENDIF

!-- big loop over j dimension
   DO j = jts,jtf
     JCOL = J

     !-- initialization
     DO i= its,itf
        ztexec   (i) = 0.0
        zqexec   (i) = 0.0
        last_ierr(i) = -999
        fixout_qv(i) = 1.0
        !
        conprr     (i,j) = 0.0
        lightn_dens(i,j) = 0.0
        !--- (i,k)
        revsu_gf_2d (i,:) = 0.0
        prfil_gf_2d (i,:) = 0.0
        Tpert_2d    (i,:) = 0.0
        !
        temp_tendqv(i,:) = 0.0
        !- tendencies (w/ maxiens)
        outt   (i,:,:)=0.0
        outu   (i,:,:)=0.0
        outv   (i,:,:)=0.0
        outq   (i,:,:)=0.0
        outqc  (i,:,:)=0.0
        outnice(i,:,:)=0.0
        outnliq(i,:,:)=0.0
        outbuoy(i,:,:)=0.0
     ENDDO

     IF(APPLY_SUB_MP == 1) THEN
       DO i= its,itf
        !- tendencies (w/ nmp and maxiens)
        outmpqi(:,i,:,:)=0.0
        outmpql(:,i,:,:)=0.0
        outmpcf(:,i,:,:)=0.0
     ENDDO
     ENDIF
     DO i= its,itf
         omeg (i,:,:)=0.0
     ENDDO
     !-- for the moisture adv trigger (Ma and Tan, AR 2009)
     IF(ADV_TRIGGER == 2) THEN
       do k=kts,kte
          do i=its,itf
             Tpert_2d(i,k)= Tpert_h(k,i,j) + Tpert_v(k,i,j) !- don't use "kr" here
       enddo;enddo
     ENDIF

     IF(USE_TRACER_TRANSP==1) THEN
         out_chem = 0.0
     ENDIF
     !
     IF(autoconv == 2) THEN
      DO I= its,itf
             ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
      ENDDO
     ELSE
      DO I= its,itf
             ccn(i) = 100.
      ENDDO
     ENDIF

     DO i=its,itf

         xlandi(i) = xland(i,j)!flag < 1 para land
                               !flag  =1 para water
         psur  (i) = sfc_press(i,j)*1.e-2 ! mbar
         tsur  (i) = temp2m(i,j)
         ter11 (i) = max(0.,topt(i,j))
         kpbli (i) = kpbl(i,j)
         xlons (i) = lons(i,j)*180./3.14159
         xlats (i) = lats(i,j)*180./3.14159
     ENDDO

     DO k=kts,ktf
       DO i=its,itf
         kr=k   !+1   !<<<< only kr=k
!
         !- heigths, current pressure, temp and water vapor mix ratio
         zo      (i,k)  = zt(kr,i,j)+topt(i,j)
         po      (i,k)  = press(kr,i,j)*1.e-2 !mbar
         temp_old(i,k)  = temp(kr,i,j)

         qv_old  (i,k)  = rvap     (kr,i,j) ! @ begin of the timestep
         qv_curr (i,k)  = curr_rvap(kr,i,j) ! current (after dynamics + physical processes called before GF)

         !- air density, TKE and cloud liq water mixing ratio
         rhoi    (i,k)  = 1.e2*po (i,k)/( 287.04*temp_old(i,k)*(1.+0.608*qv_old(i,k)))
         tkeg    (i,k)  = tkmin
         rcpg    (i,k)  = 0.

         !- wind velocities
         us      (i,k)  =  u (kr,i,j)
         vs      (i,k)  =  v (kr,i,j)
         ws      (i,k)  =  w (kr,i,j)
         dm2d    (i,k)  =  dm(kr,i,j)
         omeg    (i,k,:)= -g*rhoi(i,k)*w(kr,i,j)
         !-buoyancy excess
         buoy_exc2d(i,k)= buoy_exc(kr,i,j)
         !- temp/water vapor modified only by advection
         temp_new_ADV(i,k)= temp_old(i,k)  +  (rth_advten(kr,i,j) )*dt
           qv_new_ADV(i,k)=   qv_old(i,k)  +  (rqvften   (kr,i,j) )*dt
       ENDDO
     ENDDO

     IF(APPLY_SUB_MP == 1) THEN
       DO k=kts,ktf
          DO i=its,itf
             kr=k   !+1   !<<<< only kr=k
             !- microphysics ice and liq mixing ratio, and cloud fraction of the host model
             !- (only subsidence is applied)
             mpqi   (:,i,k) = mp_ice  (:,kr,i,j) ! kg/kg
             mpql   (:,i,k) = mp_liq  (:,kr,i,j) ! kg/kg
             mpcf   (:,i,k) = mp_cf   (:,kr,i,j) ! 1
          ENDDO
       ENDDO
     ENDIF
     IF(USE_TRACER_TRANSP==1) THEN
      DO k=kts,kte
        DO i=its,itf
          kr=k !+1
          !- atmos composition
          DO ispc=1,mtp
             se_chem(ispc,i,k) = max(CNV_Tracers(ispc)%Q(i,j,flip(kr)),mintracer)
          ENDDO
        ENDDO
      ENDDO
     ENDIF
     !- pbl  (i) = depth of pbl layer (m)
     !- kpbli(i) = index of zo(i,k)
     !call get_zi_gf(j,its,ite,kts,kte,its,itf,ktf,ierrs,kpbli,pbl,&
     !             tkeg,rcpg,zo,ter11,tkmin)
     DO i=its,itf
         pbl  (i)  = zo(i,kpbli(i)) - topt(i,j)
         !print*,"PBL=",kpbli(i),zo(i,kpbli(i)),topt(i,j),pbl  (i)
     ENDDO

     !- get execess T and Q for source air parcels
     do i=its,itf
       pten = temp_old(i,1)
       pqen = qv_old  (i,1)
       paph = 100.*psur(i)
       zrho = paph/(287.04*(temp_old(i,1)*(1.+0.608*qv_old(i,1))))
       !- sensible and latent sfc fluxes for the heat-engine closure
       h_sfc_flux (i)=zrho*cp *sflux_t(i,j)!W/m^2
       le_sfc_flux(i)=zrho*xlv*sflux_r(i,j)!W/m^2
       !
       !- local le and h fluxes for W*
       pahfs=-sflux_t(i,j) *zrho*1004.64  !W/m^2
       pqhfl=-sflux_r(i,j)                !kg/m^2/s
       !- buoyancy flux (h+le)
       zkhvfl= (pahfs/1004.64+0.608*pten*pqhfl)/zrho ! K m s-1
       !- depth of 1st model layer
       !- (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
       pgeoh =  2.*( zo(i,1)-topt(i,j) )*g ! m+2 s-2
       !-convective-scale velocity w*
       !- in the future, change 0.001 by ustar^3
       zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten) ! m+3 s-3

       if(zws(i) > tiny(pgeoh)) then
         !-convective-scale velocity w*
         zws(i) = 1.2*zws(i)**.3333
         !- temperature excess
         ztexec(i)     = max(0.,-1.5*pahfs/(zrho*zws(i)*1004.64)) ! K
         !print*,"exce1=",pahfs,zrho,ztexec(i),zws(i),pgeoh,zo(i,1),topt(i,j);call flush(6)
         !- moisture  excess
         zqexec(i)     = max(0.,-1.5*pqhfl/(zrho*zws(i)))        !kg kg-1
       endif   ! zws > 0
       !
       !- zws for shallow convection closure (Grant 2001)
       !- depth of the pbl
       pgeoh = pbl(i)*g
       !-convective-scale velocity W* (m/s)
       zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)
       zws(i) = 1.2*zws(i)**.3333
     enddo

!
!------ CALL CUMULUS PARAMETERIZATION
!

     DO ii_plume = 1, maxiens
       if (SH_MD_DP) then
         if(ii_plume == 1) plume = shal
         if(ii_plume == 2) plume = mid
         if(ii_plume == 3) plume = deep
       else
         if(ii_plume == 1) plume = shal
         if(ii_plume == 2) plume = deep
         if(ii_plume == 3) plume = mid
       endif

       hei_down_land  =  cum_hei_down_land  (plume)
       hei_down_ocean =  cum_hei_down_ocean (plume)
       hei_updf_land  =  cum_hei_updf_land  (plume)
       hei_updf_ocean =  cum_hei_updf_ocean (plume)
       min_edt_land   =  cum_min_edt_land   (plume)
       min_edt_ocean  =  cum_min_edt_ocean  (plume)
       max_edt_land   =  cum_max_edt_land   (plume)
       max_edt_ocean  =  cum_max_edt_ocean  (plume)
       fadj_massflx   =  cum_fadj_massflx   (plume)
       use_excess     =  cum_use_excess     (plume)
       ave_layer      =  cum_ave_layer      (plume)
       !print*,"plume=",plume,shal,mid,deep

       IF(icumulus_gf(plume) /= ON ) cycle

       !-- set minimum/max for excess of T and Q
       if(use_excess == 0) then
         cum_ztexec(:)= 0.        ; cum_zqexec(:)= 0.
       elseif (use_excess == 1) then
         cum_ztexec(:)= ztexec(:) ; cum_zqexec(:)= zqexec(:)
       elseif (use_excess == 2) then
         do i=its,itf
           cum_zqexec(i)=min(5.e-4, max(1.e-4,zqexec(i)))! kg kg^-1
           cum_ztexec(i)=min(0.5,   max(0.2  ,ztexec(i)))! Kelvin
         enddo
       else
         do i=its,itf
           if(xlandi(i) > 0.98) then ! ocean
             cum_zqexec(i)=min(5.e-4, max(1.e-4,zqexec(i)))! kg kg^-1
             cum_ztexec(i)=min(0.5,   max(0.2  ,ztexec(i)))! Kelvin
           else                      ! land
             cum_ztexec(:)= ztexec(:) ; cum_zqexec(:)= zqexec(:)
           endif
         enddo
       endif

       do i=its,itf
         do k=kts,ktf
           kr=k!+1 <<<<

           temp_new(i,k)=temp_old(i,k) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
           qv_new  (i,k)=  qv_old(i,k) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
           qv_new  (i,k)= max(smallerqv,qv_new  (i,k))

           !- moist static energy
           dhdt(i,k)= cp *(rthblten(kr,i,j)+rthften(kr,i,j)) + &
                      xlv*(rqvblten(kr,i,j)+rqvften(kr,i,j))

           !- temp/water vapor modified only by bl processes
           temp_new_BL(i,k)= temp_old(i,k)  +  (rthblten(kr,i,j) )*dt
           qv_new_BL  (i,k)= qv_old  (i,k)  +  (rqvblten(kr,i,j) )*dt

           !- temp/water vapor modified only by advection
           !temp_new_ADV(i,k)= temp_old(i,k)  +  (rth_advten(kr,i,j) )*dt
           !  qv_new_ADV(i,k)=   qv_old(i,k)  +  (rqvften   (kr,i,j) )*dt

         enddo
       enddo

       CALL CUP_gf(its,ite,kts,kte, itf,ktf, mtp, nmp &
                  ,cumulus_type  (plume)            &
                  ,closure_choice(plume)            &
                  ,cum_entr_rate (plume)            &
                  ,cum_use_excess(plume)            &
                  !- input data
                  ,dx2d          (:,j)              &
                  ,stochastic_sig(:,j)              &
                  ,col_sat       (:,j)              &
                  ,dt                               &
                  ,kpbli                            &
                  ,cum_ztexec                       &
                  ,cum_zqexec                       &
                  ,ccn                              &
                  ,rhoi                             &
                  ,omeg                             &
                  ,temp_old                         &
                  ,qv_old                           &
                  ,ter11                            &
                  , h_sfc_flux                      &
                  ,le_sfc_flux                      &
                  ,xlons                            &
                  ,xlats                            &
                  ,xlandi                           &
                  ,cnvfrc(:,j)                      &
                  ,srftype(:,j)                     &
                  ,temp_new                         &
                  ,qv_new                           &
                  ,temp_new_BL                      &
                  ,qv_new_BL                        &
                  ,temp_new_ADV                     &
                  ,qv_new_ADV                       &
                  ,zo                               &
                  ,po                               &
                  ,tsur                             &
                  ,psur                             &
                  ,us                               &
                  ,vs                               &
                  ,ws                               &
                  ,entr_c(:,:,j)                    &
                  ,dm2d                             &
                  ,se_chem                          &
                  ,zws                              &
                  ,dhdt                             &
                  ,buoy_exc2d                       &
                  ,mpqi                             &
                  ,mpql                             &
                  ,mpcf                             &
                  ,last_ierr            (:)         &
                  !output data
                  ,outt                 (:,:,plume) &
                  ,outq                 (:,:,plume) &
                  ,outqc                (:,:,plume) &
                  ,outu                 (:,:,plume) &
                  ,outv                 (:,:,plume) &
                  ,outnliq              (:,:,plume) &
                  ,outnice              (:,:,plume) &
                  ,outbuoy              (:,:,plume) &
                  ,outmpqi            (:,:,:,plume) &
                  ,outmpql            (:,:,:,plume) &
                  ,outmpcf            (:,:,:,plume) &
                  ,out_chem           (:,:,:,plume) &
                  !- for convective transport
                  ,ierr4d               (:,j,plume) &
                  ,jmin4d               (:,j,plume) &
                  ,klcl4d               (:,j,plume) &
                  ,k224d                (:,j,plume) &
                  ,kbcon4d              (:,j,plume) &
                  ,ktop4d               (:,j,plume) &
                  ,kstabi4d             (:,j,plume) &
                  ,kstabm4d             (:,j,plume) &
                  ,cprr4d               (:,j,plume) &
                  ,xmb4d                (:,j,plume) &
                  ,edt4d                (:,j,plume) &
                  ,pwav4d               (:,j,plume) &
                  ,sigma4d              (:,j,plume) &
                  ,pcup5d             (:,:,j,plume) &
                  ,entr5d             (:,:,j,plume) & 
                  ,up_massentr5d      (:,:,j,plume) &
                  ,up_massdetr5d      (:,:,j,plume) &
                  ,dd_massentr5d      (:,:,j,plume) &
                  ,dd_massdetr5d      (:,:,j,plume) &
                  ,zup5d              (:,:,j,plume) &
                  ,zdn5d              (:,:,j,plume) &
                  ,prup5d             (:,:,j,plume) &
                  ,prdn5d             (:,:,j,plume) &
                  ,clwup5d            (:,:,j,plume) &
                  ,tup5d              (:,:,j,plume) &
                  ,conv_cld_fr5d      (:,:,j,plume) &
                  !-- for debug/diag
                  ,AA0(:,j),AA1(:,j),AA2(:,j),AA3(:,j),AA1_BL(:,j),AA1_CIN(:,j),TAU_BL(:,j),TAU_EC(:,j) &
                  !-- for diag
                  ,lightn_dens  (:,j)               &
                  ,revsu_gf_2d                      &
                  ,prfil_gf_2d                      &
                  ,Tpert_2d                         &
                  )
        !--- accumulate precip for each plume
        DO i=its,itf
           CONPRR(i,j)= CONPRR(i,j)+cprr4d(i,j,plume)
        ENDDO
        ! Save ierr from this plume
        ! if (plume /= SHAL) then
        ! DO i=its,itf
        !   last_ierr(i) = max(last_ierr(i), ierr4d(i,j,plume))
        ! ENDDO
        ! endif
     ENDDO !- plume

     !--- reset ierr4d to value different of zero in case the correspondent
     !--- plume (shalllow, congestus, deep) was not actually used
     DO n=1,maxiens
       if(icumulus_gf(n) == OFF ) ierr4d (:,j,n) = -99
     ENDDO

     DO i=its,itf
        do_this_column(i,j) = 0
loop1:  do n=1,maxiens
          if(ierr4d (i,j,n) == 0 ) then
            do_this_column(i,j) = 1
            exit loop1
          endif
        enddo loop1
     ENDDO
     !----------- check for negative water vapor mix ratio
     DO i=its,itf
          if(do_this_column(i,j) == 0) CYCLE
          do k = kts,ktf
            temp_tendqv(i,k)= outq (i,k,shal) + outq (i,k,deep) + outq (i,k,mid )
          enddo

          do k = kts,ktf
            distance(k)= qv_curr(i,k) + temp_tendqv(i,k) * dt
          enddo

          if(minval(distance(kts:ktf)) < 0.) then
                zmax   =  MINLOC(distance(kts:ktf),1)

                if( abs(temp_tendqv(i,zmax) * dt) <  mintracer) then
                  fixout_qv(i)= 0.999999
                 !fixout_qv(i)= 0.
                else
                  fixout_qv(i)= ( (smallerQV - qv_curr(i,zmax))) / (temp_tendqv(i,zmax) *dt)
                endif
                fixout_qv(i)=max(0.,min(fixout_qv(i),1.))
          endif
          !--- apply to convective precip
          CONPRR(i,j)= CONPRR(i,j) * fixout_qv(i)
     ENDDO

     !------------ feedback
     !-- deep + shallow + mid convection
     DO i = its,itf
      if(do_this_column(i,j) == 0) CYCLE
      DO k = kts,kte
             kr=k!+1
             !- feedback the tendencies from convection
             RTHCUTEN (kr,i,j)= (outt (i,k,shal) + outt (i,k,deep) + outt (i,k,mid )) *fixout_qv(i)

             RQVCUTEN (kr,i,j)= (outq (i,k,shal) + outq (i,k,deep) + outq (i,k,mid )) *fixout_qv(i)

             RQCCUTEN (kr,i,j)= (outqc(i,k,shal) + outqc(i,k,deep) + outqc(i,k,mid )) *fixout_qv(i)

             REVSU_GF (kr,i,j)= revsu_gf_2d(i,k)*fixout_qv(i) !-- already contains deep and mid amounts.

            !---these arrays are only for the deep plume mode
             PRFIL_GF (kr,i,j)= prfil_gf_2d (i,k)*fixout_qv(i) !-- ice/liq prec flux of the deep plume
      ENDDO
     ENDDO
     IF(USE_MOMENTUM_TRANSP > 0) THEN
      DO i = its,itf
       if(do_this_column(i,j) == 0) CYCLE
       DO k = kts,kte
             kr=k!+1
             RUCUTEN (kr,i,j)= (outU(i,k,deep)+outU(i,k,mid)+outU(i,k,shal)) *fixout_qv(i)
             RVCUTEN (kr,i,j)= (outV(i,k,deep)+outV(i,k,mid)+outV(i,k,shal)) *fixout_qv(i)
       ENDDO
      ENDDO
     ENDIF


     IF(APPLY_SUB_MP == 1) THEN
      DO i = its,itf
       if(do_this_column(i,j) == 0) CYCLE
       DO k = kts,kte
             kr=k!+1
             SUB_MPQL (:,kr,i,j)= (outmpql(:,i,k,deep)+outmpql(:,i,k,mid)+outmpql(:,i,k,shal)) *fixout_qv(i)
             SUB_MPQI (:,kr,i,j)= (outmpqi(:,i,k,deep)+outmpqi(:,i,k,mid)+outmpqi(:,i,k,shal)) *fixout_qv(i)
             SUB_MPCF (:,kr,i,j)= (outmpcf(:,i,k,deep)+outmpcf(:,i,k,mid)+outmpcf(:,i,k,shal)) *fixout_qv(i)
       ENDDO
      ENDDO
     ENDIF

     IF(USE_TRACER_TRANSP==1) THEN
      DO i = its,itf
       if(do_this_column(i,j) == 0) CYCLE
       DO k = kts,kte
         kr=k!+1
         RCHEMCUTEN (:,kr,i,j)= (out_CHEM(:,i,k,deep) +out_CHEM(:,i,k,mid)+out_CHEM(:,i,k,shal)) *fixout_qv(i)
       ENDDO
      ENDDO

      !- constrain positivity for tracers
      DO i = its,itf
         if(do_this_column(i,j) == 0) CYCLE

         do ispc=1,mtp

           do k=kts,ktf
              distance(k)= se_chem(ispc,i,k) + RCHEMCUTEN(ispc,k,i,j)* dt
           enddo

           !-- fixer for mass of tracer
           IF(minval(distance(kts:ktf)) < 0.) THEN
                zmax   =  MINLOC(distance(kts:ktf),1)

                if( abs(RCHEMCUTEN(ispc,zmax,i,j)*dt) <  mintracer) then
                  fixouts= 0.999999
                 !fixouts= 0.
                else
                  fixouts=  ( (mintracer - se_chem(ispc,i,zmax))) / (RCHEMCUTEN(ispc,zmax,i,j)*dt)
                endif
                if(fixouts > 1. .or. fixouts <0.)fixouts=0.

                RCHEMCUTEN(ispc,kts:ktf,i,j)=fixouts*RCHEMCUTEN(ispc,kts:ktf,i,j)
           ENDIF
         enddo
      ENDDO
     ENDIF

     IF(CONVECTION_TRACER==1) THEN
      DO i = its,itf
       if(do_this_column(i,j) == 0) CYCLE
       DO k = kts,kte
             kr=k!+1
             RBUOYCUTEN (kr,i,j)= (outbuoy(i,k,deep)+outbuoy(i,k,mid)+outbuoy(i,k,shal)) *fixout_qv(i)
       ENDDO
      ENDDO
     ENDIF

!-----memory
     AA3(:,j)=cprr4d(:,j,deep) *fixout_qv(:)
     AA2(:,j)=cprr4d(:,j,mid)  *fixout_qv(:)
     !IF(USE_MEMORY > -1 ) THEN
     ! DO i=its,itf
     !    AA1_CIN(i,j)=cprr4d(i,j,deep) * fixout_qv(i)
     !    AA2(i,j)=cprr4d(i,j,deep) *fixout_qv(i)
     !    IF(AA2(i,j) > 0. .and. AA3(i,j)>0.) then
     !      AA1_CIN(i,j)=AA1_CIN(i,j)+1.
     !    ELSE
     !      AA1_CIN(i,j)=0.
     !    ENDIF
     ! ENDDO
     !ENDIF
!-----memory


     ENDDO

    END SUBROUTINE GF2020_DRV
!---------------------------------------------------------------------------------------------------

    SUBROUTINE CUP_gf(its,ite,kts,kte ,itf,ktf, mtp, nmp &
                     ,cumulus           &
                     ,ichoice           &
                     ,entr_rate_plume   &
                     ,use_excess        &
                     !input data
                     ,dx                &
                     ,stochastic_sig    &
                     ,col_sat           &
                     ,dtime             &
                     ,kpbl              &
                     ,ztexec            &
                     ,zqexec            &
                     ,ccn               &
                     ,rho               &
                     ,omeg              &
                     ,t                 &
                     ,q                 &
                     ,z1                &
                     , h_sfc_flux       &
                     ,le_sfc_flux       &
                     ,xlons             &
                     ,xlats             &
                     ,xland             &
                     ,cnvfrc            &
                     ,srftype           &
                     ,tn                &
                     ,qo                &
                     ,tn_bl             &
                     ,qo_bl             &
                     ,tn_adv            &
                     ,qo_adv            &
                     ,zo                &
                     ,po                &
                     ,tsur              &
                     ,psur              &
                     ,us                &
                     ,vs                &
                     ,ws                &
                     ,entr_c            &
                     ,dm2d              &
                     ,se_chem           &
                     ,zws               &
                     ,dhdt              &
                     ,buoy_exc          &
                     ,mpqi              &
                     ,mpql              &
                     ,mpcf              &
                     ,last_ierr         &
                     !output data
                     ,outt              &
                     ,outq              &
                     ,outqc             &
                     ,outu              &
                     ,outv              &
                     ,outnliq           &
                     ,outnice           &
                     ,outbuoy           &
                     ,outmpqi           &
                     ,outmpql           &
                     ,outmpcf           &
                     ,out_chem          &
                     !- for convective transport
                     ,ierr              &
                     ,jmin              &
                     ,klcl              &
                     ,k22               &
                     ,kbcon             &
                     ,ktop              &
                     ,kstabi            &
                     ,kstabm            &
                     ,pre               &
                     ,xmb               &
                     ,edto              &
                     ,pwavo             &
                     ,sig               &
                     ,po_cup            &
                     ,entr_rate         &
                     ,up_massentro      &
                     ,up_massdetro      &
                     ,dd_massentro      &
                     ,dd_massdetro      &
                     ,zuo               &
                     ,zdo               &
                     ,pwo               &
                     ,pwdo              &
                     ,qrco              &
                     ,tup               &
                     ,clfrac            &
                     !- for convective transport-end
                     !- for debug/diag
                     ,AA0_,AA1_,AA2_,AA3_,AA1_BL_,AA1_CIN_,TAU_BL_,TAU_EC_   &
                     ,lightn_dens        &
                      ,revsu_gf          &
                     ,prfil_gf          &
                     ,Tpert             &
                     )
     IMPLICIT NONE

     !-local settings
     LOGICAL, PARAMETER:: USE_LCL       =.FALSE.
     LOGICAL, PARAMETER:: USE_INV_LAYERS=.TRUE.

     CHARACTER*(*),INTENT(IN) :: cumulus
     INTEGER      ,INTENT(IN) :: itf,ktf,its,ite,kts,kte,ichoice,use_excess,mtp, nmp
     INTEGER      ,INTENT(IN),     DIMENSION (its:ite) ::   last_ierr
     INTEGER      ,INTENT(INOUT),  DIMENSION (its:ite) ::   kpbl
  !
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     REAL,    DIMENSION (its:ite,kts:kte) ,INTENT (INOUT)   ::       &
        outu,outv,outt,outq,outqc,outbuoy,revsu_gf,prfil_gf          &
       ,outnliq ,outnice

     REAL,    DIMENSION (its:ite)         ,INTENT (OUT  )   ::       &
        pre,sig,lightn_dens

  !
  ! basic environmental input includes
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     REAL,    DIMENSION (kts:kte,its:ite)       ,INTENT (IN      )    ::   &
        entr_c
     REAL,    DIMENSION (its:ite,kts:kte)       ,INTENT (INOUT   )    ::   &
        dhdt,rho,t,po,us,vs,ws,tn,dm2d,buoy_exc,tn_bl,tn_adv
     REAL,    DIMENSION (its:ite,kts:kte,1:ens4),INTENT (INOUT)    ::   &
        omeg
     REAL,    DIMENSION (its:ite,kts:kte)       ,INTENT (INOUT)    ::   &
         q,qo,Tpert,qo_bl,qo_adv
     REAL,    DIMENSION (its:ite)               ,INTENT (INOUT   )    ::   &
        ccn,Z1,PSUR,xland,xlons,xlats, h_sfc_flux,le_sfc_flux,tsur,dx

     REAL,    DIMENSION (its:ite)               ,INTENT (IN   )    ::   &
        col_sat,&
        cnvfrc,srftype,&
        stochastic_sig
     REAL,    DIMENSION (its:ite)               ,INTENT (INOUT)    ::   &
        zws,ztexec,zqexec
     REAL                                       ,INTENT (IN   )    ::   &
        dtime,entr_rate_plume
     REAL,    DIMENSION (nmp,its:ite,kts:kte)   ,INTENT (INOUT)   ::   &
         mpqi,mpql,mpcf
     REAL,    DIMENSION (nmp,its:ite,kts:kte)   ,INTENT (INOUT)   ::   &
         outmpqi,outmpql,outmpcf

  !
  ! local ensemble dependent variables in this routine
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,1:ensdim)        ::                   &
        xf_ens,pr_ens
  !
  !*******the following are your basic environmental
  !          variables. They carry a "_cup" if they are
  !          on model cloud levels (staggered). They carry
  !          an "o"-ending (z becomes zo), if they are the forced
  !          variables. They are preceded by x (z becomes xz)
  !          to indicate modification by some typ of cloud
  !
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
  !
  !
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur    = surface pressure
  ! z1      = terrain elevation
  ! pr_ens  = precipitation ensemble
  ! xf_ens  = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg    = omega from large scale model
  ! mconv   = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! ichoice     = flag if only want one closure (usually set to zero!)
  ! dby  = buoancy term
  ! ktop = cloud top (output)
  ! xmb  = total base mass flux
  ! hc   = cloud moist static energy
  ! hkb  = moist static energy at originating level
     logical :: keep_going
     character*128 :: ierrc(its:ite)

     real,    dimension (its:ite,kts:kte) ::                           &
        mentrd_rate                                                    &
       , he, hes, qes, z,heo,heso,qeso,zo, zu,zd                       &
       ,xhe,xhes,xqes,xz,xt,xq                                         &
       , qes_cup, q_cup, he_cup, hes_cup, z_cup, p_cup, gamma_cup, t_cup &
       ,qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,       gammao_cup,tn_cup &
       ,xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup                        &
       ,xt_cup,hcot,evap_bcb                                           &
       ,dby,hc,clw_all                                                 &
       ,dbyo,qco,qrcdo,hcdo,qcdo,dbydo,hco                             &
       ,xdby,xzu,xzd,xhc,cupclw,pwo_eff,                               &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,dellah,dellaq,dellat,dellaqc,dsubq,dsubh,dellabuoy,      &
        u_cup,v_cup,uc,vc,ucd,vcd,dellu,dellv,                          &
        up_massentr,up_massdetr,dd_massentr,dd_massdetr,                &
        subten_H,subten_Q,subten_T

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects
  ! edt     = epsilon

     real,    dimension (its:ite) ::                            &
       edt,edtx,aa1,aa0,xaa0,hkb,                               &
       hkbo,xhkb,qkb, pwevo,bu,bud,cap_max,xland1,              &
       cap_max_increment,psum,psumh,sigd,mconv

     integer,    dimension (its:ite) ::                         &
       kzdown,kdet,kb, ierr2,ierr3,kbmax
     integer :: iloop,nall,iedt,nens,nens3,ki,I,K,KK,iresult,nvar,nvarbegin
     integer :: jprnt,k1,k2,kbegzu,kdefi,kfinalzu,kstart,jmini,imid,k_free_trop

     real :: day,dz,dzo,radius,entrd_rate,                         &
             zcutdown,depth_min,zkbmax,z_detr,zktop,               &
             massfld,dh,trash,frh,xlamdd,radiusd,frhd,effec_entrain
     real :: detdo1,detdo2,entdo,dp,subin,detdo,entup,             &
             detup,subdown,entdoj,entupk,detupk,totmas
     real :: tot_time_hr,beta,wmeanx,env_mf,env_mf_p,env_mf_m
     real :: dts,denom,denomU

     real,    dimension (its:ite,1:maxens3) ::  xff_mid
     real,    dimension (kts:kte)   :: dummy1,dummy2
     integer :: iversion,bl=1,fa=2,step
     real :: umean,T_star

     real, dimension (its:ite)         :: aa0_bl,aa1_bl,tau_bl,tau_ecmwf,wmean,aa1_fa,aa1_tmp,hkbo_x &
                                         ,aa2,aa3,cin0,cin1,edtmax,edtmin,aa1_lift,aa_tmp,aa_ini,aa_adv,daa_adv_dt
     real, dimension (its:ite,kts:kte) :: tn_x, qo_x, qeso_x, heo_x, heso_x,zo_cup_x &
                                         ,qeso_cup_x,qo_cup_x, heo_cup_x,heso_cup_x,po_cup_x&
                                         ,gammao_cup_x,tn_cup_x,hco_x,DBYo_x,u_cup_x,v_cup_x

     real, dimension (its:ite,kts:kte) :: xhe_x,xhes_x,xt_x,xq_x,xqes_x, &
                                          xqes_cup_x,xq_cup_x,xhe_cup_x,xhes_cup_x,gamma_cup_x,xt_cup_x
     real, dimension (its:ite)         ::xaa0_x,xk_x

     real, dimension(its:ite) :: xf_dicycle,mbdt
     real :: C_up, E_dn,G_rain,trash2, pgc,bl2dp,trash3
     character(len=2) :: cty
!
     real   :: dsubh_aver,dellah_aver,x_add,cap_max_inc,tlll,plll,rlll,tlcl,plcl,dzlcl,zlll
     integer:: start_k22,start_level(its:ite)
     real,    dimension (its:ite,kts:kte) ::  dtempdz
     integer, dimension (its:ite,kts:kte) ::  k_inv_layers
     integer :: ipr=0,jpr=0,fase

     real,    dimension (its:ite,kts:kte) ::  vvel2d,tempco,tempcdo
     real,    dimension (its:ite        ) ::  vvel1d, x_add_buoy

     real,    dimension (its:ite,kts:kte) ::  p_liq_ice,melting_layer,melting
     real,    dimension (its:ite,kts:kte) ::  c1d

     real,    dimension (its:ite) :: lambau_dn,lambau_dp
     real,    dimension (its:ite,kts:kte) ::  up_massentru,up_massdetru,&
                                              dd_massentru,dd_massdetru

     real,    dimension (its:ite) :: q_wetbulb,t_wetbulb
     integer :: i_wb=0

     real,    dimension (its:ite) :: p_cwv_ave,cape, depth_neg_buoy,frh_bcon &
                                    ,check_sig,random

     real,    dimension (its:ite,kts:kte) :: prec_flx,evap_flx,qrr

     real,    dimension (its:ite,9) :: xff_shal

    !- atmos composition arrays
     real, dimension (mtp,its:ite,kts:kte),intent (inout) ::   se_chem
     real, dimension (mtp,its:ite,kts:kte),intent (inout) ::   out_chem

     !-locals
     integer :: ispc,kmp,istep,lstep
     real, dimension (mtp,its:ite,kts:kte) ::  se_cup_chem,sc_up_chem,sc_dn_chem,pw_up_chem,pw_dn_chem
     real, dimension (mtp,its:ite)         ::  tot_pw_up_chem,tot_pw_dn_chem
     real, dimension (mtp,kts:kte)         ::  trcflx_in,sub_tend,ddtr,ddtr_upd,zenv_diff,fp_mtp,fm_mtp
     real, dimension (mtp)                 ::  min_tend_chem,dummy_chem,delS_up,delS_dn,env_sub,outchem1
     real, dimension (nmp,its:ite,kts:kte) ::  dellampqi,dellampql,dellampcf
     real, dimension (kts:kte)             ::  aa,bb,cc,ddu,ddv,ddh,ddq,fp,fm
     real, dimension (nmp,kts:kte)         ::  dd
     real, dimension (its:ite,kts:kte)     ::  massflx,zenv,rho_hydr,alpha_H,alpha_Q
     real                                  ::  evap_(mtp),wetdep_(mtp),trash_(mtp),trash2_(mtp) &
                                              ,massi,massf,dtime_max,evap,wetdep,residu_(mtp)

    !----------------------------------------------------------------------
     integer, dimension (its:ite)      ,intent (inout)  :: &
              ierr              &
             ,jmin              &
             ,klcl              &
             ,k22               &
             ,kbcon             &
             ,ktop              &
             ,kstabi            &
             ,kstabm

     real,  dimension (its:ite)        ,intent (inout)  :: &
              xmb               &
             ,edto              &
             ,pwavo
     real,  dimension (its:ite,kts:kte),intent (inout)  :: &
              entr_rate         &
             ,po_cup            &
             ,up_massentro      &
             ,up_massdetro      &
             ,dd_massentro      &
             ,dd_massdetro      &
             ,zuo               &
             ,zdo               &
             ,pwo               &
             ,pwdo              &
             ,qrco              &
             ,tup               &
             ,clfrac
    !----------------------------------------------------------------------
    !-- debug/diag
     real,  dimension (its:ite)        ,intent (inout)  :: &
             aa0_,aa1_,aa2_,aa3_,aa1_bl_,aa1_cin_,tau_bl_,tau_ec_
     real,  dimension (its:ite,kts:kte) :: dtdt,dqdt
     real :: s1,s2,q1,q2,rzenv,factor,CWV,entr_threshold,resten_H,resten_Q,resten_T
     integer :: status
     real :: alp0,beta1,beta2,dp_p,dp_m,delt1,delt2,delt_Tvv,wkf,ckf,wkflcl,rcount
     real,    dimension (kts:kte,8)       ::  tend2d
     real,    dimension (8)               ::  tend1d
    !----------------------------------------------------------------------

     pre = 0.0
     sig = 0.0
     sigd = 0.0
     lightn_dens = 0.0

!
!--- maximum depth (mb) of capping inversion (larger cap = no convection)
!
      IF(ZERO_DIFF==1 .or. MOIST_TRIGGER==0) THEN
       if(trim(cumulus) == 'deep'   ) then ; cap_max_inc=20. ; endif ! cap_maxs=50.
       if(trim(cumulus) == 'mid'    ) then ; cap_max_inc=10. ; endif ! cap_maxs=50.
       if(trim(cumulus) == 'shallow') then ; cap_max_inc=25. ; endif ! cap_maxs=50.
      ELSE
       if(trim(cumulus) == 'deep'   ) then ; cap_max_inc=90. ; endif !--- test cap_maxs=10.  ; cap_max_inc=50
       if(trim(cumulus) == 'mid'    ) then ; cap_max_inc=90. ; endif !--- test cap_maxs=10.  ; cap_max_inc=50
       if(trim(cumulus) == 'shallow') then ; cap_max_inc=10. ; endif !--- test cap_maxs=25.  ; cap_max_inc=50
      ENDIF
!
!--- lambda_U parameter for momentum transport
!
      if(trim(cumulus) == 'deep'   ) then;lambau_dp (:) = lambau_deep;  lambau_dn (:) = lambau_shdn ;endif
      if(trim(cumulus) == 'mid'    ) then;lambau_dp (:) = lambau_shdn;  lambau_dn (:) = lambau_shdn ;endif
      if(trim(cumulus) == 'shallow') then;lambau_dp (:) = lambau_shdn;  lambau_dn (:) = lambau_shdn ;endif

      if(pgcon .ne. 0.) then
       lambau_dp (:) = 0.; lambau_dn (:) = 0.
      endif

      do i=its,itf
        kbmax  (i) = 1
        kstabm (i) = ktf-1
        ierr2  (i) = 0
        ierr3  (i) = 0
        xland1 (i) = xland(i) ! 1.
        cap_max(i) = cap_maxs
        ierrc  (i) = "ierrtxt"
        aa0    (i) = 0.0
        aa1    (i) = 0.0
        aa2    (i) = 0.0
        aa3    (i) = 0.0
        aa1_bl (i) = 0.0
        aa1_fa (i) = 0.0
        aa0_bl (i) = 0.0
        cin1   (i) = 0.0
        xk_x   (i) = 0.0
        edt    (i) = 0.0
        edto   (i) = 0.0
        tau_bl (i) = 0.0
        q_wetbulb (i) = 0.0
        t_wetbulb (i) = 0.0
        tau_ecmwf (i) = 0.0
        xf_dicycle(i) = 0.0
        x_add_buoy(i) = 0.0
        z     (i,:) = zo(i,:)
        xz    (i,:) = zo(i,:)
        hcdo  (i,:) = 0.0
        cupclw(i,:) = 0.0
        qrcdo (i,:) = 0.0
        hcot  (i,:) = 0.0
        c1d   (i,:) = 0.0
        xf_ens(i,:) = 0.0
        pr_ens(i,:) = 0.0
        evap_bcb(i,:) = 0.0
        cap_max_increment(i)=cap_max_inc
       enddo
!
!--- SCALE DEPENDENCE FACTOR (SIG), version new
!
      if( USE_SCALE_DEP == 0 ) sig(:)=1.
      if( USE_SCALE_DEP == 1 ) then
         if( trim(cumulus) == 'shallow') then
            sig(:)=1.
         else
            do i=its,itf
             sig(i) = 0.
             if (stochastic_sig(i) < 0.0) then
               ierr(i)=1
               ierrc(i)='scale_dep renders convection insignificant'
             endif
             if(ierr(i) /= 0) cycle
             sig(i)= sigma(dx(i))
             if (stochastic_sig(i) /= 1.0) then
               sig(i) = sig(i)**(stochastic_sig(i)*MAX(1.0,sig(i)))
             endif
             sig(i)= max(0.1,min(sig(i),1.))
             if(sig(i).le.0.1)then
               ierr(i)=1
               ierrc(i)='scale_dep renders convection insignificant'
             endif
             if(ierr(i) /= 0) cycle
            enddo
         endif
      endif
!
!---  create a real random number in the interval [-use_random_num, +use_random_num]
!
      if( trim(cumulus) == 'deep' .and. use_random_num > 1.e-6) then
         call gen_random(its,ite,use_random_num,random)
      else
         random = 0.0
      endif
!
!--- define entrainment/detrainment profiles for updrafts
!
      !- initial entrainment/detrainment
      do i=its,itf
          do k=kts,ktf
      entr_rate(i,k) = entr_c(k,i)*entr_rate_plume
      cd       (i,k) = entr_c(k,i)*entr_rate_plume
          enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!--  note : to make the evaporation stronger => increase "edtmin"
      DO i=its,itf
         if(xland(i) > 0.99 ) then !- over water
           edtmin(i)=MIN_EDT_OCEAN;  edtmax(i)=MAX_EDT_OCEAN
         else!- over land
           edtmin(i)=MIN_EDT_LAND;  edtmax(i)=MAX_EDT_LAND
         endif
      ENDDO
!
!--- minimum depth (m), clouds must have
!
      if(trim(cumulus) == 'deep'                               ) depth_min=1000.
      if(trim(cumulus) == 'mid'    ) depth_min=1000.
      if(trim(cumulus) == 'shallow') depth_min=500.
!
!--- max height(m) above ground where updraft air can originate
!
      if(trim(cumulus) == 'deep'                               ) zkbmax=4000.
      if(trim(cumulus) == 'mid'    ) zkbmax=3000.
      if(trim(cumulus) == 'shallow') zkbmax=2000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      if(trim(cumulus) == 'deep'   ) z_detr= depth_min*0.5
      if(trim(cumulus) == 'mid'    ) z_detr= depth_min*0.5
      if(trim(cumulus) == 'shallow') z_detr= depth_min*0.5

!
      !- mbdt ~ xmb * timescale
      do i=its,itf
          mbdt(i)= 0.1!*dtime*xmb_nm1(i)
         !mbdt(i)= 100.*(p_cup(i,kbcon(i))-p_cup(i,kbcon(i)+1))/(g*dtime)
         !mbdt(i)= 0.1*mbdt(i)
      enddo
!
!--- environmental conditions, FIRST HEIGHTS
!--- calculate moist static energy, heights, qes
!
      call cup_env(z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,ierr,-1,itf,ktf,its,ite, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, psur,ierr,-1,itf,ktf,its,ite, kts,kte)

!
!--- outputs a model sounding for the stand-alone code (part 1)
!
      IF(OUTPUT_SOUND == 1) THEN
       call SOUND(1,cumulus,int_time,dtime,ens4,itf,ktf,its,ite, kts,kte,xlats,xlons,jcol,whoami_all  &
                 ,z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,zo,qeso,heo,heso,tn,qo,us,vs ,omeg,xz     &
                 ,h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec,zqexec, xland      &
                 ,kpbl,k22,klcl,kbcon,ktop,aa0,aa1,sig,xaa0,hkb,xmb,pre,edto &
                 ,zo_cup,dhdt,rho,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv)
      ENDIF
!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,po,qes_cup,q_cup,he_cup, &
           us,vs,u_cup,v_cup,                                     &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,         &
           ierr,z1,itf,ktf,its,ite, kts,kte)

      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, heo_cup,   &
           us,vs,u_cup,v_cup,                                           &
           heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,tsur,          &
           ierr,z1,itf,ktf,its,ite, kts,kte)
!
!--- get air density at full layer (model levels) by hydrostatic balance (kg/m3)
!
      do i=its,itf
          rho_hydr(i,:) = 0.0
          if(ierr(i) /= 0)cycle
          do k=kts,ktf
             rho_hydr(i,k)=100.*(po_cup(i,k)-po_cup(i,k+1))/(zo_cup(i,k+1)-zo_cup(i,k))/g
             !print*,"rhohidr=",k,rho_hydr(i,k),po_cup(i,k+1),zo_cup(i,k+1)
          enddo
      enddo
!
!--- partition between liq/ice cloud contents
!
      call get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup,p_liq_ice,melting_layer,&
                                 itf,ktf,its,ite,kts,kte,cnvfrc,srftype,cumulus)
!
      do i=its,itf
        if(ierr(i).eq.0)then
          do k=kts,ktf
           if(zo_cup(i,k).gt.zkbmax+z1(i))then
             kbmax(i)=k ; go to 25
           endif
          enddo
 25       continue
!
!--- level where detrainment for downdraft starts
!
          do k=kts,ktf
            if(zo_cup(i,k).gt.z_detr+z1(i))then
              kdet(i)=k; go to 26
            endif
          enddo
 26       continue
        endif
      enddo
!
!--- determine level with highest moist static energy content - k22
!
      if(trim(cumulus) == 'shallow') then
        start_k22 = 1
      else
        start_k22 = 2
      endif
      k22(:)=kts
      DO i=its,itf
         if(ierr(i) /= 0) cycle
         k22(i)=maxloc(HEO_CUP(i,start_k22:kbmax(i)+1),1)+start_k22-1
         k22(i)=max(k22(i),start_k22)
            if(trim(cumulus) == 'shallow') then
                      k22(i)=min(2,k22(i))

                      if(K22(I).GT.KBMAX(i))then
                         ierr(i)=2
                         ierrc(i)="could not find k22"
                      endif

            else

                      if(k22(i) > kbmax(i))then
                        !- let's try k22=start_k22 for the cases k22>kbmax
                        k22(i)= start_k22
                        cycle
                      endif

            endif
      ENDDO
!
!
!-- get the pickup of ensemble ave prec, following Neelin et al 2009.
!
     call precip_cwv_factor(itf,ktf,its,ite,kts,kte,ierr,t,po,qo,po_cup,cumulus,p_cwv_ave)
!
!-- cold pool parameterization and convective memory
!
     IF(CONVECTION_TRACER == 1) THEN
        if(trim(cumulus) == 'deep') then
          if(USE_MEMORY == 0) then
             do i=its,itf
                if(ierr(i) /= 0) cycle
                x_add_buoy(i) = min(maxval(buoy_exc(i,kts:k22(i))),mx_buoy)
                !----------------------
                AA3_(i)=x_add_buoy(i) !--- temporary for output
                !----------------------
             enddo
          endif
          if(USE_MEMORY > 0 .and. USE_MEMORY <= 10) then
             do i=its,itf

               factor = float(USE_MEMORY)

               entr_rate(i,:) = entr_rate(i,:)*factor

               if(ierr(i) /= 0) cycle

               x_add_buoy(i) = min(maxval(buoy_exc(i,kts:k22(i))),mx_buoy)

               !-v1
               entr_threshold = 0.6  ! range 0-1, smaller => smaller entrainment
                                     ! =0.6 => buoy_exc/mx_buoy > 0.6 => smaller
                                     !         entrainment
               factor =  min(maxval(buoy_exc(i,kts:k22(i))),mx_buoy)/mx_buoy
               factor=   0.5*(1.-atan( (factor-entr_threshold)*5.)/1.37)
               factor=   max(0.1,min(factor,0.9))


               !-v2
               !CWV=max(0.1, min(col_sat(i),0.9))
               !factor = 1.-exp(15.6*(col_sat(i)-0.603))
               !factor = max(0.1,min(factor,0.9))

               !-v3
               !factor = min(max(4.875e-21*CWV**11.092,0.1),0.9) !CW=TPW [mm]
               !factor = 1.-factor
               !


               !-v4
               !factor = sqrt (factor_v1 * factor_v2)

               !-v5 avoid extra-buoyancy where rained before
               !- rain threshold is 1.e-7 mm/s => above this amount, x_add_buoy(i) -> 0.
               !- AA1_CIN is a dummy array that stores the conv precipitation of the previous timestep.
               x_add_buoy(i) =  x_add_buoy(i) * max(0.0,min(1.-atan(( (AA1_CIN_(i)-1.e-7)*1.e+7)/1.569),1.0))



               !- modulates the gross entraiment depending on the cold pool presence.
               entr_rate(i,:) = entr_rate(i,:)*factor



               !----------------------
               AA0_(i)= max(0.0,min(1.-atan(( (AA1_CIN_(i)-1.e-7)*1.e+7)/1.569),1.0))       !--- temporary for output
               AA2_(i)=factor        !--- temporary for output
               AA3_(i)=x_add_buoy(i) !--- temporary for output
               !----------------------
             enddo
          endif
          if(USE_MEMORY == 20) then
             do i=its,itf
               cap_max_inc = 0.

               cap_max(i)           = cap_maxs
               cap_max_increment(i) = cap_max_inc

               if(ierr(i) /= 0) cycle

               x_add_buoy(i) = min(maxval(buoy_exc(i,kts:k22(i))),mx_buoy)

               entr_threshold= 0.6 !v5
               factor=   0.5*(1.-atan( (factor-entr_threshold)*5)/1.37)
               factor=   max(0.1,min(factor,0.9))


               !!!factor          = min(1.,5*x_add_buoy(i)/mx_buoy)*200.
               cap_max(i) = cap_maxs + (1.-factor) *200.

               !----------------------
               AA2_(i)=(1.-factor) !--- temporaryfor output
               AA3_(i)=x_add_buoy(i) !--- temporary for output
               !----------------------
             enddo
          endif
        endif

     ENDIF !- convection_tracer
!
!
!------- determine LCL for the air parcels around K22
!
       DO i=its,itf
        klcl(i) = k22(i) ! default value
        if(ierr(i) == 0)then
           !tlll, rlll,plll - temp, water vapor and pressure of the source air parcel
           x_add = zqexec(i)
           call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),q_cup (i,kts:kte),rlll,k22(i),x_add)
           x_add = ztexec(i)
           call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),t_cup (i,kts:kte),tlll,k22(i),x_add)
           call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),p_cup (i,kts:kte),plll,k22(i))

           call get_lcl(tlll,100.*plll,rlll,tlcl,plcl,dzlcl)
           !print*,"MID",tlll,100.*plll,rlll,tlcl,plcl,dzlcl; call flush(6)

           !-get LCL
           if(dzlcl >= 0.) then ! LCL found (if dzlcl<0 => not found)
             call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),z_cup (i,kts:kte),zlll,k22(i))
loop0:       do k=kts,ktf
                if(z_cup(i,k).gt.zlll+dzlcl)then
                  klcl(i)=max(k,k22(i))
                  exit loop0
                endif
              enddo loop0
              klcl(i)=min(klcl(i),ktf-4)
           endif
        endif
        !write(12,111)'MDlcl',tlcl,plcl,dzlcl,klcl(i),ierr(i)
        111 format(1x,A5,3F10.2,2i4)
      enddo
      !-- check if LCL is below PBL height for shallow convection
      if(USE_LCL .and. trim(cumulus) == 'mid') then
       do i=its,itf
        if(ierr(i).eq.0)then
          if(klcl(i) > max(1,kpbl(i)+1)) then
              ierr(i)=21
              ierrc(i)='for mid convection:  LCL height > PBL height'
          endif
        endif
       enddo
      endif

      IF(ADV_TRIGGER==1 .and. trim(cumulus) /= 'shallow') THEN
         wkf = 0.02 ! m/s
         do i=its,itf
              if(ierr(i) /= 0) CYCLE

              k   = klcl(i); dzlcl = z_cup(i,k)-z1(i)

              ckf = wkf
              if(dzlcl .le. 2.e+3)ckf = wkf* dzlcl/2000.

              wkflcl =-(omeg(i,max(kts,k-1),1)/rho(i,max(k-1,kts)) + &
                        omeg(i,k           ,1)/rho(i,k                 ) + &
                        omeg(i,k+1         ,1)/rho(i,k+1         ) )/(3.*g)

              !...check to see if cloud is buoyant using fritsch-chappell trigger
              !...function described in kain and fritsch (1992)...w0avg is an
              !...aproximate value for the running-mean grid-scale vertical
              !...velocity, which gives smoother fields of convective initiation
              !...than the instantaneous value...formula relating temperature
              !...perturbation to vertical velocity has been used with the most
              !...success at grid lengths near 25 km.  for different grid-lengths,
              !...adjust vertical velocity to equivalent value for 25 km grid
              !...length, assuming linear dependence of w on grid length...
              if(dx(i) >= 25.e+3) then
                       wkflcl = wkflcl*dx(i)/25.E3 - ckf
              else
                       wkflcl = wkflcl - ckf
              endif
             !--- think about letting wkflcl <0 => Tpert<0 =>prevent convection in subsidence areas
             ! wkflcl = max(wkflcl,0.) ! -- only positive.
             !
             !-- Kain (2004) Eq. 1
              Tpert(i,:) = 4.64*wkflcl**(1./3.)

             !AA1_cin_(i) = Tpert(i,kts)

              if(maxval(Tpert(i,:))>+2.) Tpert(i,:)= 2.
              if(minval(Tpert(i,:))<-2.) Tpert(i,:)=-2.

             !DTLCL=4.64*WKL**0.33,  WSIGNE = signal of GDT
             !GDT=G*DTLCL*(ZLCL-Z0(LC))/(TV0(LC)+TVEN)
             !WLCL=1.+.5*WSIGNE*SQRT(ABS(GDT)+1.E-10)    !<< velocity at LCL
             !if (MAPL_AM_I_ROOT()) then
             ! if(maxval(Tpert)>+3.) print*,"MAX TPERT=",maxval(Tpert),minval(Tpert)
             ! if(minval(Tpert)<-3.) print*,"MIN TPERT=",maxval(Tpert),minval(Tpert)
             !endif
         enddo
      ENDIF
!
!--- determine the moist static energy of air parcels at source level
!
      do i=its,itf
       if(ierr(i) /= 0)cycle
       x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
       call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add,Tpert(i,kts:kte))
       call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add,Tpert(i,kts:kte))
      enddo
!
!--- determine the vertical entrainment/detrainment rates, the level of convective cloud base -kbcon-
!--- and the scale dependence factor (sig).
      DO i=its,itf
         if(ierr(i) /= 0) CYCLE
         do k=kts,kte
            frh = min(qo_cup(i,k)/max(qeso_cup(i,k),smallerQV),1.)
            if(k >= klcl(i)) then
               entr_rate(i,k)=entr_rate(i,k)*(1.3-frh)*(qeso_cup(i,k)/qeso_cup(i,klcl(i)))**3
            else
               entr_rate(i,k)=entr_rate(i,k)*(1.3-frh)
            endif
            if (ZERO_DIFF==1) then
               cd(i,k)=0.75e-4*(1.6-frh)
            else
               entr_rate(i,k) = max(entr_rate(i,k),min_entr_rate)
               if(trim(cumulus) == 'deep'   ) cd(i,k)=0.10*entr_rate(i,k)
               if(trim(cumulus) == 'mid'    ) cd(i,k)=0.50*entr_rate(i,k)
               if(trim(cumulus) == 'shallow') cd(i,k)=0.75*entr_rate(i,k)
            endif
         enddo
      ENDDO
!
!--- start_level
!
      start_level(:)=  KLCL(:)
!
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
        CALL cup_cloud_limits(cumulus,ierrc,ierr,cap_max_increment,cap_max,heo_cup,heso_cup,qo_cup   &
                             ,qeso_cup,po,po_cup,zo_cup,heo,hkbo,qo,qeso,entr_rate,hcot,k22,kbmax &
                             ,klcl,kbcon,ktop,depth_neg_buoy,frh_bcon,Tpert,start_level              &
                             ,use_excess,zqexec,ztexec,x_add_buoy,xland,itf,ktf,its,ite, kts,kte)

!--- define entrainment/detrainment profiles for downdrafts
     do i=its,itf
        do k=kts,ktf
           mentrd_rate(i,k) = entr_c(k,i)*entr_rate_plume*0.3
        enddo
     enddo
     cdd = mentrd_rate
     !- scale dependence factor
                         sigd(:) = 1.0
     if( DOWNDRAFT /= 0) sigd(:) = 0.0

!--- update hkb/hkbo in case of k22 is redefined in 'cup_kbon'
     do i=its,itf
        IF(ierr(i) /= 0) CYCLE
        x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
        call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),he_cup (i,kts:kte),hkb (i),k22(i),x_add,Tpert(i,kts:kte))
        call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup(i,kts:kte),hkbo(i),k22(i),x_add,Tpert(i,kts:kte))
     enddo
!
!--- increase detrainment in stable layers
!
      CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,itf,ktf,its,ite, kts,kte)
!
!--- use KTOP for plumes
!
      CALL rates_up_pdf(cumulus,ktop,ierr,po_cup,entr_rate,hkbo,heo,heso_cup,zo_cup, &
                        kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,klcl,hcot)
!
      IF(trim(cumulus) == 'deep') THEN
        !
        !-- check if ktop is too low for deep convection
        !
        do i=its,itf 
           if(ierr(i) /= 0)cycle
           if(po_cup(i,ktop(i)) > 400) then
               ierr(i)=22
               ierrc(i)='deep convection with cloud top too low'
           endif
        enddo
      ENDIF
!
      IF(trim(cumulus) == 'mid') THEN
        !
        !-- check if ktop is too high for mid convection
        !
        do i=its,itf
           if(ierr(i) /= 0)cycle
           if(po_cup(i,ktop(i)) < 400) then
               ierr(i)=22
               ierrc(i)='mid convection with cloud top too high'
           endif
        enddo
      ENDIF
!
      IF(trim(cumulus) == 'shallow') THEN
        !
        !-- check if ktop is too high for shallow convection
        !
        do i=its,itf
           if(ierr(i) /= 0)cycle
           if(po_cup(i,ktop(i)) < 700) then
               ierr(i)=23
               ierrc(i)='shallow convection with cloud top too high'
           endif
        enddo
      ENDIF
!
!-- avoid double-counting plumes
!
      DO i=its,itf   
        if(ierr(i) /= 0)cycle           
        if(last_ierr(i) == 0) then      
           ierr(i)=27
           ierrc(i)='prevented double-counting plumes'
        endif 
      ENDDO
!
!-- last checks for ktop
!
      DO i=its,itf
        if(ierr(i) /= 0)cycle
        if( (z_cup(i,ktop(i))-z_cup(i,kbcon(i))) < depth_min ) then
            ierr(i)=5
            ierrc(i)='cloud depth too small'
        endif
      ENDDO
      DO i=its,itf
        if(ktop(i) <= kbcon(i))then
            ierr(i)=5
            ierrc(i)='ktop too small'
        endif
      ENDDO
!
!--- determine the normalized mass flux profile for updraft
!
      do i=its,itf
         zuo(i,:)=0.
         if(ierr(i) /= 0) cycle
         CALL get_zu_zd_pdf(trim(cumulus),trim(cumulus)//"_up",ierr(i),k22(i),ktop(i),zuo(i,kts:kte),kts,kte,ktf  &
                           ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i),random(i))
      enddo

      do i=its,itf
         if(ierr(i) /= 0) cycle
         xzu(i,:)= zuo(i,:)
         zu (i,:)= zuo(i,:)
      enddo
!
! calculate mass entrainment and detrainment
!
      CALL get_lateral_massflux(itf,ktf, its,ite, kts,kte                            &
                               ,ierr,ktop,zo_cup,zuo,cd,entr_rate,po_cup          &
                               ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
                               ,cumulus,kbcon,k22,kpbl,up_massentru,up_massdetru,lambau_dp)
      uc  =0.
      vc  =0.
      hc  =0.
      hco =0.
      do i=its,itf
       IF(ierr(i).eq.0)THEN
         do k=kts,start_level(i)
            hc (i,k) =hkb (i)
            hco(i,k) =hkbo(i)
            !-get uc and vc as average between layers below k22
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),u_cup(i,kts:kte),uc(i,k),k22(i))
            call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),v_cup(i,kts:kte),vc(i,k),k22(i))
         enddo
       ENDIF
      enddo
!
!--- 1st guess for moist static energy and dbyo (not including ice phase)
!
     do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=start_level(i)  +1,ktop(i) + 1  ! mass cons option
          denom=(zu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
          if(denom > 0.0)then
            hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                                               up_massentro(i,k-1)*heo(i,k-1))  / denom
            if(k==start_level(i)+1) then
               x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
               hco(i,k)= hco(i,k) + x_add*up_massentro(i,k-1)/denom
            endif
           else
            hco(i,k)= hco(i,k-1)
          endif
         enddo
         do k=ktop(i)+2,ktf
           hco (i,k)=heso_cup(i,k)!=heo_cup(i,k)
         enddo
     enddo
!
!--- Get buoyancy of updrafts
!
     call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hco,heo_cup,heso_cup,dbyo,zo_cup)

!--- get "c1d" profile ----------------------------------------
     if(trim(cumulus) == 'deep' .and. USE_C1D) then
        Do i=its,itf
            if(ierr(i) /= 0) cycle
            c1d(i,kbcon(i)+1:ktop(i)-1)=abs(c1)
        ENDDO
     endif

     if(FIRST_GUESS_W .or. AUTOCONV == 4) then
        call cup_up_moisture_light(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
                                  ,cnvfrc,srftype,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,GAMMAo_cup,zuo          &
                                  ,qeso_cup,k22,qo_cup,ZQEXEC,use_excess,rho,up_massentr,up_massdetr    &
                                  ,psum,psumh,c1d,x_add_buoy,1,itf,ktf,ipr,jpr,its,ite, kts,kte         )

       call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
                       ,tempco,qco,qrco,qo,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte       )

     endif

!
!--- calculate moisture properties of updraft
     call cup_up_moisture(cumulus,start_level,klcl,ierr,ierrc,zo_cup,qco,qrco,pwo,pwavo,hco,tempco,xland   &
                         ,cnvfrc,srftype,po,p_cup,kbcon,ktop,cd,dbyo,clw_all,t_cup,qo,GAMMAo_cup,zuo,qeso_cup &
                         ,k22,qo_cup,ZQEXEC,use_excess,ccn,rho,up_massentr,up_massdetr,psum    &
                         ,psumh,c1d,x_add_buoy,vvel2d,vvel1d,zws,entr_rate                  &
                         ,1,itf,ktf,ipr,jpr,its,ite, kts,kte                                   )

     do i=its,itf
        if(ierr(i) /= 0) cycle
        cupclw(i,kts:ktop(i)+1)=qrco(i,kts:ktop(i)+1)
     enddo


!
!--- get melting profile
!
     call get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
                             ,pwo,edto,pwdo,melting                               &
                             ,itf,ktf,its,ite, kts,kte, cumulus                   )

!
!--- updraft moist static energy + momentum budget
!
!--- option to produce linear fluxes in the sub-cloud layer.
     if(trim(cumulus) == 'shallow' .and. use_linear_subcl_mf == 1) then
       do i=its,itf
          if(ierr(i) /= 0) cycle
           call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
                          ,he_cup (i,kts:kte), hc (i,kts:kte))
           call get_delmix(cumulus,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
                          ,heo_cup(i,kts:kte), hco(i,kts:kte))
       enddo
      endif
!
      do i=its,itf
         if(ierr(i) /= 0) cycle

         do k=start_level(i)+1 , ktop(i)+1  ! mass cons option
          denom =(zu(i,k-1)-.5*up_massdetr (i,k-1)+up_massentr (i,k-1))
          denomU=(zu(i,k-1)-.5*up_massdetru(i,k-1)+up_massentru(i,k-1))
          if(denom > 0.0 .and. denomU >0.0)then

            hc (i,k)=(hc (i,k-1)*zu (i,k-1)-.5*up_massdetr (i,k-1)*hc (i,k-1) + &
                                               up_massentr (i,k-1)*he (i,k-1))/ denom

            hco(i,k)=(hco(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco(i,k-1)+ &
                                               up_massentro(i,k-1)*heo(i,k-1))/ denom
            if(k==start_level(i)+1) then
             x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
             hco(i,k)= hco(i,k) + x_add*up_massentro(i,k-1)/denom
             hc (i,k)= hc (i,k) + x_add*up_massentr (i,k-1)/denom
            endif
                         !assuming zuo=zu,up_massdetro=up_massdetr, ...
                         !(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))

            uc(i,k)=(uc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*uc(i,k-1)+     &
                                            up_massentru(i,k-1)*us(i,k-1)      &
                   -pgcon*.5*(zu(i,k)+zu(i,k-1))*(u_cup(i,k)-u_cup(i,k-1))) /  denomU

            vc(i,k)=(vc(i,k-1)*zu(i,k-1)-.5*up_massdetru(i,k-1)*vc(i,k-1)+     &
                                            up_massentru(i,k-1)*vs(i,k-1)      &
                   -pgcon*.5*(zu(i,k)+zu(i,k-1))*(v_cup(i,k)-v_cup(i,k-1))) /  denomU


          else
            hc (i,k)= hc (i,k-1)
            hco(i,k)= hco(i,k-1)
            uc (i,k)= uc (i,k-1)
            vc (i,k)= vc (i,k-1)
          endif
          !---meltglac-------------------------------------------------
          !- includes glaciation effects on HC,HCO
          !                    ------ ice content --------
          !print*,"H=",hc (i,k),(1.-p_liq_ice(i,k))*qrco(i,k)*xlf,hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
           hc (i,k)= hc (i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
           hco(i,k)= hco(i,k)+(1.-p_liq_ice(i,k))*qrco(i,k)*xlf

         enddo
         !
         do k=ktop(i)+2,ktf
           hc  (i,k)= hes_cup(i,k)!= he_cup(i,k)
           uc  (i,k)=   u_cup(i,k)
           vc  (i,k)=   v_cup(i,k)
           hco (i,k)=heso_cup(i,k)!=heo_cup(i,k)
           zu  (i,k)=0.
           zuo (i,k)=0.
          enddo
     enddo
!
!--- Get buoyancy of updrafts
!
     call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hc, he_cup, hes_cup,  dby, z_cup)
     call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,hco,heo_cup,heso_cup,dbyo,zo_cup)

!
     if(.not. FIRST_GUESS_W) then
     !--- calculate in-cloud/updraft air temperature for vertical velocity
     !
        do i=its,itf
          if(ierr(i) == 0)then
            do k=kts,ktf
               tempco (i,k) = (1./cp)*(hco (i,k)-g*zo_cup(i,k)-xlv*qco (i,k))
            enddo
            tempco (i,kte)=tn_cup(i,kte)
          else
            tempco (i,:)  =tn_cup(i,:)
          endif
        enddo
    !
    !--- vertical velocity
    !
        call cup_up_vvel(vvel2d,vvel1d,zws,entr_rate,cd,zo,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup &
                        ,tempco,qco,qrco,qo,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
     endif

!---- new rain
!
!--- calculate rain mixing ratio in updrafts
!
!       call cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland,cnvfrc,srftype&
!                       ,zo_cup,qco,qrco,pwo,pwavo,po,p_cup,t_cup,tempco&
!                       ,zuo,up_massentr,up_massdetr,vvel2d,rho         &
!                       ,qrr                                            &
!                       ,itf,ktf,its,ite, kts,kte)

!---- new rain
!
!
!--- DOWNDRAFT section
!
      DO i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,ktf
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
            enddo
         endif
         37   CONTINUE
      ENDDO
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(heso_cup,k22,kzdown,jmin,ierr,itf,ktf,its,ite, kts,kte)

      call get_jmin(cumulus,itf,ktf,its,ite, kts,kte,ierr,kdet,ktop,kbcon,jmin,ierrc  &
                   ,beta,depth_min,heso_cup,zo_cup,melting_layer)
!
!--- this calls routine to get downdrafts normalized mass flux
!
       do i=its,itf
         zd(i,:)=0.
         if(ierr(i) /= 0) cycle
         call get_zu_zd_pdf(trim(cumulus),"DOWN",ierr(i),kdet(i),jmin(i),zdo(i,:),kts,kte,ktf&
                           ,kpbl(i),k22(i),kbcon(i),klcl(i),po_cup(i,kts:kte),psur(i),xland(i),random(i))
       enddo
!
!---  calls routine to get lateral mass fluxes associated with downdrafts
!
       call get_lateral_massflux_down(trim(cumulus),itf,ktf, its,ite, kts,kte &
                                     ,ierr,jmin,zo_cup,zdo,xzd,zd,cdd,mentrd_rate      &
                                     ,dd_massentro,dd_massdetro ,dd_massentr, dd_massdetr &
                                     ,cumulus,dd_massentru,dd_massdetru,lambau_dn)

!---  calls routine to get wet bulb temperature and moisture at jmin
!
       IF(USE_WETBULB == 1 .and. trim(cumulus) /= 'shallow' ) THEN
        do i=its,itf
          if(ierr(i) /= 0) cycle
          k = jmin(i)
          call get_wetbulb(jmin(i),qo_cup(i,k),t_cup(i,k),po_cup(i,k),q_wetbulb(i),t_wetbulb(i))
         !print*,"wb       =",jmin,qo_cup(i,k),t_cup(i,k),q_wetbulb(i),t_wetbulb(i)
         !print*,"evap/cool=",q_wetbulb(i)-qo_cup(i,k),t_wetbulb(i)-t_cup(i,k)
        enddo
       ENDIF
!
!
!--- downdraft moist static energy + moisture budget
!
       do i=its,itf
         hcdo (i,:)= heso_cup(i,:)
         ucd  (i,:)=    u_cup(i,:)
         vcd  (i,:)=    v_cup(i,:)
         dbydo(i,:)= 0.
       enddo

       do i=its,itf
            bud(i)=0.
            if(ierr(i)/= 0 .or. trim(cumulus) == 'shallow')cycle
            i_wb=0
            !--for future test
            if(use_wetbulb==1) then
             !--option 1
             !hcdo(i,jmin(i))=cp*t_wetbulb(i)+xlv*q_wetbulb(i)+zo_cup(i,jmin(i))*g
             !--option 2
             hcdo(i,jmin(i))=0.5*(cp*t_wetbulb(i)+xlv*q_wetbulb(i)+zo_cup(i,jmin(i))*g + hc(i,jmin(i)))
             i_wb=1
            endif

            dbydo(i,jmin(i))=hcdo(i,jmin(i))-heso_cup(i,jmin(i))
            bud(i)=dbydo(i,jmin(i))*(zo_cup(i,jmin(i)+1)-zo_cup(i,jmin(i)))

            do ki=jmin(i) - i_wb ,kts,-1!do ki=jmin(i)-1,1,-1
                denom = zdo(i,ki+1)-0.5*dd_massdetro(i,ki)+dd_massentro(i,ki)
                denomU= zdo(i,ki+1)-0.5*dd_massdetru(i,ki)+dd_massentru(i,ki)
                 !-tmp fix for denominator being zero
                if(denom > 0.0 .and. denomU >0.0)then
                    dzo=zo_cup(i,ki+1)-zo_cup(i,ki)

                    ucd(i,ki)=(ucd(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetru(i,ki)*ucd(i,ki+1)+ &
                                                          dd_massentru(i,ki)*us (i,ki)    &
                             -pgcon*zdo(i,ki+1)*(us(i,ki+1)-us(i,ki)))   /  denomU
                    vcd(i,ki)=(vcd(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetru(i,ki)*vcd(i,ki+1)+        &
                                                          dd_massentru(i,ki)*vs (i,ki)           &
                             -pgcon*zdo(i,ki+1)*(vs(i,ki+1)-vs(i,ki)))   /  denomU

                    hcdo(i,ki)=(hcdo(i,ki+1)*zdo(i,ki+1)-.5*dd_massdetro(i,ki)*hcdo(i,ki+1)+     &
                                                            dd_massentro(i,ki)*heo(i,ki))  /denom

                    dbydo(i,ki)=hcdo(i,ki)-heso_cup(i,ki)
                    !if(i.eq.ipr)write(0,*)'ki,bud = ',ki,bud(i),hcdo(i,ki)
                    bud(i)=bud(i)+dbydo(i,ki)*dzo
                else
                    ucd (i,ki)= ucd(i,ki+1)
                    vcd (i,ki)= vcd(i,ki+1)
                    hcdo(i,ki)=hcdo(i,ki+1)
                endif
            enddo
            if(bud(i).gt.0)then
              ierr(i)=7
              ierrc(i)='downdraft is not negatively buoyant '
            endif
       enddo
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture(cumulus,ierrc,zdo,hcdo,heso_cup,qcdo,qeso_cup,     &
           pwdo,qo_cup,zo_cup,dd_massentro,dd_massdetro,jmin,ierr,gammao_cup, &
!--test    pwevo,bu,qrcdo,qo,heo,t_cup,1,t_wetbulb,q_wetbulb,qco,pwavo,       &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,t_wetbulb,q_wetbulb,qco,pwavo,      &
           itf,ktf,its,ite, kts,kte)
!
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z_cup ,zu ,dby  ,GAMMA_CUP    ,t_cup   ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
      call cup_up_aa0(aa1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

      do i=its,itf
         if(ierr(i) /= 0) cycle
         if(aa1(i).eq.0.)then
               ierr(i)=17
               ierrc(i)="cloud work function zero"
         endif
      enddo
!
!--- calculate CIN for updrafts
!
      call cup_up_aa0(cin0,z_cup ,zu ,dby  ,GAMMA_CUP    ,t_cup   ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
      call cup_up_aa0(cin1,zo_cup,zuo,dbyo ,GAMMAo_CUP   ,tn_cup  ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,'CIN')
!
!----trigger function for KE+CIN < 0 => no convection
!
      IF( DICYCLE>1) THEN
        do i=its,itf
          if(ierr(i) /= 0) cycle
          print*,"cin=",cin0(i),0.5*zws(i)**2, omeg(i,kpbl(i),1)/(-g*rho(i,kpbl(i))) ;call flush(6)
          if(cin0(i) + 0.5*zws(i)**2 < 0.)then !think about including the grid scale vertical velocity at KE calculation
               ierr(i)=19
               ierrc(i)="CIN negat"
          endif
         enddo
       ENDIF
!
!
!--- calculate in-cloud/updraft and downdraft air temperature for vertical velocity
!
      do i=its,itf
          if(ierr(i) == 0)then
            do k=kts,ktf
               tempcdo(i,k) = (1./cp)*(hcdo(i,k)-g*zo_cup(i,k)-xlv*qcdo(i,k))
            enddo
          else
               tempcdo(i,:)=tn_cup(i,:)
          endif
      enddo
!
!
!--- diurnal cycle section
!
      if(trim(cumulus)=='deep') then
           T_star=5.   ! T_star = temp scale in original paper = 1 K
      else
           T_star=40.
      endif
!
!--- Bechtold et al 2008 time-scale of cape removal
!
      IF(SGS_W_TIMESCALE == 0) THEN
         DO i=its,itf
            if(ierr(i) /= 0) cycle
            !- time-scale cape removal
            if(trim(cumulus)=='deep') tau_ecmwf(i)=tau_deep * (1.0 + (1.0-sig(i)))
            if(trim(cumulus)=='mid' ) tau_ecmwf(i)=tau_mid  * (1.0 + (1.0-sig(i)))
         ENDDO
      ELSE
         DO i=its,itf
            if(ierr(i) /= 0) cycle
            !- time-scale cape removal from Bechtold et al. 2008
            dz = max(z_cup(i,ktop(i)+1)-z_cup(i,kbcon(i)),1.e-16) ! cloud depth (H)
            tau_ecmwf(i)=(dz / vvel1d(i)) * (1.0 + sig(i)) ! resolution dependent scale factor
            tau_ecmwf(i)= max(dtime,tau_ecmwf(i)*real(SGS_W_TIMESCALE))
         ENDDO
      ENDIF
      DO i=its,itf
         if(ierr(i) /= 0) cycle
         if(xland(i) > 0.99 ) then !- over water
            umean= 2.0+sqrt(0.5*(US(i,1)**2+VS(i,1)**2+US(i,kbcon(i))**2+VS(i,kbcon(i))**2))
            tau_bl(i) = (zo_cup(i,kbcon(i))- z1(i)) /umean
         else !- over land
            tau_bl(i) =( zo_cup(i,ktop(i))- zo_cup(i,kbcon(i)) ) / 3.0 ! 3.0 m/s is estimated wmean
         endif
      ENDDO
      tau_ec_ = tau_ecmwf
      tau_bl_ = tau_bl

      IF( (DICYCLE==1 .or. DICYCLE == 6) .or. (DICYCLE==0 .and. trim(cumulus)=='mid') ) THEN
        iversion=0
        !-- calculate pcape from BL forcing only
        call cup_up_aa1bl(iversion,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup, &
                          rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,                          &
                          xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus,tn_bl,qo_bl  )
        if(DICYCLE==6) then
           call cup_up_aa1bl(iversion,aa0_bl,aa1_fa,aa1,t,tn_bl,q,qo_bl,dtime,po_cup,zo_cup,zuo,dbyo,GAMMAo_CUP,tn_cup, &
                             rho,klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,&
                             xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus,tn_bl,qo_bl   )
           where(ierr == 0) aa1_bl=0.50*aa1_bl+0.50*aa0_bl !mix1
        endif
        DO i=its,itf
            if(ierr(i) /= 0) cycle
            aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i)
!           aa1_bl(i) = (aa1_bl(i)/T_star) * tau_bl(i) - cin1(i)
        ENDDO
      !==========================================================
      ELSEIF( DICYCLE==4) then!zhang(2002)

         !- T and Q profiles modified only by RAD+ADV tendencies
         DO i=its,itf
                if ( ierr(i) /= 0 ) CYCLE
                tn_x(i,kts:ktf) = tn(i,kts:ktf)-tn_bl(i,kts:ktf)+t(i,kts:ktf)
                qo_x(i,kts:ktf) = qo(i,kts:ktf)-qo_bl(i,kts:ktf)+q(i,kts:ktf)
         ENDDO

         !--- calculate moist static energy, heights, qes, ... only by free troposphere tendencies
         call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1, &
                    psur,ierr,-1,itf,ktf, its,ite, kts,kte)
         !--- environmental values on cloud levels only by FT tendencies
         call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,           &
                           us,vs,u_cup,v_cup,                                    &
                           heso_cup_x,zo_cup,po_cup,gammao_cup_x,tn_cup_x,psur,tsur,  &
                           ierr,z1,itf,ktf,its,ite, kts,kte)
         !--- this is (DT_ve/Dt)_adv+rad
         DO i=its,itf
             IF(ierr(i) /= 0) cycle
             aa3(i)=0.
             DO k=max(kbcon(i),kts+1),ktop(i)
                 dp=- (log(100.*po(i,k))-log(100.*po(i,k-1))) !no units
                 aa3(i)=aa3(i)- (tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)) - &
                                 t_cup   (i,k)*(1.+0.608*q_cup   (i,k)) ) *dp/dtime !units = K/s
                !print*,"tve=",k,aa3(i),tn_cup_x(i,k)*(1.+0.608*qo_cup_x(i,k)),&
                !               t_cup (i,k)*(1.+0.608*q_cup   (i,k)),dp
             ENDDO
         ENDDO
         DO i=its,itf
            if(ierr(i)/= 0)CYCLE
            !- this is (DCAPE_env/Dt)_adv+rad
            !aa1_bl(i) = -aa3(i)
            !- Zhang threshold:  65 J/kg/hour => 65/(Rd *3600)= 63 10^-6 K/s
            aa1_bl(i) = aa3(i)-(63.e-6)!*1.5
            !print*,"dcape_env=",aa3(i),aa1_bl(i)
            if(xland(i) > 0.90 ) aa1_bl(i)=1.4*aa1_bl(i) !- over water
         ENDDO
         !--- this is (DT_ve/Dt)_cu
         DO i=its,itf
             dtdt(i,:)=0.
             dqdt(i,:)=0.
             IF(ierr(i) /= 0) cycle
             DO k=max(kbcon(i),kts+1),ktop(i)
                 dp=100.*(po_cup(i,k+1)-po_cup(i,k))
                 RZenv=0.5*(zuo(i,k+1)+zuo(i,k) - (zdo(i,k+1)+zdo(i,k))*edto(i))
                 S2= cp*tn_cup_x(i,k+1) + g*zo_cup(i,k+1)
                 S1= cp*tn_cup_x(i,k  ) + g*zo_cup(i,k  )
                 Q2= qo_cup_x(i,k+1)
                 Q1= qo_cup_x(i,k  )

                 dqdt(i,k)=         -RZenv*(Q2-Q1)*g/dp
                 dtdt(i,k)= -(1./cp)*RZenv*(S2-S1)*g/dp

                 dqdt(i,k)= dqdt(i,k)+(up_massdetro(i,k)*0.5*(qco (i,k+1)+qco (i,k)-(Q2+Q1)) &
                             + edto(i)*dd_massdetro(i,k)*0.5*(qcdo(i,k+1)+qcdo(i,k)-(Q2+Q1)))*g/dp

                 dtdt(i,k)= dtdt(i,k)+(up_massdetro(i,k)*0.5*(tempco  (i,k+1) + tempco  (i,k)-  &
                                                             (tn_cup_x(i,k+1) + tn_cup_x(i,k))) &
                             + edto(i)*dd_massdetro(i,k)*0.5*(tempcdo (i,k+1) + tempcdo (i,k)-  &
                                                             (tn_cup_x(i,k+1) + tn_cup_x(i,k))) &
                             )*g/dp
                !print*,"dtdt=",k, dtdt(i,k),zuo(i,k+1),zdo(i,k+1),dqdt(i,k)
             ENDDO
             xk_x(i)=0.
             DO k=max(kbcon(i),kts+1),ktop(i)
                dp=-(log(100.*po_cup(i,k+1))-log(100.*po_cup(i,k)))      ! no units here
                xk_x(i)=xk_x(i)+   ( (1.+0.608*qo_x(i,k))*dtdt(i,k) + &
                                         0.608*tn_x(i,k) *dqdt(i,k) )*dp !  units=K m/Pa s2
                !=> aa3/xk_x will have units of kg/m2/s for the mass flux at cloud base.
                !print*,"xk_x=",k, xk_x(i),dtdt(i,k),dqdt(i,k)
             ENDDO
         ENDDO
      ENDIF
      !==========================================================
      IF(DICYCLE==5 .or. DICYCLE==2) then

PHY2:  DO step=fa,fa

        IF(step==bl) then
           !- version for real cloud-work function by getting the profiles modified only by bl tendencies
           DO i=its,itf
                if ( ierr(i) /= 0 ) CYCLE
                k_free_trop=kpbl(i)
                !below kbcon -> modify profiles
                tn_x(i,1:k_free_trop) = tn(i,1:k_free_trop)
                qo_x(i,1:k_free_trop) = qo(i,1:k_free_trop)
                 !above kbcon -> keep profiles
                tn_x(i,k_free_trop+1:kte) = t(i,k_free_trop+1:kte)
                qo_x(i,k_free_trop+1:kte) = q(i,k_free_trop+1:kte)
           ENDDO
        ELSEIF(step==FA) then
           !- version for real cloud-work function by getting the profiles modified only by free atmosphere tendencies
           DO i=its,itf
                if ( ierr(i) /= 0 ) CYCLE
                k_free_trop=kpbl(i)
                !THESE ARE NOT MODIFIED
                !below kbcon -> keep profiles
                tn_x(i,1:k_free_trop) = t(i,1:k_free_trop)
                qo_x(i,1:k_free_trop) = q(i,1:k_free_trop)
                !THESE MUST BE ONLY FROM ADV + RAD
                !above kbcon -> modify profiles
                tn_x(i,k_free_trop+1:kte) = tn(i,k_free_trop+1:kte)-tn_bl(i,k_free_trop+1:kte)+t(i,k_free_trop+1:kte)
                qo_x(i,k_free_trop+1:kte) = qo(i,k_free_trop+1:kte)-qo_bl(i,k_free_trop+1:kte)+q(i,k_free_trop+1:kte)
           ENDDO
        ENDIF

        !--- calculate moist static energy, heights, qes, ... only by bl tendencies
        call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1, &
                    psur,ierr,-1,itf,ktf, its,ite, kts,kte)
        !--- environmental values on cloud levels only by bl tendencies
        call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,             &
                          us,vs,u_cup,v_cup,                                                   &
                          heso_cup_x,zo_cup,po_cup,gammao_cup_x,tn_cup_x,psur,tsur,  &
                          ierr,z1,itf,ktf,its,ite, kts,kte)
        do i=its,itf
             if(ierr(I) /= 0)cycle
             if(step==bl) then
                x_add = (xlv*zqexec(i)+cp*ztexec(i))
             else
                x_add = 0.
             endif
             call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup_x(i,kts:kte),hkbo_x(i),k22(i),x_add,Tpert(i,kts:kte))
        enddo
        do i=its,itf
             hco_x (i,:)=0.
             IF(ierr(i) /= 0)cycle
             do k=kts,start_level(i)
               hco_x (i,k) = hkbo_x(i)
             enddo
        enddo!
        DO i=its,itf
            if(ierr(i) /= 0)cycle
            do k=start_level(i)  +1,ktop(i)+1  ! mass cons option
                denom=(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                if(denom > 0.0)then
                   hco_x (i,k)=(hco_x(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_x(i,k-1) + &
                                up_massentro(i,k-1)*heo_x(i,k-1))/ denom
                else
                   hco_x (i,k)=hco_x(i,k-1)
                endif
                hco_x (i,k)=hco_x (i,k) +(1.-p_liq_ice(i,k))*qrco(i,k)*xlf
             enddo
             do k=ktop(i)+2,ktf
                  hco_x (i,k)=heso_cup_x(i,k)
             enddo
        ENDDO
        call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
                        ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup)

        !--- calculate workfunctions for updrafts
        IF(step==bl) &
           call cup_up_aa0(aa1_bl,zo_cup,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x &
                          ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

        IF(step==fa) &
           call cup_up_aa0(aa1_fa,zo_cup,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x &
                          ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)
       ENDDO  PHY2
       DO i=its,itf
            if(ierr(i)/= 0)CYCLE
            aa1_bl(i) = aa1_fa(i)
       ENDDO
ENDIF
!
!--- Trigger function based on Xie et al 2019
!
     IF(ADV_TRIGGER == 3 .and. trim(cumulus) /= 'shallow') THEN
!    IF(ADV_TRIGGER == 3 ) THEN

       daa_adv_dt=0.

       do step=1,2
              !--- calculate moist static energy, heights, qes, ... only by ADV tendencies
              if(step==1) then
                tn_x = t     ;  qo_x = q
              else
                tn_x = tn_adv;  qo_x = qo_adv
              endif

              call cup_env(zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1,psur,ierr,-1,itf,ktf, its,ite, kts,kte)

              call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,us,vs   &
                               ,u_cup_x,v_cup_x,heso_cup_x,zo_cup_x,po_cup_x,gammao_cup_x,tn_cup_x,psur,tsur  &
                               ,ierr,z1,itf,ktf,its,ite, kts,kte)

             !--- get MSE
             do i=its,itf
                if(ierr(i) /= 0) CYCLE

                call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup_x(i,kts:kte),hkbo_x(i),k22(i))
                hco_x (i,kts:start_level(i)) = hkbo_x(i)

                do k = start_level(i) +1,ktop(i)+1

                   denom=(zuo(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
                   if(denom > 0.0)then
                     hco_x (i,k)=(hco_x(i,k-1)*zuo(i,k-1)-.5*up_massdetro(i,k-1)*hco_x(i,k-1) + &
                                  up_massentro(i,k-1)*heo_x(i,k-1))/ denom
                   else
                     hco_x (i,k)=hco_x(i,k-1)
                   endif
                enddo
                hco_x (i,ktop(i)+2:ktf)=heso_cup_x(i,ktop(i)+2:ktf)
             enddo
             call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
                              ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup_x)
             !--- get cloud work function
             aa_tmp=0.
             call cup_up_aa0(aa_tmp,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
                            ,k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

             if(step==1) aa_ini=aa_tmp ! cloud work function initial
             if(step==2) aa_adv=aa_tmp ! cloud work function modified by advection tendencies
       enddo
       !--- trigger function based on Xie et al 2019
       do i=its,itf
          if(ierr(i) /= 0) cycle

          daa_adv_dt(i)=(aa_adv(i)-aa_ini(i))/dtime

!         if( daa_adv_dt(i) > 0. .and. aa_ini(i) > 0.) cycle !
          if( daa_adv_dt(i) > dcape_threshold/3600. .and. aa_ini(i) > 0.) cycle !
!         if( daa_adv_dt(i) > 0.) cycle

          ierr(i)=90
          ierrc(i) = "dcape trigger not satisfied"

       enddo
       !--- only for output
      if(trim(cumulus) == 'deep') AA0_(:) = daa_adv_dt(:)*3600. ! J/kg/hour

     ENDIF
!
!---
!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(cumulus,ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
                      pwo,ccn,pwevo,edtmax,edtmin,maxens2,edtc,psum,psumh, &
                      rho,aeroevap,itf,ktf,ipr,jpr,its,ite, kts,kte)

     do iedt=1,maxens2
        do i=its,itf
         if(ierr(i).eq.0)then
            edto(i)=sigd(i)*edtc(i,iedt)
            edt (i)=edto(i)
         endif
        enddo
!
!--- get the environmental mass flux
!
        do i=its,itf
         zenv(i,:) = 0.0
         if(ierr(i) /= 0) cycle
         zenv(i,:) = zuo(i,:)-edto(i)*zdo(i,:)
        enddo
!
!--- check mass conservation
!
      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktop(i)
            ! these three are only used at or near mass detrainment and/or entrainment levels
            entupk=0. ; detupk=0. ;entdoj=0.
            ! detrainment and entrainment for downdrafts
            detdo=edto(i)*dd_massdetro(i,k)
            entdo=edto(i)*dd_massentro(i,k)
            ! entrainment/detrainment for updraft
            entup=up_massentro(i,k)
            detup=up_massdetro(i,k)
            ! subsidence by downdrafts only
            subin=-zdo(i,k+1)*edto(i)
            subdown=-zdo(i,k)*edto(i)
            if(k.eq.ktop(i))then
               detupk=zuo(i,ktop(i))
               subin=0.; subdown=0.
               detdo=0.; entdo=0.
               entup=0.; detup=0.
            endif
            totmas=subin-subdown+detup-entup-entdo+ &
                   detdo-entupk-entdoj+detupk+zuo(i,k+1)-zuo(i,k)
            if(abs(totmas).gt.1.e-6)then
              write(6,*)'**mass cons: k,ktop,zo(ktop),totmas,subin,subdown,detup,entup,detdo,entdo,entupk,detupk'
              write(6,123)'mass*1.e+6',k,ktop(i),zo(i,ktop(i)),totmas*1.e+6,subin*1.e+6,subdown*1.e+6,detup*1.e+6,entup*1.e+6&
                         ,detdo*1.e+6,entdo*1.e+6,entupk*1.e+6,detupk*1.e+6
              123     format(1X,A11,2i5,10E12.5)
              ! call error_fatal ( 'totmas .gt.1.e-6' )
            endif
         enddo   ! k
      enddo


!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
        dellu     =0.
        dellv     =0.
        dellah    =0.
        dellat    =0.
        dellaq    =0.
        dellaqc   =0.
        dellabuoy =0.
        subten_H  =0.
        subten_Q  =0.
        subten_T  =0.

!
!----------------------------------------------  cloud level ktop
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level ktop-1
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level k+2
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k+1
!
!----------------------------------------------  cloud level k+1
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level k
!
!----------------------------------------------  cloud level k
!
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!      .               .                 .
!
!----------------------------------------------  cloud level 3  _cup
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 2
!
!----------------------------------------------  cloud level 2  _cup
!
!- - - - - - - - - - - - - - - - - - - - - - - - model level 1

!
!------------------------------------------------------------------------------------------
IF(VERT_DISCR == 0) THEN

      do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktop(i)

            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) -             &
                          zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp        &
                        +(zdo(i,k+1)*(ucd(i,k+1)-u_cup(i,k+1)) -              &
                          zdo(i,k  )*(ucd(i,k  )-u_cup(i,k  )) )*g/dp*edto(i)

            dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) -             &
                          zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp        &
                        +(zdo(i,k+1)*(vcd(i,k+1)-v_cup(i,k+1) ) -             &
                          zdo(i,k  )*(vcd(i,k  )-v_cup(i,k  ) ) )*g/dp*edto(i)

         enddo   ! k
      enddo

      do i=its,itf
        trash  = 0.0
        trash2 = 0.0
        if(ierr(i).eq.0)then
         do k=kts,ktop(i)

            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)

            !---meltglac-------------------------------------------------
            dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                                              - melting(i,k))*g/dp

            !-- for output only
            subten_H(i,k) = -(zuo(i,k+1)*(-heo_cup(i,k+1)) - zuo(i,k)*(-heo_cup(i,k)))*g/dp       &
                            +(zdo(i,k+1)*(-heo_cup(i,k+1)) - zdo(i,k)*(-heo_cup(i,k)))*g/dp*edto(i)

            !- check H conservation
            trash2 = trash2+ (dellah(i,k))*dp/g

             !-- take out cloud liquid/ice water for detrainment
            detup=up_massdetro(i,k)
            if( trim(cumulus) == 'mid'  .or. trim(cumulus) == 'shallow') then

               dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

            elseif(trim(cumulus) == 'deep') then

               if(.not. USE_C1D) then

                  dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

               elseif(c1>0.0) then
                  if(k == ktop(i)) then
                     dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                  else
                     dz=zo_cup(i,k+1)-zo_cup(i,k)
                     dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g
                  endif
               else
                  if(k == ktop(i)) then
                     dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                  else
                     dz=zo_cup(i,k+1)-zo_cup(i,k)
                      dellaqc(i,k) = ( zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g  + &
                                       detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp )*0.5
                  endif
               endif
            endif
            !
            !---
            G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
            E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
            !
            !print*,"eva=",k,pwdo(i,k),E_dn,zdo(i,k  ),G_rain
            !
            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid/ice water (qrco) + converted to rain

            C_up = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) -       &
                                 zuo(i,k  )* qrco(i,k  )  )*g/dp + G_rain

            !-- water vapor budget
            !-- = flux divergence z*(Q_c - Q_env)_up_and_down &
            !--   - condensation term + evaporation
            dellaq(i,k) =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                           zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                           zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)    &
                         - C_up + E_dn

            !-- for output only
            subten_Q(i,k) =-(zuo(i,k+1)*(-qo_cup(i,k+1)) - zuo(i,k)*(-qo_cup(i,k)))*g/dp       &
                           +(zdo(i,k+1)*(-qo_cup(i,k+1)) - zdo(i,k)*(-qo_cup(i,k)))*g/dp*edto(i)

            !- check water conservation liq+condensed (including rainfall)
            trash= trash+ (dellaq(i,k)+dellaqc(i,k)+ G_rain-E_dn)*dp/g

            !---
            dellabuoy(i,k) = edto(i)*dd_massdetro(i,k)*0.5*(dbydo(i,k+1)+dbydo(i,k))*g/dp
            !---

            !write(3,*)'=>H= ',k,real(trash2,4),real(dellah(i,k),4)
            !write(4,*)'=>W= ',k,real(trash,4),real(dellaq(i,k),4)
         enddo   ! k
         !--- test only with double precision:
         !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
         !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
         !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
         !    !stop 33
         !endif

        endif

      enddo

ELSEIF(VERT_DISCR == 1) THEN

     !---- convective transport of momentum
     if(alp1 == 0.) then !-- fully time explicit
        do i=its,itf
          if(ierr(i) /= 0) cycle
          do k=kts,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            dellu(i,k) =-(zuo(i,k+1)*(uc (i,k+1)-u_cup(i,k+1) ) -             &
                          zuo(i,k  )*(uc (i,k  )-u_cup(i,k  ) ) )*g/dp        &
                        +(zdo(i,k+1)*(ucd(i,k+1)-u_cup(i,k+1)) -              &
                          zdo(i,k  )*(ucd(i,k  )-u_cup(i,k  )) )*g/dp*edto(i)

            dellv(i,k) =-(zuo(i,k+1)*(vc (i,k+1)-v_cup(i,k+1) ) -             &
                          zuo(i,k  )*(vc (i,k  )-v_cup(i,k  ) ) )*g/dp        &
                        +(zdo(i,k+1)*(vcd(i,k+1)-v_cup(i,k+1) ) -             &
                          zdo(i,k  )*(vcd(i,k  )-v_cup(i,k  ) ) )*g/dp*edto(i)
          enddo   ! k
        enddo

     elseif (alp1 > 0.) then              !-- time alp0*explict + alp1*implicit + upstream

       alp0=1.-alp1
       do i=its,itf
          if(ierr(i) /= 0) cycle
          do k=kts,ktop(i)+1
             fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
             fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
          enddo

          do k=kts,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            beta1 = dtime*g/dp
            aa(k) =    alp1*beta1*fm(k)
            bb(k) = 1.+alp1*beta1*(fp(k)-fm(k+1))
            cc(k) =   -alp1*beta1*fp(k+1)

            ddu(k) = us(i,k)-( zuo(i,k+1)*uc (i,k+1)-zuo(i,k)*uc (i,k) )*beta1 + &
                             ( zdo(i,k+1)*ucd(i,k+1)-zdo(i,k)*ucd(i,k) )*beta1*edto(i)

            ddu(k) = ddu(k) + alp0*beta1*(-fm(k)*us(i,max(kts,k-1)) +(fm(k+1)-fp(k))*us(i,k) +fp(k+1)*us(i,k+1))


            ddv(k) = vs(i,k)-( zuo(i,k+1)*vc (i,k+1)-zuo(i,k)*vc (i,k) )*beta1 + &
                             ( zdo(i,k+1)*vcd(i,k+1)-zdo(i,k)*vcd(i,k) )*beta1*edto(i)

            ddv(k) = ddv(k) + alp0*beta1*(-fm(k)*vs(i,max(kts,k-1)) +(fm(k+1)-fp(k))*vs(i,k) +fp(k+1)*vs(i,k+1))

           !print*,"trX4=",k,aa(k),bb(k),cc(k)!, 1.+alp1*beta1*zenv(i,k  ), -alp1*beta1*zenv(i,k+1)
         enddo
         call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddu (kts:ktop(i)))
         dellu(i,kts:ktop(i))=(ddu(kts:ktop(i))-us(i,kts:ktop(i)))/dtime

         call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddv (kts:ktop(i)))
         dellv(i,kts:ktop(i))=(ddv(kts:ktop(i))-vs(i,kts:ktop(i)))/dtime
       enddo
     endif


     !--- convective transport of MSE and Q/Qc
!     if(USE_FLUX_FORM == 1) then
      do i=its,itf
        if(ierr(i) /= 0) cycle

        !--- moist static energy : flux form + source/sink terms + time explicit
        !
!        if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
        if(use_fct == 0 ) then

          do k=kts,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            dellah(i,k) =-(zuo(i,k+1)*(hco (i,k+1)-heo_cup(i,k+1) ) -                    &
                           zuo(i,k  )*(hco (i,k  )-heo_cup(i,k  ) ) )*g/dp            &
                         +(zdo(i,k+1)*(hcdo(i,k+1)-heo_cup(i,k+1) ) -                    &
                           zdo(i,k  )*(hcdo(i,k  )-heo_cup(i,k  ) ) )*g/dp*edto(i)

            dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))* &
                                       0.5*(qrco(i,k+1)+qrco(i,k)) - melting(i,k))*g/dp

            !--- for output only
            subten_H(i,k) = -(zuo(i,k+1)*(-heo_cup(i,k+1)) - zuo(i,k)*(-heo_cup(i,k)))*g/dp       &
                            +(zdo(i,k+1)*(-heo_cup(i,k+1)) - zdo(i,k)*(-heo_cup(i,k)))*g/dp*edto(i)
          enddo   ! k

         else

         !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
          sub_tend (1,:) = 0. ! dummy array
          trcflx_in(1,:) = 0. ! dummy array
          massflx  (i,:) = 0.
          dtime_max      = dtime

          do k=kts,ktop(i)
               dp=100.*(po_cup(i,k)-po_cup(i,k+1))
               trcflx_in (1,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*heo_cup(i,k) !* xmb(i)
               massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))                    !* xmb(i)
               dtime_max=min(dtime_max,.5*dp)
          enddo
          call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),heo(i,:),massflx(i,:),trcflx_in(1,:),sub_tend(1,:))

          do k=kts,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            dellah(i,k) =-( zuo(i,k+1)*hco (i,k+1) - zuo(i,k)*hco (i,k) )*g/dp           &
                         +( zdo(i,k+1)*hcdo(i,k+1) - zdo(i,k)*hcdo(i,k) )*g/dp*edto(i)

            dellah(i,k) = dellah(i,k) + xlf*((1.-p_liq_ice(i,k))* &
                                        0.5*(qrco(i,k+1)+qrco(i,k)) - melting(i,k))*g/dp
            !- update with subsidence term from the FCT scheme
            dellah(i,k) = dellah(i,k) + sub_tend(1,k)
            !--- for output only
            subten_H(i,k) = sub_tend(1,k)
          enddo   ! k
         endif
      enddo

!     elseif(USE_FLUX_FORM == 2) THEN
!
!        !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
!        alp0=1.
!        do i=its,itf
!          if(ierr(i) /= 0) cycle
!          do istep=1,-1, -2
!
!            if(istep == 1) then
!               ddu(:) = heo(i,:)
!                do k=kts,ktop(i)+1
!                  fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
!                  fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
!                enddo
!             else
!               ddu(kts:ktop(i)+1) = heo(i,kts:ktop(i)+1) + dellah(i,kts:ktop(i)+1)*dtime
!               zenv_diff(1,kts) = 0.
!               do k=kts,ktop(i)+1
!                   dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
!                   zenv_diff (1,k+1) = 1.06* ( dp*abs(zenv(i,k+1))/g - dtime*zenv(i,k+1)**2 )/dp/g &
!                                     * (ddu(k+1) - ddu(k)) /(ddu(k+1) + ddu(k) + 1.e-16)
!               enddo
!               do k=kts,ktop(i)+1
!                   fp(k) = 0.5*(zenv_diff(1,k)+abs(zenv_diff(1,k)))
!                   fm(k) = 0.5*(zenv_diff(1,k)-abs(zenv_diff(1,k)))
!               enddo
!             endif
!             do k=kts,ktop(i)
!               dp=100.*(po_cup(i,k)-po_cup(i,k+1))
!               beta1 = dtime*g/dp
!               ddh(k) = ddu(k) + alp0*beta1*( -fm(k)*ddu(max(kts,k-1)) + (fm(k+1)-fp(k))*ddu(k) + fp(k+1)*ddu(k+1) )
!             enddo
!
!             dellah(i,kts:ktop(i)+1)=(ddh(kts:ktop(i)+1)-heo(i,kts:ktop(i)+1))/dtime
!
!           enddo
!
!           do k=kts,ktop(i)
!             dp=100.*(po_cup(i,k)-po_cup(i,k+1))
!             beta1 = g/dp
!
!             ddh(k) =  -( zuo(i,k+1)*hco (i,k+1) - zuo(i,k)*hco (i,k) )*beta1            &
!                       +( zdo(i,k+1)*hcdo(i,k+1) - zdo(i,k)*hcdo(i,k) )*beta1*edto(i)
!
!             ddh(k) = ddh(k) + xlf*((1.-p_liq_ice(i,k))* &
!                              0.5*(qrco(i,k+1)+qrco(i,k)) - melting(i,k))*beta1
!
!             dellah(i,k) =  dellah(i,k) + ddh(k)
!
!           enddo
!       enddo
!     endif
      !-------------------------
      !--- water vapor + condensates : flux form + source/sink terms + time explicit
      do i=its,itf

        if(ierr(i) /= 0) cycle
         !
         do k=kts,ktop(i)
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            !-- take out cloud liquid/ice water for detrainment
            detup=up_massdetro(i,k)
            if( trim(cumulus) == 'mid'  .or. trim(cumulus) == 'shallow') then

               dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

            elseif(trim(cumulus) == 'deep') then

               if(.not. USE_C1D) then

                  dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp

               elseif(c1>0.0) then
                  if(k == ktop(i)) then
                     dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                  else
                     dz=zo_cup(i,k+1)-zo_cup(i,k)
                     dellaqc(i,k) = zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g
                  endif
               else
                  if(k == ktop(i)) then
                     dellaqc(i,k) = detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp
                  else
                     dz=zo_cup(i,k+1)-zo_cup(i,k)
                      dellaqc(i,k) = ( zuo(i,k)*c1d(i,k)*qrco(i,k)*dz/dp*g  + &
                                       detup*0.5*(qrco(i,k+1)+qrco(i,k)) *g/dp )*0.5
                  endif
               endif
            endif
            !
            !---
            G_rain=  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
            E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i) ! pwdo < 0 and E_dn must > 0
            !
            !-- condensation source term = detrained + flux divergence of
            !-- cloud liquid/ice water (qrco) + converted to rain

            C_up = dellaqc(i,k)+(zuo(i,k+1)*qrco(i,k+1) - zuo(i,k)* qrco(i,k))*g/dp + G_rain

            !-- water vapor budget
            !-- = flux divergence z*(Q_c - Q_env)_up_and_down  - condensation term + evaporation
            dellaq(i,k) =-(zuo(i,k+1)*qco (i,k+1) - zuo(i,k)*qco (i,k))*g/dp                &
                         +(zdo(i,k+1)*qcdo(i,k+1) - zdo(i,k)*qcdo(i,k))*g/dp*edto(i)    &
                         - C_up + E_dn

            !--- source of cold pools
            dellabuoy(i,k)=edto(i)*dd_massdetro(i,k)*0.5*(dbydo(i,k+1)+dbydo(i,k))*g/dp

         enddo
!        if(use_fct == 0 .or. adjustl(cumulus) == 'shallow') then
         if(use_fct == 0 ) then
              do k=kts,ktop(i)
                 dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                 sub_tend(1,k) =-(zuo(i,k+1)*(-qo_cup(i,k+1)) - zuo(i,k)*(-qo_cup(i,k)))*g/dp       &
                                +(zdo(i,k+1)*(-qo_cup(i,k+1)) - zdo(i,k)*(-qo_cup(i,k)))*g/dp*edto(i)
              enddo
         else
              !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
              sub_tend (1,:) = 0. ! dummy array
              trcflx_in(1,:) = 0. ! dummy array
              massflx  (i,:) = 0.
              dtime_max      = dtime

              do k=kts,ktop(i)
                 dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                 trcflx_in (1,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*qo_cup(i,k) !* xmb(i)
                 massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))             !* xmb(i)
                 dtime_max=min(dtime_max,.5*dp)
              enddo
              call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),qo(i,:),massflx(i,:),trcflx_in(1,:),sub_tend(1,:))
          endif

         !--- add the contribuition from the environ subsidence
         dellaq(i,kts:ktop(i)) = dellaq(i,kts:ktop(i)) + sub_tend(1,kts:ktop(i))

         !--- for output only
         subten_Q (i,kts:ktop(i)) = sub_tend(1,kts:ktop(i))

         !     do k=kts,ktop(i)
          !     print*,"delq=",use_fct,k,dellaq(i,k) , sub_tend(1,k)
         !     enddo

         !- check H and water conservation liq+condensed (including rainfall)
         trash  = 0. ; trash2 = 0.0
         do k=kts,ktop(i)
                 dp     = 100.*(po_cup(i,k)-po_cup(i,k+1))
                 G_rain =  0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
                 E_dn   = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i)
                 trash  = trash + (dellaq(i,k) + dellaqc(i,k)+ G_rain-E_dn)*dp/g
                 trash2 = trash2+  dellah(i,k)*g/dp + xlf*((1.-p_liq_ice(i,k))*0.5*(qrco(i,k+1)+qrco(i,k)) &
                                   - melting(i,k))*g/dp
         enddo   ! k
         !--- test only with double precision:
         !write(0,*)'=>H/W-FINAL= ',real(trash2,4),real(trash,4),k22(i),kbcon(i),ktop(i)
         !if(abs(trash)>1.e-6 .or. abs(trash2) > 1.e-6) then
         !    write(0,*)'=> not water mass or H cons for deep= ',i,trash,trash2
         !    !stop 33
         !endif

      enddo


ENDIF ! vertical discretization formulation

!
!--- apply environmental subsidence on grid-scale ice and liq water contents, and cloud fraction (Upwind scheme)
!
      IF(APPLY_SUB_MP == 1) THEN
       dellampqi =0.
       dellampql =0.
       dellampcf =0.

       do i=its,itf
         if(ierr(i) /= 0) cycle
         do k=kts,ktop(i)
            dp=100.*(po_cup(i,k  )-po_cup(i,k+1))

           !--- apply environmental subsidence on grid-scale/anvil ice and liq water contents (Upwind scheme)
           !
            env_mf   = - 0.5* (zenv(i,k+1) + zenv(i,k))
            env_mf_m = min(env_mf,0.)*g/dp
            env_mf_p = max(env_mf,0.)*g/dp

            dellampqi(:,i,k) = - (  env_mf_m*(mpqi(:,i,k+1)-mpqi(:,i,k))  +            &
                                    env_mf_p*(mpqi(:,i,k  )-mpqi(:,i,max(k-1,kts))))

            dellampql(:,i,k) = - (  env_mf_m*(mpql(:,i,k+1)-mpql(:,i,k))  +               &
                                    env_mf_p*(mpql(:,i,k  )-mpql(:,i,max(k-1,kts))))

           !--- apply environmental subsidence on grid-scale/anvil cloud fraction
           !
            dellampcf(:,i,k) = - (  env_mf_m*(mpcf(:,i,k+1)-mpcf(:,i,k))  +               &
                                    env_mf_p*(mpcf(:,i,k  )-mpcf(:,i,max(k-1,kts))))
         enddo

         !--- apply environmental subsidence on grid-scale and anvil cloud fraction using time implicit/explict method
         if(alp1 > 0.) then
           alp0=1.0-alp1
           do k=kts,ktop(i)
                   dp=100.*(po_cup(i,k  )-po_cup(i,k+1))
                   env_mf   = - 0.5* (zenv(i,k+1) + zenv(i,k))
                   env_mf_m = min(env_mf,0.)*g/dp
                   env_mf_p = max(env_mf,0.)*g/dp

                   beta1 = -env_mf_m
                   beta2 = -env_mf_p

                   aa(k) =    alp1*beta2             ! coef of f(k-1,t+1),
                   bb(k) = 1.+alp1*beta1-alp1*beta2  ! coef of f(k  ,t+1),
                   cc(k) =   -alp1*beta1             ! coef of f(k+1,t+1),

                   !-- this is the rhs of the discretization
                   dd(:,k) = (1.-alp0*beta1+alp0*beta2)*mpcf(:,i,k  ) +&       ! coef of  f(k  ,t),
                                            alp0*beta1 *mpcf(:,i,k+1) -&       ! coef of  f(k+1,t),
                                            alp0*beta2 *mpcf(:,i,max(kts,k-1)) ! coef of  f(k-1,t),
           enddo
           do kmp =1,nmp
             !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
             call tridiag (ktop(i),aa(kts:ktop(i)), bb(kts:ktop(i)), cc(kts:ktop(i)), dd(kmp,kts:ktop(i)))

             dellampcf(kmp,i,kts:ktop(i)) = dd(kmp,kts:ktop(i))-mpcf(kmp,i,kts:ktop(i))
           enddo
         endif
       enddo
      ENDIF
!
!
!--- make the smoothness procedure
!
      IF(USE_SMOOTH_TEND >= 1) THEN
        do i=its,itf
           if(ierr(i) /= 0) cycle
           tend2d=0.

           do k=kts,ktop(i)
             rcount = 1.e-8 ; tend1d=0.
             do kk= max(kts,k-USE_SMOOTH_TEND),min(ktop(i),k+USE_SMOOTH_TEND)
                  dp=(po_cup(i,kk)-po_cup(i,kk+1))
                  rcount     = rcount     +  dp
                  tend1d(1)  = tend1d(1)  +  dp* DELLAH  (i,kk)
                  tend1d(2)  = tend1d(2)  +  dp* DELLAQ  (i,kk)
                  tend1d(3)  = tend1d(3)  +  dp* DELLAQC (i,kk)
                  tend1d(4)  = tend1d(4)  +  dp* DELLU   (i,kk)
                  tend1d(5)  = tend1d(5)  +  dp* DELLV   (i,kk)

             enddo
             tend2d(k,1:5)  = tend1d(1:5) /rcount
          enddo
          !--- get the final/smoother tendencies
          do k=kts,ktop(i)
            DELLAH  (i,k) = tend2d(k,1)
            DELLAQ  (i,k) = tend2d(k,2)
            DELLAQC (i,k) = tend2d(k,3)
            DELLU   (i,k) = tend2d(k,4)
            DELLV   (i,k) = tend2d(k,5)
          enddo
        enddo
      ENDIF
!
!--- using dellas, calculate changed environmental profiles
!
      do k=kts,ktf
       do i=its,itf
         dellat(i,k)=0.
         if(ierr(i) /= 0) cycle
         !
         XHE(I,K)=(DELLAH(I,K)             )*MBDT(i)+HEO(I,K)
         XQ (I,K)=(DELLAQ(I,K)+DELLAQC(i,k))*MBDT(i)+QO(I,K)
         if(XQ(I,K).LE.0.)XQ(I,K)=1.E-08

         !- do not feed dellat with dellaqc if the detrainment of liquid water
         !- will be used as a source for cloud microphysics
         if(COUPL_MPHYSICS) then
           DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xlv*DELLAQ(I,K))
         else
         !---meltglac-------------------------------------------------
           DELLAT (I,K)=(1./cp)*( DELLAH(I,K) - xlv*(DELLAQ(I,K) + DELLAQC(i,k))*(1.+(xlf/xlv)*(1.-p_liq_ice(i,k))))
          !DELLAT (I,K)=(1./cp)*( DELLAH(I,K)  -xlv*(DELLAQ(I,K) + DELLAQC(i,k)))

           !-adding dellaqc to dellaq:
           DELLAQ (I,K)= DELLAQ(I,K)+DELLAQC(I,K)
           DELLAQC(I,K)= 0.0
         endif
         !---meltglac-------------------------------------------------
         XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)*(1.+(xlf/xlv)*(1.-p_liq_ice(i,k)))))*MBDT(i) &
                + TN(I,K)
        !XT(I,K)=((1./cp)*DELLAH(i,k)-(xlv/cp)*(DELLAQ(i,k)+DELLAQC(i,k)))*MBDT(i)+TN(I,K)

         !--- temp tendency due to the environmental subsidence
         subten_T(i,k)=(1./cp)*(subten_H(i,k)-xlv*subten_Q(i,k))

       enddo
      enddo
      do i=its,itf
         if(ierr(i) /= 0) cycle
         !XHKB(I)=(dsubh(i,k22(i))+DELLAH(I,K22(i)))*MBDT+HKBO(I)
         XHE(I,ktf)=HEO(I,ktf)
         XQ (I,ktf)=QO(I,ktf)
         XT (I,ktf)=TN(I,ktf)
         IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      enddo
      !- new way for defining XHKB
      do i=its,itf
         if(ierr(i) /= 0)CYCLE
         !XHKB(I)= DELLAH(I,K22(i))*MBDT+HKBO(I)
         !-note that HKBO already contains the contribuition from
         !-ztexec and zqexec
         call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),DELLAH (i,kts:kte),DELLAH_aver,k22(i))
         XHKB(I)= DELLAH_aver*MBDT(i) + HKBO(I)
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1,psur,ierr,-1,itf,ktf,its,ite,kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, xhe_cup,   &
                        us,vs,u_cup,v_cup,                                    &
                        xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,tsur,    &
                        ierr,z1,itf,ktf,its,ite, kts,kte)
!
!--- static control
!
!--- moist static energy inside cloud
!
      do i=its,itf
        xhc(i,:)=0.
        if(ierr(i)/= 0)CYCLE
        do k=kts,start_level(i) !k22(i)
            xhc(i,k)=xhkb(i)
        enddo
      enddo
!
      do i=its,itf
       if(ierr(i)/=0)cycle
       do k=start_level(i)  +1,ktop(i)+1  ! mass cons option
         denom= (xzu(i,k-1)-.5*up_massdetro(i,k-1)+up_massentro(i,k-1))
         if(denom==0.0) then
           xhc(i,k)= xhc(i,k-1)
         else
           xhc(i,k)=(xhc(i,k-1)*xzu(i,k-1)-.5*up_massdetro(i,k-1)*xhc(i,k-1)+ &
                                              up_massentro(i,k-1)*xhe(i,k-1)) / denom
          if(k==start_level(i)+1)  then
             x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
             xhc(i,k)= xhc(i,k) + x_add*up_massentro(i,k-1)/denom
          endif
         endif
!
!- include glaciation effects on XHC
!                                   ------ ice content --------
          xhc (i,k)= xhc (i,k)+ xlf*(1.-p_liq_ice(i,k))*qrco(i,k)
        enddo
        do k=ktop(i)+2,ktf
           xHC (i,k)=xhes_cup(i,k)
           xzu (i,k)=0.
        enddo
      enddo
      call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop,xhc,xhe_cup,xhes_cup,xdby,xz_cup)
!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz_cup,xzu,xdby,GAMMA_CUP,xt_cup, k22,klcl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte)

      do nens=1,maxens
       do i=its,itf
         if(ierr(i) /= 0) cycle
         !~ xaa0_ens(i,nens)=xaa0(i)
         do k=kts,ktop(i)
             do nens3=1,maxens3
               if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
!--- b=beta
               else if(nens3.eq.8)then
                 pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
!--- b=beta/2
               else if(nens3.eq.9)then
                 pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
               else
                 pr_ens(i,nens3)=pr_ens(i,nens3) +pwo(i,k)+edto(i)*pwdo(i,k)
               endif
             enddo
         enddo
         if(pr_ens(i,7).lt.1.e-6 .and. c0_mid > 0. .and.  trim(cumulus) /= 'shallow' )then
            ierr(i)=18
            ierrc(i)="total normalized condensate too small"
            do nens3=1,maxens3
               pr_ens(i,nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,nens3).lt.1.e-5)then
            pr_ens(i,nens3)=0.
           endif
         enddo
       enddo
      enddo
!
!--- LARGE SCALE FORCING
!
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
      enddo
!
!--- calculate cloud base mass flux
!
      if(trim(cumulus) == 'deep') &
      call cup_forcing_ens_3d(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3 &
                             ,ierr,ierr2,ierr3,k22,kbcon,ktop         &
                             ,xland1,aa0,aa1,xaa0,mbdt,dtime          &
                             ,xf_ens,mconv,qo                         &
                             ,po_cup,omeg,zdo,zuo,pr_ens,edto         &
                             ,tau_ecmwf,aa1_bl,xf_dicycle, xk_x)
!
      if(trim(cumulus) == 'mid') &
      call cup_forcing_ens_3d_mid(aa0,aa1,xaa0,mbdt,dtime,ierr          &
                                 ,po_cup,ktop,k22,kbcon,kpbl,ichoice,maxens,maxens3   &
                                 ,itf,ktf,its,ite, kts,kte                      &
                                 ,tau_ecmwf,aa1_bl,xf_dicycle                   &
                                 ,dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup)



      if(trim(cumulus) == 'shallow') then
       call cup_up_cape(cape,z,zu,dby,gamma_cup,t_cup,k22,kbcon,ktop,ierr    &
                            ,tempco,qco,qrco,qo_cup,itf,ktf,its,ite, kts,kte )

       call cup_forcing_ens_3d_shal(itf,ktf,its,ite,kts,kte,dtime,ichoice    &
                                  ,ierrc,ierr,klcl,kpbl,kbcon,k22,ktop              &
                                  ,xmb,tsur,cape,h_sfc_flux,le_sfc_flux,zws  &
                                  ,po, hco, heo_cup,po_cup,t_cup,dhdt,rho    &
                                  ,xff_shal,xf_dicycle)
      endif


      do k=kts,ktf
        do i=its,itf
          if(ierr(i) /= 0) cycle
          pwo_eff(i,k)=pwo(i,k)+edto(i)*pwdo(i,k)
        enddo
      enddo

     ENDDO
!
!--- Include kinetic energy dissipation converted to heating
!
       call ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr &
                         ,po_cup,us,vs,dellu,dellv,dellat)
!
!--- FEEDBACK
!
       call cup_output_ens_3d(cumulus,xff_shal,xff_mid,xf_ens,ierr,dellat,dellaq,          &
                              dellaqc,outt, outq,outqc,zuo,pre,pwo_eff,xmb,ktop,           &
                              maxens2,maxens,ierr2,ierr3,                                  &
                              pr_ens,maxens3,ensdim,sig,xland1,                            &
                              ichoice,ipr,jpr,itf,ktf,its,ite, kts,kte,                    &
                              xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon,         &
                              dellabuoy,outbuoy,                                           &
                              dellampqi,outmpqi,dellampql,outmpql,dellampcf,outmpcf,nmp    )

!
!
!--- get the net precipitation flux (after downdraft evaporation)
       call get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb  &
                             ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco          &
                             ,prec_flx,evap_flx                               &
                             ,itf,ktf,its,ite, kts,kte)

!
!--- rainfall evap below cloud base
!
       if(USE_REBCB == 1)                                                              &
       call rain_evap_below_cloudbase(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop&
                                      ,xmb,psur,xland,qo_cup,t_cup                     &
                                      ,po_cup,qes_cup,pwavo,edto,pwevo,pwo,pwdo        &
                                      ,pre,prec_flx,evap_flx,outt,outq,outbuoy,evap_bcb)


!
!--- includes effects of the remained cloud dissipation into the enviroment
!
       if(use_cloud_dissipation >= 0.)                                             &
       call cloud_dissipation(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop    &
                             ,dtime,xmb,xland,qo_cup,qeso_cup,po_cup,outt,outq     &
                             ,outqc,zuo,vvel2d,rho_hydr,qrco,sig,tempco,qco,tn_cup &
                             ,heso_cup,zo)

!
!--- get the total (deep+congestus) evaporation flux for output (units kg/kg/s)
!
       do i=its,itf
           if(ierr(i) /= 0) cycle
           do k=kts,ktop(i)
              dp=100.*(po_cup(i,k)-po_cup(i,k+1))
              !--- add congestus and deep plumes, and convert to kg/kg/s
              revsu_gf(i,k) = revsu_gf(i,k) + evap_flx(i,k)*g/dp
           enddo
       enddo
!
!
!--- get lightning flashes density (parameterization from Lopez 2016, MWR)
!
       if( LIGHTNING_DIAG == 1 .and. trim(cumulus) == 'deep') then
            call cup_up_cape(cape,z,zu,dby,gamma_cup,t_cup,k22,kbcon,ktop,ierr &
                            ,tempco,qco,qrco,qo_cup,itf,ktf,its,ite, kts,kte   )

            call cup_up_lightning(itf,ktf,its,ite, kts,kte, ierr, kbcon,ktop,xland,cape &
                                 ,cnvfrc,srftype,zo,zo_cup,t_cup,t,tempco,qrco,po_cup,rho,prec_flx     &
                                 ,lightn_dens                                        )
       endif
!
!--- for outputs (only deep plume)
!
       if( trim(cumulus) == 'deep') then
          do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=kts,ktop(i)+1
               prfil_gf  (i,k) = prec_flx(i,k)
            enddo
          enddo
       endif
!
!--- for tracer convective transport / outputs
!
       do i=its,itf
          if(ierr(i) /= 0) cycle
          do k=kts,ktf
          !tup          (i,k) = (1./cp)*(hco(i,k)-g*zo_cup(i,k)-xlv*qco(i,k))!in-updraft temp
           tup          (i,k) = tempco(i,k) !in-updraft temp
          enddo
          tup           (i,kte) = t_cup(i,kte)
       enddo

!--- convert mass fluxes, etc...
   do i=its,itf
        if(ierr(i) /= 0) cycle
        pwavo       (i)   = xmb(i)*pwavo       (i)
        pwevo       (i)   = xmb(i)*pwevo       (i)
        zuo         (i,:) = xmb(i)*zuo         (i,:)
        zdo         (i,:) = xmb(i)*zdo         (i,:)
        pwo         (i,:) = xmb(i)*pwo         (i,:)
        pwdo        (i,:) = xmb(i)*pwdo        (i,:)
        up_massentro(i,:) = xmb(i)*up_massentro(i,:)
        up_massdetro(i,:) = xmb(i)*up_massdetro(i,:)
        dd_massentro(i,:) = xmb(i)*dd_massentro(i,:)
        dd_massdetro(i,:) = xmb(i)*dd_massdetro(i,:)
        zenv        (i,:) = xmb(i)*zenv        (i,:)
   enddo

!--for output only.
   do i=its,itf
        subten_Q    (i,:) = xmb(i)*subten_Q    (i,:)
        subten_H    (i,:) = xmb(i)*subten_H    (i,:)
        subten_T    (i,:) = xmb(i)*subten_T    (i,:)
   enddo

!
!--- outputs a model sounding for the stand-alone code (part 2)
!
   IF(OUTPUT_SOUND == 1) THEN
       call SOUND(2,cumulus,int_time,dtime,ens4,itf,ktf,its,ite, kts,kte,xlats,xlons,jcol,whoami_all          &
                 ,z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,zo,qeso,heo,heso,tn,qo,us,vs ,omeg,xz     &
                 ,h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec,zqexec, xland      &
                 ,kpbl,k22,klcl,kbcon,ktop,aa0,aa1,sig,xaa0,hkb,xmb,pre,edto                   &
                 ,zo_cup,dhdt,rho,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv)
   ENDIF

   if( trim(cumulus) == 'deep') then
    do i=its,itf
            AA1_(i) =0.
            AA0_(i) =0.
            if(ierr(i) /= 0) cycle
!           AA1_(i) = depth_neg_buoy   (i)
!           AA0_(i) = frh_bcon(i)
            AA1_(i) = AA1(i)
            AA0_(i) = AA0(i)
    enddo
   endif

   IF(LIQ_ICE_NUMBER_CONC == 1) THEN
     call get_liq_ice_number_conc(itf,ktf,its,ite, kts,kte,ierr,ktop &
                                 ,cnvfrc,srftype,dtime,rho,outqc,tempco,outnliq,outnice)
   ENDIF
!
!
!--------------------------------------------------------------------------------------------!
!- section for atmospheric composition
!--------------------------------------------------------------------------------------------!
IF(USE_TRACER_TRANSP==1)  THEN

!-1) get mass mixing ratios at the cloud levels

   call cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

!-2) determine in-cloud tracer mixing ratios
!
! a) chem - updraft
   !- note: here "sc_up_chem" stores the total in-cloud tracer mixing ratio (i.e., including the portion
   !        embedded in the condensates).
   call get_incloud_sc_chem_up(cumulus,mtp,se_chem,se_cup_chem,sc_up_chem,pw_up_chem,tot_pw_up_chem      &
                             ,zo_cup,rho,po,po_cup,qco,qrco,tempco,pwo,zuo,up_massentro,up_massdetro               &
                             ,vvel2d,vvel1d,start_level,k22,kbcon,ktop,klcl,ierr,xland,itf,ktf,its,ite, kts,kte)

! b) chem - downdraft
   call get_incloud_sc_chem_dd(cumulus,mtp,se_chem,se_cup_chem,sc_dn_chem,pw_dn_chem ,pw_up_chem,sc_up_chem &
                             ,tot_pw_up_chem,tot_pw_dn_chem                                                       &
                             ,zo_cup,rho,po_cup,qrcdo,pwdo,pwevo,edto,zdo,dd_massentro,dd_massdetro ,pwavo,pwo     &
                             ,jmin,ierr,itf,ktf,its,ite, kts,kte)
!
!-3) determine the vertical transport including mixing, scavenging and evaporation
!
!---a) change per unit mass that a model cloud would modify the environment
   do i=its,itf
        if(ierr(i) /= 0) cycle

        !- flux form + source/sink terms + time explicit + FCT
        IF(USE_FLUX_FORM == 1 .and. alp1 == 0. ) THEN

          if(use_fct == 0 ) then
            do k=kts,ktop(i)
                dp=100.*(po_cup(i,k)-po_cup(i,k+1))

                out_chem(:,i,k) =-(zuo(i,k+1)*(sc_up_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                                   zuo(i,k  )*(sc_up_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp             &
                                 +(zdo(i,k+1)*(sc_dn_chem(:,i,k+1)-se_cup_chem(:,i,k+1) ) -                 &
                                   zdo(i,k  )*(sc_dn_chem(:,i,k  )-se_cup_chem(:,i,k  ) ))*g/dp*edto(i)
            enddo

          else

            !-- FCT scheme for the subsidence transport: d(M_env*S_env)/dz
            sub_tend  = 0.
            trcflx_in = 0.
            dtime_max = dtime
            massflx (i,:)=0.

            do k=kts+1,ktop(i)+1
               dp                 = 100.*(po_cup(i,k)-po_cup(i,k+1))
               trcflx_in (:,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))*se_cup_chem(:,i,k) !* xmb(i)
               massflx   (i,k) =-(zuo(i,k)  -edto(i)*zdo(i,k))                       !* xmb(i)
               dtime_max=min(dtime_max,.5*dp)
            enddo
            !- if dtime_max<dtime => needs a loop to update from t to t+dtime (check this!)
            !if( dtime_max < dtime ) stop "dtime_max < dtime in GF scheme"

            do ispc=1,mtp
               call fct1d3 (ktop(i),kte,dtime_max,po_cup(i,:),se_chem(ispc,i,:),massflx(i,:),trcflx_in(ispc,:),sub_tend(ispc,:))
            enddo

            do k=kts,ktop(i)
              dp=100.*(po_cup(i,k)-po_cup(i,k+1))
              out_chem(:,i,k) = -(zuo(i,k+1)*(sc_up_chem(:,i,k+1)) - zuo(i,k)*(sc_up_chem(:,i,k)))*g/dp          &
                                +(zdo(i,k+1)*(sc_dn_chem(:,i,k+1)) - zdo(i,k)*(sc_dn_chem(:,i,k)))*g/dp*edto(i)

              !- update with the subsidence term from FCT scheme
               out_chem(:,i,k) = out_chem(:,i,k) + sub_tend(:,k)

            enddo
          endif

          !- include evaporation
          if(USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
            do k=kts,ktop(i)
                 dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                 out_chem(:,i,k) = out_chem(:,i,k)    &
                                 - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*g/dp !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
            enddo
          endif

          !- include scavenging
          if(USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
            do k=kts,ktop(i)
                 dp=100.*(po_cup(i,k)-po_cup(i,k+1))
                 out_chem(:,i,k) = out_chem(:,i,k) &
                                 - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*g/dp  ! incorporated in rainfall (<0)
            enddo
          endif
        ENDIF

        !- flux form + source/sink terms + time explicit/implicit + upstream
        IF(USE_FLUX_FORM == 1 .and. alp1 > 0. ) THEN

            alp0=1.-alp1
            do k=kts,ktop(i)+1
                fp(k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                fm(k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
            enddo

            do k=kts,ktop(i)
              dp=100.*(po_cup(i,k)-po_cup(i,k+1))
              beta1 = dtime*g/dp
              aa(k) =    alp1*beta1*fm(k)
              bb(k) = 1.+alp1*beta1*(fp(k)-fm(k+1))
              cc(k) =   -alp1*beta1*fp(k+1)

              ddtr(:,k) = se_chem(:,i,k) - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1      &
                                         + (zdo(i,k+1)*sc_dn_chem(:,i,k+1) - zdo(i,k)*sc_dn_chem(:,i,k))*beta1*edto(i)

              !- include evaporation
              if(USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                 out_chem(:,i,k) = out_chem(:,i,k)    &
                                 - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
              endif

              !- include scavenging
              if(USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                 out_chem(:,i,k) = out_chem(:,i,k) &
                                 - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*beta1  ! incorporated in rainfall (<0)
              endif

              ddtr(:,k) = ddtr(:,k) + out_chem(:,i,k) + &
                          alp0*beta1*(-fm(k)*se_chem(:,i,max(kts,k-1)) +(fm(k+1)-fp(k))*se_chem(:,i,k) +fp(k+1)*se_chem(:,i,k+1))

            enddo
            do ispc = 1, mtp
                call tridiag (ktop(i),aa (kts:ktop(i)),bb (kts:ktop(i)),cc (kts:ktop(i)),ddtr (ispc,kts:ktop(i)))
                out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
            enddo
        ENDIF

        !- flux form + source/sink terms + time explicit + upstream with anti-diffusion step (Smolarkiewicz 1983)
        IF(USE_FLUX_FORM == 2 .or. USE_FLUX_FORM == 3) THEN
           if(USE_FLUX_FORM == 2)  lstep = -1 ! upstream + anti-diffusion step
           if(USE_FLUX_FORM == 3)  lstep =  1 ! only upstream
           alp0=1.

           if(ierr(i) /= 0) cycle
           !--- Zenv here have the following reference:  < 0 => downward motion
           zenv(i,:) = -(zuo(i,:)-edto(i)*zdo(i,:))

           do istep=1,lstep,-2

            if(istep == 1) then
               ddtr_upd(:,:) = se_chem(:,i,:)
               do k=kts,ktop(i)+1
                  fp_mtp(:,k) = 0.5*(zenv(i,k)+abs(zenv(i,k)))
                  fm_mtp(:,k) = 0.5*(zenv(i,k)-abs(zenv(i,k)))
               enddo
            else
               ddtr_upd(:,kts:ktop(i)+1) = se_chem(:,i,kts:ktop(i)+1) + out_chem(:,i,kts:ktop(i)+1)*dtime
               zenv_diff(:,kts) = 0.
               do k=kts,ktop(i)+1
                  dz =  zo_cup(i,k+1)-zo_cup(i,k)
                  zenv_diff (:,k+1) = 1.08*( dz*abs(zenv(i,k+1)) - dtime * zenv(i,k+1)**2 ) &
                                    * (ddtr_upd(:,k+1) - ddtr_upd(:,k))               &
                                    /((ddtr_upd(:,k+1) + ddtr_upd(:,k) + 1.e-16)*dz)
                enddo
               do k=kts,ktop(i)+1
                  fp_mtp(:,k) = 0.5*(zenv_diff(:,k)+abs(zenv_diff(:,k)))
                  fm_mtp(:,k) = 0.5*(zenv_diff(:,k)-abs(zenv_diff(:,k)))
               enddo
            endif

            do k=kts,ktop(i)
              dp=-100.*(po_cup(i,k)-po_cup(i,k+1))
              beta1 = dtime*g/dp
              ddtr(:,k) = ddtr_upd(:,k) + alp0*beta1*( &
                                           (fp_mtp(:,k+1)*ddtr_upd(:,k)           +fm_mtp(:,k+1)*ddtr_upd(:,k+1)) &
                                          -(fp_mtp(:,k  )*ddtr_upd(:,max(kts,k-1))+fm_mtp(:,k  )*ddtr_upd(:,k  )) )
            enddo
            do ispc = 1, mtp
                out_chem(ispc,i,kts:ktop(i))=(ddtr(ispc,kts:ktop(i))-se_chem(ispc,i,kts:ktop(i)))/dtime
            enddo

           enddo ! anti-diff steps

           do k=kts,ktop(i)
              dp=100.*(po_cup(i,k)-po_cup(i,k+1))
              beta1 = g/dp

              out_chem(:,i,k)  =   out_chem(:,i,k)                                                           &
                                  - (zuo(i,k+1)*sc_up_chem(:,i,k+1) - zuo(i,k)*sc_up_chem(:,i,k))*beta1      &
                                  + (zdo(i,k+1)*sc_dn_chem(:,i,k+1) - zdo(i,k)*sc_dn_chem(:,i,k))*beta1*edto(i)

              !- include evaporation
              if(USE_TRACER_EVAP == 1 .and. trim(cumulus) /= 'shallow') then
                 out_chem(:,i,k) = out_chem(:,i,k)    &
                                 - 0.5*edto(i)*(zdo(i,k)*pw_dn_chem(:,i,k)+zdo(i,k+1)*pw_dn_chem(:,i,k+1))*beta1 !&  ! evaporated ( pw_dn < 0 => E_dn > 0)
              endif

              !- include scavenging
              if(USE_TRACER_SCAVEN > 0 .and. trim(cumulus) /= 'shallow') then
                 out_chem(:,i,k) = out_chem(:,i,k) &
                                 - 0.5*(zuo(i,k)*pw_up_chem(:,i,k)+zuo(i,k+1)*pw_up_chem(:,i,k+1))*beta1  ! incorporated in rainfall (<0)
              endif

            enddo

        ENDIF

        !--- check mass conservation for tracers
        do ispc = 1, mtp
          trash_ (:) = 0.
          trash2_(:) = 0.
          evap_  (:) = 0.
          wetdep_(:) = 0.
          residu_(:) = 0.
          do k=kts,ktop(i)
             dp=100.*(po_cup(i,k)-po_cup(i,k+1))
             evap   =  -0.5*(zdo(i,k)*pw_dn_chem(ispc,i,k)+zdo(i,k+1)*pw_dn_chem(ispc,i,k+1))*g/dp*edto(i)
             wetdep =   0.5*(zuo(i,k)*pw_up_chem(ispc,i,k)+zuo(i,k+1)*pw_up_chem(ispc,i,k+1))*g/dp

             evap_  (ispc) =   evap_  (ispc) + evap  *dp/g
             wetdep_(ispc) =   wetdep_(ispc) + wetdep*dp/g
             residu_(ispc) =   residu_(ispc) + (wetdep - evap)*dp/g

!            trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,i,k) - evap + wetdep)*dp/g
             trash_ (ispc) =   trash_ (ispc) + (out_chem (ispc,i,k)                )*dp/g

             trash2_(ispc) =   trash2_(ispc) + se_chem(ispc,i,k)*dp/g
           enddo
           if(residu_(ispc) < 0.) then
             beta1 = g/(po_cup(i,kts)-po_cup(i,ktop(i)+1))
             do k=kts,ktop(i)
                out_chem(ispc,i,k)=out_chem(ispc,i,k)+residu_(ispc)*beta1
             enddo
           endif
        enddo


   enddo ! loop 'i'

!--------------------------------------------------------------------------------------------!
ENDIF !- end of section for atmospheric composition
!--------------------------------------------------------------------------------------------!

!
!- begin: for GATE soundings-------------------------------------------
   if(wrtgrads) then
    if(trim(cumulus) == 'deep'   ) then ; cty='1'; nvarbegin =  0 ;endif
    if(trim(cumulus) == 'shallow') then ; cty='2'; nvarbegin =101 ;endif
    if(trim(cumulus) == 'mid'    ) then ; cty='3'; nvarbegin =201 ;endif
    do i=its,itf
     !if(ierr(i).eq.0) then
      !- 2-d section
      do k=kts,ktf  !max(1,ktop(i))
       nvar=nvarbegin

       if(trim(cumulus) == 'deep'   ) &
         call set_grads_var(jl,k,nvar,zo(i,k),"zo"//cty ,' height','3d')
!        call set_grads_var(jl,k,nvar,po(i,k),"po"//cty ,' press','3d')

       dp=100.*(po_cup(i,k)-po_cup(i,k+1))
       E_dn  = -0.5*(pwdo(i,k)+pwdo(i,k+1))*g/dp*edto(i)*86400.*xlv/cp*xmb(i) ! pwdo < 0 and E_dn must > 0
       C_up  = dellaqc(i,k)+(zuo(i,k+1)* qrco(i,k+1) - zuo(i,k  )* qrco(i,k  )  )*g/dp &
                       +0.5*(pwo (i,k)+pwo (i,k+1))*g/dp
       C_up = - C_up*86400.*xlv/cp*xmb(i)

       trash =-(zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) -                 &
                zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) )*g/dp
       trash2=+(zdo(i,k+1)*(qcdo(i,k+1)-qo_cup(i,k+1) ) -                 &
                zdo(i,k  )*(qcdo(i,k  )-qo_cup(i,k  ) ) )*g/dp*edto(i)

       trash  = trash *86400.*xlv/cp*xmb(i); trash2 = trash2*86400.*xlv/cp*xmb(i)

       env_mf = 0.5* ((zuo(i,k+1)-zdo(i,k+1)*edto(i)) + (zuo(i,k)-zdo(i,k)*edto(i)))
       resten_H = dellah(i,k) - subten_H(i,k)
       resten_Q = dellaQ(i,k) - subten_Q(i,k)
       resten_T =(1./cp)*(resten_H-xlv*resten_Q)
      !trash2 = qco   (i,k  )! zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) ) !*g/dp
      !trash  = qo_cup(i,k  )! zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) ) !*g/dp
       trash2 = zuo(i,k+1)*(qco (i,k+1)-qo_cup(i,k+1) )*1000 !*g/dp
       trash  = zuo(i,k  )*(qco (i,k  )-qo_cup(i,k  ) )*1000  !*g/dp

       call set_grads_var(jl,k,nvar,out_chem(1,i,k)*86400,"outchem"//cty ,' outchem','3d')
       call set_grads_var(jl,k,nvar,sc_up_chem(1,i,k),"scup"//cty ,' sc_chem','3d')
       call set_grads_var(jl,k,nvar,sc_dn_chem(1,i,k),"scdn"//cty ,' sc_chem','3d')
       call set_grads_var(jl,k,nvar,massi,"mi"//cty ,' initial mass','2d')
       call set_grads_var(jl,k,nvar,massf,"mf"//cty ,' final mass','2d')
       call set_grads_var(jl,k,nvar,se_chem(1,i,k),"se"//cty ,' se_chem','3d')
       call set_grads_var(jl,k,nvar,se_cup_chem(1,i,k),"secup"//cty ,' se_cup_chem','3d')
       IF(APPLY_SUB_MP == 1) THEN
       kmp=lsmp
       call set_grads_var(jl,k,nvar,outmpqi(kmp,i,k)*86400*1000,"outqi"//cty ,' outmpqi','3d')
       call set_grads_var(jl,k,nvar,outmpql(kmp,i,k)*86400*1000,"outql"//cty ,' outmpql','3d')
       call set_grads_var(jl,k,nvar,outmpcf(kmp,i,k)*86400,"outcf"//cty ,' outmpcf','3d')

       call set_grads_var(jl,k,nvar,mpqi(kmp,i,k),"mpqi"//cty ,' mpqi','3d')
       call set_grads_var(jl,k,nvar,mpql(kmp,i,k),"mpql"//cty ,' mpql','3d')
       call set_grads_var(jl,k,nvar,mpcf(kmp,i,k),"mpcf"//cty ,' mpcf','3d')
       ENDIF
       call set_grads_var(jl,k,nvar,env_mf,"sub"//cty ,' sub','3d')


       IF(LIQ_ICE_NUMBER_CONC == 1) THEN
         call set_grads_var(jl,k,nvar,outnice(i,k)*86400.,"outnice"//cty ,'out # ice1/day','3d')
         call set_grads_var(jl,k,nvar,outnliq(i,k)*86400.,"outnliq"//cty ,'out # liq /day','3d')
       ENDIF
       call set_grads_var(jl,k,nvar,zuo(i,k)/xmb(i),"zup"//cty,'norm m flux up ','3d')
       call set_grads_var(jl,k,nvar,zdo(i,k)/xmb(i),"zdn"//cty,'norm m flux dn ','3d')
       call set_grads_var(jl,k,nvar,zenv(i,k),"zenv"//cty,'norm m flux env ','3d')

       call set_grads_var(jl,k,nvar,-edto(i)*xmb(i)*zdo(i,k),"mdn"//cty ,'m flux down (kg/s/m^2)','3d')
       call set_grads_var(jl,k,nvar,up_massentro(i,k),"upent"//cty ,'up_massentr(kg/s/m^2)','3d')
       call set_grads_var(jl,k,nvar,xmb(i)*up_massdetro(i,k),"updet"//cty ,'up_massdetr(kg/s/m^2)','3d')
       call set_grads_var(jl,k,nvar,outt(i,k)*86400.,"outt"//cty ,'outt K/day','3d')

       call set_grads_var(jl,k,nvar,resten_T*86400.,           "rest"//cty ,'residuo T K/day','3d')
       call set_grads_var(jl,k,nvar,resten_H*86400./cp,    "resh"//cty ,'residuo H J/kg/day','3d')
       call set_grads_var(jl,k,nvar,resten_Q*86400.*xlv/cp,"resq"//cty ,'residuo q K/day   ','3d')
       call set_grads_var(jl,k,nvar,subten_T(i,k)*86400.,       "subt"//cty ,'subT K/day','3d')
       call set_grads_var(jl,k,nvar,subten_H(i,k)*86400./cp,    "subh"//cty ,'subH J/kg/day','3d')
       call set_grads_var(jl,k,nvar,subten_Q(i,k)*86400.*xlv/cp,"subq"//cty ,'subq K/day   ','3d')

       call set_grads_var(jl,k,nvar,outq(i,k)*86400.*xlv/cp,"outq"//cty ,'outq K/s','3d')
       call set_grads_var(jl,k,nvar,outqc(i,k)*86400.*xlv/cp,"outqc"//cty ,'outqc K/day','3d')
       call set_grads_var(jl,k,nvar,pre(i)*3600.,"precip"//cty ,'precip mm','2d')
       call set_grads_var(jl,k,nvar,prec_flx(i,k)*3600.,"precflx"//cty ,'prec flx mm','3d')
       call set_grads_var(jl,k,nvar,pwo(i,k),"pwo"//cty ,' xx','3d')
       call set_grads_var(jl,k,nvar,outu(i,k)*86400.,"outu"//cty ,'out_U m/s/day','3d')
       call set_grads_var(jl,k,nvar,outv(i,k)*86400.,"outv"//cty ,'out_V m/s/day','3d')
       call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty ,'xmb kg/m2/s','2d')
       call set_grads_var(jl,k,nvar,vvel2d(i,k),"W2d"//cty ,'W /m/s','3d')
       call set_grads_var(jl,k,nvar,vvel1d(i),"W1d"//cty ,'W1s /m/s','2d')
       call set_grads_var(jl,k,nvar,us(i,k),"us"//cty ,'U /m/s','3d')
       call set_grads_var(jl,k,nvar,outu(i,k)*86400./(1.e-16+xmb(i)),"delu"//cty ,'dellu','3d')
       call set_grads_var(jl,k,nvar,evap_bcb(i,k)*1000.,"evcb"//cty ,'g/kg','3d')

       call set_grads_var(jl,k,nvar,tot_pw_up_chem(1,i),"pwup"//cty ,'pwup','2d')
       call set_grads_var(jl,k,nvar,tot_pw_dn_chem(1,i),"pwdn"//cty ,'pwdn','2d')
!----
!----
       call set_grads_var(jl,k,nvar,xmb(i)*dellah(i,k)*86400./cp,"delh"//cty ,'dellah K/day','3d')
       call set_grads_var(jl,k,nvar,xmb(i)*dellaq(i,k)*86400.*xlv/cp, "dellq"//cty ,'dellaq K/day','3d')
       call set_grads_var(jl,k,nvar,xmb(i)*dellaqc(i,k)*86400.*xlv/cp,"dellqc"//cty ,'dellaqc K/day','3d')
       call set_grads_var(jl,k,nvar,xmb(i),"xmb"//cty,'m flux up (kg/s/m^2)','2d')
       call set_grads_var(jl,k,nvar,aa1(i),"aa1"//cty,'AA1 J/kg3)','2d')
       call set_grads_var(jl,k,nvar,float(ierr(i)),"ierr"//cty ,'ierr #','2d')
       call set_grads_var(jl,k,nvar,xmb(i)*dd_massentro(i,k),"ddent"//cty ,'dd_massentr(kg/s/m^2)','3d')
       call set_grads_var(jl,k,nvar,xmb(i)*dd_massdetro(i,k),"dddet"//cty ,'dd_massdetr(kg/s/m^2)','3d')
!!      go to 333
       call set_grads_var(jl,k,nvar,hc(i,k),"hc"//cty ,' hc','3d')
       call set_grads_var(jl,k,nvar,hco(i,k),"hco"//cty ,' hco','3d')
       call set_grads_var(jl,k,nvar,dby(i,k),"dby"//cty ,' dbuo','3d')
       !call set_grads_var(jl,k,nvar,QCUP(i,k),"qcup"//cty ,'C_UP','3d')
       call set_grads_var(jl,k,nvar,t_cup(i,k)-273.15,"te"//cty ,' K','3d')
       call set_grads_var(jl,k,nvar,1000.*q_cup(i,k),"qe"//cty ,' kg kg-1','3d')
       call set_grads_var(jl,k,nvar,he_cup(i,k),"he"//cty ,' he','3d')
       call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
       call set_grads_var(jl,k,nvar,HKB(i),"hkb"//cty ,' H','2d')
       call set_grads_var(jl,k,nvar,z_cup(i,max(1,k22  (i))),"zs"//cty ,' m','2d')
       call set_grads_var(jl,k,nvar,z_cup(i,max(1,kbcon(i))),"zbcon"//cty ,' m','2d')
       call set_grads_var(jl,k,nvar,z_cup(i,max(1,ktop (i))),"ztop"//cty ,' m','2d')
       call set_grads_var(jl,k,nvar,z_cup(i,max(1,klcl (i))),"zlcl"//cty ,' m','2d')
       call set_grads_var(jl,k,nvar,z_cup(i,max(1,jmin (i))),"zjmin"//cty ,' m','2d')
       call set_grads_var(jl,k,nvar,zws(i),"ws"//cty ,' m/s','2d')
       call set_grads_var(jl,k,nvar,clfrac(i,k),"clfrac"//cty ,'shcf #','3d')
       call set_grads_var(jl,k,nvar,entr_rate(i,k),"entr"//cty ,' m-1','3d')
       call set_grads_var(jl,k,nvar,cd(i,k),"detr"//cty ,' m-1','3d')
       call set_grads_var(jl,k,nvar,pwdo(i,k),"pwd"//cty ,' xx','3d')
       call set_grads_var(jl,k,nvar,edto(i),"edt"//cty ,'edt kg/m2/s','2d')
       call set_grads_var(jl,k,nvar,E_DN,"EVAP"//cty ,' xx','3d')
       call set_grads_var(jl,k,nvar,C_UP,"CUP"//cty ,' xx','3d')
!       call set_grads_var(jl,k,nvar,trash,"TUP"//cty ,' xx','3d')
!       call set_grads_var(jl,k,nvar,trash2,"TDN"//cty ,' xx','3d')
       call set_grads_var(jl,k,nvar,trash,"F1"//cty ,' F1','3d')
       call set_grads_var(jl,k,nvar,trash2,"F2"//cty ,' F2','3d')

       call set_grads_var(jl,k,nvar,p_liq_ice(i,k),"pli"//cty ,'#','3d')
       call set_grads_var(jl,k,nvar,melting_layer(i,k),"cpli"//cty ,'#','3d')
       call set_grads_var(jl,k,nvar,t(i,k),"t"//cty ,'temp K','3d')
       call set_grads_var(jl,k,nvar,tn(i,k),"tn"//cty ,'temp K','3d')
       call set_grads_var(jl,k,nvar,1000.*q(i,k),"q"//cty ,'q g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*qo(i,k),"qn"//cty ,'q g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*qrco(i,k),"qrc"//cty ,'q g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*(q(i,k)+outq(i,k)*dtime),"qnc"//cty ,'q upd conv g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*(qo(i,k)+outq(i,k)*dtime),"qnall"//cty ,'q upd all g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*qrr(i,k),"qrr"//cty ,'qrr g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*qco(i,k),"qc"//cty ,'qc g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*qo_cup(i,k),"qcup"//cty ,'qc g/kg','3d')
       call set_grads_var(jl,k,nvar,1000.*qeso_cup(i,k),"qescup"//cty ,'qc g/kg','3d')

       !~ call set_grads_var(jl,k,nvar,aa0(i),"a0"//cty,'aa0','2d')
       !~ call set_grads_var(jl,k,nvar,aa1_fa(i),"aa1fa"//cty,'aa1fa','2d')
       !~ call set_grads_var(jl,k,nvar,aa1_bl(i),"aa1bl"//cty,'aa1bl','2d')
       !~ call set_grads_var(jl,k,nvar,aa0_bl(i),"aa0bl"//cty,'aa0bl','2d')
       !~ call set_grads_var(jl,k,nvar,aa1(i),"a1"//cty,'aa1','2d')
       !~ call set_grads_var(jl,k,nvar,aa1(i)/(1.e-6+tau_ecmwf(i)),"mb13"//cty,'aa0','2d')
       !~ call set_grads_var(jl,k,nvar,xaa0(i),"xa0"//cty,'xaa0','2d')
       !~ call set_grads_var(jl,k,nvar,(XAA0(I)-AA1(I))/MBDT(I),"xk"//cty,'xk','2d')
  333 continue
      enddo
      if(wrtgrads) then
          call wrt_bin_ctl(1,kte,po(1,1:kte),cumulus)
      endif
    enddo
   endif
!- end  : for GATE soundings-------------------------------------------
!
!
!-------------------------- not in use ------------------------------------------------------!
!--- get cloud fraction
!
! do i=its,itf
!    clfrac(i,:)=0.
!    if(ierr(i) /= 0) cycle
!    dummy1(kts:ktf) = xmb(i)* zuo(i,kts:ktf)
!    dummy2(kts:ktf) = 100.*po_cup(i,kts:ktf)
!    call get_cloud_fraction(ktf,kts,ktf                                                   &
!     ,dummy2(kts:ktf),zo_cup(i,kts:ktf),tn_cup(i,kts:ktf),qo_cup(i,kts:ktf) &
!     ,qco (i,kts:ktf),  qrco(i,kts:ktf),  dummy1(kts:ktf),clfrac(i,kts:ktf) )
! enddo
!--------------------------------------------------------------------------------------------!
!

   END SUBROUTINE CUP_gf
!------------------------------------------------------------------------------------
   SUBROUTINE cup_dd_edt(cumulus,ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pw,ccn,pwev,edtmax,edtmin,maxens2,edtc,psum2,psumh,    &
              rho,aeroevap,itf,ktf,ipr,jpr,its,ite, kts,kte )

     IMPLICIT NONE

     character *(*),        intent (in )  :: cumulus
     integer  ,intent (in  )             ::                           &
        ipr,jpr,aeroevap,itf,ktf,its,ite, kts,kte
     integer, intent (in   )              :: maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        rho,us,vs,z,p,pw
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev,ccn,psum2,psumh,edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
!
!  local variables in this routine
!
     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::   vshear,sdp,vws
     real :: pefc,aeroadd,rhoc,dp,prop_c
     real, parameter ::  alpha3 = 1.9 ,beta3  = -1.13

!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       edt   =0.
       vws   =0.
       sdp   =0.
       vshear=0.
       edtc  =0.

       if(trim(cumulus)=='shallow') return

       do i=its,itf
         if(ierr(i) /= 0)cycle
         do kk = kbcon(i),ktop(i)
             dp = p(i,kk) - p(i,kk+1)
             vws(i) = vws(i) + (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk)))    &
                             +  abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * dp
             sdp(i) = sdp(i) + dp
         enddo
         vshear(i) = 1.e3 * vws(i) / sdp(i)
       end do

       do i=its,itf
         if(ierr(i) /= 0)cycle
            pef = (1.591-0.639*vshear(i)+0.0953*(vshear(i)**2) -0.00496*(vshear(i)**3))
            pef = min(pef,0.9)
            pef = max(pef,0.1)

!
!--- cloud base precip efficiency
!
            zkbc  = z(i,kbcon(i))*3.281e-3
            prezk = 0.02
            if(zkbc > 3.0) prezk=0.96729352+zkbc*(-0.70034167+zkbc* &
                                (0.162179896+zkbc*(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            if(zkbc > 25.) prezk=2.4

            pefb = 1./(1.+prezk)
            pefb = min(pefb,0.9)
            pefb = max(pefb,0.1)

            edt(i) = 1.-0.5*(pefb+pef)

            if(aeroevap.gt.1)then
               aeroadd=(ccnclean**beta3)*((psumh(i))**(alpha3-1)) !*1.e6
              !if(i.eq.ipr)write(0,*)'edt',ccnclean,psumh(i),aeroadd
              !prop_c=.9/aeroadd
               prop_c=.5*(pefb+pef)/aeroadd
               aeroadd=(ccn(i)**beta3)*((psum2(i))**(alpha3-1)) !*1.e6
              !if(i.eq.ipr)write(0,*)'edt',ccn(i),psum2(i),aeroadd,prop_c
               aeroadd=prop_c*aeroadd
               pefc=aeroadd
               if(pefc.gt.0.9)pefc=0.9
               if(pefc.lt.0.1)pefc=0.1
               edt(I)=1.-pefc
               if(aeroevap.eq.2)edt(I)=1.-.25*(pefb+pef+2.*pefc)
            endif

!--- edt here is 1-precipeff!
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
       enddo
       do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               edtc(I,K)=-edtc(I,K)*pwav(I)/pwev(I)
               IF(edtc(I,K).GT.edtmax(i))edtc(I,K)=edtmax(i)
               IF(edtc(I,K).LT.edtmin(i))edtc(I,K)=edtmin(i)
            enddo
         endif
       enddo

   END SUBROUTINE cup_dd_edt
!------------------------------------------------------------------------------------
   SUBROUTINE cup_dd_moisture(cumulus,ierrc,zd,hcd,hes_cup,qcd,qes_cup,                &
              pwd,q_cup,z_cup,dd_massentr,dd_massdetr,jmin,ierr,                       &
              gamma_cup,pwev,bu,qrcd, q,he,t_cup,iloop,t_wetbulb,q_wetbulb,qco,pwavo,  &
              itf,ktf,its,ite, kts,kte                              )

     IMPLICIT NONE

     character(len=*), intent(in) :: cumulus
     integer         , intent(in) :: itf,ktf,its,ite, kts,kte
  ! q       = environmental q on model levels
  ! q_cup   = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd  = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd  = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr = entrainment rate
  ! cdd  = detrainment function
  !
     real,    dimension (its:ite) ,intent (in  )          ::           &
        t_wetbulb,q_wetbulb, pwavo
     real,    dimension (its:ite,kts:kte) ,intent (in   ) ::           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,                      &
        dd_massentr,dd_massdetr,gamma_cup,q,he,qco
     integer  ,intent (in   )                             ::           &
        iloop
     integer, dimension (its:ite) ,intent (in   )         ::           &
        jmin
     integer, dimension (its:ite) ,intent (inout)         ::           &
        ierr
     real,    dimension (its:ite,kts:kte) ,intent (out  ) ::           &
        qcd,qrcd,pwd
     real,    dimension (its:ite) ,intent (out  )         ::           &
        pwev,bu
     character*128 :: ierrc(its:ite)
 !
 !  local variables in this routine
 !
     integer   ::     i,k
     real      ::     dh,dz,dq_eva,denom,fix_evap
 !
     bu  =0.  !-- buoyancy
     qcd =0.  !-- in-downdradt water vapor mixing ratio
     qrcd=0.  !-- saturation water vapor mixing ratio
     pwev=0.  !-- column integrated rain evaporation (normalized)
     pwd =0.  !-- rain evaporation at layer k

     if(trim(cumulus) == 'shallow') return
!
     DO i=its,itf
        if(ierr(i) /= 0) cycle

        !-- boundary condition in jmin ('level of free sinking')
        k=jmin(i)
        dz=z_cup(i,k+1)-z_cup(i,k)

        qcd(i,k)=q_cup(i,k)

        if(use_wetbulb==1) then
          !--option 1
          !qcd(i,k)=q_wetbulb(i)
          !--option 2
          qcd(i,k)=0.5*(q_wetbulb(i)+qco(i,k)) ! mixture 50% env air + updraft
        endif

        dh=hcd(i,k)-hes_cup(i,k)

        if(dh.lt.0)then
          qrcd(i,k)=(qes_cup(i,k)+(1./xlv)*(gamma_cup(i,k) /(1.+gamma_cup(i,k)))*dh)
        else
          qrcd(i,k)=qes_cup(i,k)
        endif

        pwd (i,k) = zd(i,k)*min(0.,qcd(i,k)-qrcd(i,k))
        qcd (i,k) = qrcd(i,k)
        pwev(i)   = pwev(i)+pwd(i,k)
        bu  (i)   = dz*dh

        do k=jmin(i)-1,kts,-1

           dz=z_cup(i,k+1)-z_cup(i,k)

          !-- downward transport + mixing
           denom = (zd(i,k+1)-0.5*dd_massdetr(i,k)+dd_massentr(i,k))
           if( denom == 0.0 )then
             qcd(i,k)= qcd(i,k+1)
           else
             qcd(i,k)=(qcd(i,k+1)*zd(i,k+1) -0.5*dd_massdetr(i,k)*qcd(i,k+1)+ &
                                                 dd_massentr(i,k)*q  (i,k)    )/ denom
           endif
          !
          !--- to be negatively buoyant, hcd should be smaller than hes!
          !--- ideally, dh should be negative till dd hits ground, but that is not always
          !--- the case
          !
           dh    = hcd(i,k)-hes_cup(i,k)
           bu(i) = bu(i)+dz*dh
           qrcd(i,k)=qes_cup(i,k)+(1./xlv)*(gamma_cup(i,k) /(1.+gamma_cup(i,k)))*dh

          !-- rain water evaporation amount at layer k
           dq_eva=qcd(i,k)-qrcd(i,k)

           if(dq_eva.gt.0.)then
            dq_eva=0.
            qrcd(i,k)=qcd(i,k)
           endif
          !-- amount of the evaporated rain water
           pwd(i,k)=zd(i,k)*dq_eva  ! kg[water vapor]/kg[air]

          !-- source term for in-downdraft water vapor mixing ratio
           qcd(i,k)=qrcd(i,k)     ! => equiv to qcd = qcd - dq_eva !( -dq_eva >0 => source term for qcd)

          !-- total evaporated rain water
           pwev(i)=pwev(i)+pwd(i,k)

          !-- for GEOS diagnostic
          ! evap(i,k) = - edt * xmb * zd * dq_eva = - edt * xmb * pwd (i,k)
          ! downdfrat temp = (hcd(i,k)-qcd(i,k)*xlv-g*z_cup(i,k))/cp - 273.15

        enddo

        if(pwev(i).ge.0.and.iloop.eq.1)then
         ierr(i)=70
         ierrc(i)="problem with buoy in cup_dd_moisture"
        endif
        if(bu(i).ge.0.and.iloop.eq.1)then
         ierr(i)=73
         ierrc(i)="problem2 with buoy in cup_dd_moisture"
        endif

        if(ZERO_DIFF==0 .and. EVAP_FIX==1) then
          if(abs(pwev(i)) > pwavo(i) .and. ierr(i) == 0)then
             fix_evap = pwavo(i)/(1.e-16+abs(pwev(i)))
             pwev(i)  = 0.

             do k=jmin(i),kts,-1
                pwd(i,k) = pwd (i,k)*fix_evap
                pwev(i)  = pwev(i) + pwd(i,k)
                dq_eva   = pwd (i,k)/(1.e-16+zd(i,k))
                qcd(i,k) = qrcd(i,k) + dq_eva
             enddo
             if(pwev(i) .ge. 0.)then
               ierr(i)=70; ierrc(i)="problem with buoy in cup_dd_moisture"
             endif
           endif
         endif
     ENDDO!--- end loop over i

   END SUBROUTINE cup_dd_moisture

!------------------------------------------------------------------------------------

   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,psur,ierr,itest,                     &
                      itf,ktf,its,ite, kts,kte             )

     implicit none

     integer ,intent (in)                ::                         &
        itf,ktf,its,ite, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t,q
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

       integer                              ::                           &
       i,k,iph
!      real, dimension (1:2) :: AE,BE,HT
       real, dimension (its:ite,kts:kte) :: tv
       real :: e,tvbar
!      real, external :: satvap
!      real :: satvap

       he  =0.0
       hes =0.0
       qes =0.0

      if(SATUR_CALC == 0) then
       do k=kts,ktf
        do i=its,itf
        if(ierr(i).eq.0)then

        e=satvap(t(i,k))
        qes(i,k)=0.622*e/max(1.e-8,(p(i,k)-e))

        IF(QES(I,K).LE.1.E-08  ) QES(I,K)=1.E-08
        IF(QES(I,K).GT.max_qsat) QES(I,K)=max_qsat
        IF(QES(I,K).LT.Q(I,K)  ) QES(I,K)=Q(I,K)
!       IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
        enddo
       enddo

      else

       !--- better formulation for the mixed phase regime
       do k=kts,ktf
        do i=its,itf
         if(ierr(i).eq.0)then
           qes(i,k) = satur_spec_hum(t(i,k),p(i,k))
           qes(i,k) = min(max_qsat, max(1.E-08,qes(i,k)))
           qes(i,k) = max(qes(i,k), q(i,k))
           tv (i,k) = t(i,k)+.608*q(i,k)*t(i,k)
         endif
        enddo
      enddo
     endif

!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.eq.1 .or. itest.eq.0)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/g
           endif
         enddo

! --- calculate heights
         do K=kts+1,ktf
          do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/g
           endif
          enddo
         enddo
      else if(itest.eq.2)then
         do k=kts,ktf
          do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/g
             z(i,k)=max(1.e-3,z(i,k))
           endif
          enddo
         enddo
      else if(itest.eq.-1)then
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       do k=kts,ktf
        do i=its,itf
            if(ierr(i) /= 0) cycle
            if(itest.le.0) he (i,k)=g*z(i,k)+cp*t(i,k)+xlv*q(i,k)

            hes(i,k)=g*z(i,k)+cp*t(i,k)+xlv*qes(i,k)

            if(he(i,k).ge.hes(i,k)) he(i,k)=hes(i,k)
       enddo
      enddo

   END SUBROUTINE cup_env
!------------------------------------------------------------------------------------
   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
              us, vs,u_cup,v_cup,                                   &
              hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur,tsur,        &
              ierr,z1,itf,ktf,its,ite, kts,kte                      )
     implicit none
     integer ,intent (in   )              :: itf,ktf, its,ite, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  !
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t,us, vs
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,u_cup,v_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1,tsur
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!
     integer                              :: i,k
     real                                 :: p1,p2,ct1,ct2,rho
     integer, save                        :: irun = 0

      qes_cup  =0.
      q_cup    =0.
      hes_cup  =0.
      he_cup   =0.
      z_cup    =0.
      p_cup    =0.
      t_cup    =0.
      gamma_cup=0.
      u_cup    =0.
      v_cup    =0.

      IF( clev_grid == 2 ) THEN
      !--original formulation
       do k=kts+1,ktf
        do i=its,itf
        if(ierr(i) /= 0) cycle
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup  (i,k)=.5*(q  (i,k-1)+q  (i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup (i,k)=.5*(he (i,k-1)+he (i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup  (i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup  (i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup  (i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xlv/cp)*(xlv/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        u_cup  (i,k)=.5*(us(i,k-1)+us(i,k))
        v_cup  (i,k)=.5*(vs(i,k-1)+vs(i,k))

        enddo
       enddo
       do i=its,itf
        if(ierr(i) /= 0) cycle
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
       !hes_cup(i,1)=hes(i,1)
       !he_cup(i,1)=he(i,1)
        hes_cup(i,1)=g*z1(i)+cp*t(i,1)+xlv*qes(i,1)
        he_cup (i,1)=g*z1(i)+cp*t(i,1)+xlv*q  (i,1)
       !z_cup(i,1)=.5*(z(i,1)+z1(i))
       !p_cup(i,1)=.5*(p(i,1)+psur(i))
        z_cup(i,1)=z1(i)
        p_cup(i,1)=psur(i)
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        u_cup(i,1)=us(i,1)
        v_cup(i,1)=vs(i,1)
       enddo

      ELSEIF( clev_grid == 0) THEN
      !--- weigthed mean
         do i=its,itf
            if(ierr(i) /= 0) cycle
            p_cup  (i,1)=psur(i)
            z_cup  (i,1)=z1(i)
            do k=kts,ktf-1
                     p_cup (i,k+1) = 2.0*p(i,k) - p_cup(i,k)
                     z_cup (i,k+1) = 2.0*z(i,k) - z_cup(i,k)
            enddo

            ! ----------- p,T          k+1
            !p1
            ! ----------- p_cup,T_cup  k+1
            !p2
            ! ----------- p,T          k
            !
            ! ----------- p_cup,T_cup  k

            do k=kts,ktf-1
               p1=abs((p    (i,k+1) - p_cup(i,k+1))/(p(i,k+1)-p(i,k)))
               p2=abs((p_cup(i,k+1) - p    (i,k  ))/(p(i,k+1)-p(i,k)))

               t_cup  (i,k+1) = p1*t  (i,k) + p2*t  (i,k+1)

               u_cup  (i,k+1) = p1*us (i,k) + p2*us (i,k+1)
               v_cup  (i,k+1) = p1*vs (i,k) + p2*vs (i,k+1)
               q_cup  (i,k+1) = p1*q  (i,k) + p2*q  (i,k+1)
               he_cup (i,k+1) = p1*he (i,k) + p2*he (i,k+1)

               qes_cup(i,k+1) = p1*qes(i,k) + p2*qes(i,k+1)
               hes_cup(i,k+1) = p1*hes(i,k) + p2*hes(i,k+1)

               if(he_cup(i,k+1).gt.hes_cup(i,k+1))he_cup(i,k+1)=hes_cup(i,k+1)

               gamma_cup(i,k+1)=(xlv/cp)*(xlv/(rv*t_cup(i,k+1) &
                                *t_cup(i,k+1)))*qes_cup(i,k+1)

            enddo
            !--- surface level from X(kts) and X_cup(kts+1) determine X_cup(kts)
            k=kts
            p1=abs(p    (i,k  )-p_cup(i,k))
            p2=abs(p_cup(i,k+1)-p_cup(i,k))

            ct1=(p1+p2)/p2
            ct2=p1/p2

            t_cup  (i,k) = ct1*t  (i,k) - ct2*t_cup(i,k+1)
            q_cup  (i,k) = ct1*q  (i,k) - ct2*q_cup(i,k+1)

            u_cup  (i,k) = ct1*us (i,k) - ct2*u_cup(i,k+1)
            v_cup  (i,k) = ct1*vs (i,k) - ct2*v_cup(i,k+1)
            qes_cup(i,k) = ct1*qes(i,k) - ct2*qes_cup(i,k+1)

            hes_cup(i,k)=g*z_cup(i,k)+cp*t_cup(i,k)+xlv*qes_cup(i,k)
            he_cup (i,k)=g*z_cup(i,k)+cp*t_cup(i,k)+xlv*q_cup  (i,k)

            if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)

            gamma_cup(i,k)=xlv/cp*(xlv/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
         enddo

      ELSEIF(clev_grid == 1) THEN
      !--- based on Tiedke (1989)
          do i=its,itf
            if(ierr(i) /= 0) cycle
            do k=ktf, kts+1,-1

              qes_cup(i,k) = qes(i,k)
              q_cup  (i,k) = q  (i,k)
              p_cup  (i,k) = 0.5*(p(i,k-1)+p(i,k))
              z_cup  (i,k) = 0.5*(z(i,k-1)+z(i,k))
              t_cup  (i,k) = (max(cp*t(i,k-1)+g*z(i,k-1),cp*t(i,k)+g*z(i,k)) - g*z_cup(i,k))/cp

              if(qes(i,k) < max_qsat) &
              call get_interp(qes_cup(i,k),t_cup(i,k),p_cup(i,k),qes_cup(i,k),t_cup(i,k))

              q_cup  (i,k) = min(q(i,k),qes(i,k)) + qes_cup(i,k) - qes(i,k)
              q_cup  (i,k) = max(q_cup  (i,k) ,0.0)

            enddo
            !---level kts
            qes_cup(i,1)= qes (i,1)
            q_cup  (i,1)= q   (i,1)
            z_cup  (i,1)= z1  (i)
            p_cup  (i,1)= psur(i)

            t_cup  (i,1)= (cp*t(i,1)+g*z(i,1) - g*z_cup(i,1))/cp

            hes_cup(i,1)=g*z_cup(i,1)+cp*t_cup(i,1)+xlv*qes_cup(i,1)
            he_cup (i,1)=g*z_cup(i,1)+cp*t_cup(i,1)+xlv*q_cup  (i,1)

            gamma_cup(i,1)=xlv/cp*(xlv/(rv*t_cup(i,1)*t_cup(i,1)))*qes_cup(i,1)
            u_cup(i,1)=us(i,1)
            v_cup(i,1)=vs(i,1)

            do k=ktf,kts+1,-1
                p1=max(cp*t_cup(i,k)+g*z_cup(i,k), cp*t_cup(i,k-1)+g*z_cup(i,k-1))
                t_cup(i,k) = (p1-g*z_cup(i,k))/cp

                hes_cup(i,k)=cp*t_cup(i,k)+xlv*qes_cup(i,k)+g*z_cup  (i,k)
                he_cup (i,k)=cp*t_cup(i,k)+xlv*q_cup  (i,k)+g*z_cup  (i,k)
                if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)

                gamma_cup(i,k)=(xlv/cp)*(xlv/(rv*t_cup(i,k)*t_cup(i,k)))*qes_cup(i,k)
                u_cup    (i,k)=us(i,k)
                v_cup    (i,k)=vs(i,k)
            enddo
          enddo
      ELSE
            STOP "cup_env_clev"
      ENDIF

      RETURN
     !IF( MAPL_AM_I_ROOT() .and. irun == 0) then
       irun = 1
       DO i=its,itf
           IF(ierr(i) == 0) then
            do k=kts,kte-1
              rho=100*(p_cup(i,k)-p_cup(i,k+1))/(z_cup(i,k+1)-z_cup(i,k))/g ! air dens by hidrostatic balance (kg/m3)
              write(23,101) i,k,z_cup(i,k),p_cup(i,k),t_cup(i,k),q_cup(i,k)*1000.,he_cup(i,k),u_cup(i,k),v_cup(i,k),rho

              rho=100*(p(i,k)-p(i,k+1))/(z(i,k+1)-z(i,k))/g
              write(25,101) i,k,z    (i,k),p    (i,k),t    (i,k),q    (i,k)*1000.,he        (i,k),us   (i,k),vs   (i,k),rho

              101 format(2i3,8F15.5)
            ENDDO
            GOTO 400
           ENDIF
       ENDDO
     !ENDIF
     400 continue

   END SUBROUTINE cup_env_clev
!------------------------------------------------------------------------------------

   SUBROUTINE cup_forcing_ens_3d_mid(aa0,aa1,xaa0,mbdt,dtime,ierr,&
                                     po_cup,ktop,k22,kbcon,kpbl,ichoice, maxens,maxens3, &
                                     itf,ktf,its,ite, kts,kte, &
                                     tau_ecmwf,aa1_bl,xf_dicycle, &
                                     dhdt,xff_mid,zws,hc,hco,he_cup,heo_cup)

     IMPLICIT NONE
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! kbcon   = LFC of parcel from k22
  ! k22     = updraft originating level
  ! ichoice = flag if only want one closure
  ! name    = deep,mid or shallow convection flag
  !
     integer  ,intent (in   )    :: itf,ktf,its,ite, kts,kte,maxens,maxens3
     integer, dimension (its:ite)          ,intent (in) ::      &
        k22,kbcon,ktop,kpbl
     real,    dimension (its:ite,kts:kte)  ,intent (in) ::      &
        po_cup,dhdt,hc,hco,he_cup,heo_cup
     real,    dimension (its:ite)          ,intent (in) ::      &
        xaa0
     real,    dimension (its:ite)          ,intent (in) ::      &
        aa1,zws,mbdt,  aa0
     real                                  ,intent (in) :: dtime
     integer                               ,intent (in) :: ichoice
     real,    dimension (its:ite)          ,intent (IN) :: aa1_bl,tau_ecmwf
     integer, dimension (its:ite)          ,intent (inout) :: ierr
     real,    dimension (its:ite)          ,intent (inout) :: xf_dicycle
     real,    dimension (its:ite,1:maxens3),intent (out)   :: xff_mid
!
!  local variables in this routine
!
     real,    dimension (1:maxens) :: xk
     integer                       :: i,k
     real                          :: xff_dicycle, trash, blqe,xff_ens1,mf_ens1
       DO i=its,itf
          !-initialization
          xff_mid   (i,:)= 0.
          xf_dicycle(i)  = 0.

          if(ierr(i) /= 0)cycle

          !- Think about this:
          !xff0= (AA1(I)-AA0(I))/DTIME
          !if(xff0.lt.0.) xff_dicycle = 0.

          XK(1)=(XAA0(I)-(AA1(I)))/MBDT(i)

          if(xk(1).le.0.and.xk(1).gt.-0.1*mbdt(i)) xk(1)=-0.1*mbdt(i)
          if(xk(1).gt.0.and.xk(1).lt.1.e-2       ) xk(1)=1.e-2

          !- diurnal cycle mass flux
          if(DICYCLE==1 .or. DICYCLE==6 .or. DICYCLE==0 ) then
              !----  Betchold et al (2014)
              xff_dicycle = AA1_BL(i)/tau_ecmwf(i)
              xf_dicycle(i) = max(0.,-xff_dicycle /xk(1))
          endif
          if(DICYCLE==5 ) then
              xff_ens1 = max(0.,(AA1(I)-AA0(I))/dtime)
              mf_ens1  = 0.0
              if(XK(1).lt.0.) mf_ens1=max(0.,-xff_ens1/xk(1))
              xf_dicycle(i)=-( max(0., -(AA1_BL(I)-AA0(I))/dtime/xk(1)) - mf_ens1 )
             !xf_dicycle(i)=-( max(0., -(AA1_BL(I)-AA0(I))/dtime/xk(1)) - max(0.,-max(0.,(AA1(I)-AA0(I))/dtime)/xk(1) ))
          endif

          !- closures 3 and 4 for mid
          if(xk(1) < 0.) then
             xff_mid(i,3)=max(0., -(AA1(i)/tau_ecmwf(i))/xk(1))
             xff_mid(i,4)=xf_dicycle(i)
          endif

          !IF( ichoice == 3) xf_dicycle(i) =0. ! No dicycle for congestus
       ENDDO

       do i=its,itf
         if(ierr(i) /= 0) cycle
         !- Boundary layer quasi-equilibrium (Raymond 1995)
         if(k22(i).lt.kpbl(i)+1)then
           blqe=0.
           do k=kts,kbcon(i) !- orig formulation
           !do k=kts,kpbl(i)
             blqe = blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
           enddo
          !trash = max((hc (i,kbcon(i))-he_cup (i,kbcon(i))),1.e1)!- orig formulation
           trash = max((hco(i,kbcon(i))-heo_cup(i,kbcon(i))),1.e1)
           xff_mid(i,2) = max(0.,blqe/trash)
         endif

         !- W* closure (Grant,2001)
         xff_mid(i,1)=0.03*zws(i)
       enddo

       !-set xf_dicycle(i)=0 in case the closure is 4 and DICYCLE closure will be applied
       IF((DICYCLE>0) .and. ichoice == 4) xf_dicycle(:)=0.

       IF( ichoice == 5) then
         if(DICYCLE>0) then
          do i=its,itf
                if(ierr(i) /= 0) cycle
                !xff_mid (i,5) = 0.5*xff_mid(i,4)-xff_mid(i,2)
                xff_mid (i,5) = 0.05*xff_mid(i,4)!-xff_mid(i,2)
                xf_dicycle(i) = 0.
          enddo
         elseif(DICYCLE==0) then
          do i=its,itf
                if(ierr(i) /= 0) cycle
                !xff_mid (i,5) = 0.25*xff_mid(i,4)
                xff_mid (i,5) = 0.05*xff_mid(i,4)
                xf_dicycle(i) = 0.
          enddo
         endif
       ENDIF

 END SUBROUTINE cup_forcing_ens_3d_mid
!------------------------------------------------------------------------------------

   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              itf,ktf,its,ite, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         itf,ktf,                                    &
         its,ite, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

       DO i=its,itf
         KT(I)=KS(I)
         if(ierr(i).eq.0)then
           X(I)=ARRAY(I,KS(I))
           KSTOP=MAX(KS(I)+1,KEND(I))
!
           DO K=KS(I)+1,KSTOP
              IF(ARRAY(I,K).LT.X(I)) THEN
                X(I)=ARRAY(I,K)
                KT(I)=K
              ENDIF
           ENDDO
         endif
       ENDDO

   END SUBROUTINE cup_MINIMI
!------------------------------------------------------------------------------------
   SUBROUTINE cup_up_aa0(aa0,z_cup,zu,dby,GAMMA_CUP,t_cup,   &
              k22,klcl,kbcon,ktop,ierr,                      &
              itf,ktf,its,ite, kts,kte,integ_interval        )
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     IMPLICIT NONE
  !
  ! on input
     integer,intent (in   )              ::  itf,ktf,its,ite, kts,kte
     real,    dimension (its:ite,kts:kte) ,intent (in   )  ::  &
        z_cup,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)         ,intent (in   )  ::  &
        k22,klcl,kbcon,ktop
     character(LEN=*),OPTIONAL, intent(in) :: integ_interval
  ! input and output
     integer, dimension (its:ite) ,intent (inout)  :: ierr
     real,    dimension (its:ite) ,intent (out  )  :: aa0
  !
  !  local variables in this routine
     integer                             ::  i,k
     real                                ::  dz,da,aa_2,aa_1
     integer, dimension (its:ite)        ::  kbeg,kend
  !
  !
  !  initialize array to zero.
      aa0(:)=0.
  !  set domain of integration
      IF(present(integ_interval)) then
       if(integ_interval == 'BL') then
        kbeg(:) = kts
        kend(:) = kbcon(:)-1
       elseif(integ_interval == 'CIN') then
        kbeg(:) = k22(:) !klcl (:) ! kts
        kend(:) = kbcon(:)-1
       else
        STOP "unknown range in cup_up_aa0"
       endif
      ELSE
        kbeg(:) = kbcon(:)
        kend(:) = ktop (:)
      ENDIF

      do i=its,itf
        if(ierr(i) /= 0)cycle
        do k= kbeg(i),kend(i)

          dz=z_cup(i,k+1)-z_cup(i,k)
          aa_1=zu(i,k  )*(g/(cp*t_cup(i,k  )))*dby(i,k  )/(1.+gamma_cup(i,k  ))
          aa_2=zu(i,k+1)*(g/(cp*t_cup(i,k+1)))*dby(i,k+1)/(1.+gamma_cup(i,k+1))
          da=0.5*(aa_1+aa_2)*dz

          aa0(i)=aa0(i)+da

          !aa0(i)=aa0(i)+max(0.,da)
        enddo
      enddo

   END SUBROUTINE cup_up_aa0
!------------------------------------------------------------------------------------

   SUBROUTINE cup_up_moisture(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc,xland,&
                                  cnvfrc,srftype,po,p_cup,kbcon,ktop,cd,dby,clw_all,                  &
                                  t_cup,q,gamma_cup,zu,qes_cup,k22,qe_cup,            &
                                  zqexec,use_excess,ccn,rho,                          &
                                  up_massentr,up_massdetr,psum,psumh,c1d,x_add_buoy,  &
                                  vvel2d,vvel1d,zws,entr_rate,                          &
                                  itest,itf,ktf,ipr,jpr,its,ite, kts,kte                  )

    implicit none
    real, parameter :: bdispm = 0.366       !berry--size dispersion (maritime)
    real, parameter :: bdispc = 0.146       !berry--size dispersion (continental)
    real, parameter :: T_BF   = 268.16 , T_ice_BF = 235.16
    real, parameter :: rk = 3 ! or 2
    real, parameter :: xexp = 2.
!
!  on input
     integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
                                 ,its,ite,ipr,jpr, kts,kte
  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     character *(*)                    ,intent (in) ::  name
     integer, dimension (its:ite)      ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
     real,  dimension (its:ite,kts:kte),intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
                                                       ,qe_cup,hc,po,up_massentr,up_massdetr &
                                                       ,dby,qes_cup,z_cup,cd,c1d

     real,  dimension (its:ite)        ,intent (in) ::  cnvfrc,srftype
     real,  dimension (its:ite)        ,intent (in) ::  zqexec,xland,x_add_buoy
     real,  dimension (its:ite)        ,intent (in) ::  zws,ccn
     real,  dimension (its:ite,kts:kte),intent (in) ::  entr_rate
     real,  dimension (its:ite,kts:kte),intent (in) ::  vvel2d
     real,  dimension (its:ite        ),intent (in) ::  vvel1d
!
! input and output
!
   ! ierr error value, maybe modified in this routine
     integer, dimension (its:ite)  ,intent (inout)                   ::  ierr
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,   dimension (its:ite,kts:kte),intent (out)   :: qc,qrc,pw,clw_all,tempc
     real,   dimension (its:ite)        ,intent (out)   :: pwav,psum,psumh
     character*128                      ,intent (inout) :: ierrc(its:ite)
!
!  local variables in this routine
!
     integer                              ::                           &
        iounit,iprop,i,k,k1,k2,n,nsteps
     real                                 ::                           &
        dp,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,qavail
     real delt,tem1,cup

        !--- no precip for small clouds
        if(name.eq.'shallow')  c0 = c0_shal
        if(name.eq.'mid'    )  c0 = c0_mid
        if(name.eq.'deep'   )  c0 = c0_deep
        do i=its,itf
          pwav (i)=0.
          psum (i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
         do i=its,itf
          pw      (i,k)=0.
          clw_all (i,k)=0.
          tempc   (i,k)=t_cup (i,k)
          qrc     (i,k)=0.          !--- liq/ice water
          qc      (i,k)=qe_cup(i,k) !--- total water: liq/ice = vapor water
         enddo
        enddo

        !--- get boundary condition for qc
        do i=its,itf
            if(ierr(i)  /= 0) cycle
            call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),qe_cup (i,kts:kte),qaver,k22(i))
            qc  (i,kts:start_level(i)) = qaver + zqexec(i) + 0.5*x_add_buoy(i)/xlv
            qrc (i,kts:start_level(i)) = 0.
        enddo

        !--- option to produce linear fluxes in the sub-cloud layer.
        if(name == 'shallow' .and. use_linear_subcl_mf == 1) then
             do i=its,itf
                if(ierr(i) /= 0) cycle
                call get_delmix(name,kts,kte,ktf,xland(i),start_level(i),po(i,kts:kte) &
                               ,qe_cup(i,kts:kte), qc(i,kts:kte))
            enddo
        endif
        DO i=its,itf
          IF (ierr(i) /= 0) cycle

          DO k=start_level(i) + 1,ktop(i) + 1

            DZ=Z_cup(i,K)-Z_cup(i,K-1)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            QRCH = qes_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom =  (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
            if(denom > 0.) then

                qc (i,k)=  ( qc (i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1) +   &
                                                     up_massentr(i,k-1)* q (i,k-1)     &
                           )/ denom

                if(k==start_level(i)+1) qc(i,k)= qc(i,k) + zqexec(i) &
                                                         * up_massentr(i,k-1)/denom
                !--- assuming no liq/ice water in the environment
                qrc(i,k)=  ( qrc(i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qrc(i,k-1)   &
                           )/ denom

            else
                qc (i,k)=    qc (i,k-1)
                qrc(i,k)=    qrc(i,k-1)
            endif

            !-- updraft temp
            tempc(i,k) = (1./cp)*(hc(i,k)-g*z_cup(i,k)-xlv*QRCH)

            !--- total condensed water before rainout
            clw_all(i,k)= max(0.,qc(i,k)-qrch)

            qrc   (i,k) = min(clw_all(i,k),qrc(i,k))

            !--- production term => condensation/diffusional growth
            cup         = max(0.,qc(i,k)-qrch-qrc(i,k))/dz

            if(c0 < 1.e-6)then
                qrc (i,k) = clw_all(i,k)
                qc  (i,k) = qrc(i,k)+min(qc(i,k),qrch)
                pwav(i)   = 0.
                psum(i)   = psum(i)+clw_all(i,k)*zu(i,k) *dz
                cycle
            endif

            IF (AUTOCONV == 1 ) then
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd )
                cx0     = (c1d(i,k)+c0)*DZ
                qrc(i,k)= clw_all(i,k)/(1.+cx0)
                pw (i,k)= cx0*max(0.,qrc(i,k) - min_liq)! units kg[rain]/kg[air]
                !--- convert pw to normalized pw
                pw (i,k)=pw(i,k)*zu(i,k)

            ELSEIF (AUTOCONV == 2 ) then
              ! this is similar to AUTOCONV == 1 with temperature dependence
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd )
                cx0     = (c1d(i,k)+c0)*DZ*fract_liq_f(tempc(i,k),cnvfrc(i),srftype(i))
                qrc(i,k)= clw_all(i,k)/(1.+cx0) 
                pw (i,k)= cx0*max(0.,qrc(i,k) - min_liq)! units kg[rain]/kg[air]
                !--- convert PW to normalized PW
                pw (i,k)=pw(i,k)*zu(i,k)

            ELSEIF (AUTOCONV == 3 ) then
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd )
                if(clw_all(i,k) <= min_liq) then
                   qrc(i,k)= clw_all(i,k)
                   pw(i,k) = 0.
                else
                   cx0     = c0*fract_liq_f(tempc(i,k),cnvfrc(i),srftype(i))
                   cx0     = max(cx0,0.50*c0)
                   qrc(i,k)= qrc(i,k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                   qrc(i,k)= min(clw_all(i,k), qrc(i,k))
                   pw (i,k)= clw_all(i,k) - qrc(i,k)
                  !--- convert pw to normalized pw
                   pw (i,k)= pw(i,k)*zu(i,k)
                endif

            ELSEIF (AUTOCONV == 4 ) then
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd )
                if(clw_all(i,k) <= min_liq) then
                   qrc(i,k)= clw_all(i,k)
                   pw(i,k) = 0.
                else
                  tem1 = fract_liq_f(tempc(i,k),cnvfrc(i),srftype(i))
                  cbf  = 1.
                  if(tempc(i,k) < T_BF) cbf=1.+0.5*sqrt(min(max(T_BF-tempc(i,k),0.),T_BF-T_ice_BF))
                  qrc_crit_BF = ccn(i)/cbf
                  cx0 = c0*cbf*(tem1*1.3+(1.-tem1))/(0.75*min(15.,max(vvel2d(i,k),1.)))
                  !---analytical solution
                  cx0     = cx0* (1.- exp(- (qrc(i,k)/qrc_crit_BF)**2))
                  qrc(i,k)= qrc(i,k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                  !---
                  pw (i,k)= max(clw_all(i,k) - qrc(i,k),0.)
                  !--- convert PW to normalized PW
                  pw (i,k)= pw(i,k)*zu(i,k)
                endif

            ELSEIF (AUTOCONV == 5 ) then
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd )
                if(clw_all(i,k) <= min_liq) then
                   qrc(i,k)= clw_all(i,k)
                   pw(i,k) = 0.
                else
                   cx0     = (c1d(i,k)+c0)*(1.+ 0.33*fract_liq_f(tempc(i,k),cnvfrc(i),srftype(i)))
                   qrc(i,k)= qrc(i,k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                   pw (i,k)= max(0.,clw_all(i,k)-qrc(i,k)) ! units kg[rain]/kg[air]
                   qrc(i,k)= clw_all(i,k)-pw(i,k)
                  !--- convert pw to normalized pw
                   pw (i,k)= pw(i,k)*zu(i,k)
                endif

            ELSEIF (AUTOCONV == 6 ) then
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd ) 
                if(clw_all(i,k) <= min_liq) then
                   qrc(i,k)= clw_all(i,k)
                   pw(i,k) = 0.
                else
                   cx0     = (c1d(i,k)+c0)*dz
                   qrc(i,k)= (clw_all(i,k))*exp(-cx0)
                   pw (i,k)= clw_all(i,k) - qrc(i,k)
                  !--- convert pw to normalized pw
                   pw (i,k)= pw(i,k)*zu(i,k)
                endif

            ELSEIF (AUTOCONV == 7 ) then
                min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd ) 
                if(clw_all(i,k) <= min_liq) then
                   qrc(i,k)= clw_all(i,k)
                   pw(i,k) = 0.
                else
                   cx0     = c1d(i,k)+c0
                   qrc(i,k)= qrc(i,k)*exp(-cx0*dz) + (cup/cx0)*(1.-exp(-cx0*dz))
                   pw (i,k)= max(clw_all(i,k) - qrc(i,k),0.)
                   qrc(i,k)= clw_all(i,k) - pw (i,k)
                  !--- convert pw to normalized pw
                   pw (i,k)= pw(i,k)*zu(i,k)
                endif

            ENDIF
            !- total water (vapor + condensed) in updraft after the rainout
            qc(i,k)=qrc(i,k)+min(qc(i,k),qrch)

            !--- integrated normalized condensates
            pwav(i)=pwav(i)+pw(i,k)
            psum(i)=psum(i)+clw_all(i,k)*zu(i,k) *dz

          ENDDO

          IF(ZERO_DIFF==0) THEN
            if(pwav(i) < 0.) then
              ierr(i)=66
              ierrc(i)="pwav negative"
            endif
          ENDIF

        ENDDO

        !--- get back water vapor qc
         do i=its,itf
             if (ierr(i)  /= 0) cycle
             do k=kts,ktop(i)+1
               qc(i,k)= qc(i,k)-qrc(i,k)
               !if(qc(i,k) < 0.)stop " qc negative"
             enddo
         enddo

   END SUBROUTINE cup_up_moisture

!------------------------------------------------------------------------------------
   SUBROUTINE cup_up_moisture_light(name,start_level,klcl,ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc,xland &
                                   ,cnvfrc,srftype,po,p_cup,kbcon,ktop,cd,dby,clw_all,t_cup,q,gamma_cup,zu  &
                                   ,qes_cup,k22,qe_cup,zqexec,use_excess,rho                 &
                                   ,up_massentr,up_massdetr,psum,psumh,c1d,x_add_buoy        &
                                   ,itest,itf,ktf,ipr,jpr,its,ite, kts,kte                   )

    implicit none
!
!  on input
     integer  ,intent (in   ) ::  use_excess,itest,itf,ktf    &
                                 ,its,ite,ipr,jpr, kts,kte
  ! cd= detrainment function
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  !
     character *(*)                    ,intent (in) ::  name
     integer, dimension (its:ite)      ,intent (in) ::  kbcon,ktop,k22,klcl,start_level
     real,  dimension (its:ite,kts:kte),intent (in) ::  t_cup,p_cup,rho,q,zu,gamma_cup       &
                                                       ,qe_cup,hc,po,up_massentr,up_massdetr &
                                                       ,dby,qes_cup,z_cup,cd,c1d

     real,  dimension (its:ite)        ,intent (in) ::  zqexec,xland,x_add_buoy,cnvfrc,srftype
!
! input and output
!
   ! ierr error value, maybe modified in this routine
     integer, dimension (its:ite)  ,intent (inout)                   ::  ierr
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,   dimension (its:ite,kts:kte),intent (out)   :: qc,qrc,pw,clw_all,tempc
     real,   dimension (its:ite)        ,intent (out)   :: pwav,psum,psumh
     character*128                      ,intent (inout) :: ierrc(its:ite)
!
!  local variables in this routine
!
     integer                              ::                           &
        iounit,iprop,i,k,k1,k2,n,nsteps
     real                                 ::                           &
        dp,rhoc,dh,qrch,c0,dz,radius,berryc0,q1,berryc
     real :: qaver,denom,aux,cx0,qrci,step,cbf,qrc_crit_BF,min_liq,qavail,delt_hc_glac
     real delt,tem1

        !--- no precip for small clouds
        if(name.eq.'shallow')  c0 = c0_shal
        if(name.eq.'mid'    )  c0 = c0_mid
        if(name.eq.'deep'   )  c0 = c0_deep
        do i=its,itf
          pwav (i)=0.
          psum (i)=0.
          psumh(i)=0.
        enddo
        do k=kts,ktf
         do i=its,itf
          pw      (i,k)=0.
          qrc     (i,k)=0.
          clw_all (i,k)=0.
          tempc   (i,k)=t_cup (i,k)
          qc      (i,k)=qe_cup(i,k)
         enddo
        enddo

        !--- get boundary condition for qc
        do i=its,itf
            if(ierr(i)  /= 0) cycle
            call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),qe_cup (i,kts:kte),qaver,k22(i))
            qc (i,kts:start_level(i)) = qaver + zqexec(i) + 0.5*x_add_buoy(i)/xlv
        enddo

        DO i=its,itf
          IF (ierr(i)  /= 0) cycle

          DO k=start_level(i)+1,ktop(i) + 1

            DZ=Z_cup(i,K)-Z_cup(i,K-1)
            !
            !--- saturation  in cloud, this is what is allowed to be in it
            !
            QRCH = qes_cup(I,K)+(1./XLV)*(GAMMA_cup(i,k)/(1.+GAMMA_cup(i,k)))*DBY(I,K)

            !-    1. steady state plume equation, for what could
            !-       be in cloud without condensation
            denom =  (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
            if(denom > 0.) then
                qc (i,k)=  ( qc (i,k-1)*zu(i,k-1)-.5*up_massdetr(i,k-1)* qc(i,k-1) +   &
                                                     up_massentr(i,k-1)* q (i,k-1)     &
                           )/ denom
                if(k==start_level(i)+1) qc(i,k)= qc(i,k) + zqexec(i) &
                                                         * up_massentr(i,k-1)/denom
            else
                qc (i,k)=    qc (i,k-1)
            endif

            !--- total condensed water before rainout
            clw_all(i,k)=max(0.,QC(I,K)-QRCH)
            !--- updraft temp
            tempc(i,k) = (1./cp)*(hc(i,k) - g*z_cup(i,k)-xlv*QRCH)

            !--add glaciation effect on the MSE
            if(MELT_GLAC) then
               delt_hc_glac = clw_all(i,k)*(1.- fract_liq_f(tempc(i,k),cnvfrc(i),srftype(i)))*xlf

               tempc(i,k) = tempc(i,k)+(1./cp)*delt_hc_glac
            endif

            cx0     = (c1d(i,k)+c0)*DZ
            if(c0 < 1.e-6) cx0 = 0.

            qrc(i,k)= clw_all(i,k)/(1.+cx0)
            min_liq  = ( xland(i)*qrc_crit_ocn + (1.-xland(i))*qrc_crit_lnd ) 
            pw (i,k)= cx0*max(0.,qrc(i,k) - min_liq)! units kg[rain]/kg[air]
            !--- convert pw to normalized pw
            pw (i,k)= pw(i,k)*zu(i,k)

            !- total water (vapor + condensed) in updraft after the rainout
            qc(i,k)=qrc(i,k)+min(qc(i,k),qrch)

          ENDDO
         ENDDO

        !- get back water vapor qc
         do i=its,itf
             if (ierr(i)  /= 0) cycle
             do k=kts,ktop(i)+1
               qc(i,k)= qc(i,k)-qrc(i,k)
             enddo
         enddo

   END SUBROUTINE cup_up_moisture_light

!------------------------------------------------------------------------------------

   REAL FUNCTION satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
              ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
        toot = 273.16 / temp2
        toto = 1 / toot
        eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
          log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
        satvap = 10 ** eilog
      else
        tsot = 373.16 / temp2
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
                   (log(tsot) / log(10.))
        ewlog2 = ewlog - 1.3816e-07 * &
                   (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + .0081328 * &
                   (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.))
        satvap = 10 ** ewlog4
      endif

   END FUNCTION
!------------------------------------------------------------------------------------

   SUBROUTINE cup_up_aa1bl(version,aa1_bl,aa1_fa,aa1,t,tn,q,qo,dtime,po_cup, &
              z_cup,zu,dby,GAMMA_CUP,t_cup,rho,                  &
              klcl,kpbl,kbcon,ktop,ierr,itf,ktf,its,ite, kts,kte,&
              xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux,tau_bl,tau_ecmwf,t_star,cumulus ,tn_bl,qo_bl )

   IMPLICIT NONE
   character*(*), intent(in)                         :: cumulus
   integer      , intent (in   )                     :: &
        itf,ktf,its,ite, kts,kte,version
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte) ,intent (in  ) ::    &
        z_cup,zu,gamma_cup,t_cup,dby,t,tn,q,qo,po_cup,rho,tn_bl,qo_bl
     integer, dimension (its:ite)         ,intent (in  ) ::    &
        klcl,kbcon,ktop,kpbl
     real                                 ,intent(in   ) :: dtime, t_star

     real,    dimension (its:ite)         ,intent (IN  ) ::    &
        xland,ztexec,xlons,xlats, h_sfc_flux,le_sfc_flux       &
       ,aa1,tau_bl,tau_ecmwf
!
! input and output
     integer, dimension (its:ite),intent (inout) ::  ierr
     real,    dimension (its:ite),intent (out  ) ::  aa1_bl,aa1_fa
!  local variables in this routine
!
     integer                              ::    i,k,iprloc
     real                                 ::    dz,da,aa_1,aa_2,tcup,da_bl,a1_bl

!==================
     aa1_bl (:)=0.
     aa1_fa (:)=0.
     IF(version == 0 ) then
        do i=its,itf

           if(ierr(i) /= 0 ) cycle
!***       do k=kts,kbcon(i)
           do k=kts,kpbl(i)
              dz = (z_cup (i,k+1)-z_cup (i,k))*g
              da = dz*(tn(i,k)*(1.+0.608*qo(i,k))-t(i,k)*(1.+0.608*q(i,k)))/dtime

!--
!             tcup=0.5*(t_cup(i,k+1)+t_cup(i,k))
!             da=(da/tcup)*dtime !UNIT J/kg
!--
              aa1_bl(i)=aa1_bl(i)+da ! Units : J K / (kg seg)
            enddo
        enddo
     ELSEIF(version==1) then
        do i=its,itf
           if(ierr(i) /= 0 ) cycle
           do k=kts,klcl(i)
             dz = (z_cup (i,k+1)-z_cup (i,k))
             aa_1=(g/(cp*t_cup(i,k  )))*dby(i,k  )*zu(i,k  )
             aa_2=(g/(cp*t_cup(i,k+1)))*dby(i,k+1)*zu(i,k+1)
             da=0.5*(aa_1+aa_2)*dz! Units : J / kg
             aa1_bl(i)=aa1_bl(i)+da
             write(10,*) k,dby(i,k  ),da,aa1_bl(i),zu(i,k)
           enddo
           do k=klcl(i)+1,kbcon(i)-1
             dz = (z_cup (i,k+1)-z_cup (i,k))
             aa_1=(g/(cp*t_cup(i,k  )))*dby(i,k  )/(1.+gamma_cup(i,k  ))*zu(i,k  )!
             aa_2=(g/(cp*t_cup(i,k+1)))*dby(i,k+1)/(1.+gamma_cup(i,k+1))*zu(i,k+1)!
             da=0.5*(aa_1+aa_2)*dz! Units : J / kg
             aa1_bl(i)=aa1_bl(i)+da
             write(10,*)k,dby(i,k  ),da,aa1_bl(i),zu(i,k)
           enddo

        enddo
     ELSE
        stop "unknown version option in routine: cup_up_aa1bl"
     ENDIF

     return

     do i=its,itf
        if(ierr(i) /= 0)cycle
        do k= kbcon(i),ktop(i)

          dz=z_cup(i,k+1)-z_cup(i,k)
          aa_1=(g/(cp*((t_cup(i,k  )))))*dby(i,k  )/(1.+gamma_cup(i,k  ))*zu(i,k)
          aa_2=(g/(cp*((t_cup(i,k+1)))))*dby(i,k+1)/(1.+gamma_cup(i,k+1))*zu(i,k+1)
          da=0.5*(aa_1+aa_2)*dz

          aa1_fa(i)=aa1_fa(i)+da
        enddo
      enddo

  END SUBROUTINE cup_up_aa1bl
!------------------------------------------------------------------------------------
  SUBROUTINE get_lateral_massflux(itf,ktf, its,ite, kts,kte                             &
                                  ,ierr,ktop,zo_cup,zuo,cd,entr_rate,po_cup          &
                                  ,up_massentro, up_massdetro ,up_massentr, up_massdetr &
                                  ,draft,kbcon,k22,kpbl,up_massentru,up_massdetru,lambau)
     implicit none
     character *(*), intent (in) :: draft
     integer, intent(in) :: itf,ktf, its,ite, kts,kte
     integer, intent(in)   , dimension(its:ite)            :: ierr,ktop,kbcon,k22,kpbl
     real,    intent(in)   , dimension(its:ite), OPTIONAL  :: lambau
     real,    intent(in)   , dimension(its:ite,kts:kte) :: zo_cup,zuo,po_cup
     real,    intent(in)   , dimension(its:ite,kts:kte) :: entr_rate
     real,    intent(inout), dimension(its:ite,kts:kte) :: cd
     real,    intent(  out), dimension(its:ite,kts:kte) :: up_massentro, up_massdetro&
                                                          ,up_massentr,  up_massdetr
     real,    intent(  out), dimension(its:ite,kts:kte),  &
                                               OPTIONAL :: up_massentru,up_massdetru
     !
     !-- local vars
     integer :: i,k, turn, ismooth1,ismooth2
     real :: dz, mass1,mass2,dp,rho,dz_zuo_ave
     integer, parameter :: MASS_U_OPTION = 1
     integer :: k_ent
     !---

     up_massentro(:,:)=0.
     up_massdetro(:,:)=0.
     if(present(up_massentru) .and. present(up_massdetru))then
        up_massentru(:,:)=0.
        up_massdetru(:,:)=0.
     endif

     do i=its,itf
         if(ierr(i)/= 0)cycle

         k_ent=maxloc(zuo(i,:),1)
         !-will not allow detrainment below cloud base or in the PBL
         if(draft=='shallow') then
              cd(i,1:max(kbcon(i),kpbl(i))+1)=0.0
         else
              cd(i,1:k_ent+1)=0.0
         endif

        !- mass entrainment and detrainment are defined on model levels
         do k=kts,k_ent-1
           !-- below location of maximum value zu -> change entrainment

           dz_zuo_ave = (zo_cup(i,k+1)-zo_cup(i,k))*zuo(i,k)

           up_massdetro(i,k)=cd(i,k)*dz_zuo_ave

           up_massentro(i,k)=zuo(i,k+1)-zuo(i,k)+up_massdetro(i,k)
           up_massentro(i,k)=max(up_massentro(i,k),0.0)

           !-- check dd_massdetro in case of dd_massentro has been changed above
           up_massdetro(i,k)=-zuo(i,k+1)+zuo(i,k)+up_massentro(i,k)
         enddo

         do k=k_ent,ktop(i)
           !-- above location of maximum value zu -> change detrainment
           dz_zuo_ave = (zo_cup(i,k+1)-zo_cup(i,k))*zuo(i,k)

           up_massentro(i,k)=entr_rate(i,k)*dz_zuo_ave

           up_massdetro(i,k)=zuo(i,k)+up_massentro(i,k)-zuo(i,k+1)
           up_massdetro(i,k)=max(up_massdetro(i,k),0.0)
           !-- check up_massentro in case of dd_up_massdetro has been changed above
           up_massentro(i,k)=-zuo(i,k)+up_massdetro(i,k)+zuo(i,k+1)
         enddo

         do k=kts,kte
           up_massentr(i,k)=up_massentro(i,k)
           up_massdetr(i,k)=up_massdetro(i,k)
         enddo
         IF(present(up_massentru) .and. present(up_massdetru))THEN
           if(mass_U_option==1) then
              do k=kts+1,kte
              !--       for weaker mixing
                up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
              !--       for stronger mixing
              ! up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
              ! up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
              enddo
           else
              do k=kts+1,k_ent
                up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massentro(i,k-1)
                up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massentro(i,k-1)
              enddo
              do k=k_ent+1,kte
                up_massentru(i,k-1)=up_massentro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
                up_massdetru(i,k-1)=up_massdetro(i,k-1)+lambau(i)*up_massdetro(i,k-1)
              enddo
           endif
         ENDIF

    ENDDO
   !!---- check mass conservation
   !DO i=its,itf
   !     if(ierr(i)/= 0)cycle
   !     do k=kts+1,kte
   !        dz =      zo_cup(i,k)-zo_cup(i,k-1)
   !        dp = 100*(po_cup(i,k)-po_cup(i,k-1))
   !        rho= -dp/dz/g
   !        mass1= (zuo(i,k)-zuo(i,k-1)) - up_massentro(i,k-1)+up_massdetro(i,k-1)
   !        !print*,"masscons=",mass1!,-rho*g*(zuo(i,k)-zuo(i,k-1))/dp, (zuo(i,k)-zuo(i,k-1))/dz,( up_massentro(i,k-1)-up_massdetro(i,k-1))/dz,rho
   !        mass2= (zuo(i,k)-zuo(i,k-1)) - up_massentru(i,k-1)+up_massdetru(i,k-1)
   !     enddo
   !ENDDO
 END SUBROUTINE get_lateral_massflux

!------------------------------------------------------------------------------------
 SUBROUTINE get_lateral_massflux_down(cumulus,itf,ktf, its,ite, kts,kte                    &
                                ,ierr,jmin,zo_cup,zdo,xzd,zd,cdd,mentrd_rate         &
                                ,dd_massentro,dd_massdetro ,dd_massentr, dd_massdetr &
                                ,draft,dd_massentru,dd_massdetru,lambau)

     implicit none
     character *(*), intent (in) :: draft,cumulus
     integer, intent(in):: itf,ktf, its,ite, kts,kte
     integer, intent(in)   , dimension(its:ite)         :: ierr,jmin
     real,    intent(in)   , dimension(its:ite        ) :: lambau
     real,    intent(in)   , dimension(its:ite,kts:kte) :: zo_cup,zdo
     real,    intent(inout), dimension(its:ite,kts:kte) :: cdd,mentrd_rate,xzd,zd
     real,    intent(  out), dimension(its:ite,kts:kte) :: dd_massentro, dd_massdetro&
                                                          ,dd_massentr,  dd_massdetr
     real,    intent(  out), dimension(its:ite,kts:kte), OPTIONAL &
                                                        :: dd_massentru, dd_massdetru
     integer ::i,ki
     real :: dzo

     cdd          = 0.
     dd_massentr  = 0.
     dd_massdetr  = 0.
     dd_massentro = 0.
     dd_massdetro = 0.
     if(present(dd_massentru).and.present(dd_massdetru))then
       dd_massentru = 0.
       dd_massdetru = 0.
     endif
     if(trim(cumulus) == 'shallow') return

     do i=its,itf
        if(ierr(i) /= 0) cycle

        cdd           (i,1:jmin(i)-1) =mentrd_rate(i,1:jmin(i)-1)
        mentrd_rate(i,1)=0.

        do ki=jmin(i)   ,maxloc(zdo(i,:),1),-1

          !=> from jmin to maximum value zd -> change entrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massdetro(i,ki)=cdd(i,ki)*dzo*zdo(i,ki+1)
!XXX
          dd_massentro(i,ki)= zdo(i,ki)-zdo(i,ki+1)+dd_massdetro(i,ki)
          dd_massentro(i,ki)= MAX(0.,dd_massentro(i,ki))
          !-- check dd_massdetro in case of dd_massentro has been changed above
          dd_massdetro(i,ki)= dd_massentro(i,ki)-zdo(i,ki)+zdo(i,ki+1)

          !~ if(dd_massentro(i,ki).lt.0.)then
             !~ dd_massentro(i,ki)=0.
             !~ dd_massdetro(i,ki)=zdo(i,ki+1)-zdo(i,ki)
             !~ if(zdo(i,ki+1) > 0.0)&
               !~ cdd(i,ki)=dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
          !~ endif
          !~ if(zdo(i,ki+1) > 0.0)&
            !~ mentrd_rate(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
        enddo

        do ki=maxloc(zdo(i,:),1)-1,kts,-1
          !=> from maximum value zd to surface -> change detrainment
          dzo=zo_cup(i,ki+1)-zo_cup(i,ki)
          dd_massentro(i,ki)=mentrd_rate(i,ki)*dzo*zdo(i,ki+1)
!XXX
          dd_massdetro(i,ki) = zdo(i,ki+1)+dd_massentro(i,ki)-zdo(i,ki)
          dd_massdetro(i,ki) = MAX(0.0,dd_massdetro(i,ki))
          !-- check dd_massentro in case of dd_massdetro has been changed above
          dd_massentro(i,ki) = dd_massdetro(i,ki)+zdo(i,ki)-zdo(i,ki+1)


          !~ if(dd_massdetro(i,ki).lt.0.)then
            !~ dd_massdetro(i,ki)=0.
            !~ dd_massentro(i,ki)=zdo(i,ki)-zdo(i,ki+1)
            !~ if(zdo(i,ki+1) > 0.0)&
              !~ mentrd_rate(i,ki)=dd_massentro(i,ki)/(dzo*zdo(i,ki+1))
          !~ endif
          !~ if(zdo(i,ki+1) > 0.0)&
            !~ cdd(i,ki)= dd_massdetro(i,ki)/(dzo*zdo(i,ki+1))
        enddo

        do ki=jmin(i),kts,-1
          xzd(i,ki)= zdo(i,ki)
          zd (i,ki)= zdo(i,ki)
          dd_massentr(i,ki)=dd_massentro(i,ki)
          dd_massdetr(i,ki)=dd_massdetro(i,ki)
        enddo
        if(present(dd_massentru).and.present(dd_massdetru))then
          do ki=jmin(i),kts,-1
            dd_massentru(i,ki)=dd_massentro(i,ki)+lambau(i)*dd_massdetro(i,ki)
            dd_massdetru(i,ki)=dd_massdetro(i,ki)+lambau(i)*dd_massdetro(i,ki)
          enddo
        endif
    enddo

 END SUBROUTINE get_lateral_massflux_down
!------------------------------------------------------------------------------------

 SUBROUTINE get_zi_gf(j,its,ite,kts,kte,istart,iend,ktf,ierr,kzi,pbl,tkeg, &
                  rcpg,z,ztop,tkmin)

  implicit none
  integer,intent(in):: j,its,ite,kts,kte, ktf,istart,iend
  integer :: kzimax,ktke_max,i,k
  real tkmin,tke_tmp
  real,    dimension(its:ite,kts:kte) :: tkeg,rcpg,z
  real,    dimension(its:ite)          :: ztop,pbl
  integer, dimension(its:ite)          :: kzi,ierr

  real, parameter :: rcpmin=1.e-5 , pblhmax=3000.

  do i=istart,iend
     kzi(i)  = 1
!    if(ierr(i).eq.0)then
!         tke_tmp = 0.
         ktke_max= 1
         kzimax  =ktf-1
         !---  max level for kzi
         DO K=kts,ktf
           if(z(i,k).ge. pblhmax+ztop(i)) then
              kzimax = min(k,ktf-1)
              !if(j==8 .and. i==10) print*,"1",z(i,k), pblhmax,ztop(i),kzimax
              exit
           endif
         enddo
         !---

!         !level of max tke  below kzimax and w/out clouds
!         do  k=kts,kzimax
!           !print*,k,tkeg(i,k), tke_tmp,ktke_max,kzimax
!           if(rcpg(i,k) .lt. rcpmin) then
!             if( tkeg(i,k) .ge. tke_tmp) then
!               tke_tmp = tkeg(i,k)
!               cycle
!             else
!               ktke_max= max(1,k-1)
!               exit
!             endif
!           endif
!         enddo
!    201         continue
!             print*,ktke_max

         do k=ktke_max,kzimax+1

!           if(tkeg(i,k) .gt. 1.1*tkmin .and. rcpg(i,k) .lt. rcpmin )  then
            if(tkeg(i,k) .gt. 1.1*tkmin )  then
              kzi(i) = k
              !if(j==8 .and. i==10) print*,"I",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
              cycle

            else
               kzi(i) = max(1,k-1)
               !if(j==8 .and. i==10) print*,"F",k,rcpg(i,k),tkeg(i,k),kzi(i),z(i,k)-ztop(i)
               exit
            endif


         enddo
         kzi(i) = max(1     ,kzi(i))
         !print*,"1",kzi(i),i
         kzi(i) = min(kzimax,kzi(i))
         !print*,"2",kzi(i),i;call flush(6)
         pbl(i) = max( z(i,kzi(i))-ztop(i), z(i,1)-ztop(i) )
 enddo

 END SUBROUTINE get_zi_gf
!------------------------------------------------------------------------------------

 SUBROUTINE rates_up_pdf(name,ktop,ierr,p_cup,entr_rate,hkbo,heo,heso_cup,z_cup, &
                         kstabi,k22,kbcon,its,ite,itf,kts,kte,ktf,zuo,kpbl,klcl,hcot)
     implicit none
     character *(*)                      ,intent (in) :: name
     integer                             ,intent (in) :: its,ite,itf,kts,kte,ktf
     integer, dimension (its:ite)        ,intent (in) :: kstabi,k22,kbcon,kpbl,klcl
     real   , dimension (its:ite,kts:kte),intent (in) :: entr_rate,zuo
     real   , dimension (its:ite,kts:kte),intent (in) :: p_cup, heo,heso_cup,z_cup
     real   , dimension (its:ite)        ,intent (in) :: hkbo
     integer, dimension (its:ite)        ,intent (inout) :: ierr,ktop
     !--local vars
     real   , dimension (its:ite,kts:kte) :: hcot
     real :: dz,dh,Z_overshoot
     integer :: i,k,ipr,kdefi,kstart,kbegzu,kfinalzu
     integer, dimension (its:ite) :: start_level
     real  :: delz_oversh !--- height of cloud overshoot is 10% higher than the LNB.
                                          !--- Typically it can 2 - 2.5km higher, but it depends on
                                          !--- the severity of the thunderstorm.

     delz_oversh = OVERSHOOT
     hcot        = 0.0
     if(name == 'shallow') return

     DO i=its,itf
       ktop(i)=ktf-2
       if(ierr(i) /= 0)cycle

       start_level(i)=kbcon(i)
       !-- hcot below kbcon
       hcot(i,kts:start_level(i))=hkbo(i)

       do k=start_level(i)+1,ktf-2
           dz=z_cup(i,k)-z_cup(i,k-1)

           hcot(i,k)=( (1.-0.5*entr_rate(i,k-1)*dz)*hcot(i,k-1)    &
                              +entr_rate(i,k-1)*dz *heo (i,k-1) )/ &
                       (1.+0.5*entr_rate(i,k-1)*dz)
       enddo
       do k=start_level(i)+1,ktf-2
          if(hcot(i,k) < heso_cup(i,k) )then
              ktop(i) = k - 1
              exit
          endif
       enddo
       if(ktop(i).le.kbcon(i)+1) ierr(i)=41
       !----------------
       if(OVERSHOOT > 0. .and. ierr(i) == 0) then
           Z_overshoot = (1. + delz_oversh) * z_cup(i,ktop(i))
           do k=ktop(i),ktf-2
              if(Z_overshoot < z_cup(i,k)) then
                !print*,"top=",ktop(i), min(k-1, ktf-2),z_cup(i,ktop(i)), z_cup(i,min(k-1, ktf-2))
                !call flush(6)

                ktop(i) = min(k-1, ktf-2)
                exit
              endif
           enddo
       endif
     ENDDO
  END SUBROUTINE rates_up_pdf
!------------------------------------------------------------------------------------
  SUBROUTINE get_zu_zd_pdf(cumulus, draft,ierr,kb,kt,zu,kts,kte,ktf,kpbli,k22,kbcon,klcl,po_cup,psur&
                          ,xland,random)

  implicit none
  integer, intent(in   ) :: kts,kte,ktf,kpbli,k22,kbcon,kt,kb,klcl
  integer, intent(inout) :: ierr
  real   , intent(in   ) :: po_cup(kts:kte),psur,xland,random
  real   , intent(inout) :: zu(kts:kte)
  character*(*), intent(in) ::draft,cumulus
  !- local var
  integer :: kk,add,i,nrec=0,k,kb_adj,kpbli_adj,level_max_zu,ktarget
  real :: zumax,ztop_adj,a2,beta, alpha,kratio,tunning,FZU,krmax,dzudk,hei_updf,hei_down
  real :: zuh(kts:kte),zul(kts:kte),  pmaxzu ! pressure height of max zu for deep
  real,   parameter :: px =45./120. ! px sets the pressure level of max zu. its range is from 1 to 120.
  real,   parameter :: px2=45./120. ! px sets the pressure level of max zu. its range is from 1 to 120.
  integer:: itest                   ! 5=gamma+beta, 4=gamma, 1=beta

  integer:: minzu,maxzul,maxzuh,kstart
  logical :: do_smooth

  !-------- gama pdf
  real, parameter :: beta_deep=1.25,g_beta_deep=0.8974707
  INTEGER :: k1
  real :: lev_start,g_alpha2,g_a,y1,x1,g_b,a,b,alpha2,csum,zubeg,wgty,dp_layer,slope
  real, dimension(30) :: x_alpha,g_alpha
  data  (x_alpha(k),k=1,30)/                                    &
                3.699999,3.699999,3.699999,3.699999,            &
                  3.024999,2.559999,2.249999,2.028571,1.862500, &
                  1.733333,1.630000,1.545454,1.475000,1.415385, &
                  1.364286,1.320000,1.281250,1.247059,1.216667, &
                  1.189474,1.165000,1.142857,1.122727,1.104348, &
                  1.087500,1.075000,1.075000,1.075000,1.075000, &
                1.075000/
  data (g_alpha(k),k=1,30)/                                         &
                4.1706450,4.1706450,4.1706450,4.1706450,            &
                  2.0469250,1.3878370,1.1330030,1.012418,0.9494680, &
                  0.9153771,0.8972442,0.8885444,0.8856795,0.8865333,&
                  0.8897996,0.8946404,0.9005030,0.9070138,0.9139161,&
                  0.9210315,0.9282347,0.9354376,0.9425780,0.9496124,&
                  0.9565111,0.9619183,0.9619183,0.9619183,0.9619183,&
                0.9619183/
  !-------- gama pdf
  DO_SMOOTH = .FALSE.

  !-- fill zu with zeros
  itest=-999
  zu =0.0
  zuh=0.0
  zul=0.0
  IF(zero_diff==1) then
   if(draft == "deep_up" .and. xland >  0.90) itest=11 !ocean
   if(draft == "deep_up" .and. xland <= 0.90) itest=12 !land
   if(draft == "mid_up"                     ) itest= 5
  ELSE
!  if(draft == "deep_up"                    ) itest=21  !ocean/land
   if(draft == "deep_up"                    ) itest=20  !ocean/land
   if(draft == "mid_up"                     ) itest=20
  ENDIF

  !---------------------------------------------------------
  IF(itest==5 .and. draft == "mid_up" ) then

      !--- part 1 GAMMA format
      csum =0.
      zubeg=0.
      lev_start=min(.9,.1+csum*.013)
      kb_adj=max(kb,2)
      kb_adj=min(kb_adj,kt-1)
      if(kb_adj==kt) stop "kb_adj==kt"

      tunning = 0.30
      alpha2  = (tunning*(beta_deep -2.)+1.)/(1.-tunning)

      do k=27,3,-1
            if(x_alpha(k) >= alpha2)exit
      enddo
      k1=k+1
      if(x_alpha(k1) .ne. x_alpha(k1-1))then
        a=x_alpha(k1)-x_alpha(k1-1)
        b=x_alpha(k1-1)*(k1) -(k1-1)*x_alpha(k1)
        x1= (alpha2-b)/a
        y1=a*x1+b
        g_a=g_alpha(k1)-g_alpha(k1-1)
        g_b=g_alpha(k1-1)*k1 - (k1-1)*g_alpha(k1)
        g_alpha2=g_a*x1+g_b
      else
        g_alpha2=g_alpha(k1)
      endif

      fzu = gammaBrams(alpha2 + beta_deep)/(g_alpha2*g_beta_deep)
      fzu=0.01*fzu
      do k=kb_adj,min(kte,kt)
         kratio= (po_cup(k)-po_cup(kb_adj))/(po_cup(kt)-po_cup(kb_adj))
         zu(k) = zubeg+FZU*kratio**(alpha2-1.0) * (1.0-kratio)**(beta_deep-1.0)

      enddo
      !- normalize ZU
      zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-12+maxval(zu(kts:min(kte,kt+1))))

      !--- part 2: BETA format
      pmaxzu=psur-px*(psur-po_cup(kt))
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)
      !beta=4.  !=> must be larger than 1
                !=> higher makes the profile sharper
                !=> around the maximum zu
      !- 2nd approach for beta and alpha parameters
      !- the tunning parameter must be between 0.5 (low  level max zu)
      !-                                   and 1.5 (high level max zu)
      !tunning= 1.0
      tunning = 0.6
      !
      beta    = 2.0/tunning
      alpha   = tunning*beta
      !
      !- this alpha constrains the location of the maximun ZU to be at
      !- "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      !
      ! imposing zu(ktop) = 0
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      !- normalize ZU
      zuh(kts:min(kte,kt+1))= zuh(kts:min(kte,kt+1))/ (1.e-12+maxval(zuh(kts:min(kte,kt+1))))

      !--- part 3: BETA format from sfc to max zuh, then GAMMA format
      Do k=kts,max(kts,maxloc(zuh(:),1)-2)
          zu(k)=zuh(k)
      ENDDO
      Do k=max(kts,maxloc(zuh(:),1)-1),min(maxloc(zuh(:),1)+1,kt)
          zu(k)=0.5*(zu(k)+zuh(k))
      ENDDO

      !-- special treatment below k22/klcl
      DO k=klcl,kts+1,-1
        zu(k)=zu(k+1)*0.5
      enddo
     !-- smooth section
     IF(do_smooth) then
          !--from surface
          zul(kts+1)=zu(kts+1)*0.25
          do k=kts+2,maxloc(zu,1)
             zul(k)=(zu(k-1)+zu(k))*0.5
          enddo
          do k=kts+1,maxloc(zu,1)
             zu(k)=(zul(k)+zu(k))*0.5
          enddo

          !--from ktop
          zul(kt-1)=zu(kt-1)*0.1
          !print*,"ZUMD=",kt,zu(kt),zul(kt)
          do k=kt-2,max( kt-min(maxloc(zu,1),5), kts ),-1
             zul(k)=(zul(k+1)+zu(k))*0.5
          enddo
          wgty=0.
          do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
             wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
             zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
             !print*,"zuMD=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
          enddo
      ENDIF
      zu(kts)=0.
  !---------------------------------------------------------
 !ELSEIF(itest==20 .and. draft == "deep_up") then       !--- land/ocean
  ELSEIF(itest==20                         ) then       !--- land/ocean

      hei_updf=(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN
      !- add a randomic perturbation
      hei_updf = hei_updf + random

      !- for gate soundings
      !hei_updf = max(0.1, min(1.,float(JL)/100.))
      !beta =1.0+float(JL)/100. * 5.

      !--hei_updf parameter goes from 0 to 1 = rainfall decreases with hei_updf
      pmaxzu  =  (psur-100.) * (1.- 0.5*hei_updf) + 0.6*( po_cup(kt) ) * 0.5*hei_updf

      !- beta parameter: must be larger than 1, higher makes the profile sharper around the maximum zu
      beta    = max(1.1, 2.1 - 0.5*hei_updf)

!--- tmp IF(trim(cumulus) == 'deep') beta =beta_sh
      !print*,"hei=",jl,pmaxzu,hei_updf,beta!(pmaxzu-(psur-100.))/( -(psur-100.) +  0.5*( 0.25*(psur-100.) + 0.75*po_cup(kt) ))

      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
      kb_adj=max(kb,kb_adj) ; kb_adj=min(kb_adj,kt)

      !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

      !
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo

      !-- special treatment below k22/klcl
      do k=klcl,kts+1,-1
        zu(k)=zu(k+1)*0.5
      enddo
      !-- smooth section
      IF(do_smooth) then
          !--from surface
          zul(kts+1)=zu(kts+1)*0.25
          do k=kts+2,maxloc(zu,1)
             zul(k)=(zu(k-1)+zu(k))*0.5
          enddo
          do k=kts+1,maxloc(zu,1)
             zu(k)=(zul(k)+zu(k))*0.5
          enddo

          !--from ktop
          zul(kt)=zu(kt)*0.1
          do k=kt-1,max( kt-min(maxloc(zu,1),5),kts) ,-1
             zul(k)=(zul(k+1)+zu(k))*0.5
          enddo
          wgty=0.0
          do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
             wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
             zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
          enddo
      ENDIF

      zu(kts)=0.
  !---------------------------------------------------------
  ELSEIF(itest==21) then

      hei_updf=(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN

      !- for gate soundings
      !if(gate) hei_updf = max(0.1, min(1.,float(JL)/160.))
      !print*,"JL=",jl,hei_updf

      pmaxzu=850.
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"1=",kb_adj,po_cup(kb_adj)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)
      tunning = 0.6
      beta    = 2.0/tunning
      !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zul(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      zul(kts:min(kte,kt))= zul(kts:min(kte,kt))/ (1.e-9+maxval(zul(kts:min(kte,kt)),1))

      !-----------
      pmaxzu=po_cup(kt)+150.
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)
      tunning = 0.8
      beta    = 1.0/tunning

      !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      zuh(kts:min(kte,kt))= zuh(kts:min(kte,kt))/ (1.e-9+maxval(zuh(kts:min(kte,kt)),1))

     !increasing contribuition of zuh => more heating at upper levels/less precip
      zu(:)=(1.-hei_updf)*zul(:) + hei_updf*zuh(:)

      !-- special treatment below k22/klcl
      DO k=klcl,kts+1,-1
        zu(k)=zu(k+1)*0.5
      enddo
     !-- smooth section
     IF(do_smooth) then
          !--from surface
          zul(kts+1)=zu(kts+1)*0.25
          do k=kts+2,maxloc(zu,1)
             zul(k)=(zu(k-1)+zu(k))*0.5
          enddo
          do k=kts+1,maxloc(zu,1)
             zu(k)=(zul(k)+zu(k))*0.5
          enddo

          !--from ktop
          zul(kt)=zu(kt)*0.1
          do k=kt-1,max( kt-min(maxloc(zu,1),5), kts ),-1
             zul(k)=(zul(k+1)+zu(k))*0.5
          enddo

          wgty=0.
          do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
             wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
             zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
             !print*,"zu=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
          enddo
      ENDIF
      zu(kts)=0.
  !---------------------------------------------------------
  ELSEIF(itest==12 .and. draft == "deep_up") then
      !- kb cannot be at 1st level

      if(xland  < 0.90 ) then !- over land
        hei_updf= hei_updf_LAND
      else
        hei_updf= hei_updf_OCEAN
      endif

      !- for gate soundings
      !if(gate) hei_updf = max(0.1, min(1.,float(JL)/100.)) ! for gate soundings

!---non-zero-diff-APR-08-2020
      IF( zero_diff==1) then
        pmaxzu=psur-px*(psur-po_cup(kt))
      ELSE
        pmaxzu=psur-hei_updf*(psur-po_cup(kt))
      ENDIF
!---non-zero-diff-APR-08-2020
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)

      !beta=4.  !=> must be larger than 1
                !=> higher makes the profile sharper
                !=> around the maximum zu
      !- 2nd approach for beta and alpha parameters
      !- the tunning parameter must be between 0.5 (low  level max zu)
      !-                                   and 1.5 (high level max zu)
      tunning = 1.15
      beta    = 2.0/tunning
      alpha   = tunning*beta
      !
      !- this alpha constrains the location of the maximun ZU to be at
      !- "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      !
      ! imposing zu(ktop) = 0
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo

      !-- zu = zero below k22
      !zu(1:max(k22-1,1)) = 0.0
      !-- special treatment below k22/klcl
      DO k=klcl,kts+1,-1
        zu(k)=zu(k+1)*0.5
      enddo
     !-- smooth section
     IF(do_smooth) then
          !--from surface
          zul(kts+1)=zu(kts+1)*0.25
          do k=kts+2,maxloc(zu,1)
             zul(k)=(zu(k-1)+zu(k))*0.5
          enddo
          do k=kts+1,maxloc(zu,1)
             zu(k)=(zul(k)+zu(k))*0.5
          enddo

          !--from ktop
          zul(kt)=zu(kt)*0.1
          do k=kt-1,max( kt-min(maxloc(zu,1),5), kts ),-1
             zul(k)=(zul(k+1)+zu(k))*0.5
          enddo
          wgty=0.0
          do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
             wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
             zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
          enddo
      ENDIF

      zu(kts)=0.
  !---------------------------------------------------------
  ELSEIF(itest==11 .and. draft == "deep_up") then

      if(xland  < 0.90 ) then !- over land
        hei_updf= hei_updf_LAND
      else
        hei_updf= hei_updf_OCEAN
      endif

     !- for gate soundings
      !if(gate) hei_updf = max(0.1, min(1.,float(JL)/100.))
      !print*,"JL=",jl,hei_updf

      pmaxzu=850.
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)!;print*,"1=",kb_adj,po_cup(kb_adj)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)
      tunning = 0.6
      beta    = 2.0/tunning
      !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zul(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
          !print*,"1",k,zul(k),kb_adj,pmaxzu
      enddo
      zul(kts:min(kte,kt))= zul(kts:min(kte,kt))/ (1.e-9+maxval(zul(kts:min(kte,kt)),1))

      !-----------
      pmaxzu=po_cup(kt)+150.
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)
      kb_adj=max(kb,kb_adj)
      kb_adj=min(kb_adj,kt)
      tunning = 0.8
      beta    = 1.0/tunning

      !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))
      do k=klcl-1,min(kte,kt)
          kratio= float(k)/float(kt+1)
          zuh(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      zuh(kts:min(kte,kt))= zuh(kts:min(kte,kt))/ (1.e-9+maxval(zuh(kts:min(kte,kt)),1))

     !increasing contribuition of zuh => more heating at upper levels/less precip
!---non-zero-diff-APR-08-2020
      IF(zero_diff==1) then
         zu(:)= 0.65        *zul(:)+ 0.35     *zuh(:)
      ELSE
         zu(:)=(1.-hei_updf)*zul(:) + hei_updf*zuh(:)
      ENDIF
!---non-zero-diff-APR-08-2020

      !-- special treatment below k22/klcl
      DO k=klcl,kts+1,-1
        zu(k)=zu(k+1)*0.5
      enddo
     !-- smooth section
     IF(do_smooth) then
          !--from surface
          zul(kts+1)=zu(kts+1)*0.25
          do k=kts+2,maxloc(zu,1)
             zul(k)=(zu(k-1)+zu(k))*0.5
          enddo
          do k=kts+1,maxloc(zu,1)
             zu(k)=(zul(k)+zu(k))*0.5
          enddo

          !--from ktop
          zul(kt)=zu(kt)*0.1
          do k=kt-1,max( kt-min(maxloc(zu,1),5), kts ),-1
             zul(k)=(zul(k+1)+zu(k))*0.5
          enddo

          wgty=0.
          do k=kt,max( kt-min(maxloc(zu,1),5), kts ),-1
             wgty=wgty+1./(float(min(maxloc(zu,1),5))+1)
             zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
             !print*,"zu=",k,zu(k),zul(k),(zul(k)+zu(k))*0.5,min(maxloc(zu,1),5),wgty
          enddo
      ENDIF
      zu(kts)=0.
 !---------------------------------------------------------
 ELSEIF(draft == "shallow_up") then
      kb_adj   =kts     ! level where mass flux starts
      kpbli_adj=kpbli
      if(kpbli_adj < kb_adj .or. kpbli_adj >= kt ) then
         kpbli_adj = kb_adj + 1
      endif

      !- location of the maximum Zu: dp_layer mbar above PBL height
      !dp_layer     = 10. !mbar
      !level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(kpbli_adj)-dp_layer)),1)
      !

      k1           = max(kbcon,kpbli_adj)
      !- location of the maximum Zu: dp_layer mbar above k1 height
      hei_updf     =(1.-xland)*hei_updf_LAND+xland*hei_updf_OCEAN

     !hei_updf = (float(JL)-20)/40. ; print*,"JL=",jl,hei_updf

      dp_layer     = hei_updf*(po_cup(k1)-po_cup(kt))

      level_max_zu = minloc(abs(po_cup(kts:kt+1)-(po_cup(k1)-dp_layer)),1)
      level_max_zu = min(level_max_zu,kt -1)
      level_max_zu = max(level_max_zu,kts+1)

      krmax        = float(level_max_zu)/float(kt+1)
      krmax        = min(krmax,0.99)

      beta= beta_sh!smaller => sharper detrainment layer
     !beta= ((1.-xland)*0.43 +xland)*beta_sh

      !beta= 3.0!smaller => sharper detrainment layer
     !beta = 1.+4.*(float(JL))/40.

      !- this alpha imposes the maximum zu at kpbli
      alpha=1.+krmax*(beta-1.)/(1.-krmax)
      !alpha=min(6.,alpha)

      !- to check if dZu/dk = 0 at k=kpbli_adj
      !kratio=krmax
      !dzudk=(alpha-1.)*(kratio**(alpha-2.)) * (1.-kratio)**(beta-1.) - &
      !          (kratio**(alpha-1.))*((1.-kratio)**(beta-2.))*(beta-1.)

      !- Beta PDF
      do k=kts+1,min(kte,kt)
         kratio= float(k)/float(kt+1)
         zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      zu(kts)=0.
      !
      !-- special treatment below kbcon - linear Zu
      if(use_linear_subcl_mf == 1) then
        kstart=kbcon
        slope=(zu(kstart)-zu(kts))/(po_cup(kstart)-po_cup(kts))
        do k=kstart-1,kts+1,-1
           zu(k) = zu(kstart)-slope*(po_cup(kstart)-po_cup(k))
          !print*,"k=",zu(kstart),zu(k),zu(kts)
        enddo
      endif
     !-- special treatment below kclcl
     !do k=(klcl-1),kts+1,-1
     !  zu(k)=zu(k+1)*0.5
     !enddo
     !
     !-- smooth section
     !IF( .not. do_smooth) then
     ! zul(kts+1)=zu(kts+1)*0.1
     ! do k=kts+2,maxloc(zu,1)
     !    zul(k)=(zu(k-1)+zu(k))*0.5
     ! enddo
     ! do k=kts+1,maxloc(zu,1)
     !    zu(k)=(zul(k)+zu(k))*0.5
     ! enddo
     !ENDIF
     !zu(kts)=0.

  !---------------------------------------------------------
  ELSEIF(draft == "DOWN" ) then
      IF(trim(cumulus) == 'shallow') return
      IF(trim(cumulus) == 'mid' ) beta =2.5
      IF(trim(cumulus) == 'deep') beta =2.5

      hei_down=(1.-xland)*hei_down_LAND+xland*hei_down_OCEAN

!---non-zero-diff-APR-08-2020
      IF(zero_diff==1) hei_down= 0.5
!---non-zero-diff-APR-08-2020

      pmaxzu= hei_down * po_cup(kt) + (1.-hei_down)*psur
      kb_adj=minloc(abs(po_cup(kts:kt)-pmaxzu),1)

      !- this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
      alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

      do k=kts+1,min(kt+1,ktf)
              kratio= float(k)/float(kt+1)
              zu(k)= kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
      enddo
      !-- smooth section
      IF(do_smooth) then
        zul(kts+1)=zu(kts+1)*0.1
        wgty=0.
        do k=kts+2,maxloc(zu,1)
           wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
           wgty=0.5
           !print*,"zD1=",k,zu(k),zul(k-1),wgty,zul(k-1)*(1.-wgty)+ zu(k)*wgty
           zul(k)=zul(k-1)*(1.-wgty)+ zu(k)*wgty
        enddo
        wgty=0.
        do k=kts+1,maxloc(zu,1)
           wgty=wgty+1./(float(max(2,maxloc(zu,1)))-1.)
           wgty=0.5
           !print*,"zD2=",k,zu(k),zul(k),wgty,zul(k)*(1.-wgty)+ zu(k)*wgty
           zu(k)=zul(k)*(1.-wgty)+ zu(k)*wgty
        enddo
      ENDIF
      zu(kts)=0.

  ENDIF

  !---------------------------------------------------------

  IF( maxval(zu(kts:min(kte,kt+1)),1) <= 0.0) then
     zu=0.0
     ierr=51 !ierr(i)=51
  ELSE
    !- normalize ZU
    zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ (1.e-9+maxval(zu(kts:min(kte,kt+1)),1))
  ENDIF

 end SUBROUTINE get_zu_zd_pdf
!------------------------------------------------------------------------------------

 SUBROUTINE get_zu_zd_pdf_orig(draft,ierr,kb,kt,zs,zuf,ztop,zu,kts,kte,ktf)

 implicit none
 integer, intent(in) ::kb,kt,kts,kte,ktf
 real, intent(in) :: Zs,Zuf,Ztop
 real, intent(inout) :: zu(kts:kte)
 integer, intent(inout) :: ierr
 character*(*), intent(in) ::draft

 !- local var
 integer :: add,i,nrec=0,k,kb_adj
 real ::zumax,ztop_adj
 real ::beta, alpha,kratio,tunning

 !- kb cannot be at 1st level
 kb_adj=max(kb,2)

 !-- fill zu with zeros
 zu=0.0

 IF(draft == "UP" .or. draft == "up" ) then
  if(kt<=kb_adj) then
    !stop "ktop must be larger than kbcon"
    ierr=99
    return
  endif
  !beta=4.  !=> must larger than 1
            !=> higher makes the profile sharper
            !=> around the maximum zu
  add=0     !=> additional levels above kbcon, where
            !=> the maximum zu will resides
  kb_adj=kb_adj+add
  kb_adj=max(10,kb_adj)

  !- this alpha constrains the location of the maximun ZU to be at
  !- "kb_adj" vertical level
  !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
  tunning = 0.6
  beta    = 2.0/tunning
  alpha   = tunning*beta

   !- Beta PDF
  do k=kts,min(kte,kt+1)
      kratio= float(k)/float(kt+1)

      zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
  enddo

ELSEIF(draft == "DOWN" .or. draft == "DOWNM") then
  add=0    !=> additional levels above kbcon, where
           !=> the maximum zu will resides
  beta=4.  !=> must larger than 1
           !=> higher makes the profile sharper
           !=> around the maximum zu
  alpha= 0.25*beta

  !- for downdrafts kb = jmin(i)-levadj
  kb_adj=kb_adj+add

  !- 2nd approach for beta and alpha parameters
  !- the tunning parameter must be between 0.5 (low  level max zu)
  !-                                   and 1.5 (high level max zu)
  tunning = 1.
  beta    = 2.0/tunning
  alpha   = tunning*beta

  !- Beta PDF
  do k=kts,min(kte,kt)
      kratio= float(k)/float(kt)

      zu(k+1) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
  enddo

ELSEIF(draft == "shallow" .or. draft == "SHALLOW") then

 alpha= 3.
 beta = 2.*alpha
 kb_adj=1 ! level where mass flux starts

 !- Beta PDF
 do k=kts+kb_adj-1,min(kte,kt+1)
    kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)

    zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
 enddo

ELSE
 print*, "unknown type of flow" ,draft
 stop "routine get_zu_zd"

ENDIF

   !- normalize ZU
   zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))

   !--- Sanity checks
   if(beta <= 1) stop "beta must be larger than 1"

   if(minval(zu(:)) < 0.0 ) then
     print*," zu has negative values for ", draft
     stop   " zu < zero"
   endif
   if(maxval(zu(:)) > 1.0 ) then
     print*," zu has values greater than 1 for ", draft
     stop   " zu  >  one"
   endif

   return

!OPEN(19,FILE= 'zu.gra', FORM='unformatted',ACCESS='direct'&
!       ,STATUS='unknown',RECL=4)
! DO k = kts,kte
!    nrec=nrec+1
!     WRITE(19,REC=nrec) zu(k)
! END DO
!close (19)

  END SUBROUTINE get_zu_zd_pdf_orig
!------------------------------------------------------------------------------------
  SUBROUTINE cup_up_cape(aa0,z,zu,dby,GAMMA_CUP,t_cup,  &
              k22,kbcon,ktop,ierr,tempco,qco,qrco, qo_cup,   &
              itf,ktf,its,ite, kts,kte      )

   IMPLICIT NONE
   integer ,intent (in   )                   ::        &
        itf,ktf, its,ite, kts,kte

  ! aa0 = dummy array for CAPE (total cape)
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels
  ! ierr = error value, maybe modified in this routine
  ! tempco = in-cloud temperature (Kelvin) on model cloud levels
  ! qco    = in-cloud water vapor mixing ratio on model cloud levels
  ! qo_cup = environ water vapor mixing ratio on model cloud levels
  ! qrco   = in-cloud liquid water mixing ratio on model cloud levels

   real,    dimension (its:ite,kts:kte) ,intent (in   )    ::       &
        z,zu,gamma_cup,t_cup,dby,tempco,qco,qrco, qo_cup
   integer, dimension (its:ite)         ,intent (in   )    ::       &
        k22,kbcon,ktop
!
! input and output
   integer, dimension (its:ite)         ,intent (inout)    ::       &
        ierr
   real,    dimension (its:ite)         ,intent (out  )    ::       &
        aa0
!
!  local variables in this routine
   integer       ::   i,k
   real          ::   dz,daa0
!
   aa0(:)=0.
   DO i=its,itf
          IF(ierr(i) == 0) then
             DO k=kbcon(i),ktop(i)
              dz=z(i,k)-z(i,max(1,k-1))
              daa0=g*dz*(   (tempco(i,k)*(1.+0.608*qco   (i,k))  - t_cup(i,k)*(1.+0.608*qo_cup(i,k)))&
                          / (t_cup (i,k)*(1.+0.608*qo_cup(i,k))) &
                        )
              aa0(i)=aa0(i)+max(0.,daa0)
              !~ print*,"cape",k,AA0(I),tempco(i,k),t_cup(i,k), qrco  (i,k)
            ENDDO
          ENDIF
   ENDDO
  END SUBROUTINE cup_up_cape
!------------------------------------------------------------------------------------


 SUBROUTINE FLIPZ(flip,mzp)
    implicit none
    integer, intent(In) :: mzp
    integer, dimension(mzp), INtent(inout) :: flip
    integer :: m,k
    m=mzp
    do k=1,mzp
     flip(k)=m
     m=m-1
    Enddo
 END SUBROUTINE FLIPZ
!------------------------------------------------------------------------------------

 SUBROUTINE set_index_loops( ims,ime, jms,jme, kms,kme,    &
                             its,ite, jts,jte, kts,kte,    &
                             mxp,myp,mzp                   )

    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: mxp,myp,mzp
    INTEGER, INTENT(INOUT)      :: ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte


    ims=1  ;ime=mxp ;jms=1  ;jme=myp ;kms=1 ;kme=mzp
    its=1  ;ite=mxp ;jts=1  ;jte=myp ;kts=1 ;kte=mzp

 END SUBROUTINE set_index_loops
!------------------------------------------------------------------------------------

 SUBROUTINE get_vars(LM,mxp,myp,Q,T,PLE,ZLE,ZLO,PLO,PK)
     implicit none
     integer, intent(in) :: LM,mxp,myp
     real, intent(in) , Dimension(mxp,myp,0:LM) :: PLE
     real, intent(in) , Dimension(mxp,myp,LM)   :: T,Q

     real, intent(OUT), Dimension(mxp,myp,0:LM) :: ZLE
     real, intent(OUT), Dimension(mxp,myp,LM)   :: ZLO,PLO,PK

     INTEGER :: L
     real, Dimension(mxp,myp,0:LM) :: PKE

      PLO      = 0.5*(PLE(:,:,0:LM-1) +  PLE(:,:,1:LM  ) )
      PKE      = (PLE/MAPL_P00)**(MAPL_RGAS/MAPL_CP)
      PK       = (PLO/MAPL_P00)**(MAPL_RGAS/MAPL_CP)

      ZLE(:,:,LM) = 0.
      do L=LM,1,-1
         ZLE(:,:,L-1) = (T(:,:,L)/PK(:,:,L)) * (1.+MAPL_VIREPS*Q(:,:,L))
         ZLO(:,:,L  ) = ZLE(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PKE(:,:,L)-PK (:,:,L  ) ) * ZLE(:,:,L-1)
         ZLE(:,:,L-1) = ZLO(:,:,L) + (MAPL_CP/MAPL_GRAV)*( PK (:,:,L)-PKE(:,:,L-1) ) * ZLE(:,:,L-1)
      end do

 END SUBROUTINE get_vars
!------------------------------------------------------------------------------------
 SUBROUTINE get_lcl(t0,pp0,r0,tlcl,plcl,dzlcl)
    implicit none
    real,INTENT(IN ) :: t0,pp0,r0
    real,INTENT(OUT) :: tlcl,plcl,dzlcl

    real, parameter :: &
            cpg      = 102.45             &
        ,   rgas     = 287.               &
        ,   cp       = 1004.              &
        ,   p00      = 1.e5               &
        ,   g        = 9.80               &
        ,   rocp     = rgas / cp          &
        ,   p00i     = 1. / p00           &
        ,   cpor     = cp / rgas          &
        ,   cpi      = 1. / cp            &
        ,   p00k     = 26.870941          &  !  = p00 ** rocp
        ,   p00ki    = 1. / p00k

    integer :: nitt,ip
    real :: p0k,pi0i,ttth0,ttd,dz,pki,pppi,ti,rvs,e
    !real, external :: td,satvap

    !================
    !-simpler,cheaper method
    ttd=td(pp0,r0)
    tlcl=ttd-(0.001296*ttd+0.1963)*(t0-ttd)
    plcl=pp0*(tlcl/t0)**cpor
    dzlcl=127*(t0-ttd)
    if(dzlcl.le.0.)dzlcl=-999.
    !print*,"1st meth",tlcl,plcl,dzlcl;call flush(6)
    RETURN
    !================
    !-2nd method
    dzlcl=-999.
    ip=0
    11 continue

    plcl=pp0
    tlcl=t0
    p0k=pp0**rocp
    pi0i=p0k/p00k*cp
    ttth0=t0*p00k/p0k
    ttd=td(pp0,r0)
    dz=cpg*(t0-ttd)
    if(dz.le.0.)then
       dzlcl=-999.
       return
    endif
    do nitt=1,50
         pki=pi0i-g*dz/(ttth0*(1.+.61*r0))
         pppi=(pki/cp)**cpor*p00
         ti=ttth0*pki/cp
         e=100.*satvap(ti)
         rvs= ( 0.622*e )/ max(1.e-8,(pppi-e))
         !print*,'1',nitt,rvs,r0,ttd,ti,dz
         if(abs(rvs-r0).lt..00003)go to 110
         ttd=td(pppi,r0)
         dz=dz+cp/g*(ti-ttd)
         !print*,'2',nitt,rvs-r0,ttd,ti,dz
    enddo
    print*, 'no converge in LCL:',t0,pp0,r0
    ip=ip+1
    if(ip==1)go to 11
    RETURN

    110 continue
    !- solution for LCL
    plcl=pppi
    tlcl=ti
    dzlcl=dz !displacement
    !print*,"2nd meth",tlcl,plcl,dz
 END SUBROUTINE get_lcl
!------------------------------------------------------------------------------------

 real function td(p,rs)
  implicit none
  real :: rr,rs,es,esln,p
  rr=rs+1e-8
  es=p*rr/(.622+rr)
  esln=log(es)
  td=(35.86*esln-4947.2325)/(esln-23.6837)
  return
 end function td
!------------------------------------------------------------------------------------
SUBROUTINE get_inversion_layers(cumulus,ierr,psur,po_cup,to_cup,zo_cup,k_inv_layers,&
                                 dtempdz,itf,ktf,its,ite, kts,kte)

        IMPLICIT NONE
        INTEGER,                              INTENT (IN )  :: itf,ktf,its,ite,kts,kte
        CHARACTER *(*),                       INTENT (IN )  :: cumulus
        INTEGER, DIMENSION (its:ite),         INTENT (INOUT):: ierr
        REAL,    DIMENSION (its:ite),         INTENT (IN )  :: psur
        REAL,    DIMENSION (its:ite,kts:kte), INTENT (IN )  :: po_cup,to_cup,zo_cup
        REAL,    DIMENSION (its:ite,kts:kte), INTENT (OUT)  :: dtempdz
        INTEGER, DIMENSION (its:ite,kts:kte), INTENT (OUT)  :: k_inv_layers
        REAL :: dzm,delp, first_deriv(kts:kte),sec_deriv(kts:kte),distance(kts:kte)
        INTEGER:: i,k,ilev,kk,k1,ix,k800,k550,ist
        INTEGER, PARAMETER :: extralayer = 0 !- makes plume top higher
        INTEGER, DIMENSION (its:ite,kts:kte)  :: local_k_inv_layers
        !
        !-initialize k_inv_layers as 1 (non-existent layer)_
        k_inv_layers= 1 !integer
        dtempdz     = 0.0
        first_deriv = 0.0
        sec_deriv   = 0.0
        distance    = 0.0
        local_k_inv_layers=1
        ist=3

        DO i = its,itf
         if(ierr(i) /= 0) cycle
         !- displacement from local surface pressure level
         delp=1000.-psur(i)

        !- 2nd method
        ! DO k = kts+1,ktf-2
           !dtempdz(i,k)=   ( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 1,ierr(i)))
           !!! sec_deriv(k)=abs( deriv3(zo_cup(i,k), zo_cup(i,kts:ktf), to_cup(i,kts:ktf), ktf-kts+1, 2))
           !print*,"2=",k,z_cup(i,k),dtempdz(i,k),
        ! ENDDO
        ! if(ierr(i) /= 0) cycle

         !-1st method
         !-  get the 1st derivative
         DO k = kts+ist,ktf-ist
           first_deriv(k)  = (to_cup(i,k+1)-to_cup(i,k-1))/(zo_cup(i,k+1)-zo_cup(i,k-1))
         ENDDO
         first_deriv(kts      :kts+ist-1)  =first_deriv(kts+ist)
         first_deriv(ktf-ist+1:kte      )  =first_deriv(ktf-ist)

         dtempdz  (i,:)  = first_deriv(:)

         !-  get the abs of the 2nd derivative
         DO k = kts+ist+1,ktf-ist-1
           sec_deriv(k)= abs((first_deriv(k+1)-first_deriv(k-1))/(zo_cup(i,k+1)-zo_cup(i,k-1)))
         ENDDO
         sec_deriv(kts    :kts+ist)=sec_deriv(kts+ist+1)
         sec_deriv(ktf-ist:kte    )=sec_deriv(ktf-ist-1)

         ix=1
         DO kk=kts+ist+2,ktf-ist-2
             if(sec_deriv(kk) < sec_deriv(kk+1) .and. sec_deriv(kk) < sec_deriv(kk-1)) then
                local_k_inv_layers(i,ix)=kk
                ix  =ix+1
             endif
         ENDDO

         !- 2nd criteria
         DO k=kts+ist+2,ktf-ist-2
           kk=local_k_inv_layers(i,k)
           if(kk == 1) cycle
           if( dtempdz(i,kk) < dtempdz(i,kk-1) .and. dtempdz(i,kk) < dtempdz(i,kk+1) ) then ! the layer is not a local maximum
               local_k_inv_layers(i,k) = 1
           endif
         ENDDO

        ENDDO


        !- find the locations of inversions around 800 and 550 hPa
        DO i = its,itf
!----------------
!k_inv_layers(i,mid)=1
!----------------
         if(ierr(i) /= 0) cycle
         !- displacement from local surface pressure level
         delp=1000.-psur(i)
!----------------
!k_inv_layers(i,mid)=21
!cycle
!----------------
         IF( trim(cumulus)=='shallow') THEN
            !- now find the closest layers of 800 and 550 hPa.
            !- this is for shallow convection k800
            DO k=kts,ktf
              distance(k)=abs(po_cup(i,local_k_inv_layers(i,k))-(750.-delp))
            ENDDO
            k800=minloc(abs(distance(kts:ktf)),1)

            if( k800 <= kts .or. k800 >= ktf - 4) then
               k_inv_layers(i,shal)= ktf
               !ierr(i)=8
            else
              !-save k800 in the k_inv_layers array
              k_inv_layers(i,shal)=local_k_inv_layers(i,k800) +extralayer
            endif
            !if(  k_inv_layers(i,shal) <= kts .or. k_inv_layers(i,shal) >= ktf-4) then
             !print*,"SHAL_k_inv_layers=",k_inv_layers(i,shal),ierr(i)
             !ierr(i)=11
            !endif

         ELSEIF( trim(cumulus)=='mid') THEN
            !- this is for mid/congestus convection k500
            DO k=kts,ktf
              distance(k)=abs(po_cup(i,local_k_inv_layers(i,k))-(550.-delp))
            ENDDO
            k550=minloc(abs(distance(kts:ktf)),1)

            if( k550 <= kts .or. k550 >= ktf - 4) then
               k_inv_layers(i,mid) = 1
               ierr(i)=8
            else
               !-save k550 in the k_inv_layers array
               k_inv_layers(i,mid )=local_k_inv_layers(i,k550) +extralayer
            endif
            if(  k_inv_layers(i,mid) <= kts .or. k_inv_layers(i,mid) >= ktf-4) then
             !print*,"MID_k_inv_layers=",k_inv_layers(i,MID),ierr(i)
             ierr(i)=12
            endif
         ELSE
             k_inv_layers(i,:)=1
             ierr(i)=88
         ENDIF

        ENDDO

   contains
   real function deriv3(xx, xi, yi, ni, m,ierr)
    !====================================================================
    ! Evaluate first- or second-order derivatives
    ! using three-point Lagrange interpolation
    ! written by: Alex Godunov (October 2009)
    !--------------------------------------------------------------------
    ! input ...
    ! xx    - the abscissa at which the interpolation is to be evaluated
    ! xi()  - the arrays of data abscissas
    ! yi()  - the arrays of data ordinates
    ! ni - size of the arrays xi() and yi()
    ! m  - order of a derivative (1 or 2)
    ! output ...
    ! deriv3  - interpolated value
    !============================================================================*/

    implicit none
    integer, parameter :: n=3
    real   , intent(in):: xx
    integer, intent(in):: ni, m
    real   , intent(in) :: xi(ni), yi(ni)
    real:: x(n), f(n)
    integer i, j, k, ix
    integer, intent(inout) :: ierr

    ! exit if too high-order derivative was needed,
    if (m > 2) then
      deriv3 = 0.0
      return
    end if

    ! if x is ouside the xi(1)-xi(ni) interval set deriv3=0.0
    if (xx < xi(1) .or. xx > xi(ni)) then
      deriv3 = 0.0
      ierr=8
      !stop "problem with 2nd derivative-deriv3 routine"
      return
    endif

    ! a binary (bisectional) search to find i so that xi(i-1) < x < xi(i)
    i = 1
    j = ni
    do while (j > i+1)
      k = (i+j)/2
      if (xx < xi(k)) then
        j = k
      else
        i = k
      endif
    enddo

    ! shift i that will correspond to n-th order of interpolation
    ! the search point will be in the middle in x_i, x_i+1, x_i+2 ...
      i = i + 1 - n/2

    ! check boundaries: if i is ouside of the range [1, ... n] -> shift i
    if (i < 1) i=1
    if (i + n > ni) i=ni-n+1

    !  old output to test i
    !  write(*,100) xx, i
    !  100 format (f10.5, I5)

    ! just wanted to use index i
    ix = i

    ! initialization of f(n) and x(n)
    do i=1,n
      f(i) = yi(ix+i-1)
      x(i) = xi(ix+i-1)
    end do

    ! calculate the first-order derivative using Lagrange interpolation
    if (m == 1) then
        deriv3 =          (2.0*xx - (x(2)+x(3)))*f(1)/((x(1)-x(2))*(x(1)-x(3)))
        deriv3 = deriv3 + (2.0*xx - (x(1)+x(3)))*f(2)/((x(2)-x(1))*(x(2)-x(3)))
        deriv3 = deriv3 + (2.0*xx - (x(1)+x(2)))*f(3)/((x(3)-x(1))*(x(3)-x(2)))
    ! calculate the second-order derivative using Lagrange interpolation
    else
        deriv3 =          2.0*f(1)/((x(1)-x(2))*(x(1)-x(3)))
        deriv3 = deriv3 + 2.0*f(2)/((x(2)-x(1))*(x(2)-x(3)))
        deriv3 = deriv3 + 2.0*f(3)/((x(3)-x(1))*(x(3)-x(2)))
    endif
  end function deriv3
 END SUBROUTINE get_inversion_layers
!------------------------------------------------------------------------------------

 SUBROUTINE alloc_grads_arr(n,mzp,task,jl)
     implicit none
     integer, intent(in)    :: n,mzp,task
     integer, intent(inout) :: jl
     integer :: nvar

     if(task == 1) then
        jl = n
        allocate (cupout(nvar_grads))
        do nvar=1,nvar_grads
         allocate(cupout(nvar)%varp(n,mzp))
         allocate(cupout(nvar)%varn(3))
         cupout(nvar)%varp(:,:)=0.0
         cupout(nvar)%varn(:)  ="xxxx"
        enddo
     else
        do nvar=1,nvar_grads
         deallocate(cupout(nvar)%varp)
         deallocate(cupout(nvar)%varn)
        enddo
        deallocate(cupout)
     endif

 END SUBROUTINE alloc_grads_arr
!------------------------------------------------------------------------------------

 SUBROUTINE set_grads_var(i,k,nvar,f,name1,name2,name3)
     implicit none
     integer, intent(in)    :: i,k
     integer, intent(inOUT) :: nvar
     real, intent(in) :: f
     character*(*), intent(in) :: name1,name2,name3

     cupout(nvar)%varp(i,k)= f
     cupout(nvar)%varn(1)=name1
     cupout(nvar)%varn(2)=name2
     cupout(nvar)%varn(3)=name3
     nvar=nvar+1
     if(nvar>nvar_grads) stop 'nvar>nvar_grads'

 END SUBROUTINE set_grads_var
!------------------------------------------------------------------------------------

 SUBROUTINE wrt_bin_ctl(n,mzp,p2d,cumulus)
   implicit none
   integer, intent(IN):: n,mzp
   character*(*), intent(in) :: cumulus
   real, dimension(mzp),intent(in):: p2d
   integer:: nvartotal,klevgrads(200),jk,int_byte_size,nvar,maxklevgrads
   real   :: real_byte_size
   real, parameter :: undef=-9.99e33
   integer :: nrec=0
   integer :: recSize

   maxklevgrads=min(60,mzp)
   runname='15geos5_'//cumulus
   runlabel=runname

   print*,"writing grads control file:',trim(runname)//'.ctl",ntimes;call flush(6)
   !
   !number of variables to be written
   nvartotal=0
   do nvar=1,nvar_grads
     if(cupout(nvar)%varn(1) .ne. "xxxx") nvartotal=nvartotal+1
     if(cupout(nvar)%varn(3)  ==  "3d"  ) klevgrads(nvar)=maxklevgrads
     if(cupout(nvar)%varn(3)  ==  "2d"  ) klevgrads(nvar)=1
   enddo

   !- binary file
   inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list

   print*, 'opening grads file:',trim(runname)//'.gra'
   recSize=size(cupout(nvar)%varp,1)*real_byte_size
   if(ntimes == 1) then
    open(19,file= trim(runname)//'.gra',form='unformatted',&
            access='direct',status='replace',recl=recSize)
   else
    open(19,file= trim(runname)//'.gra',form='unformatted',&
            access='direct',status='old', recl=recSize)
   endif

   do nvar=1,nvar_grads
       if(cupout(nvar)%varn(1) .ne. "xxxx") then
        do jk=1,klevgrads(nvar)
          nrec=nrec+1
          !write(19)          real((cupout(nvar)%varp(:,jk)),4)
          write(19,rec=nrec)  real((cupout(nvar)%varp(:,jk)),4)
        enddo
       endif
   enddo
   close (19)
   !-setting vertical dimension '0' for 2d var
   where(klevgrads==1)klevgrads=0
   !- ctl file
   open(20,file=trim(runname)//'.ctl',status='unknown')
   write(20,2001) '^'//trim(runname)//'.gra'
   write(20,2002) 'undef -9.99e33'
   write(20,2002) 'options sequential byteswapped' ! zrev'
   write(20,2002) 'title '//trim(runlabel)
   write(20,2003) 1,0.,1. ! units m/km
   write(20,2004) n,1.,1.
   write(20,2005) maxklevgrads,(p2d(jk),jk=1,maxklevgrads)
   write(20,2006) ntimes,'00:00Z24JAN1999','10mn'
   write(20,2007) nvartotal
   do nvar=1,nvar_grads
    if(cupout(nvar)%varn(1) .ne. "xxxx") then
     !
     write(20,2008) cupout(nvar)%varn(1)(1:len_trim(cupout(nvar)%varn(1))),klevgrads(nvar)&
                   ,cupout(nvar)%varn(2)(1:len_trim(cupout(nvar)%varn(2)))
    endif
   enddo
   write(20,2002) 'endvars'
   close(20)

  2001 format('dset ',a)
  2002 format(a)
  2003 format('xdef ',i4,' linear ',2f15.3)
  2004 format('ydef ',i4,' linear ',2f15.3)
  2005 format('zdef ',i4,' levels ',60f8.3)
  2006 format('tdef ',i4,' linear ',2a15)
  2007 format('vars ',i4)
  2008 format(a10,i4,' 99 ',a40)!'[',a8,']')
  2055 format(60f7.0)
   133 format (1x,F7.0)

 END SUBROUTINE wrt_bin_ctl
!------------------------------------------------------------------------------------

 SUBROUTINE writetxt(mzp,t,ple,th1,pk,q1,u1,zle,zlo      &
                      ,DYNF_Q ,DYNF_T ,DYNF_PLE, DYNF_UA  &
                      !
                      ,theta,pp,rv,up,zm3d,zt3d,vp,omega&
                          ,sflux_r,sflux_t,topt,xland,sfc_press,dx,kpbl,temp2m,dt_moist&
                      )

    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: mzp
    real, dimension(mzp)   :: t ,th1 ,pk ,q1 ,u1 ,zlo &
                             ,DYNF_Q ,DYNF_T , DYNF_UA

    real, dimension(mzp)   :: theta ,pp ,rv ,up ,zm3d ,zt3d,vp,omega
    real, dimension(0:mzp) :: ple ,zle ,DYNF_PLE
    REAL :: sflux_r,sflux_t,topt,xland,sfc_press,dx,temp2m,dt_moist

    integer :: k,kpbl

    write(8,*) "================================================"
    write(7,*) "================================================"
    write(7,*) kpbl,sflux_r,sflux_t,topt,xland,sfc_press,dx,temp2m,dt_moist
    do k=1,mzp
     write(8,10) k,ple(k),t(k),th1(k),pk(k),1000.*q1(k),u1(k),zle(k),zlo(k),86400.*DYNF_Q(k)
     write(7,11) k,theta(k),pp(k),1000.*rv(k),up(k),zm3d(k),zt3d(k),vp(k),omega(k)
    enddo
    call flush(7)
    call flush(8)
 10 FORMAT(1x,i4,9F11.3)
 11 Format(1x,i4,8F11.3)
 END SUBROUTINE writetxt
!------------------------------------------------------------------------------------

 SUBROUTINE get_cloud_fraction(  mzp, kts, ktf, &
                PPABS, PZZ, PT, PRV, QCO, QRCO, PMFLX, PCLDFR ,PRC, PRI )
    !!    PURPOSE
    !!    -------
    !!**  Routine to diagnose cloud fraction and liquid and ice condensate mixing ratios
    !!**  METHOD
    !!    ------
    !!    Based on the large-scale fields of temperature, water vapor, and possibly
    !!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed
    !!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
    !!
    !!    The total variance is parameterized as the sum of  stratiform/turbulent variance
    !!    and a convective variance.
    !!    The turbulent variance is parameterized as a function of first-order moments, and
    !!    the convective variance is modelled as a function of the convective mass flux (units kg/s m^2)
    !!    as provided by a mass flux convection scheme.
    !!
    !!    Nota: if the host model does not use prognostic values for liquid and solid condensate
    !!    or does not provide a convective mass flux, put all these values to zero.
    !!    Also, it is supposed that vertical model levels are numbered from
    !!    1 to MZP, where 1 is the first model level above the surface
    !!
    !!    ------------------
    !!    REFERENCE
    !!    ---------
    !!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
    !!      Chaboureau J.P. and P. Bechtold (JGR/AGU 2005)
    !!
    !!    AUTHOR
    !!    ------
    !!      P. BECHTOLD       * Laboratoire d'Aerologie *
    !!
    !!    MODIFICATIONS
    !!    -------------
    !!      Original    13/06/2001
    !!      modified    20/03/2002 : add convective Sigma_s and improve turbulent
    !!                               length-scale in boundary-layer and near tropopause
    !!      adapted     09/12/2016 : adapted to GEOS-5 by Saulo Freitas
    !-------------------------------------------------------------------------------
    !*       0.    DECLARATIONS
    !              ------------
    IMPLICIT NONE
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    Set the fundamental thermodynamical constants
    !              these have the same values (not names) as in ARPEGE IFS
    !              -------------------------------------------------------
    real, parameter :: XP00   = 1.E5        ! reference pressure
    real, parameter :: XPI    = 3.141592654 ! Pi
    real, parameter ::  XG    = 9.80665     ! gravity constant
    real, parameter :: XMD    = 28.9644E-3  ! molecular weight of dry air
    real, parameter :: XMV    = 18.0153E-3  ! molecular weight of water vapor
    real, parameter :: XRD    = 287.05967   ! gaz constant for dry air
    real, parameter :: XRV    = 461.524993  ! gaz constant for water vapor
    real, parameter :: XCPD   = 1004.708845 ! specific heat of dry air
    real, parameter :: XCPV   = 1846.1      ! specific heat of water vapor
    real, parameter :: XRHOLW = 1000.       ! density of liquid water
    real, parameter :: XCL    = 4218.       ! specific heat of liquid water
    real, parameter :: XCI    = 2106.       ! specific heat of ice
    real, parameter :: XTT    = 273.16      ! triple point temperature
    real, parameter :: XLVTT  = 2.5008E6    ! latent heat of vaporisation at XTT
    real, parameter :: XLSTT  = 2.8345E6    ! latent heat of sublimation at XTT
    real, parameter :: XLMTT  = 0.3337E6    ! latent heat of melting at XTT
    real, parameter :: XESTT  = 611.14      ! saturation pressure at XTT
    real, parameter :: XALPW  = 60.22416    ! constants in saturation pressure over liquid water
    real, parameter :: XBETAW = 6822.459384
    real, parameter :: XGAMW  = 5.13948
    real, parameter :: XALPI  = 32.62116    ! constants in saturation pressure over ice
    real, parameter :: XBETAI = 6295.421
    real, parameter :: XGAMI  = 0.56313
    LOGICAL, PARAMETER :: LUSERI = .TRUE. ! logical switch to compute both
                                              ! liquid and solid condensate (LUSERI=.TRUE.)
                                              ! or only liquid condensate (LUSERI=.FALSE.)
    !
    !*       0.1   Declarations of dummy arguments :
    !
    !
    INTEGER,              INTENT(IN)   :: mzp     ! vertical dimension
    INTEGER,              INTENT(IN)   :: kts     ! vertical  computations start at
    !                                             ! KTS that is at least 1
    INTEGER,              INTENT(IN)   :: ktf     ! vertical computations can be
                                                  ! limited to MZP + 1 - KTF
                                                  ! default=1
    REAL, DIMENSION(mzp), INTENT(IN)    :: PPABS  ! pressure (Pa)
    REAL, DIMENSION(mzp), INTENT(IN)    :: PZZ    ! height of model levels (m)
    REAL, DIMENSION(mzp), INTENT(IN)    :: PT     ! grid scale T  (K)
    REAL, DIMENSION(mzp), INTENT(IN)    :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
    REAL, DIMENSION(mzp), INTENT(IN)    :: PMFLX  ! convective mass flux (kg/(s m^2))
    REAL, DIMENSION(mzp), INTENT(IN)    :: QRCO   ! sub-grid scale liq water mixing ratio (kg/kg)
    REAL, DIMENSION(mzp), INTENT(IN)    :: QCO    ! in-cloud water mixing ratio (kg/kg)
    REAL, DIMENSION(mzp), INTENT(INOUT),OPTIONAL :: PRC    ! grid scale r_c mixing ratio (kg/kg)
    REAL, DIMENSION(mzp), INTENT(INOUT),OPTIONAL :: PRI    ! grid scale r_i (kg/kg)
    REAL, DIMENSION(mzp), INTENT(OUT)   :: PCLDFR ! fractional cloudiness (between 0 and 1)
    !
    !
    !*       0.2   Declarations of local variables :
    !
    INTEGER  ::  JKT, JKP, JKM,K     ! loop index
    REAL, DIMENSION(mzp) :: ZTLK, ZRT       ! work arrays for T_l, r_t
    REAL, DIMENSION(mzp) :: ZL              ! length-scale
    INTEGER   :: ITPL    ! top levels of tropopause/highest inversion
    REAL      :: ZTMIN   ! min Temp. related to ITPL
    REAL, DIMENSION(mzp) :: LOC_PRC,LOC_PRI
    !
    REAL :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
    REAL :: ZLL, DZZ, ZZZ ! length scales
    REAL :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
    REAL :: ZSIG_CONV,  ZSIGMA_NOCONV,  ZQ1_NOCONV      ! convective part of Sig_s
    !
    !*       0.3  Definition of constants :
    !
    !-------------------------------------------------------------------------------
    !
    REAL :: ZL0     = 600.        ! tropospheric length scale
                                  ! changed to 600 m instead of 900 m to give a consistent
                                  ! value (linear increase) in general 500 m deep oceanic
                                  ! mixed layer - but could be put back to 900 m if wished
    REAL :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
    REAL :: ZCSIG_CONV = 0.30E-2  ! scaling factor for ZSIG_CONV as function of mass flux
    !
    !
    LOGICAL :: ONLY_CONVECTIVE_CLOUD_FRACTION=.TRUE. ! set .false. for the total cloud fraction
    !-------------------------------------------------------------------------------
    !RETURN
    !
    IF(PRESENT(PRC)) THEN
       LOC_PRC(:)=PRC(:)
    ELSE
       LOC_PRC(:)=0.0
    ENDIF
    IF(PRESENT(PRI)) THEN
       LOC_PRI(:)=PRI(:)
    ELSE
       LOC_PRI(:)=0.0
    ENDIF

    PCLDFR(:) = 0. ! Initialize values
    !

    JKT = MZP+1-KTS
    !-will limit the model vertical column to 60 hPa
    DO K=KTF,KTS,-1
       if(PPABS(k) > 60.*100.) then
          JKT = k
          !PRINT*,"JKT=",K,MZP+1-KTS ;CALL FLUSH(6)
          exit
       endif
    ENDDo

    DO K=KTS,JKT
       ZTEMP  = PT(k)
        !latent heat of vaporisation/sublimation
       ZLV    = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
       ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

       !store temperature at saturation and total water mixing ratio
       ZRT(k)   = PRV(k) + LOC_PRC(k) + LOC_PRI(k)
       ZCPD     = XCPD  + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
       ZTLK(k)  = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD
    END DO

    !-------------------------------------------------------------------------------
    ! Determine tropopause/inversion  height from minimum temperature

    ITPL  = KTS+1
    ZTMIN = 400.
    DO k = KTS+1,JKT-1
             IF ( PT(k) < ZTMIN ) THEN
                  ZTMIN = PT(k)
                  ITPL  = K
             ENDIF
    END DO

    ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

    ZL(:) = 20.
    DO k = KTS+1,JKT

       ! free troposphere
       ZL(k) = ZL0
       JKP   = ITPL
       ZZZ   = PZZ(k) -  PZZ(KTS)
          ! approximate length for boundary-layer : linear increase
       IF ( ZL0 > ZZZ )  ZL(k) = ZZZ
          ! gradual decrease of length-scale near and above tropopause/top inversion
       IF ( ZZZ > 0.9*(PZZ(JKP)-PZZ(KTS)) ) &
            ZL(k) = .6 * ZL(K-1)
    END DO
    !-------------------------------------------------------------------------------

    DO k=KTS+1,JKT-1
       JKP=k+1
       JKM=k-1
       ZTEMP  = PT(k)

       !latent heat of vaporisation/sublimation
       ZLV    = XLVTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
       ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

       ZCPD   = XCPD + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
       !temperature at saturation
       ZTL    = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD

       !saturated water vapor mixing ratio over liquid water
       ZPV    = MIN(EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) ),0.99*PPABS(k))
       ZQSL   = XRD / XRV * ZPV / ( PPABS(k) - ZPV )

        !saturated water vapor mixing ratio over ice
       ZPIV   = MIN(EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) ),0.99*PPABS(k))
       ZQSI   = XRD / XRV * ZPIV / ( PPABS(k) - ZPIV )

       !interpolate between liquid and solid as function of temperature
       ! glaciation interval is specified here to 20 K
       ZFRAC = ( ZTL  - 250.16 ) / ( XTT - 250.16 )  ! liquid/solid fraction
       ZFRAC = MAX( 0., MIN(1., ZFRAC ) )

       IF(.NOT. LUSERI) ZFRAC=1.
       ZQSL = ( 1. - ZFRAC ) * ZQSI + ZFRAC * ZQSL
       ZLV  = ( 1. - ZFRAC ) * ZLS  + ZFRAC * ZLV

       !coefficients a and b
       ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 ) * (XRV * ZQSL / XRD + 1.)
       !orig  ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )

       ZA   = 1. / ( 1. + ZLV/ZCPD * ZAH )
       ZB   = ZAH * ZA

       !- parameterize Sigma_s with first_order closure
       DZZ    =  PZZ (JKP)  - PZZ(JKM)
       ZDRW   =  ZRT (JKP)  - ZRT(JKM)
       ZDTL   =  ZTLK(JKP) - ZTLK(JKM) + XG/ZCPD * DZZ
       ZLL    =  ZL(k)

       !- standard deviation due to convection
       ZSIG_CONV = ZCSIG_CONV * PMFLX(k) / ZA

       !- turb + conv
       ZSIGMA = SQRT( MAX( 1.E-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
                                   ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
                                 + ZB*ZB*ZDTL*ZDTL                      ) &
                                 + ZSIG_CONV * ZSIG_CONV ) )

       !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
       ZSIGMA = MAX( ZSIGMA, 1.E-10 )

       !- normalized saturation deficit
       ZSBAR = ZA * ( ZRT (k) - ZQSL )
       !- "Q1" parameter
       ZQ1   = ZSBAR / ZSIGMA

       !- total cloud fraction
       PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

       IF(ONLY_CONVECTIVE_CLOUD_FRACTION) THEN
          !- get cloud fraction associated with ONLY the sub-grid scale convective part
          !- this sigma does not include the sub-grid scale convective part
          ZSIGMA_NOCONV = SQRT( MAX( 1.E-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
                                             ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
                                           + ZB*ZB*ZDTL*ZDTL  )))
          !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
          ZSIGMA_NOCONV = MAX( ZSIGMA_NOCONV, 1.E-10 )
          ZQ1_NOCONV = ZSBAR / ZSIGMA_NOCONV

          !- cloud fraction associated with ONLY convective part ("total-turb")
          PCLDFR(k) = 0.36*(ATAN(1.55*ZQ1)-ATAN(1.55*ZQ1_NOCONV))

          PCLDFR(k) = MAX( 0., MIN(1.,PCLDFR(k)) )

       ENDIF
       !- newer formulation, see GMD 2015
       !PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.34*ATAN(1.85*ZQ1+2.33)) )
       !- this is area fraction of cloud cores
       !PCLDFR(k) = MAX( 0., MIN(1.,0.292/ZQ1**2) )

       CYCLE
       !total condensate diagnostic (not being used)
       IF (ZQ1 > 0. .AND. ZQ1 <= 2. ) THEN
        !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
        ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
       ELSE IF (ZQ1 > 2.) THEN
          ZCOND = ZQ1
       ELSE
          ZCOND = EXP( 1.2*ZQ1-1. )
       END IF
       ZCOND = ZCOND * ZSIGMA

       if ( zcond < 1.e-12) then
           zcond = 0.
           pcldfr(k) = 0.
       end if
        if ( pcldfr(k) == 0.) then
           zcond = 0.
       end if

       LOC_PRC(k) = ZFRAC * ZCOND ! liquid condensate
       IF (LUSERI) THEN
          LOC_PRI(k) = (1.-ZFRAC) * ZCOND   ! solid condensate
       END IF

!---
! compute s'rl'/Sigs^2
! used in w'rl'= w's' * s'rl'/Sigs^2
!  PSIGRC(k) = PCLDFR(k)   ! Gaussian relation
!
! s r_c/ sig_s^2
!    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
!
!    multiply PSRCS by the lambda3 coefficient
!
!      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
! in the 3D case lambda_3 = 1.
!      INQ1 = MIN( MAX(-22,FLOOR(2*ZQ1) ), 10)
!      ZINC = 2.*ZQ1 - INQ1
!
!      PSIGRC(k) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
!
!      PSIGRC(k) = PSIGRC(k)* MIN( 3. , MAX(1.,1.-ZQ1) )
!---
    END DO
    !
  END SUBROUTINE get_cloud_fraction
!------------------------------------------------------------------------------------
  SUBROUTINE cup_cloud_limits(name,ierrc,ierr,cap_inc,cap_max_in                                   &
                             ,heo_cup,heso_cup,qo_cup,qeso_cup,po,po_cup,z_cup,heo,hkbo,qo,qeso    &
                             ,entr_rate,hcot,k22,kbmax,klcl,kbcon,ktop,depth_neg_buoy,frh,Tpert &
                             ,start_level_,use_excess,zqexec,ztexec, x_add_buoy,xland              &
                             ,itf,ktf,its,ite, kts,kte)

     IMPLICIT NONE
     character *(*), intent (in)         ::     name

     integer ,intent (in   )             ::    &
        itf,ktf,its,ite, kts,kte,use_excess

     real,    dimension (its:ite,kts:kte)  ,intent (in   )     ::    &
        heo_cup,heso_cup,po_cup,z_cup,heo,qo_cup,qeso_cup,po,qo,qeso,Tpert
     real,    dimension (its:ite)          ,intent (in   )     ::    &
        cap_max_in,cap_inc,xland
     real,    dimension (its:ite)          ,intent (in   )     ::    &
        zqexec,ztexec,x_add_buoy
     integer, dimension (its:ite)          ,intent (in   )     ::    &
        kbmax,start_level_
     integer, dimension (its:ite)          ,intent (inout)     ::    &
        kbcon,ierr,ktop,klcl,k22
     character*128                    ,intent (inout) :: ierrc(its:ite)
     real, dimension (its:ite)        ,intent (inout) :: hkbo,depth_neg_buoy,frh
     real, dimension (its:ite,kts:kte),intent (in)    :: entr_rate
     real, dimension (its:ite,kts:kte),intent (inout) :: hcot

!  local variables in this routine

     real, parameter              :: frh_crit_O=0.7
     real, parameter              :: frh_crit_L=0.7  !--- test 0.5
     real                         :: delz_oversh !--- height of cloud overshoot is 10% higher than the LNB.
                                                       !--- Typically it can 2 - 2.5km higher, but it depends on
                                                       !--- the severity of the thunderstorm.

     real,    dimension (its:ite) ::   cap_max
     integer                      ::   i,k,k1,k2,kfinalzu
     real                         ::   plus,hetest,dz,dbythresh,denom &
                                      ,dzh,del_cap_max,fx,x_add,Z_overshoot,frh_crit
     real   , dimension (kts:kte) ::   dby
     integer, dimension (its:ite) ::   start_level

     delz_oversh = OVERSHOOT
     hcot = 0.0
     dby  = 0.0
     start_level = 0
     cap_max(:)  = cap_max_in(:)

      DO i=its,itf
        if(ierr(i) /= 0) cycle
        if(ZERO_DIFF==1) then
           start_level(i) = klcl(i)
        else
           start_level(i) = start_level_(i)
        endif

        do k=kts,start_level(i)
           hcot(i,k) = hkbo(i) ! assumed no entraiment between these layers
        enddo
      ENDDO
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
loop0: DO i=its,itf
        !-default value
        kbcon         (i)=kbmax(i)+3
        depth_neg_buoy(i)=0.
        frh           (i)=0.
        if(ierr(i) /= 0) cycle


loop1:  DO WHILE(ierr(i) == 0)

            kbcon(i)=start_level(i)
            do k=start_level(i)+1,KBMAX(i)+3
                dz=z_cup(i,k)-z_cup(i,k-1)
                hcot(i,k)= ( (1.-0.5*entr_rate(i,k-1)*dz)*hcot(i,k-1)     &
                                   + entr_rate(i,k-1)*dz *heo (i,k-1) ) / &
                             (1.+0.5*entr_rate(i,k-1)*dz)
                if(k==start_level(i)+1) then
                  x_add    = (xlv*zqexec(i)+cp*ztexec(i)) + x_add_buoy(i)
                  hcot(i,k)= hcot(i,k) +  x_add
                endif
            enddo

loop2:      do while (hcot(i,kbcon(i)) < HESO_cup(i,kbcon(i)))
                kbcon(i)=kbcon(i)+1
                if(kbcon(i).gt.kbmax(i)+2) then
                    ierr(i)=3
                    ierrc(i)="could not find reasonable kbcon in cup_kbcon : above kbmax+2 "
                    exit loop2
                endif
                !print*,"kbcon=",kbcon(i);call flush(6)
            enddo loop2

            IF(ierr(i) /= 0) cycle loop0

            !---     cloud base pressure and max moist static energy pressure
            !---     i.e., the depth (in mb) of the layer of negative buoyancy
            depth_neg_buoy(i) = - (po_cup(i,kbcon(i))-po_cup(i,start_level(i)))

            IF(MOIST_TRIGGER == 1) THEN
                frh(i)=0. ; dzh = 0
                do k=k22(i),kbcon(i)
                   dz     = z_cup(i,k)-z_cup(i,max(k-1,kts))
                   frh(i) = frh(i) + dz*(qo(i,k)/qeso(i,k))
                   dzh    = dzh + dz
                   !print*,"frh=", k,dz,qo(i,k)/qeso(i,k)
                enddo
                frh(i) = frh(i)/(dzh+1.e-16)
                frh_crit =frh_crit_O*xland(i) + frh_crit_L*(1.-xland(i))

               !fx     = 2.*(frh(i) - frh_crit) !- linear
               !fx     = 4.*(frh(i) - frh_crit)* abs(frh(i) - frh_crit) !-quadratic
                fx     = ((2./0.78)*exp(-(frh(i) - frh_crit)**2)*(frh(i) - frh_crit)) !- exponential
                fx     = max(-1.,min(1.,fx))

                del_cap_max = fx* cap_inc(i)
                cap_max(i)  = min(max(cap_max_in(i) + del_cap_max, 10.),150.)
               !print*,"frh=", frh(i),kbcon(i),del_cap_max, cap_max(i),  cap_max_in(i)
            ENDIF

            !- test if the air parcel has enough energy to reach the positive buoyant region
            if(cap_max(i) > depth_neg_buoy(i)) cycle loop0


!--- use this for just one search (original k22)
!            if(cap_max(i) < depth_neg_buoy(i)) then
!                    ierr(i)=3
!                    ierrc(i)="could not find reasonable kbcon in cup_cloud_limits"
!            endif
!            cycle loop0
!---

            !- if am here -> kbcon not found for air parcels from k22 level
            k22(i)=k22(i)+1
            !--- increase capmax
            IF(USE_MEMORY == 20) cap_max(i)=cap_max(i)+cap_inc(i)

            !- get new hkbo
            x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
            call get_cloud_bc(name,kts,kte,ktf,xland(i),po(i,kts:kte),heo_cup (i,kts:kte),hkbo (i),k22(i),x_add,Tpert(i,kts:kte))
            !
            start_level(i)=start_level(i)+1
            !
            hcot(i,start_level(i))=hkbo (i)
        ENDDO loop1
        !--- last check for kbcon
        if(kbcon(i) == kts) then
            ierr(i)=33
            ierrc(i)="could not find reasonable kbcon in cup_kbcon = kts"
        endif
     ENDDO loop0

!
!      if(NAME /= 'shallow') return
!
!--- DETERMINE THE LEVEL OF NEUTRAL BUOYANCY - KTOP
!
     DO i=its,itf
         ktop(i) = ktf-1
         IF(ierr(i) /= 0) cycle
         !~ dby(:)=0.0

         start_level(i)=kbcon(i)

         do k=start_level(i)+1,ktf-1
           dz=z_cup(i,k)-z_cup(i,k-1)
           denom = 1.+0.5*entr_rate(i,k-1)*dz
           if(denom == 0.) then
                hcot(i,k)=hcot(i,k-1)
           else
                hcot(i,k)=( (1.-0.5*entr_rate(i,k-1)*dz)*hcot(i,k-1)    &
                                   +entr_rate(i,k-1)*dz *heo (i,k-1) )/ denom
           endif
         enddo
         do k=start_level(i)+1,ktf-1

             if(hcot(i,k) < heso_cup(i,k) )then
               ktop(i)  =  k - 1
               exit
            endif
         enddo
         if(ktop(i).le.kbcon(i)+1) ierr(i)=41

         !----------------
         if(OVERSHOOT > 1.e-6 .and. ierr(i) == 0) then
           Z_overshoot = (1. + delz_oversh) * z_cup(i,ktop(i))
           do k=ktop(i),ktf-2
              if(Z_overshoot < z_cup(i,k)) then
                ktop(i) = min(k-1, ktf-2)
                exit
              endif
           enddo
         endif
     ENDDO
  END SUBROUTINE cup_cloud_limits
!------------------------------------------------------------------------------------
     subroutine get_buoyancy(itf,ktf, its,ite, kts,kte,ierr,klcl,kbcon,ktop &
                            ,hc,he_cup,hes_cup,dby,z_cup)

        IMPLICIT NONE
        integer, intent (in   )   :: itf,ktf,its,ite, kts,kte
        integer, dimension (its:ite)           ,intent (in)     ::      &
                                  ierr,klcl,kbcon,ktop
        real,    dimension (its:ite,kts:kte)   ,intent (in   )  ::      &
                                  hc,he_cup,hes_cup,z_cup
        real,    dimension (its:ite,kts:kte)   ,intent (OUT  )  ::      &
                                  dby
        INTEGER :: i,k

        do i=its,itf
          dby (i,:)=0.
          if(ierr(i) /= 0) cycle
          do k=kts,klcl(i)
            dby (i,k)=hc (i,k)-he_cup (i,k)
          enddo
          DO  k=klcl(i)+1,ktop(i)+1
            dby (i,k)=hc (i,k)-hes_cup (i,k)
          ENDDO
        ENDDO
     end subroutine get_buoyancy

!------------------------------------------------------------------------------------
     subroutine cup_up_vvel(vvel2d,vvel1d,zws,entr_rate,cd ,z,z_cup,zu,dby,GAMMA_CUP,t_cup &
                           ,tempco,qco,qrco,qo,klcl,kbcon,ktop &
                           ,ierr,itf,ktf,its,ite, kts,kte)

     implicit none
     real, parameter :: ctea=1./3. ,cteb=2., visc=2000., eps=0.622
     integer,intent (in   )              ::  itf,ktf,its,ite, kts,kte
     real,    dimension (its:ite,kts:kte) ,intent (in   )  ::  &
        z,z_cup,zu,gamma_cup,t_cup,dby,entr_rate,cd,tempco,qco,qrco,qo

     integer, dimension (its:ite)         ,intent (in   )  ::  &
        klcl,kbcon,ktop
     real,    dimension (its:ite)         ,intent (in   )  ::  &
        zws

  ! input and output
     integer, dimension (its:ite)        ,intent (inout) :: ierr
     real,    dimension (its:ite,kts:kte),intent (out  ) ::  vvel2d
     real,    dimension (its:ite        ),intent (out  ) ::  vvel1d
  !
  !  local variables in this routine
     integer                             ::  i,k,k1,nvs
     real                                ::  dz,BU,dw2,dw1,kx,dz1m,Tv,Tve,vs,ftun1,ftun2
     real   , parameter :: f=2., C_d=0.506, gam=0.5, beta=1.875 !,ftun1=0.5, ftun2=0.8
     logical, parameter :: smooth=.true.
     integer, parameter :: n_smooth=1

     ftun1=0.25 ; ftun2=1.
     if(ZERO_DIFF==1) then
        ftun1=1. ; ftun2=0.5
     endif

     do i=its,itf
        !-- initialize arrays to zero.
        vvel1d(i  ) = 0.0
        vvel2d(i,:) = 0.0

        if(ierr(i) /= 0) cycle
        vvel2d(i,kts:kbcon(i))= max(1.,zws(i)**2)

loop0:  do k= kbcon(i),ktop(i)

          dz=z_cup(i,k+1)-z_cup(i,k)

          Tve= 0.5* ( t_cup (i,k  )*(1.+(qo (i,k  )/eps)/(1.+qo (i,k  ))) + &
                      t_cup (i,k+1)*(1.+(qo (i,k+1)/eps)/(1.+qo (i,k+1))))

          Tv = 0.5* ( tempco(i,k  )*(1.+(qco(i,k  )/eps)/(1.+qco(i,k  ))) + &
                      tempco(i,k+1)*(1.+(qco(i,k+1)/eps)/(1.+qco(i,k+1)) ))

          BU = g*( (Tv-Tve)/Tve -  ftun2*0.50*(qrco(i,k+1)+qrco(i,k) ))

          dw1 = 2./(f*(1.+gam)) * BU * dz
          if(ZERO_DIFF==1) then
             kx  =               max(entr_rate(i,k),cd(i,k))*dz
          else
             kx  = (1.+beta*C_d)*max(entr_rate(i,k),cd(i,k))*dz*ftun1
          endif

          dw2 =  (vvel2d(i,k)) -2.*kx * (vvel2d(i,k))

          vvel2d(i,k+1)=(dw1+dw2)/(1.+kx)

          if( vvel2d(i,k+1)< 0.) then
            vvel2d(i,k+1) = 0.5* vvel2d(i,k)
          endif

        enddo loop0
     enddo
     if(smooth) then
      if(ZERO_DIFF==1) then
       do i=its,itf
         if(ierr(i) /= 0)cycle
         do k=kts,ktop(i)-2
           nvs=0; vs =0.
           do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
             nvs = nvs + 1
              vs =  vs + vvel2d(i,k1)
           enddo
           vvel2d(i,k) = vs/(1.e-16+float(nvs))
         enddo
       enddo
      else
       do i=its,itf
         if(ierr(i) /= 0)cycle
         do k=kts,ktop(i)+1
           vs =0.; dz1m= 0.
           do k1 = max(k-n_smooth,kts),min(k+n_smooth,ktf)
                dz   = z_cup(i,k1+1)-z_cup(i,k1)
                vs   =  vs + dz*vvel2d(i,k1)
            dz1m = dz1m + dz
           enddo
           vvel2d(i,k) = vs/(1.e-16+dz1m)
           !if(k>ktop(i)-3)print*,"v2=",k,ktop(i),sqrt(vvel2d(i,k)),sqrt(vvel2d(i,ktop(i)))
         enddo
       enddo
      endif
     endif

     !-- convert to vertical velocity
     do i=its,itf
         if(ierr(i) /= 0)cycle
         vvel2d(i,:)= sqrt(max(0.1,vvel2d(i,:)))

         if(maxval(vvel2d(i,:)) < 1.0) then
           ierr(i)=54
         !  print*,"ierr=54",maxval(vvel2d(i,:))
         endif


         !-- sanity check
         where(vvel2d(i,:) < 1. ) vvel2d(i,:) = 1.
         where(vvel2d(i,:) > 20.) vvel2d(i,:) = 20.
         if(ZERO_DIFF==0)         vvel2d(i,ktop(i)+1:kte) = 0.1

         !-- get the column average vert velocity
         do k= kbcon(i),ktop(i)
            dz=z_cup(i,k+1)-z_cup(i,k)
            vvel1d(i)=vvel1d(i)+vvel2d(i,k)*dz
            !print*,"w=",k,z_cup(i,k),vvel2d(i,k)
         enddo
         vvel1d(i)=vvel1d(i)/(z_cup(i,ktop(i)+1)-z_cup(i,kbcon(i))+1.e-16)
         vvel1d(i)=max(1.,vvel1d(i))
     enddo

   end subroutine cup_up_vvel

!------------------------------------------------------------------------------------
   SUBROUTINE cup_output_ens_3d(name,xff_shal,xff_mid,xf_ens,ierr,dellat,dellaq,dellaqc,  &
                                outtem,outq,outqc,zu,pre,pw,xmb,ktop,                     &
                                nx,nx2,ierr2,ierr3,pr_ens, maxens3,ensdim,sig,xland1,     &
                                ichoice,ipr,jpr,itf,ktf,its,ite, kts,kte,                 &
                                xf_dicycle,outu,outv,dellu,dellv,dtime,po_cup,kbcon,       &
                                dellabuoy,outbuoy, dellampqi,outmpqi,dellampql,outmpql,   &
                                dellampcf,outmpcf ,nmp)
   IMPLICIT NONE
!
!  on input
!

   ! only local dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ichoice,ipr,jpr,itf,ktf,                                       &
        its,ite, kts,kte
     integer, intent (in   )              ::                           &
        ensdim,nx,nx2,maxens3,nmp
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     character *(*), intent (in)          ::    name
     real,    dimension (its:ite,1:ensdim)                             &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc,outu,outv,outbuoy

     real,    dimension (nmp,its:ite,kts:kte)                          &
        ,intent (out  )                   ::                           &
        outmpqi,outmpql,outmpcf
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                            &
        zu,po_cup
     real,   dimension (its:ite)                                       &
         ,intent (in  )                   ::                           &
        sig
     real,   dimension (its:ite,maxens3)                               &
         ,intent (in  )                   ::                           &
        xff_mid
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                 ::                         &
        xland1
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        dellat,dellaqc,dellaq,pw,dellu,dellv,dellabuoy
     real,    dimension (nmp,its:ite,kts:kte)                          &
        ,intent (in   )                   ::                           &
        dellampqi,dellampql,dellampcf

     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     real,    intent(in), dimension (its:ite) :: xf_dicycle
     real,    intent(in) :: dtime
     real,   dimension (its:ite,9)                                     &
         ,intent (in  )                   ::                           &
        xff_shal
!
!  local variables in this routine
!

     integer                              ::  i,k,n,ncount,zmax,kk
     real                                 ::  outtes,ddtes,dtt,dtq,dtqc&
                                             ,dtpw,prerate,fixouts,dp
     real                                 ::  dtts,dtqs,fsum, rcount
     real,    dimension (its:ite)         ::  xmb_ave,xmbmax
     real,    dimension (kts:kte,8)       ::  tend2d
     real,    dimension (8)               ::  tend1d
     real,    dimension (its:ite,8)       ::  check_cons_I,check_cons_F
!
      do k=kts,ktf
       do i=its,itf
        outtem   (i,k)=0.
        outq     (i,k)=0.
        outqc    (i,k)=0.
        outu     (i,k)=0.
        outv     (i,k)=0.
        outbuoy  (i,k)=0.
       enddo
      enddo
      do i=its,itf
        pre(i)  =0.
        xmb(i)  =0.
      enddo

      do i=its,itf
        if(ierr(i).eq.0)then
         do n=1,maxens3
           if(pr_ens(i,n).le.0.)then
             xf_ens(i,n)=0.
           endif
         enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!
      IF(NAME == 'deep') THEN
        DO i=its,itf
         if(ierr(i).eq.0)then
          k=0
          xmb_ave(i)=0.
          do n=1,maxens3
           k=k+1
           xmb_ave(i)=xmb_ave(i)+xf_ens(i,n)
          enddo
          !- 'ensemble' average mass flux
          xmb_ave(i)=xmb_ave(i)/float(k)
         endif
        ENDDO

      !- mid (congestus type) convection
      ELSEIF(NAME=='mid') then
         if(ichoice .le. 5) then
             do i=its,itf
              if(ierr(i) /= 0) cycle
              if(ichoice == 0) then
                 xmb_ave(i)=0.3333*(xff_mid(i,1)+xff_mid(i,2)+xff_mid(i,3))
              else
                 xmb_ave(i)= xff_mid(i,ichoice)
              endif
             enddo
         else
             stop 'For mid ichoice must be 0,..,5'
         endif

      !- shallow  convection
      ELSEIF(NAME=='shallow') then
        do i=its,itf
          if(ierr(i) /= 0) cycle

          if(ichoice > 0) then
            xmb_ave(i)=xff_shal(i,ichoice)
          else
            fsum=0.
            xmb_ave(i)=0.
            do k=1,9
               !- heat engine closure is providing too low values for mass fluxes.
               !- until this is checked, the ensemble closure will be calculatd
               !- only using the closures BLQE and Wstar
               !if(k.ge.4 .and. k.le.6) cycle
               xmb_ave(i)=xmb_ave(i)+xff_shal(i,k)
               fsum=fsum+1.
            enddo
            !- ensemble average of mass flux
            xmb_ave(i)=xmb_ave(i)/fsum
          endif
       enddo
      ENDIF

      !- set the updradt mass flux and do not allow negative values and apply the diurnal cycle closure
      DO i=its,itf
          if(ierr(i) /= 0) cycle
          !- mass flux of updradt at cloud base
          xmb(i) = xmb_ave(i)

          !- apply the adjust factor for tunning
          xmb(i) = FADJ_MASSFLX * xmb(i)

          !- add uplift by cold pools
          !if(name == 'deep')  xmb(i) = xmb(i) + xmbdn(i)

          !- diurnal cycle closure
          xmb(i) = xmb(i) - xf_dicycle(i)
          if(xmb(i) .le. 0.)then
               ierr(i)=13; xmb (i)=0.
          endif
      ENDDO
      !-apply the scale-dependence Arakawa's approach
      DO i=its,itf
          if(ierr(i) /= 0) cycle
          !- scale dependence
          xmb(i)=sig(i)*xmb(i)

          if(xmb(i) == 0. ) ierr(i)=14
          if(xmb(i) > 100.) ierr(i)=15
      ENDDO

!--- sanity check for mass flux
!
      DO i=its,itf
           if(ierr(i) /= 0) cycle
           xmbmax(i)=100.*(po_cup(i,kbcon(i))-po_cup(i,kbcon(i)+1))/(g*dtime)
           xmb(i) = min(xmb(i),xmbmax(i))
      ENDDO

!--- check outtem and and outq for high values
!--- criteria: if abs (dT/dt or dQ/dt) > 100 K/day => fix xmb
      IF(MAX_TQ_TEND > 0.) THEN
      DO i=its,itf
        IF(ierr(i) /= 0) CYCLE
        fixouts=xmb(i) *86400.*max(maxval(abs(dellat(i,kts:ktop(i)))),&
                          (xlv/cp)*maxval(abs(dellaq(i,kts:ktop(i)))) )

          if(fixouts > MAX_TQ_TEND) then ! K/day
            fixouts=MAX_TQ_TEND/(fixouts)
          xmb   (i)  = xmb   (i)  *fixouts
          xf_ens(i,:)= xf_ens(i,:)*fixouts
       endif
      ENDDO
      ENDIF
!
!-- now do feedback
!
      DO i=its,itf
        IF(ierr(i) /= 0) CYCLE
        DO k=kts,ktop(i)
           pre    (i)  = pre(i)+pw(i,k)*xmb(i)

           outtem   (i,k)= dellat     (i,k)*xmb(i)
           outq     (i,k)= dellaq     (i,k)*xmb(i)
           outqc    (i,k)= dellaqc    (i,k)*xmb(i)
           outu     (i,k)= dellu      (i,k)*xmb(i)
           outv     (i,k)= dellv      (i,k)*xmb(i)
           outbuoy  (i,k)= dellabuoy  (i,k)*xmb(i)
        ENDDO
        xf_ens (i,:)= sig(i)*xf_ens(i,:)

        IF(APPLY_SUB_MP == 1) THEN
         DO k=kts,ktop(i)
           outmpqi(:,i,k)= dellampqi(:,i,k)*xmb(i)
           outmpql(:,i,k)= dellampql(:,i,k)*xmb(i)
           outmpcf(:,i,k)= dellampcf(:,i,k)*xmb(i)
        ENDDO
        outmpqi(:,i,ktop(i):ktf)=0.
        outmpql(:,i,ktop(i):ktf)=0.
         outmpcf(:,i,ktop(i):ktf)=0.
        ENDIF
      ENDDO
!
!--  smooth the tendencies (future work: include outbuoy, outmpc* and tracers)
!
      IF(USE_SMOOTH_TEND < 0) THEN
        DO i=its,itf
           IF(ierr(i) /= 0) CYCLE
           tend2d=0.

           !--- get the initial integrals
           rcount = 1.e-6; tend1d=0.
           DO k=kts,ktop(i)
              dp         = (po_cup(i,k)-po_cup(i,k+1))
              rcount     = rcount + dp
              tend1d(1)  = tend1d(1)  +  dp*outtem (i,k)
              tend1d(2)  = tend1d(2)  +  dp*outq   (i,k)
              tend1d(3)  = tend1d(3)  +  dp*outqc  (i,k)
              tend1d(4)  = tend1d(4)  +  dp*outu   (i,k)
              tend1d(5)  = tend1d(5)  +  dp*outv   (i,k)

           ENDDO
           check_cons_I(i,1:5) = tend1d(1:5)/rcount
           !check_cons_I(i,2) = tend1d(2)/rcount
           !check_cons_I(i,3) = tend1d(3)/rcount
           !check_cons_I(i,4) = tend1d(4)/rcount
           !check_cons_I(i,5) = tend1d(5)/rcount
           !---

           !--- make the smoothness procedure
           DO k=kts,ktop(i)
             rcount = 1.e-6 ; tend1d=0.
             !print*,"xx=",max(kts,k-USE_SMOOTH_TEND),min(ktop(i),k+USE_SMOOTH_TEND)
             DO kk= max(kts,k-USE_SMOOTH_TEND),min(ktop(i),k+USE_SMOOTH_TEND)
                  dp=(po_cup(i,kk)-po_cup(i,kk+1))
                  rcount = rcount + dp

                  tend1d(1)  = tend1d(1)  +  dp*outtem (i,kk)
                  tend1d(2)  = tend1d(2)  +  dp*outq   (i,kk)
                  tend1d(3)  = tend1d(3)  +  dp*outqc  (i,kk)
                  tend1d(4)  = tend1d(4)  +  dp*outu   (i,kk)
                  tend1d(5)  = tend1d(5)  +  dp*outv   (i,kk)

             ENDDO
             tend2d(k,1:5)  = tend1d(1:5) /rcount
             !tend2d(k,2)  = tend1d(2) /rcount
             !tend2d(k,3)  = tend1d(3) /rcount
             !tend2d(k,4)  = tend1d(4) /rcount
             !tend2d(k,5)  = tend1d(5) /rcount
           ENDDO
           !--- get the final/smoother tendencies
           DO k=kts,ktop(i)
             outtem (i,k) = tend2d(k,1)
             outq   (i,k) = tend2d(k,2)
             outqc  (i,k) = tend2d(k,3)
             outu   (i,k) = tend2d(k,4)
             outv   (i,k) = tend2d(k,5)
           ENDDO
cycle


           !--- check the final integrals
           rcount = 1.e-6; tend1d=0.
           DO k=kts,ktop(i)

              dp         = (po_cup(i,k)-po_cup(i,k+1))
              rcount     = rcount + dp
              tend1d(1)  = tend1d(1)  +  dp*outtem (i,k)
              tend1d(2)  = tend1d(2)  +  dp*outq   (i,k)
              tend1d(3)  = tend1d(3)  +  dp*outqc  (i,k)
              tend1d(4)  = tend1d(4)  +  dp*outu   (i,k)
              tend1d(5)  = tend1d(5)  +  dp*outv   (i,k)

           ENDDO
           !--- get the ratio between initial and final integrals.
           check_cons_F(i,1:5) = tend1d(1:5)/rcount
           !check_cons_F(i,2) = tend1d(2)/rcount
           !check_cons_F(i,3) = tend1d(3)/rcount
           !check_cons_F(i,4) = tend1d(4)/rcount
           !check_cons_F(i,5) = tend1d(5)/rcount
          !
           !--- apply correction to preserve the integrals
           DO kk=1,5

            if(abs(check_cons_F(i,kk))>0.) then
                      check_cons_F(i,kk) = abs(check_cons_I(i,kk)/check_cons_F(i,kk))
            else
                   check_cons_F(i,kk) = 1.
            endif
           ENDDO

           !check_cons_F(i,1) = abs(check_cons_I(i,1)/check_cons_F(i,1))
           !check_cons_F(i,2) = abs(check_cons_I(i,2)/check_cons_F(i,2))
           !check_cons_F(i,3) = abs(check_cons_I(i,3)/check_cons_F(i,3))
           !check_cons_F(i,4) = abs(check_cons_I(i,4)/check_cons_F(i,4))
           !check_cons_F(i,5) = abs(check_cons_I(i,5)/check_cons_F(i,5))

!cycle
           DO k=kts,ktop(i)
             outtem (i,k) = outtem (i,k) * check_cons_F(i,1)
             outq   (i,k) = outq   (i,k) * check_cons_F(i,2)
             outqc  (i,k) = outqc  (i,k) * check_cons_F(i,3)
                  outu   (i,k) = outu   (i,k) * check_cons_F(i,4)
                   outv   (i,k) = outv   (i,k) * check_cons_F(i,5)
           ENDDO

!           print*,"check=",real( (check_cons_F(i,3)+check_cons_F(i,2))/(check_cons_I(i,3)+check_cons_I(i,2)),4)!&
!                          ,real( check_cons_F(i,2)/check_cons_I(i,2),4)&
!                          ,real( check_cons_F(i,1)/check_cons_I(i,1),4)
           !-----
        ENDDO

      ENDIF

   END SUBROUTINE cup_output_ens_3d
!------------------------------------------------------------------------------------

  SUBROUTINE cup_forcing_ens_3d(itf,ktf,its,ite, kts,kte,ens4,ensdim,ichoice,maxens,maxens2,maxens3&
                                ,ierr,ierr2,ierr3,k22,kbcon,ktop      &
                                ,xland,aa0,aa1,xaa0,mbdt,dtime        &
                                ,xf_ens,mconv,qo                      &
                                ,p_cup,omeg,zd,zu,pr_ens,edt          &
                                ,tau_ecmwf,aa1_bl,xf_dicycle,xk_x  )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        itf,ktf,its,ite, kts,kte,ens4
     integer, intent (in   )              ::                           &
        ensdim,maxens,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! ichoice       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (its:ite,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (its:ite,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,zu,p_cup,qo
     real,    dimension (its:ite,kts:kte,1:ens4)                       &
        ,intent (in   )                   ::                           &
        omeg
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,xland
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        mconv
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        ichoice
      real,    intent(IN)   , dimension (its:ite) :: aa1_bl,tau_ecmwf
      real,    intent(INOUT), dimension (its:ite) :: xf_dicycle,xk_x
      !- local var
      real  :: xff_dicycle
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (its:ite)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3
     real                                 ::                           &
       a1,a_ave,xff0,xomg

     real :: betajb
     integer :: kk
     real, dimension (its:ite) :: ens_adj!,xmbmax
!
       ens_adj(:)=1.

!--- LARGE SCALE FORCING
!
       DO i=its,itf
          xf_ens(i,1:16)= 0.
          IF(ierr(i) /=  0)cycle

                xff0 = (AA1(I)-AA0(I))/DTIME
                !-- default
                xff_ens3(1) = max(0.,(AA1(I)   -AA0(I))/dtime)

                xff_ens3(2) = xff_ens3(1)
                xff_ens3(3) = xff_ens3(1)
                xff_ens3(16)= xff_ens3(1)
!
!--- more like Brown (1979), or Frank-Cohen (199?)
!--- omeg is in Pa/s
                xomg=0.
                kk=0
                xff_ens3(4)=0.
                do k=max(kts,kbcon(i)-1),kbcon(i)+1
                   !-  betajb=(zu(i,k)-edt(i)*zd(i,k))
                   betajb=1.
                   !if(betajb .gt. 0.)then
                       xomg=xomg-omeg(i,k,1)/g/betajb
                       kk=kk+1
                   !endif
                enddo
                if(kk.gt.0)xff_ens3(4)=xomg/float(kk) ! kg[air]/m^3 * m/s
                xff_ens3(4) = max(0.0, xff_ens3(4))
                xff_ens3(5) = xff_ens3(4)
                xff_ens3(6) = xff_ens3(4)
                xff_ens3(14)= xff_ens3(4)
!
!--- more like Krishnamurti et al.;
!
                !mconv(i) = 0.
                !do k=k22(i),ktop(i)
                !    mconv(i)=mconv(i)+omeg(i,k,1)*(qo(i,k+1)-qo(i,k))/g
                !enddo
                !- 2nd option (assuming that omeg(ktop)*q(ktop)<< omeg(kbcon)*q(kbcon))
                mconv(i)  = -omeg(i,kbcon(i),1)*qo(i,kbcon(i))/g ! (kg[air]/m^3)*m/s*kg[water]/kg[air]

                mconv(i)  = max(0., mconv(i))
                xff_ens3(7) = mconv(i)
                xff_ens3(8) = xff_ens3(7)
                xff_ens3(9) = xff_ens3(7)
                xff_ens3(15)= xff_ens3(7)
!
!---- more like  Betchold et al (2014). Note that AA1 already includes the forcings tendencies
                xff_ens3(10)= AA1(i)/tau_ecmwf(i)

                xff_ens3(11)= xff_ens3(10)
                xff_ens3(12)= xff_ens3(10)
                xff_ens3(13)= xff_ens3(10)

!
                if(ichoice == 0)then
                   if(xff0 < 0.)then
                     xff_ens3( 1)=0.
                     xff_ens3( 2)=0.
                     xff_ens3( 3)=0.
                     xff_ens3(16)=0.

                     xff_ens3(10)=0.
                     xff_ens3(11)=0.
                     xff_ens3(12)=0.
                     xff_ens3(13)=0.
                  endif
                endif

                xk(i)=(XAA0(I)-(AA1(I)          ))/MBDT(I)
                if(xk(i).le.0. .and. xk(i).gt.-0.1*mbdt(i)) xk(i)=-0.1*mbdt(i)
                if(xk(i).gt.0. .and. xk(i).lt.1.e-2       ) xk(i)=1.e-2
!
!---  over water, enfor!e small cap for some of the closures
!
                if(xland(i).lt.0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
                      xff_ens3(1:16) = ens_adj(i)*xff_ens3(1:16)
                 endif
                endif
!
!--- special treatment for stability closures
!
                if(xk(i).lt.0.)then
                   if(xff_ens3( 1).gt.0.)xf_ens(i, 1)=max(0.,-xff_ens3( 1)/xk(i))
                   if(xff_ens3( 2).gt.0.)xf_ens(i, 2)=max(0.,-xff_ens3( 2)/xk(i))
                   if(xff_ens3( 3).gt.0.)xf_ens(i, 3)=max(0.,-xff_ens3( 3)/xk(i))
                   if(xff_ens3(16).gt.0.)xf_ens(i,16)=max(0.,-xff_ens3(16)/xk(i))
                else
                   xff_ens3(1 )=0.
                   xff_ens3(2 )=0.
                   xff_ens3(3 )=0.
                   xff_ens3(16)=0.
                endif

                xf_ens(i,4) =max(0.,xff_ens3(4) )
                xf_ens(i,5) =max(0.,xff_ens3(5) )
                xf_ens(i,6) =max(0.,xff_ens3(6) )
                xf_ens(i,14)=max(0.,xff_ens3(14))

                a1=max(1.e-3,pr_ens(i,7) ); xf_ens(i,7) =max(0.,xff_ens3(7)/a1)
                a1=max(1.e-3,pr_ens(i,8) ); xf_ens(i,8) =max(0.,xff_ens3(8)/a1)
                a1=max(1.e-3,pr_ens(i,9) ); xf_ens(i,9) =max(0.,xff_ens3(9)/a1)
                a1=max(1.e-3,pr_ens(i,15)); xf_ens(i,15)=max(0.,xff_ens3(15)/a1)
                if(xk(i).lt.0.)then
                   xf_ens(i,10)= max(0.,-xff_ens3(10)/xk(i))
                   xf_ens(i,11)= max(0.,-xff_ens3(11)/xk(i))
                   xf_ens(i,12)= max(0.,-xff_ens3(12)/xk(i))
                   xf_ens(i,13)= max(0.,-xff_ens3(13)/xk(i))
                 else
                   xf_ens(i,10)= 0.
                   xf_ens(i,11)= 0.
                   xf_ens(i,12)= 0.
                   xf_ens(i,13)= 0.
                endif


                if(ichoice.ge.1)then
                  xf_ens(i,1:16) =xf_ens(i,ichoice)
                endif

!---special combination for 'ensemble closure':
!---over the land, only applies closures 1 and 10.
!if(ichoice == 0 .and. xland(i) < 0.1)then
!  xf_ens(i,1:16) =0.5*(xf_ens(i,10)+xf_ens(i,1))
!endif

!---over the land, only applies closure 10.
if(zero_diff == 0 .and. ichoice == 0) then
  xf_ens(i,1:16)=(1.-xland(i))*xf_ens(i,10)+xland(i)*xf_ens(i,1:16)
endif

!------------------------------------


        ENDDO
!-
!- diurnal cycle mass flux
!-
IF(DICYCLE==1 .or. DICYCLE==6 )THEN

       DO i=its,itf
          xf_dicycle(i) = 0.
          IF(ierr(i) /=  0)cycle

            xff_dicycle  = (AA1(i)-AA1_BL(i))/tau_ecmwf(i)
            if(xk(i).lt.0) xf_dicycle(i)= max(0.,-xff_dicycle/xk(i))
            xf_dicycle(i)= xf_ens(i,10)-xf_dicycle(i)
       ENDDO
ELSEIF( DICYCLE==2) THEN ! for trigger function only
       DO i=its,itf
          xf_dicycle(i) = 0.
          IF(ierr(i) /=  0)cycle

               xff_ens3(1) = max(0.,-(AA1_BL(I) - AA0(I))/dtime)

               if( xff_ens3(1)   <= 0. ) then
                 xf_ens    (i,:) = 0.0
                 xf_dicycle(i)   = 0.0
                 ierr(i)=45
                 !ierrc(i)="DCclosure"
               endif
       ENDDO
ELSEIF( DICYCLE==3) THEN

       DO i=its,itf
          xf_dicycle(i) = 0.
          IF(ierr(i) /=  0)cycle
               xff_ens3(1) = 5.*max(0.,-(AA1(I) - AA0(I))/dtime)

               if(xk(i)<0.   )xf_dicycle(i)=max(0.,-xff_ens3(1)/xk(i))
               xf_ens    (i,:) = xf_dicycle(i)
               xf_dicycle(i)   = 0.0

       ENDDO
ELSEIF( DICYCLE==4) THEN

       DO i=its,itf
          xf_dicycle(i) = 0.
          IF(ierr(i) /=  0)cycle
                !the signal "-" is to convert from Pa/s to kg/m2/s
                if(xk_x(i) > 0.) xf_dicycle(i)= max(0., -AA1_BL(I))/xk_x(i)

                xf_ens    (i,:) = xf_dicycle(i)
                xf_dicycle(i)   = 0.0
       ENDDO
ELSEIF( DICYCLE==5) THEN

       DO i=its,itf
          xf_dicycle(i) = 0.
          IF(ierr(i) /=  0)cycle

               xff_ens3(1) =max(0.,-(AA1_BL(I) - AA0(I))/dtime)

               if(xk(i)<0.                     )xf_dicycle(i)=max(0.,-xff_ens3(1)/xk(i))
!              if(xk(i)<0. .and. xff_ens3(1)>0.)xf_dicycle(i)=max(0.,-xff_ens3(1)/xk(i))
               xf_ens    (i,:) = xf_dicycle(i)
               xf_dicycle(i)   = 0.0
       ENDDO
ELSE

       xf_dicycle(:)=0.0

ENDIF
!--


   END SUBROUTINE cup_forcing_ens_3d

!------------------------------------------------------------------------------------
   SUBROUTINE get_partition_liq_ice(ierr,tn,z1,zo_cup,po_cup, p_liq_ice,melting_layer         &
                                   ,itf,ktf,its,ite, kts,kte, cnvfrc, srftype, cumulus )
     IMPLICIT NONE
     CHARACTER *(*), INTENT (IN)                          :: cumulus
     INTEGER  ,INTENT (IN   )                             :: itf,ktf, its,ite, kts,kte
     INTEGER  ,INTENT (IN   ), DIMENSION(its:ite)         :: ierr
     REAL     ,INTENT (IN   ), DIMENSION(its:ite)         :: z1, cnvfrc, srftype
     REAL     ,INTENT (IN   ), DIMENSION(its:ite,kts:kte) :: tn,po_cup,zo_cup
     REAL     ,INTENT (INOUT), DIMENSION(its:ite,kts:kte) :: p_liq_ice,melting_layer
     INTEGER :: i,k
     REAL    :: dp, height
     REAL, DIMENSION(its:ite) :: norm
     REAL, PARAMETER ::  T1=276.16, Z_meltlayer1=4000.,Z_meltlayer2=6000.,delT=3.
     p_liq_ice    (:,:) = 1.
     melting_layer(:,:) = 0.
     !-- get function of T for partition of total condensate into liq and ice phases.
     IF(MELT_GLAC .and. trim(cumulus) == 'deep') then
        DO k=kts,ktf
          DO i=its,itf
             if(ierr(i) /= 0) cycle
             p_liq_ice(i,k) = fract_liq_f(tn(i,k),cnvfrc(i),srftype(i))
         ENDDO
        ENDDO
!        go to 650
!
        !-- define the melting layer (the layer will be between T_0+1 < TEMP < T_1
        !-- definition em terms of temperatura
        DO k=kts,ktf
          DO i=its,itf
             if(ierr(i) /= 0) cycle
             if    (tn(i,k) <= T_0-delT) then
                melting_layer(i,k) = 0.

             elseif(  tn(i,k) < T_0+delT .and. tn(i,k) > T_0-delT) then
                melting_layer(i,k) =  ((tn(i,k)-(T_0-delt))/(2.*delT))**2

             else
                melting_layer(i,k) = 1.
             endif
             melting_layer(i,k) = melting_layer(i,k)*(1.-melting_layer(i,k))
          ENDDO
        ENDDO
!go to 655
!650 continue
!        !-- definition em terms of height above local terrain
!        DO k=kts,ktf
!          DO i=its,itf
!             if(ierr(i) /= 0) cycle
!             height= zo_cup(i,k)+z1(i)
!             if   (height > Z_meltlayer2 ) then
!                melting_layer(i,k) = 0.
!
!             elseif(height > Z_meltlayer1  .and. height < Z_meltlayer2 ) then
!
!                melting_layer(i,k) =  ((height - Z_meltlayer1)/(Z_meltlayer2-Z_meltlayer1))**2.
!
!
!             else
!                melting_layer(i,k) = 1.
!             endif
!             melting_layer(i,k) = melting_layer(i,k)*(1.-melting_layer(i,k))
!          ENDDO
!        ENDDO
!
!
!         655 continue
        !-normalize vertical integral of melting_layer to 1
        norm(:)=0.
        DO k=kts,ktf-1
          DO i=its,itf
             if(ierr(i) /= 0) cycle
             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
             norm(i) = norm(i) + melting_layer(i,k)*dp/g
          ENDDO
        ENDDO
        DO i=its,itf
         if(ierr(i) /= 0) cycle
         melting_layer(i,:)=melting_layer(i,:)/(norm(i)+1.e-6)*(100*(po_cup(i,kts)-po_cup(i,ktf))/g)
         !print*,"i2=",i,maxval(melting_layer(i,:)),minval(melting_layer(i,:)),norm(i)
        ENDDO
        !--check
!       norm(:)=0.
!        DO k=kts,ktf-1
!          DO i=its,itf
!             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
!             norm(i) = norm(i) + melting_layer(i,k)*dp/g/(100*(po_cup(i,kts)-po_cup(i,ktf))/g)
!             !print*,"n=",i,k,norm(i)
!          ENDDO
!        ENDDO

     !~ ELSE
        !~ p_liq_ice    (:,:) = 1.
        !~ melting_layer(:,:) = 0.
     ENDIF
   END  SUBROUTINE get_partition_liq_ice

!------------------------------------------------------------------------------------
   SUBROUTINE get_melting_profile(ierr,tn_cup,po_cup, p_liq_ice,melting_layer,qrco    &
                                 ,pwo,edto,pwdo,melting                               &
                                 ,itf,ktf,its,ite, kts,kte, cumulus                   )
     IMPLICIT NONE
     CHARACTER *(*), INTENT (IN)                          :: cumulus
     INTEGER  ,INTENT (IN   )                             :: itf,ktf, its,ite, kts,kte
     INTEGER  ,INTENT (IN   ), DIMENSION(its:ite)         :: ierr
     REAL     ,INTENT (IN   ), DIMENSION(its:ite)         :: edto
     REAL     ,INTENT (IN   ), DIMENSION(its:ite,kts:kte) :: tn_cup,po_cup,qrco,pwo &
                                                            ,pwdo,p_liq_ice,melting_layer
     REAL     ,INTENT (INOUT), DIMENSION(its:ite,kts:kte) :: melting
     INTEGER :: i,k
     REAL    :: dp
     REAL, DIMENSION(its:ite)         :: norm,total_pwo_solid_phase
     REAL, DIMENSION(its:ite,kts:kte) :: pwo_solid_phase,pwo_eff

     IF(MELT_GLAC .and. trim(cumulus) == 'deep') then

        norm                  = 0.0
        pwo_solid_phase       = 0.0
        pwo_eff               = 0.0
        melting               = 0.0
        !-- set melting mixing ratio to zero for columns that do not have deep convection
        DO i=its,itf
            if(ierr(i) > 0) melting(i,:) = 0.
        ENDDO

        !-- now, get it for columns where deep convection is activated
        total_pwo_solid_phase(:)=0.

        DO k=kts,ktf-1
          DO i=its,itf
           if(ierr(i) /= 0) cycle
             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))

             !-- effective precip (after evaporation by downdraft)
             !-- pwdo is not defined yet
            !pwo_eff(i,k) = 0.5*(pwo(i,k)+pwo(i,k+1) + edto(i)*(pwdo(i,k)+pwdo(i,k+1)))
             pwo_eff(i,k) = 0.5*(pwo(i,k)+pwo(i,k+1))

             !-- precipitation at solid phase(ice/snow)
             pwo_solid_phase(i,k) = (1.-p_liq_ice(i,k))*pwo_eff(i,k)

             !-- integrated precip at solid phase(ice/snow)
             total_pwo_solid_phase(i) = total_pwo_solid_phase(i)+pwo_solid_phase(i,k)*dp/g
          ENDDO
        ENDDO

        DO k=kts,ktf
          DO i=its,itf
           if(ierr(i) /= 0) cycle
             !-- melting profile (kg/kg)
             melting(i,k) = melting_layer(i,k)*(total_pwo_solid_phase(i)/(100*(po_cup(i,kts)-po_cup(i,ktf))/g))
             !print*,"mel=",k,melting(i,k),pwo_solid_phase(i,k),po_cup(i,k)
          ENDDO
        ENDDO

!-- check conservation of total solid phase precip
!       norm(:)=0.
!        DO k=kts,ktf-1
!          DO i=its,itf
!             dp = 100.*(po_cup(i,k)-po_cup(i,k+1))
!             norm(i) = norm(i) + melting(i,k)*dp/g
!          ENDDO
!        ENDDO
!
!       DO i=its,itf
!         print*,"cons=",i,norm(i),total_pwo_solid_phase(i)
!        ENDDO
!--

     ELSE
        !-- no melting allowed in this run
        melting     (:,:) = 0.
     ENDIF
   END  SUBROUTINE get_melting_profile
!------------------------------------------------------------------------------------
   SUBROUTINE  ke_to_heating(itf,ktf,its,ite, kts,kte,ktop,ierr &
                             ,po_cup,us,vs,dellu,dellv,dellat)

   implicit none
   integer                              ,intent (in   ) :: itf,ktf,its,ite, kts,kte
   integer, dimension (its:ite)         ,intent (in   ) :: ierr,ktop
   real   , dimension (its:ite,kts:kte) ,intent (in   ) :: po_cup,us,vs,dellu,dellv
   real   , dimension (its:ite,kts:kte) ,intent (inout) :: dellat

   real :: dts,fp,dp,fpi
   integer ::i,k

! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
!
   do i=its,itf
          if(ierr(i) /= 0) cycle
          dts=0.
          fpi=0.
          do k=kts,ktop(i)
             dp=(po_cup(i,k)-po_cup(i,k+1))*100.
             !total KE dissiptaion estimate
             dts = dts - (dellu(i,k)*us(i,k)+dellv(i,k)*vs(i,k))*dp/g
             !
             ! fpi needed for calcualtion of conversion to pot. energyintegrated
             fpi = fpi + sqrt(dellu(i,k)*dellu(i,k) + dellv(i,k)*dellv(i,k))*dp
          enddo
          if(fpi.gt.0.)then
             do k=kts,ktop(i)
                fp= sqrt((dellu(i,k)*dellu(i,k)+dellv(i,k)*dellv(i,k)))/fpi
                dellat(i,k) = dellat(i,k) + fp*dts*g/cp
             enddo
          endif
   enddo

  end SUBROUTINE  ke_to_heating
!------------------------------------------------------------------------------------
  function auto_rk(n,step,aux,xexp,qrc1) RESULT(PW)
      INTEGER, INTENT(in) :: n
      real   , INTENT(in) :: step,aux,qrc1,xexp
      real                :: PW

      PW=step*qrc1*(1.0-exp(-aux**xexp))/float(n)

  end function auto_rk

!-----------------------------------------------------------------------------------------

   SUBROUTINE fct1d3 (ktop,n,dt,z,tracr,massflx,trflx_in,del_out)

   ! --- modify a 1-D array of tracer fluxes for the purpose of maintaining
   ! --- monotonicity (including positive-definiteness) in the tracer field
   ! --- during tracer transport.

   ! --- the underlying transport equation is   (d tracr/dt) = - (d trflx/dz)
   ! --- where  dz = |z(k+1)-z(k)| (k=1,...,n) and  trflx = massflx * tracr

   ! --- note: tracr is carried in grid cells while z and fluxes are carried on
   ! --- interfaces. interface variables at index k are at grid location k-1/2.
   ! --- sign convention: mass fluxes are considered positive in +k direction.

   ! --- massflx and trflx_in  must be provided independently to allow the
   ! --- algorithm to generate an auxiliary low-order (diffusive) tracer flux
   ! --- as a stepping stone toward the final product trflx_out.

   implicit none
   integer,intent(IN) :: n,ktop                        ! number of grid cells
   real   ,intent(IN) :: dt                            ! transport time step
   real   ,intent(IN) :: z(n+0)                        ! location of cell interfaces
   real   ,intent(IN) :: tracr(n)                      ! the transported variable
   real   ,intent(IN) :: massflx  (n+0)                ! mass flux across interfaces
   real   ,intent(IN) :: trflx_in (n+0)                ! original tracer flux
   real   ,intent(OUT):: del_out  (n+0)                ! modified tracr flux
   real               :: trflx_out(n+0)                ! modified tracr flux
   integer k,km1,kp1
   logical :: NaN, error=.false., vrbos=.false.
   real dtovdz(n),trmax(n),trmin(n),flx_lo(n+0),antifx(n+0),clipped(n+0),  &
        soln_hi(n),totlin(n),totlout(n),soln_lo(n),clipin(n),clipout(n),arg
   real,parameter :: epsil=1.e-22           ! prevent division by zero
   real,parameter :: damp=1.                ! damper of antidff flux (1=no damping)

   logical, parameter :: hi_order = .false.

   NaN(arg) = .not. (arg.ge.0. .or. arg.le.0.)        ! NaN detector
   soln_lo(:)=0.
   antifx (:)=0.
   clipout(:)=0.
   flx_lo (:)=0.

   do k=1,ktop
     dtovdz(k)=.01*dt/abs(z(k+1)-z(k))                ! time step / grid spacing
!     if (z(k).eq.z(k+1)) error=.true.
   end do
   if (vrbos .or. error) print '(a/(8es10.3))','(fct1d) dtovdz =',dtovdz(1:ktop)

   do k=2,ktop
     if (massflx(k) > 0.) then
       flx_lo(k)=massflx(k)*tracr(k-1)              ! low-order flux, upstream
     else
       flx_lo(k)=massflx(k)*tracr(k)                ! low-order flux, upstream
     endif
     antifx(k)=trflx_in(k)-flx_lo(k)                ! antidiffusive flux
   end do
   flx_lo(  1)   =trflx_in(  1)
   flx_lo(ktop+1)=trflx_in(ktop+1)
   antifx(  1)   =0.
   antifx(ktop+1)=0.
! --- clip low-ord fluxes to make sure they don't violate positive-definiteness
   do k=1,ktop
    totlout(k)=max(0.,flx_lo(k+1))-min(0.,flx_lo(k  ))         ! total flux out
    clipout(k)=min(1.,tracr(k)/max(epsil,totlout(k))/ (1.0001*dtovdz(k)))
   end do

   do k=2,ktop
   if (massflx(k).ge.0.)  then
     flx_lo(k)=flx_lo(k)*clipout(k-1)
   else
     flx_lo(k)=flx_lo(k)*clipout(k)
   endif
   end do
   if (massflx(1)     .lt.0.) flx_lo(1)     =flx_lo(1)     *clipout(1)
   if (massflx(ktop+1).gt.0.) flx_lo(ktop+1)=flx_lo(ktop+1)*clipout(ktop)

! --- a positive-definite low-order (diffusive) solution can now be  constructed

   do k=1,ktop
     soln_lo  (k)=tracr(k)-(flx_lo(k+1)-flx_lo(k))*dtovdz(k)        ! low-ord solutn
     del_out  (k)=-g*(flx_lo(k+1)-flx_lo(k))*dtovdz(k)/dt
   end do

   if(.not. hi_order) return

   soln_hi  (:)=0.
   clipin   (:)=0.
   trmin    (:)=0.
   trmax    (:)=0.
   clipped  (:)=0.
   trflx_out(:)=0.


   do k=1,ktop
     km1=max(1,k-1)
     kp1=min(n,k+1)
     trmax(k)=       max(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
                         tracr  (km1),tracr  (k),tracr  (kp1))        ! upper bound
     trmin(k)=max(0.,min(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
                         tracr  (km1),tracr  (k),tracr  (kp1)))       ! lower bound
   end do

   do k=1,ktop
     totlin (k)=max(0.,antifx(k  ))-min(0.,antifx(k+1))                ! total flux in
     totlout(k)=max(0.,antifx(k+1))-min(0.,antifx(k  ))                ! total flux out

     clipin (k)=min(damp,(trmax(k)-soln_lo(k))/max(epsil,totlin (k)) / (1.0001*dtovdz(k)))
     clipout(k)=min(damp,(soln_lo(k)-trmin(k))/max(epsil,totlout(k)) / (1.0001*dtovdz(k)))

     if (NaN(clipin (k))) print *,'(fct1d) error: clipin is NaN,  k=',k
     if (NaN(clipout(k))) print *,'(fct1d) error: clipout is NaN,  k=',k

     if (clipin(k).lt.0.) then
       print 100,'(fct1d) error: clipin < 0 at k =',k,                        &
       'clipin',clipin(k),'trmax',trmax(k),'soln_lo',soln_lo(k),        &
       'totlin',totlin(k),'dt/dz',dtovdz(k)
       error=.true.
     end if
     if (clipout(k).lt.0.) then
       print 100,'(fct1d) error: clipout < 0 at k =',k,                        &
       'clipout',clipout(k),'trmin',trmin(k),'soln_lo',soln_lo(k),        &
       'totlout',totlout(k),'dt/dz',dtovdz(k)
       error=.true.
     end if
 100 format (a,i3/(4(a10,"=",es9.2)))
   end do

   do k=2,ktop
     if (antifx(k).gt.0.)  then
       clipped(k)=antifx(k)*min(clipout(k-1),clipin(k))
     else
       clipped(k)=antifx(k)*min(clipout(k),clipin(k-1))
     end if
     trflx_out(k)=flx_lo(k)+clipped(k)
     if (NaN(trflx_out(k)))  then
       print *,'(fct1d) error: trflx_out is NaN,  k=',k
       error=.true.
     end if
   end do

   trflx_out(     1)=trflx_in(     1)
   trflx_out(ktop+1)=trflx_in(ktop+1)
   do k=1,ktop
     soln_hi(k)=tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
     del_out(k) =     -g*(trflx_out(k+1)-trflx_out(k))*dtovdz(k)/dt
     !write(32,*)'3',k,soln_lo(k),soln_hi(k)
   end do

   if (vrbos .or. error) then
     do k=2,ktop
       write(32,99)k,                   &
       'tracr(k)', tracr(k),            &
       'flx_in(k)', trflx_in(k),        &
       'flx_in(k+1)', trflx_in(k+1),    &
       'flx_lo(k)', flx_lo(k),          &
       'flx_lo(k+1)', flx_lo(k+1),      &
       'soln_lo(k)', soln_lo(k),        &
       'trmin(k)', trmin(k),            &
       'trmax(k)', trmax(k),            &
       'totlin(k)', totlin(k),          &
       'totlout(k)', totlout(k),        &
       'clipin(k-1)', clipin(k-1),      &
       'clipin(k)', clipin(k),          &
       'clipout(k-1)', clipout(k-1),    &
       'clipout(k)', clipout(k),        &
       'antifx(k)', antifx(k),          &
       'antifx(k+1)', antifx(k+1),      &
       'clipped(k)', clipped(k),        &
       'clipped(k+1)', clipped(k+1),    &
       'flx_out(k)', trflx_out(k),      &
       'flx_out(k+1)', trflx_out(k+1),  &
       'dt/dz(k)', dtovdz(k),           &
       'final', tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
 99    format ('(trc1d)   k =',i4/(3(a13,'=',es13.6)))
     end do
     if (error) stop '(fct1d error)'
   end if
  END SUBROUTINE fct1d3
!---------------------------------------------------------------------------------------------------
  SUBROUTINE tridiag (m,a,b,c,f)
   !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
   !-- an updated "f" at time t+1 is the output
   implicit none
   integer, intent(in) :: m
   real, dimension(m), intent(inout) :: a,b,c
   real, dimension(m), intent(inout) :: f
   !--locals
   real, dimension(m) :: q
   integer :: k
   real :: p

   c(m)=0.
   q(1)=-c(1)/b(1)
   f(1)= f(1)/b(1)
   do k=2,m
    p         = 1./( b(k)+a(k)*q(k-1) )
    q(k) = -c(k)*p
    f(k) = p*(f(k) - a(k)*f(k-1))
   enddo
   do k=m-1,1,-1
    f(k) = f(k) +q(k)*f(k+1)
   enddo

  END SUBROUTINE tridiag
!---------------------------------------------------------------------------------------------------
  SUBROUTINE bidiag (m,b,c,f)
   !-- this routine solves the problem:  bb*f(k,t+1) + cc*f(k+1,t+1) = dd
   !-- an updated "f" at time t+1 is the output
   implicit none
   integer, intent(in) :: m
   real, dimension(m), intent(inout) :: b,c
   real, dimension(m), intent(inout) :: f
   !--locals
   real, dimension(m) :: q
   integer :: k
   real :: p

   c(m)=0.
   q(1)=-c(1)/b(1)
   f(1)= f(1)/b(1)
   do k=2,m
    p         = 1./b(k)
    q(k) = -c(k)*p
    f(k) =  f(k)*p
   enddo
   do k=m-1,1,-1
    f(k) = f(k) +q(k)*f(k+1)
   enddo

  END SUBROUTINE bidiag
!---------------------------------------------------------------------------------------------------
  SUBROUTINE cup_env_clev_chem(mtp,se_chem,se_cup_chem,ierr,itf,ktf,its,ite, kts,kte)

     IMPLICIT NONE
     !-inputs
     integer ,intent (in   )                   :: itf,ktf, its,ite, kts,kte
     integer ,intent (in   )                   :: mtp
     integer, dimension (its:ite),intent (in)  :: ierr
     real, dimension (mtp,its:ite,kts:kte),intent (in)  ::   se_chem
     !-outputs
     real, dimension (mtp,its:ite,kts:kte),intent (out) ::   se_cup_chem
     !-locals
     integer ::  i,k
     integer, parameter ::  clev_option = 2 !- use option 2

     !
     IF(clev_option == 1) THEN
      !-- original version
      do i=its,itf
        if(ierr(i) /= 0) cycle
        do k=kts+1,ktf
            se_cup_chem(1:mtp,i,k)=0.5*(se_chem(1:mtp,i,k-1)+se_chem(1:mtp,i,k))
        enddo
        se_cup_chem(1:mtp,i,kts)=se_chem(1:mtp,i,kts)
        se_cup_chem(1:mtp,i,kte)=se_chem(1:mtp,i,ktf)
      enddo
     ELSE
      !-- version 2: se_cup (k+1/2) = se(k) => smoother profiles
      do i=its,itf
                if(ierr(i) /= 0) cycle
                do k=kts,ktf
                    se_cup_chem(1:mtp,i,k)=se_chem(1:mtp,i,k)
           enddo
      enddo
     ENDIF

  END SUBROUTINE cup_env_clev_chem
!---------------------------------------------------------------------------------------------------

 SUBROUTINE rain_evap_below_cloudbase(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop,xmb,psur,xland&
                                    ,qo_cup,t_cup,po_cup,qes_cup,pwavo,edto,pwevo,pwo,pwdo&
                                    ,pre,prec_flx,evap_flx,outt,outq,outbuoy,evap_bcb)

 implicit none
 real, parameter :: alpha1=5.44e-4 & !1/sec
                   ,alpha2=5.09e-3 & !unitless
                   ,alpha3=0.5777  & !unitless
                   ,c_conv=0.05      !conv fraction area, unitless

 character*(*)                      ,intent(in)    :: cumulus
 integer                            ,intent(in)    :: itf,ktf, its,ite, kts,kte
 integer, dimension(its:ite)        ,intent(in)    :: ierr,kbcon,ktop
 real,    dimension(its:ite)        ,intent(in)    :: psur,xland,pwavo,edto,pwevo,xmb
 real,    dimension(its:ite,kts:kte),intent(in)    :: po_cup,qo_cup,qes_cup,pwo,pwdo,t_cup
 real,    dimension(its:ite)        ,intent(inout) :: pre
 real,    dimension(its:ite,kts:kte),intent(inout) :: outt,outq,outbuoy,prec_flx,evap_flx

 real,    dimension(its:ite,kts:kte),intent(out)   :: evap_bcb

 !-- locals
 integer :: i,k
 real    :: RH_cr , del_t,del_q,dp,q_deficit, pqsat
 real    :: RH_cr_OCEAN,RH_cr_LAND
 real,    dimension(its:ite) :: tot_evap_bcb,eff_c_conv

 if(trim(cumulus) == 'shallow') then
     RH_cr_OCEAN   = 1.
     RH_cr_LAND    = 1.
     eff_c_conv(:) = min(0.2,max(xmb(:),c_conv))
 else
     RH_cr_OCEAN   = 0.95 !test 0.90
     RH_cr_LAND    = 0.85
     eff_c_conv(:) = c_conv
 endif

 prec_flx     = 0.0
 evap_flx     = 0.0
 tot_evap_bcb = 0.0
!! if(trim(cumulus) == 'shallow') return

 do i=its,itf

     if(ierr(i) /= 0) cycle

     !-- critical rel humidity  - check this, if the value is too small, not evapo will take place.
      RH_cr=RH_cr_OCEAN*xland(i)+RH_cr_LAND*(1.0-xland(i))

     !if(xland(i)  < 0.90 ) then !- over land
     !  RH_cr = RH_cr_LAND
     !else
     !  RH_cr = RH_cr_OCEAN
     !endif

     do k=ktop(i),kts,-1

        dp = 100.*(po_cup(i,k)-po_cup(i,k+1))

        !---rainfall evaporation below cloud base
        if(k <= kbcon(i)) then
         q_deficit = max(0.,(RH_cr*qes_cup(i,k) -qo_cup(i,k)))
        !pqsat=satur_spec_hum(t_cup(i,k),po_cup(i,k))

         !--units here: kg[water]/kg[air}/sec
         evap_bcb(i,k) = eff_c_conv(i) * alpha1 * q_deficit * &
                      ( sqrt(po_cup(i,k)/psur(i))/alpha2 * prec_flx(i,k+1)/eff_c_conv(i) )**alpha3

         !--units here: kg[water]/kg[air}/sec * kg[air]/m3 * m = kg[water]/m2/sec
         evap_bcb(i,k)= evap_bcb(i,k)*dp/g

        else

         evap_bcb(i,k)=0.0

        endif

         !-- get the net precitation flux after the local evaporation and downdraft
         prec_flx(i,k) = prec_flx(i,k+1) - evap_bcb(i,k) + xmb(i)*(pwo(i,k) + edto(i)*pwdo(i,k))
         prec_flx(i,k) = max(0.,prec_flx(i,k))

         evap_flx(i,k) = evap_flx(i,k+1) + evap_bcb(i,k) - xmb(i)*edto(i)*pwdo(i,k)
         evap_flx(i,k) = max(0.,evap_flx(i,k))

         tot_evap_bcb(i) = tot_evap_bcb(i)+evap_bcb(i,k)

         !-- feedback
         del_q =  evap_bcb(i,k)*g/dp          ! > 0., units: kg[water]/kg[air}/sec
         del_t = -evap_bcb(i,k)*g/dp*(xlv/cp) ! < 0., units: K/sec

         outq   (i,k) = outq   (i,k) + del_q
         outt   (i,k) = outt   (i,k) + del_t
         outbuoy(i,k) = outbuoy(i,k) + cp*del_t+xlv*del_q

         pre(i) = pre(i) - evap_bcb(i,k)

        !--for future use (rain and snow precipitation fluxes)
        !p_liq_ice(i,k) = fract_liq_f(tempco(i,k),,cnvfrc(i),srftype(i))
        !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
        !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)


     enddo
 enddo

 END SUBROUTINE rain_evap_below_cloudbase
!---------------------------------------------------------------------------------------------------

   SUBROUTINE get_precip_fluxes(cumulus,klcl,kbcon,ktop,k22,ierr,xland,pre,xmb  &
                      ,pwo,pwavo,edto,pwevo,pwdo,t_cup,tempco                   &
                      ,prec_flx,evap_flx                                        &
                      ,itf,ktf,its,ite, kts,kte)

     implicit none
     character *(*)            , intent (in) :: cumulus
     integer                    ,intent (in) :: itf,ktf,its,ite,kts,kte
     integer, dimension(its:ite),intent (in) :: kbcon,ktop,k22,klcl,ierr
     real,    dimension(its:ite),intent (in) :: xland,pwavo,pwevo,edto,pre,xmb
     real,    dimension(its:ite,kts:kte),intent (in)  :: pwo,pwdo,t_cup,tempco
     real,    dimension(its:ite,kts:kte),intent (out) :: prec_flx,evap_flx !-- units kg[water]/m2/s

    !-- locals
     integer :: i,k
     prec_flx = 0.0
     evap_flx = 0.0
!!     if(trim(cumulus) == 'shallow') return

     DO i=its,itf
         if (ierr(i)  /= 0) cycle

         do k=ktop(i),kts,-1

           !--- precipitation flux (at 'cup' levels), units: kg[water]/m2/s
           prec_flx(i,k) = prec_flx(i,k+1) + xmb(i)*(pwo(i,k) + edto(i)*pwdo(i,k))
           prec_flx(i,k) = max(0.,prec_flx(i,k))

           !--- evaporation flux (at 'cup' levels), units: kg[water]/m2/s
           evap_flx(i,k) = evap_flx(i,k+1) - xmb(i)*edto(i)*pwdo(i,k)
           evap_flx(i,k) = max(0.,evap_flx(i,k))


           !
           !--for future use (rain and snow precipitation fluxes)
           !p_liq_ice(i,k) = fract_liq_f(tempco(i,k),,cnvfrc(i),srftype(i))
           !prec_flx_rain(k) = prec_flx(i,k)*(1.-p_liq_ice(k))
           !prec_flx_snow(k) = prec_flx(i,k)*    p_liq_ice(k)
         enddo

         !if(prec_flx   (i,kts) .ne. pre(i)) then
         !print*,"error=",100.*(prec_flx   (i,kts) - pre(i))/(1.e-16+pre(i)),pre(i),prec_flx   (i,kts)
         !STOP 'problem with water balance'
         !endif
     ENDDO
   END   SUBROUTINE get_precip_fluxes
!------------------------------------------------------------------------------------
   REAL FUNCTION satur_spec_hum(pt,press) result(pqsat)
    implicit none
    real   , intent(in ) :: pt,press ! Kelvin, hPa
   !real   , intent(out) :: pqsat  !saturation specific humidity kg/kg

   !---locals
    real :: zew,zqs,zcor,foealfcu,foeewmcu
    real, parameter ::                 &
    RD=287.06                          &
   ,RV=461.52                               &
   ,RTT=273.16                               &
   ,RETV=RV/RD-1.0                        &
   ,R2ES=611.21*RD/RV                       &
   ,R3LES=17.502                       &
   ,R3IES=22.587                       &
   ,R4LES=32.19                               &
   ,R4IES=-0.7                               &
   ,RTWAT=RTT                               &
   ,RTICE=RTT-23.                         &
   ,RTICECU=RTT-23.                       &
   ,RTWAT_RTICE_R=1./(RTWAT-RTICE)     &
   ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)

    foealfcu = min(1.0,((max(rticecu,min(rtwat,pt))-rticecu)*rtwat_rticecu_r)**2)
    foeewmcu = r2es *(foealfcu *exp(r3les*(pt-rtt)/(pt-r4les))+&
                    (1.0-foealfcu)*exp(r3ies*(pt-rtt)/(pt-r4ies)))

    zew  = foeewmcu
    zqs  = zew/(100.*press)
    if(1.0-retv*zqs > 0. )then
       zcor = 1.0/(1.0-retv*zqs)  ! divide by zero
       pqsat= zqs*zcor
    else
       pqsat= max_qsat
    endif

   END FUNCTION satur_spec_hum
!---------------------------------------------------------------------------------------------------
  SUBROUTINE get_jmin(cumulus,itf,ktf,its,ite, kts,kte,ierr,kdet,ktop,kbcon,jmin,ierrc  &
                   ,beta,depth_min,heso_cup,zo_cup,melting_layer)

     implicit none
     character *(*)             ,intent (in)    :: cumulus
     real                       ,intent (in)    :: depth_min
     integer                    ,intent (in)    :: itf,ktf,its,ite,kts,kte
     integer, dimension(its:ite),intent (in)    :: ktop,kbcon
     real,    dimension(its:ite,kts:kte),intent (in)  ::  heso_cup,zo_cup,melting_layer

     integer, dimension(its:ite),intent (inout) :: ierr,jmin,kdet
     real                       ,intent (out)   :: beta
     character*128,              intent (out)   :: ierrc(its:ite)


    !-- locals
     integer :: i,k,jmini,ki
     real    :: dh,dz
     real,    dimension(its:ite,kts:kte)  ::  hcdo
     logical :: keep_going

     if(trim(cumulus) == 'deep') beta=0.05
     if(trim(cumulus) == 'mid' ) beta=0.02

     if(trim(cumulus) == 'shallow'  ) then
        beta    = 0.02
        jmin(:) = 0
        return
     endif

     do i=its,itf
         if(ierr(i) /= 0) cycle

         if(trim(cumulus) == 'deep' .and. melt_glac) jmin(i)=max(jmin(i),maxloc(melting_layer(i,:),1))

         !
         !--- check whether it would have buoyancy, if there where
         !--- no entrainment/detrainment
         !
         jmini = jmin(i)
         keep_going = .true.
         do while ( keep_going )
           keep_going = .false.
           if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
           if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
           ki = jmini
           hcdo(i,ki)=heso_cup(i,ki)
           dz=zo_cup(i,ki+1)-zo_cup(i,ki)
           dh=0.
           do k=ki-1,1,-1
             hcdo(i,k)=heso_cup(i,jmini)
             dz=zo_cup(i,k+1)-zo_cup(i,k)
             dh=dh+dz*(hcdo(i,k)-heso_cup(i,k))
             if(dh.gt.0.)then
               jmini=jmini-1
               if ( jmini .gt. 5 ) then
                 keep_going = .true.
               else
                 ierr(i)=9
                 ierrc(i) = "could not find jmini9"
                 exit
               endif
             endif
           enddo
         enddo
         jmin(i) = jmini
         if ( jmini .le. 5 ) then
           ierr(i)=4
           ierrc(i) = "could not find jmini4"
         endif
     enddo
     !
     ! - must have at least depth_min m between cloud convective base and cloud top.
     !
     do i=its,itf
         if(ierr(i) /= 0) cycle
         if ( jmin(i) - 1 .lt. kdet(i)) kdet(i) = jmin(i)-1
         if (-zo_cup(i,kbcon(i))+zo_cup(i,ktop(i)).lt.depth_min)then
             ierr(i)=6
             ierrc(i)="cloud depth very shallow"
         endif
     enddo

  END SUBROUTINE get_jmin
!------------------------------------------------------------------------------------
  SUBROUTINE precip_cwv_factor(itf,ktf,its,ite,kts,kte,ierr,t,po,qo,po_cup,cumulus,p_cwv_ave)
     implicit none
     character *(*), intent (in)                        :: cumulus
     integer  ,intent (in )                             :: itf,ktf, its,ite, kts,kte
     integer  ,intent (in ), dimension(its:ite)         :: ierr
     real     ,intent (in ), dimension(its:ite,kts:kte) :: t,po,qo,po_cup
     real     ,intent (out), dimension(its:ite)         :: p_cwv_ave

     !--locals
     integer :: i,k
     real    :: dp, trash
     real    ,dimension(its:ite) :: w_col,w_ccrit,t_troposph
     real, parameter :: fpkup=0.8  !-- 90% of precip occurs above 80% of critical w

     p_cwv_ave = 0.0
     if(trim(cumulus) /= 'deep') return
!
!-- get the pickup of ensemble ave prec, following Neelin et al 2009.
!
      do i=its,itf
          w_col     (i)= 0.
          w_ccrit   (i)= 0.
          t_troposph(i)= 0.
          if(ierr(i) /= 0) cycle
          trash=0.
loopN:    do k=kts,ktf
            if(po(i,k) .lt. 200.) exit loopN

            dp=100.*(po_cup(i,k)-po_cup(i,k+1))
            trash=trash+dp/g

            w_col     (i)= w_col     (i) + qo(i,k)*dp/g ! unit mm
            t_troposph(i)= t_troposph(i) + t (i,k)*dp/g
          enddo loopN
          !--average temperature
          t_troposph(i) =   t_troposph(i)/(1.e-8+trash)! unit K
          !
          !--- wcrit given by Neelin et al 2009.
          w_ccrit(i) = max(0.,56.2 + 2.1*(t_troposph(i)-268.)) ! unit mm
          !
          !--- pickup (normalized by the factor 'a')
          !-- <p>=a[(w-w_c)/w_c]**beta, a=0.15, beta=0.23
          !
          p_cwv_ave(i) = (max(0.,w_col(i)-fpkup*w_ccrit(i))/(1.e-8+fpkup*w_ccrit(i)))**0.23
          p_cwv_ave(i) = max(0., min(1., p_cwv_ave(i)))

          !print*,"NEE=",i,w_col(i),t_troposph(i),w_ccrit(i),p_cwv_ave    (i)
          !print*,"=================================================="
      enddo
  END SUBROUTINE precip_cwv_factor
!------------------------------------------------------------------------------------
SUBROUTINE get_wetbulb(jmin,qo_cup,t_cup,po_cup ,q_wetbulb,t_wetbulb)

 implicit none
 integer ,intent (in   ) :: jmin
 real    ,intent (in   ) :: qo_cup,t_cup,po_cup
 real    ,intent (inout) :: q_wetbulb,t_wetbulb

 !---locals
 real ::  zqp, zcond, zcond1, zcor, zqsat
 real :: psp, pt , pq
 real :: z3es,   z4es, z5alcp, zaldcp
 real :: ptare, evap
 real :: foedelta,foeewmcu,foealfcu,foedemcu,foeldcpmcu

 real, parameter :: &
   RD=287.06                             &
  ,RV=461.52                             &
  ,RCPD=1004.71                          &
  ,RTT=273.16                            &
  ,RHOH2O=1000.                          &
  ,RLVTT=2.5008E+6                       &
  ,RLSTT=2.8345E+6                       &
  ,RETV=RV/RD-1.0                        &
  ,RLMLT=RLSTT-RLVTT                     &
  ,RCPV=4.*RV                            &
  ,R2ES=611.21*RD/RV                     &
  ,R3LES=17.502                          &
  ,R3IES=22.587                          &
  ,R4LES=32.19                           &
  ,R4IES=-0.7                            &
  ,R5LES=R3LES*(RTT-R4LES)               &
  ,R5IES=R3IES*(RTT-R4IES)               &
  ,R5ALVCP=R5LES*RLVTT/RCPD              &
  ,R5ALSCP=R5IES*RLSTT/RCPD              &
  ,RALVDCP=RLVTT/RCPD                    &
  ,RALSDCP=RLSTT/RCPD                    &
  ,RALFDCP=RLMLT/RCPD                    &
  ,RTWAT=RTT                             &
  ,RTBER=RTT-5.                          &
  ,RTBERCU=RTT-5.0                       &
  ,RTICE=RTT-23.                         &
  ,RTICECU=RTT-23.                       &
  ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
  ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
  ,RVTMP2=RCPV/RCPD-1.                   &
  ,ZQMAX=0.5

  !-- for testing
  !              PSP                   TEMP                         Q                     ZCOND1
  ! input   85090.0000000000        289.140030372766     1.105078557441815E-002
  ! output  85090.0000000000        287.230570412846     1.181792062536557E-002 -2.761256206705639E-005
  ! PT  = 289.140030372766
  ! PQ  = 1.105078557441815E-002
  ! PSP = 85090.
  !----------------------

  !-- environmental values
   PT  = t_cup       ! K
   PQ  = qo_cup      ! kg/kg
   PSP = po_cup*100. ! hPa

   IF (PT > RTT) THEN
    Z3ES=R3LES
    Z4ES=R4LES
    Z5ALCP=R5ALVCP
    ZALDCP=RALVDCP
   ELSE
    Z3ES=R3IES
    Z4ES=R4IES
    Z5ALCP=R5ALSCP
    ZALDCP=RALSDCP
   ENDIF

   !--- get wet bulb thermo properties --------------------------
    PTARE = PT
    ZQP    =1.0/PSP

    FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
    FOEEWMCU = R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
                    (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
    ZQSAT=FOEEWMCU*ZQP

    ZQSAT=MIN(max_qsat,ZQSAT)
    ZCOR=1.0/(1.0-RETV  *ZQSAT)
    ZQSAT=ZQSAT*ZCOR

    FOEDEMCU =  FOEALFCU *R5ALVCP*(1.0/(PTARE-R4LES)**2)+&
              (1.0-FOEALFCU)*R5ALSCP*(1.0/(PTARE-R4IES)**2)

    ZCOND=(PQ-ZQSAT)/(1.0+ZQSAT*ZCOR*FOEDEMCU)

    ZCOND=MIN(ZCOND,0.0)

    FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
    PT=PT+FOELDCPMCU*ZCOND

    PQ=PQ-ZCOND

    !--update PTARE
    PTARE = PT

    FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
    FOEEWMCU = R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
                    (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
    ZQSAT=FOEEWMCU*ZQP

    ZQSAT=MIN(0.5,ZQSAT)
    ZCOR=1.0/(1.0-RETV  *ZQSAT)
    ZQSAT=ZQSAT*ZCOR

    FOEDEMCU =  FOEALFCU *R5ALVCP*(1.0/(PTARE-R4LES)**2)+&
              (1.0-FOEALFCU)*R5ALSCP*(1.0/(PTARE-R4IES)**2)
    ZCOND1=(PQ-ZQSAT)/(1.0+ZQSAT*ZCOR*FOEDEMCU)

    IF(ZCOND == 0.0)ZCOND1=MIN(ZCOND1,0.0)
    FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
    PT=PT+FOELDCPMCU*ZCOND1
    PQ=PQ-ZCOND1

   !-- set output --------------------------
    q_wetbulb =  PQ
    t_wetbulb =  PT
    evap      = -ZCOND1 != q_wetbulb-qo_cup, source for water vapor
END SUBROUTINE get_wetbulb
!------------------------------------------------------------------------------------
  SUBROUTINE cup_forcing_ens_3d_shal(itf,ktf,its,ite,kts,kte,dtime,ichoice     &
                                    ,ierrc,ierr,klcl,kpbl,kbcon,k22,ktop       &
                                    ,xmb,tsur,cape,h_sfc_flux,le_sfc_flux,zws  &
                                    ,po, hco, heo_cup,po_cup,t_cup,dhdt,rho    &
                                    ,xff_shal2d,xf_dicycle)

      implicit none
      integer                               ,intent (in) :: itf,ktf,its,ite, kts,kte,ichoice
      integer ,dimension (its:ite)          ,intent (in) :: klcl,kpbl,kbcon,k22,ktop
      real                                  ,intent (in) :: dtime
      real    ,dimension (its:ite)          ,intent (in) :: tsur,cape,h_sfc_flux,le_sfc_flux,zws
      real    ,dimension (its:ite,kts:kte)  ,intent (in) :: po,hco,heo_cup,po_cup,t_cup,dhdt,rho
      integer ,dimension (its:ite)          ,intent (inout):: ierr
      character*128,dimension (its:ite)     ,intent (inout):: ierrc
      real    ,dimension (its:ite)          ,intent (out)  :: xmb,xf_dicycle
      real    ,dimension (its:ite,9)        ,intent (out)  :: xff_shal2d

      !---local vars
      real   ,dimension (its:ite)    :: xmbmax
      integer :: i,k,kbase
      real    :: blqe,trash,tcold,fin,fsum,efic,thot,dp
      real    ,dimension (9)        :: xff_shal

      do i=its,itf
        xmb       (i)     = 0.
        xf_dicycle(i)     = 0.
        if(ierr(i) /= 0 ) cycle

        xmbmax(i)=100.*(po(i,kbcon(i))-po(i,kbcon(i)+1))/(g*dtime)

        !- limiting the mass flux at cloud base
        xmbmax(i)=min(xmbmaxshal,xmbmax(i))

        !- closure from Grant (2001)
        !- increased by ~125%
       !xff_shal(1)=.030*zws(i)*rho(i,kpbl(i))
        xff_shal(1)=.0675*zws(i)*rho(i,kpbl(i))
        xff_shal(2)=xff_shal(1)
        xff_shal(3)=xff_shal(1)
        !- cloud base
        kbase=kbcon(i)
       !kbase=klcl(i)

        !- closure from boundary layer QE (Raymond 1995)
        blqe=0.
        trash=0.
        if(k22(i).lt.kpbl(i)+1)then
             do k=kts,kbase
                blqe=blqe+100.*dhdt(i,k)*(po_cup(i,k)-po_cup(i,k+1))/g
             enddo
             trash = max((hco(i,kbase)-heo_cup(i,kbase)),1.e1)
             xff_shal(7)=max(0.,blqe/trash)
             !print*,"blqe=", xff_shal(7),blqe,trash
        else
             xff_shal(7)=0.0
        endif
        !print*,"blqe=", xff_shal(7),blqe,hco(i,kbcon(i))
        xff_shal(8)= xff_shal(7)
        xff_shal(9)= xff_shal(7)

        !- closure from the heat-engine principle
        !- Renno and Ingersoll(1996), Souza et al (1999)
        !- get the averaged environment temperature between cloud base
        !- and cloud top
        tcold=0.
        do k=kbase,ktop(i)
          dp   = po_cup(i,k)-po_cup(i,k+1)
          tcold= tcold + t_cup(i,k)*dp
        enddo
        tcold=tcold/(po_cup(i,kbase)-po_cup(i,ktop(i)+1))

        !-surface temperature
        thot=tsur(i)  ! + ztexec(i)
        !- thermodynamic eficiency
        !efic = max(0.05, (thot-tcold)/thot )
        efic = max(0.0, (thot-tcold)/thot )

        !- total heat flux from surface
        fin = max(0.0, h_sfc_flux(i)+le_sfc_flux(i))

        !- mass flux at cloud base
        !if(cape(i) > 0.0 .and. h_sfc_flux(i) >0.0 ) then
        if(cape(i) > 0.0  ) then
         xff_shal(4) = efic * fin / cape(i)
        else
         xff_shal(4) = 0.0
        endif
        xff_shal(5)=xff_shal(4)
        xff_shal(6)=xff_shal(4)


        xff_shal2d(i,:) = xff_shal(:)

       enddo
   END subroutine cup_forcing_ens_3d_shal

!------------------------------------------------------------------------------------
  SUBROUTINE cup_up_lightning(itf,ktf,its,ite, kts,kte, ierr, kbcon,ktop,xland,cape &
                             ,cnvfrc,srftype,zo,zo_cup,t_cup,t,tempco,qrco,po_cup,rho,prec_flx     &
                             ,lightn_dens)

   !=====================================================================================
   !- Lightning parameterization based on:
   !- "A Lightning Parameterization for the ECMWF Integrated Forecasting System"
   !-  P. Lopez, 2016 MWR
   !
   !- Coded/adapted to the GF scheme by Saulo Freitas (10-Aug-2019)
   !=====================================================================================
   implicit none
   integer                            ,intent(in)  :: itf,ktf, its,ite, kts,kte
   integer, dimension(its:ite)        ,intent(in)  :: ierr,kbcon,ktop
   real,    dimension(its:ite)        ,intent(in)  :: cape,xland,cnvfrc,srftype
   real,    dimension(its:ite,kts:kte),intent(in)  :: po_cup,zo_cup,t_cup,t,tempco,zo &
                                                     ,qrco,rho,prec_flx

   real,    dimension(its:ite)        ,intent(out) :: lightn_dens ! lightning flash density
                                                                  ! rate (units: 1/km2/day)

   !-- locals
   real, parameter :: V_graup     = 3.0  ! m/s
   real, parameter :: V_snow      = 0.5  ! m/s
   real, parameter :: beta_land   = 0.70 ! 1
   real, parameter :: beta_ocean  = 0.45 ! 1
   real, parameter :: alpha       = 37.5 ! 1
   real, parameter :: t_initial   =  0.0 + 273.15 ! K
   real, parameter :: t_final     = -25. + 273.15 ! K

   integer :: i, k, k_initial, k_final
   real    :: Q_R, z_base,beta,prec_flx_fr,dz
   real,    dimension(kts:kte) :: p_liq_ice, q_graup,q_snow

   do i=its,itf
      lightn_dens(i) = 0.0
      if(ierr(i) /= 0) cycle

      beta= xland(i)*beta_ocean + (1.-xland(i))*beta_land

      q_graup(:) = 0.
      q_snow (:) = 0.

      do k=kts,ktop(i)

        p_liq_ice(k) = fract_liq_f(tempco(i,k),cnvfrc(i),srftype(i))

        prec_flx_fr=   p_liq_ice(k)*prec_flx(i,k)/rho(i,k)

        q_graup(k) =      beta *prec_flx_fr/V_graup ! - graupel mixing ratio (kg/kg)
        q_snow (k) =  (1.-beta)*prec_flx_fr/V_snow  ! - snow    mixing ratio (kg/kg)

      enddo

      k_initial = minloc(abs(tempco(i,kbcon(i):ktop(i))-t_initial),1)+kbcon(i)-1
      k_final   = minloc(abs(tempco(i,kbcon(i):ktop(i))-t_final  ),1)+kbcon(i)-1

      Q_R = 0.0
      do k = k_initial, k_final
        dz  = zo(i,k)-zo(i,k-1)
        Q_R = Q_R + dz*rho(i,k)*(q_graup(k)*(qrco(i,k)+q_snow(k)))
        !print*,"qr=",q_r,tempco(i,k)-273.15,k,tempco(i,k)-t_initial
      enddo

      z_base = zo_cup(i,kbcon(i))/1000. ! km

     !---
     !--- lightning flash density (units: number of flashes/km2/day) - equation 5
     !--- (to compare with Lopez 2016's results, convert to per year: lightn_dens*365)
     !
      lightn_dens(i) = alpha * Q_R *sqrt (max(0.,cape(i))) * min(z_base,1.8)**2
     !
   enddo
  END SUBROUTINE cup_up_lightning

!------------------------------------------------------------------------------------
   SUBROUTINE cup_up_rain(cumulus,klcl,kbcon,ktop,k22,ierr,xland,cnvfrc,srftype &
                      ,zo_cup,qco,qrco,pwo,pwavo,po,p_cup,t_cup,tempco  &
                      ,zuo,up_massentr,up_massdetr,vvel2d,rho           &
                      ,qrr                                              &
                      ,itf,ktf,its,ite, kts,kte)

     implicit none
     character *(*)            , intent (in) :: cumulus
     integer                    ,intent (in) :: itf,ktf,its,ite,kts,kte
     integer, dimension(its:ite),intent (in) :: kbcon,ktop,k22,klcl,ierr
     real,    dimension(its:ite),intent (in) :: xland,cnvfrc,srftype,pwavo
     real,    dimension(its:ite,kts:kte),intent (in)  ::   &
              zo_cup,qco,qrco,pwo,po,p_cup,t_cup,zuo       &
             ,up_massentr,up_massdetr,vvel2d,tempco,rho

    !--for future use (rain water mixing ratio)
     real,    dimension(its:ite,kts:kte),intent (out) :: qrr      !-- units kg[water]/kg[air]

    !-- locals
     integer :: i,k
     real :: tmp
     integer, dimension(its:ite) :: start_level
     real :: dz,z1,zrnew,zc,zd,zint,z2,zrold,denom,fall_fact,wup, exp1,R_vr
     real,    dimension(kts:kte) :: prec_flx_rain,prec_flx_snow
     real,    dimension(kts:kte) :: pw,p_liq_ice ! - rainfall source
     real,    parameter :: rho1000mb = 1.2 , rhow = 1000., N_r = 0.1 ! cm^-3, rainfall drops number concen
     real,    parameter :: exp_KR    = 1./5. &! Kuo & Raymond 1983
                          ,exp_SB    = 2./3.  ! Seifert & Beheng 2006 eq 27

     qrr = 0.
     if(trim(cumulus) == 'shallow') return

     !--- rain water mixing ratio
     do i=its,itf
       if (ierr(i)  /= 0) cycle

       do k=ktop(i),kts,-1

           p_liq_ice(k) = fract_liq_f(tempco(i,k),cnvfrc(i),srftype(i))

           !--- transport + mixing
           denom = zuo(i,k+1)-.5*up_massdetr(i,k)+up_massentr(i,k)
           if(denom > 0.) then
               qrr(i,k) = (qrr(i,k+1)*zuo(i,k+1)-.5*up_massdetr(i,k)* qrr(i,k+1))/ denom
           else
               qrr(i,k) =  qrr(i,k+1)
           endif

          !--- rain source
           pw(k)= pwo(i,k)/(1.e-16+zuo(i,k))

          !-- rainfall sedimentation
          !-- Kuo & Raymond 1983
          !-- fallout of rain (21.18 * qrr^0.2 have m/s as units with qrr in kg/kg)
          !---                                half velocity for ice phase
           fall_fact = 21.18 *( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
           exp1      = exp_KR
          !
          !-- Seifert & Beheng 2006 eq 27, units= m/s
          ! fall_fact = 159.*sqrt(rho1000mb/rho(i,k))*( p_liq_ice(k) + 0.5*(1.-p_liq_ice(k) ))
          ! exp1      = exp_SB

          !-- Kogan 2013
          R_vr = (4.*MAPL_PI*rhow/(3.*rho(i,k)))**(-1./3.) * (qrr(i,k) + pw(k))**(1./3.) * &
                 (N_r)**(-1./3.)! cm, N_r is the rainfall drops number concentration, and not CCN or N_c

          R_vr = max(40., R_vr * 1.e-2 * 1.e+6 )   ! micrometer
          fall_fact = 1.e-2*(2.4*R_vr-62.0)        ! m/s
          exp1      = 0.

           wup      = min(15.,max(2.,vvel2d(i,k)))
           z2       = fall_fact/wup

           !--exact solution
           if(qrr(i,k) + pw(k)>0.) then
            !-- this is the sedimentation speed divided by W_up
            zd= z2* (qrr(i,k) + pw(k))**exp1
            zint= EXP(-zd)
            zrold=qrr(i,k)
            zc=pw(k)
            zrnew=zrold*zint+zc/zd*(1.-zint)
            zrnew=max(0.0,min(qrr(i,k) + pw(k),zrnew))
           else
            zrnew=0.
           endif
           qrr(i,k)=zrnew

           !--solution with rk3
           !z1= qrr(i,k)
           !z1= qrr(i,k) +  1./3.*(pw(k) - z2*(z1)**exp1 * z1)
           !z1= qrr(i,k) +  0.500*(pw(k) - z2*(z1)**exp1 * z1)
           !z1= qrr(i,k) +  1.000*(pw(k) - z2*(z1)**exp1 * z1)
           !qrr(i,k)=z1
           !---
           !--- solution 2
           !qrr(i,k)= qrr(i,k) + pw(k) - (fall_fact*(qrr(i,k) + pw(k))**exp1) / wup * qrr(i,k)

           !--- solution 3
           !tmp = 0.5*(qrr(i,k) + pw(k) + qrr(i,k+1) + pw(k+1))
           !qrr(i,k)= qrr(i,k) + pw(k) - z2*(tmp)**exp1 * 0.5*(qrr(i,k)+qrr(i,k+1))

           !print*,"rr3=",k,zo_cup(i,k), pwo(i,k)*1000.,qrr(i,k)*1000.,fall_fact*(qrr(i,k) + pw(k))**exp1
          ENDDO
     ENDDO

   END SUBROUTINE cup_up_rain

!------------------------------------------------------------------------------------
SUBROUTINE get_condensation(q_old,t_old,po_cup,q_new,t_new)

!-- calculate condensation and adjust t and q accordingly
 implicit none
 real    ,intent (in   ) :: po_cup,q_old,t_old ! before condensation
 real    ,intent (inout) ::        q_new,t_new ! after  condensation

 !---locals
 real ::  zqp, zcond, zcond1, zcor, zqsat,zi,zl,zf
 real :: psp, pt , pq
 real :: z3es,   z4es, z5alcp, zaldcp
 real :: ptare, cond
 real :: foeewmcu,foealfcu,foedemcu,foeldcpmcu

 real, parameter :: &
   RD=287.06                             &
  ,RV=461.52                             &
  ,RCPD=1004.71                          &
  ,RTT=273.16                            &
  ,RHOH2O=1000.                          &
  ,RLVTT=2.5008E+6                       &
  ,RLSTT=2.8345E+6                       &
  ,RETV=RV/RD-1.0                        &
  ,RLMLT=RLSTT-RLVTT                     &
  ,RCPV=4.*RV                            &
  ,R2ES=611.21*RD/RV                     &
  ,R3LES=17.502                          &
  ,R3IES=22.587                          &
  ,R4LES=32.19                           &
  ,R4IES=-0.7                            &
  ,R5LES=R3LES*(RTT-R4LES)               &
  ,R5IES=R3IES*(RTT-R4IES)               &
  ,R5ALVCP=R5LES*RLVTT/RCPD              &
  ,R5ALSCP=R5IES*RLSTT/RCPD              &
  ,RALVDCP=RLVTT/RCPD                    &
  ,RALSDCP=RLSTT/RCPD                    &
  ,RALFDCP=RLMLT/RCPD                    &
  ,RTWAT=RTT                             &
  ,RTBER=RTT-5.                          &
  ,RTBERCU=RTT-5.0                       &
  ,RTICE=RTT-23.                         &
  ,RTICECU=RTT-23.                       &
  ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
  ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
  ,RVTMP2=RCPV/RCPD-1.                   &
  ,ZQMAX=0.5

  !----------------------
  !-- for testing
  !----------------------
  !                          PSP (hPA)            TEMP(K)              Q(kg/kg)                ZCOND1(kg/kg)
  ! input          1   98020.0000000000        295.163188640152     1.745679200989956E-002
  ! output         1   98020.0000000000        295.513490789916     1.731605637618801E-002 -9.779916453243843E-007
  !----------------------
  ! input       157   85090.0000000000        288.089188935407     1.399404805052166E-002
  ! output      157   85090.0000000000        289.294751760460     1.350970717999820E-002 -1.146268822756454E-005
  !----------------------
  ! PT  = 288.089188935407
  ! PQ  = 1.399404805052166E-002
  ! PSP = 85090
  !----------------------

  !-- initial values
    PT  = t_old       ! K
    PQ  = q_old       ! kg/kg
    PSP = po_cup*100. ! hPa


   !--- get condensation in moist ascent --------------------------
    PTARE = PT
    ZQP    =1.0/PSP

    ZL=1.0/(PT-R4LES)
    ZI=1.0/(PT-R4IES)

    FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
    ZQSAT=R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)*ZL)+&
               (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)*ZI))

    ZQSAT=ZQSAT*ZQP
    ZQSAT=MIN(0.5,ZQSAT)
    ZCOR=1.0-RETV*ZQSAT

    ZF=FOEALFCU*R5ALVCP*ZL**2 + (1.0-FOEALFCU)*R5ALSCP*ZI**2
    ZCOND=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)

    IF(ZCOND > 0.0)THEN
      FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
      PT=PT+FOELDCPMCU*ZCOND
      PTARE = PT
      PQ=PQ-ZCOND

      ZL=1.0/(PT-R4LES)
      ZI=1.0/(PT-R4IES)


      FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
      ZQSAT=R2ES *(FOEALFCU* EXP(R3LES*(PT-RTT)*ZL)+&
                 (1.0-FOEALFCU)*EXP(R3IES*(PT-RTT)*ZI))

      ZQSAT=ZQSAT*ZQP
      ZQSAT=ZQSAT-0.5*(ABS(0.5-ZQSAT)-(0.5-ZQSAT))


      ZCOR=1.0-RETV*ZQSAT
      ZF=FOEALFCU*R5ALVCP*ZL**2 + (1.0-FOEALFCU)*R5ALSCP*ZI**2

      ZCOND1=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
      IF(ZCOND ==  0.0)ZCOND1=0.0
      FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
      PT=PT+FOELDCPMCU*ZCOND1
      PQ=PQ-ZCOND1
    ENDIF

   !-- FINAL --------------------------
    q_new =  PQ
    t_new =  PT
    cond  = -ZCOND1 != q_old-qnew, source for the liquid water
END SUBROUTINE get_condensation

!------------------------------------------------------------------------------------
SUBROUTINE get_interp(q_old,t_old,po_cup,q_new,t_new)
   implicit none
   real    ,intent (in   ) :: po_cup ! original
   real    ,intent (inout) :: q_old,t_old,q_new,t_new ! extrapolated

   !---locals
   real ::  zqp, zcond1, zcor, zqsat
   real ::  psp, pt , pq, ptare
   real ::  FOEALFCU, FOEEWMCU,FOEDEMCU,FOELDCPMCU
   real, parameter :: &
   RD=287.06                             &
  ,RV=461.52                             &
  ,RCPD=1004.71                          &
  ,RTT=273.16                            &
  ,RHOH2O=1000.                          &
  ,RLVTT=2.5008E+6                       &
  ,RLSTT=2.8345E+6                       &
  ,RETV=RV/RD-1.0                        &
  ,RLMLT=RLSTT-RLVTT                     &
  ,RCPV=4.*RV                            &
  ,R2ES=611.21*RD/RV                     &
  ,R3LES=17.502                          &
  ,R3IES=22.587                          &
  ,R4LES=32.19                           &
  ,R4IES=-0.7                            &
  ,R5LES=R3LES*(RTT-R4LES)               &
  ,R5IES=R3IES*(RTT-R4IES)               &
  ,R5ALVCP=R5LES*RLVTT/RCPD              &
  ,R5ALSCP=R5IES*RLSTT/RCPD              &
  ,RALVDCP=RLVTT/RCPD                    &
  ,RALSDCP=RLSTT/RCPD                    &
  ,RALFDCP=RLMLT/RCPD                    &
  ,RTWAT=RTT                             &
  ,RTBER=RTT-5.                          &
  ,RTBERCU=RTT-5.0                       &
  ,RTICE=RTT-23.                         &
  ,RTICECU=RTT-23.                       &
  ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
  ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
  ,RVTMP2=RCPV/RCPD-1.                   &
  ,ZQMAX=0.5

  integer :: i

    PT  = t_old       ! K
    PQ  = q_old       ! kg/kg
    PSP = po_cup*100. ! hPa

  !-- for testing
  !              PSP                   TEMP                         Q                     ZCOND1
  ! input    27940.0000000000        236.604976804749       3.220181796223121E-004
  ! output   27940.0000000000        236.361132108860       4.084506812610067E-004
  !  PT  = 236.604976804749      ! K
  !  PQ  = 3.220181796223121E-004       ! kg/kg
  !  PSP = 27940. ! hPa
  !----------------------
  !print*,"1",PSP,PT,PQ

   ZQP   =1.0/PSP
   do i=1,2
    PTARE = PT

    FOEALFCU = MIN(1.0,((MAX(RTICECU,MIN(RTWAT,PTARE))-RTICECU)*RTWAT_RTICECU_R)**2)
    FOEEWMCU = R2ES *(FOEALFCU *EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
                    (1.0-FOEALFCU)*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))
    ZQSAT=FOEEWMCU*ZQP

!    if(1.0-RETV  *ZQSAT == 0.) then
!
!      print*,"ZQSAT=",ZQP,FOEEWMCU,q_old,t_old,po_cup,q_new,t_new
!3.5491847E-02   46.36052      0.5000000       249.8219
!  0.2817549      0.5000000       249.8219
!      call flush(6)
!      stop 3333
!    endif

    ZCOR=1.0/(1.0-RETV  *ZQSAT)
    ZQSAT=ZQSAT*ZCOR

    FOEDEMCU =  FOEALFCU *R5ALVCP*(1.0/(PTARE-R4LES)**2)+&
              (1.0-FOEALFCU)*R5ALSCP*(1.0/(PTARE-R4IES)**2)


    ZCOND1=(PQ-ZQSAT)/(1.0+ZQSAT*ZCOR*FOEDEMCU)

    FOELDCPMCU= FOEALFCU*RALVDCP+(1.0-FOEALFCU)*RALSDCP
    PT=PT+FOELDCPMCU*ZCOND1
    PQ=PQ-ZCOND1
   enddo
   !-- FINAL --------------------------
   q_new =  PQ
   t_new =  PT
   !print*,"2",PSP,PT,PQ
   !print*,"E",100*(PT-236.361132108860)/236.361132108860,100*(PQ-4.084506812610067E-004)/4.084506812610067E-004
end SUBROUTINE get_interp

!------------------------------------------------------------------------------------
REAL FUNCTION fract_liq_f(temp2,cnvfrc,srftype) ! temp2 in Kelvin, fraction between 0 and 1.
   implicit none
   real,intent(in)  :: temp2 ! K
   real,intent(in)  :: cnvfrc,srftype
   SELECT CASE(FRAC_MODIS)
   CASE (1)
       fract_liq_f = 1.0 - ice_fraction(temp2,cnvfrc,srftype)
   CASE DEFAULT
       fract_liq_f =  min(1., (max(0.,(temp2-t_ice))/(t_0-t_ice))**2)
   END SELECT
 END FUNCTION
!------------------------------------------------------------------------------------

 subroutine prepare_temp_pertubations(kts,kte,ktf,its,ite,itf,jts,jte,jtf,dt,xland,topo,zm &
                                     ,temp,rqvften,rthblten,rthften,Tpert_h,Tpert_v&
                                     ,AA1_CIN,AA1_BL)!<<<<<<<<<<
      implicit none
      integer ,intent(in) :: kts,kte,ktf,its,ite,itf,jts,jte,jtf
      real    ,intent(in) :: dt

      real    ,dimension(its:ite,jts:jte)        , intent(in) :: xland   !flag < 1 para land
                                                                         !flag  =1 para water
      real    ,dimension(its:ite,jts:jte)        , intent(in) :: topo

      real    ,dimension(kts:kte,its:ite,jts:jte), intent(in) ::        &
                                                               zm       &
                                                              ,rthften  &
                                                              ,rqvften  &
                                                              ,rthblten &
                                                              ,temp

      real    ,dimension(kts:kte,its:ite,jts:jte), intent(out)::Tpert_h,Tpert_v


      REAL   ,DIMENSION(its:ite,jts:jte)  ,INTENT(INOUT)  ::AA1_CIN,AA1_BL!<<<<<<<<<<

     !-- local vars
      integer :: i,j,k,kr,i1,i2,i3,j1,j2,j3
      real    :: aveh_qmax,aveh_qmin,avev_qmax,avev_qmin,coef_h,coef_v,dz1,dz2,dz3
      real, dimension(kts:kte)                  :: avev_adv_q,avev_t
      real, dimension(kts:kte,its:ite,jts:jte)  :: ave_T,ave_Q_adv
      integer, DIMENSION(its:ite,jts:jte)       :: is_ocean

        !--avoid the trigger over the land areas
        is_ocean(:,:) = 0
        do j = jts,jtf
          do i= its,itf
           if(xland(i,j) > 0.999 .and. topo(i,j) < 10.) is_ocean(i,j) = 1
        enddo;enddo

        !-reset/initialization
        Tpert_h=0.; Tpert_v=0. ; ave_T=0. ; ave_Q_adv=0. ;avev_adv_q=0.; avev_adv_q =0. ; avev_t=0.

        !-- calculate 9-point average of moisture advection and temperature using halo (Horizontal)
        !-- in the future these lines must be moved to the dynamics, because there is not "halo" in physics.
        do j = jts,jtf
          j1 = max (j-1,jts) ; j2 = j ; j3 = min(j+1,jtf)
          do i= its,itf
             i1 = max (i-1,its) ; i2 = i ; i3 = min(i+1,itf)
             if(is_ocean(i,j) == 0) cycle  ! do it only over ocean regions

             do k = kts, ktf
               kr = k

               !--think about using temp _OR_ temp_new(kr,i,j)= temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
               ave_T     (k,i,j) = ( temp(kr,i1,j1) + temp(kr,i1,j2) + temp(kr,i1,j3) + & ! row 1
                                     temp(kr,i2,j1) + temp(kr,i2,j2) + temp(kr,i2,j3) + & ! row 2
                                     temp(kr,i3,j1) + temp(kr,i3,j2) + temp(kr,i3,j3)   & ! row 3
                                   ) / 9.
              !-- advection forcing of q
               ave_Q_adv (k,i,j) = ( rqvften(kr,i1,j1) + rqvften(kr,i1,j2) + rqvften(kr,i1,j3) + & ! row 1
                                     rqvften(kr,i2,j1) + rqvften(kr,i2,j2) + rqvften(kr,i2,j3) + & ! row 2
                                     rqvften(kr,i3,j1) + rqvften(kr,i3,j2) + rqvften(kr,i3,j3)   & ! row 3
                                   ) / 9.

        enddo;enddo;enddo

        !--search for max/min moisture advection in 9-point range, calculate horizontal T-perturbation (tpert_h)

        do j = jts,jtf
           j1 = max (j-1,jts) ; j2 = j ; j3 = min(j+1,jtf)
           do i= its,itf
              i1 = max (i-1,its) ; i2 = i ; i3 = min(i+1,itf)

              if(is_ocean(i,j) == 0) cycle

              do k = kts, ktf
                 kr = k
                 aveh_qmax = maxval( ave_Q_adv(k,i1:i3,j1:j3) )
                 aveh_qmin = minval( ave_Q_adv(k,i1:i3,j1:j3) )

                 if(aveh_qmax > aveh_qmin )then
                    coef_h = (ave_Q_adv(k,i,j) - aveh_qmin)/(aveh_qmax-aveh_qmin)
                 else
                    coef_h = 0.
                 endif
                 coef_h=min(coef_h,1.0);coef_h=max(coef_h,0.0)

                 !--think about using temp _OR_ temp_new(kr,i,j)= temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                 Tpert_h(k,i,j)=coef_h *(temp(kr,i,j)-ave_T(k,i,j))

         enddo;enddo;enddo


        !--search for max/min moisture advection in 3-point vertical range, calculate vertical T-perturbation (tpert_v)
        do j = jts,jtf
           do i= its,itf
              if(is_ocean(i,j) == 0) cycle

              do k = kts+1,ktf-2
                 kr=k
                 dz1 = zm(kr  ,i,j) -  zm(kr-1,i,j)
                 dz2 = zm(kr+1,i,j) -  zm(kr  ,i,j)
                 dz3 = zm(kr+2,i,j) -  zm(kr+1,i,j)

                 !--think about using temp _OR_ temp_ne(kr,i,j)w=temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                 avev_t     (k) = ( dz1*temp   (kr-1,i,j) + dz2*temp   (kr,i,j)+ dz3*temp   (kr+1,i,j) )/ (dz1+dz2+dz3)
                 avev_adv_q (k) = ( dz1*rqvften(kr-1,i,j) + dz2*rqvften(kr,i,j)+ dz3*rqvften(kr+1,i,j) )/ (dz1+dz2+dz3)

              enddo
              avev_t     (kts)       =  avev_t     (kts+1)
              avev_adv_q (kts)       =  avev_adv_q (kts+1)
              avev_t     (ktf-1:ktf) =  avev_t     (ktf-2)
              avev_adv_q (ktf-1:ktf) =  avev_adv_q (ktf-2)


              do k = kts+1,ktf-2
                 kr=k
                 avev_qmax = maxval( avev_adv_q (k-1:k+1) )
                 avev_qmin = minval( avev_adv_q (k-1:k+1) )
                 if(avev_qmax > avev_qmin) then
                     coef_v = (avev_adv_q(k)-avev_qmin)/(avev_qmax-avev_qmin)
                 else
                     coef_v = 0.
                 endif

                 !--think about using temp _OR_ temp_new(kr,i,j)=temp(kr,i,j) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
                 Tpert_v(k,i,j)=coef_v*( temp(kr,i,j)-avev_t(k) )

              enddo
              Tpert_v(kts,i,j)=  Tpert_v(kts+1,i,j)
              Tpert_v(ktf,i,j)=  Tpert_v(ktf-1,i,j)

        enddo;enddo

        !--avoid the trigger over the land areas
        do j = jts,jtf
           do i= its,itf
                if(is_ocean(i,j) == 1) cycle  ! ocean areas
                 do k=kts,ktf
                  Tpert_v(k,i,j)= 0.! Tpert_v(k,i,j) * xland(i,j)
                  Tpert_h(k,i,j)= 0.! Tpert_h(k,i,j) * xland(i,j)
                enddo

              !print*,"TperH=",i,j,whoami_all,maxval(Tpert_h(:,i,j)),minval(Tpert_h(:,i,j))
              !print*,"TperV=",i,j,whoami_all,maxval(Tpert_V(:,i,j)),minval(Tpert_V(:,i,j))
        enddo;enddo

!-----
        return
!-----

      !---- check balance
        do j = jts,jtf
           do i= its,itf
              dz2=0.
              AA1_BL (i,j)=0.
              AA1_CIN(i,j)=0.

              do k = kts,13 !-- 13 ~ 900 hPa
                dz1 = zm(k+1,i,j) -  zm(k,i,j)
                dz2 = dz2 + dz1
                AA1_BL  (i,j)=AA1_BL  (i,j)+dz1*Tpert_v(k,i,j)
                AA1_cin (i,j)=AA1_cin (i,j)+dz1*Tpert_h(k,i,j)
              enddo
                AA1_BL  (i,j)=AA1_BL  (i,j)/(dz2+1.e-16)
                AA1_cin (i,j)=AA1_cin (i,j)/(dz2+1.e-16)
        enddo;enddo


     ! New trigger function
     ! Na vertical, fa\E7a isto apenas dentro da camada de 60 hPa entorno do k22.
     ! IF(trigger.eq.2) then
     !         DTLCL=amax1(tpart_h1D(KLCL)+tpart_v1D(KLCL),0.0)
     !  ENDIF

 end subroutine prepare_temp_pertubations
!------------------------------------------------------------------------------------
 SUBROUTINE SOUND(part,cumulus,int_time,dtime,ens4,itf,ktf,its,ite, kts,kte,xlats,xlons,jcol,whoami_all  &
                ,z ,qes ,he ,hes ,t ,q ,po,z1 ,psur,zo,qeso,heo,heso,tn,qo,us,vs ,omeg,xz &
                ,h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec,zqexec, xland &
                ,kpbl,k22,klcl,kbcon,ktop,aa0,aa1,sig,xaa0,hkb,xmb,pre,edto              &
                ,zo_cup,dhdt,rho,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv)

       implicit none
       character*(*) ,intent(in)    :: cumulus
       integer, intent(in) ::ens4, itf,ktf,its,ite, kts,kte,jcol,whoami_all,part
       real   , intent(in) ::int_time,dtime
       real,    dimension (its:ite,kts:kte) :: z ,qes ,he ,hes ,t ,q ,po,zo,qeso,heo,heso,tn,qo &
                                              ,us,vs ,xz,zo_cup,dhdt,rho &
                                              ,zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv
       real,    dimension (its:ite,kts:kte,1:ens4) :: omeg

       integer, dimension (its:ite) ::kpbl,k22,klcl,kbcon,ktop
       real,    dimension (its:ite) ::h_sfc_flux,le_sfc_flux,tsur, dx,stochastic_sig,zws,ztexec &
                                     ,zqexec,xlats,xlons,xland,z1,psur
       real,    dimension (its:ite) ::aa0,aa1,xaa0,hkb,xmb,pre,edto,sig

       !---locals
       integer :: i,k,X_kte,X_i,X_jcol,X_k,X_WHOAMI_ALL
       real :: X_time
       character(200) :: lixo
       real,    dimension (its:ite) :: X_stochastic_sig,X_xland
       real, parameter :: LATSND=-10., LONSND= 301., DELTX=0.2
!       real, parameter :: LATSND= -8.72, LONSND= 186.6, DELTX=0.2

        IF(trim(rundata)=="NONE") THEN

         IF( mod(int_time,3600.) < dtime ) THEN
         open(15, file="dataLXXX.dat_"//trim(cumulus),status='unknown',position="APPEND")

         IF(part == 1) THEN

          do i=its,itf

             if( xlats(i) > LATSND-DELTX .and. xlats(i) < LATSND+DELTX ) then
              if(xlons(i) > LONSND-DELTX .and. xlons(i) < LONSND+DELTX ) then

               print*,"==============================================="
               print*,"00>", i,jcol,xlats(i),xlons(i),whoami_all,int_time/3600.
               call flush(6)

               write(15,*) "====begin====="
               write(15,*) "i,jcol,xlats(i),xlons(i),int_time/3600."
               write(15,*)  i,jcol,xlats(i),xlons(i),int_time/3600.

               write(15,*) "kte,z1(i),psur(i),tsur(i),xland(i)"
               write(15,*)  kte,z1(i),psur(i),tsur(i),xland(i)

               write(15,*) "h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)"
               write(15,*)  h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)

               write(15,*) "stochastic_sig(i), dx(i),zws(i),kpbl(i)"
               write(15,*)  stochastic_sig(i), dx(i),zws(i),kpbl(i)

               write(15,*) "=>k zo po t tn-t q qo-q us vs qes he hes qeso-qes heo-he heso-hes dhdt omeg"
               do k = kts,kte
                write(15,100)    k,zo(i,k), po(i,k) ,t  (i,k) ,tn(i,k)-t(i,k), q(i,k) ,qo(i,k)-q(i,k)  &
                                 ,us (i,k), vs(i,k) ,qes(i,k) ,he(i,k) ,hes(i,k) &
                                 ,qeso(i,k)- qes(i,k) &
                                 ,heo (i,k)- he (i,k) &
                                 ,heso(i,k)- hes(i,k) &
                                 ,dhdt(i,k), omeg(i,k,1:ens4)
               enddo


             endif;endif
          enddo

         ELSE

           do i=its,itf
             if(xlats(i)  > LATSND-DELTX .and. xlats(i) < LATSND+DELTX ) then
              if(xlons(i) > LONSND-DELTX .and. xlons(i) < LONSND+DELTX ) then

                write(15,*) "====outputs======="
                write(15,*) "L=",i,jcol,xlats(i),xlons(i),whoami_all
                write(15,*) "A=",aa0(i),aa1(i),xaa0(i),sig(i)
                write(15,*) "K=",k22(i),klcl(i),kpbl(i),kbcon(i),ktop(i)
                write(15,*) "Z=",zo_cup(i,k22(i))-z1(i),zo_cup(i,klcl(i))-z1(i),zo_cup(i,kpbl(i))-z1(i)&
                                ,zo_cup(i,kbcon(i))-z1(i),zo_cup(i,ktop(i))-z1(i)
                write(15,*) "H=",hkb(i)/cp,edto(i)
                write(15,*) "T=",maxval(outt(i,1:ktop(i)))*86400.,maxval(outq(i,1:ktop(I)))*86400.*1000.,&
                                 minval(outt(i,1:ktop(i)))*86400.,minval(outq(i,1:ktop(I)))*86400.*1000.
                write(15,*) "P=",xmb(i)*1000.,'g/m2/s',3600*pre(i),'mm/h'
                if(xmb(i)>0.0) then
                  write(15,*) "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                  do k = kts,kte
                    write(15,101) k,zo(i,k), po(i,k) &
                               ,zuo(i,k),zdo(i,k),up_massentro(i,k),up_massdetro(i,k),outt(i,k)*86400. &
                               ,outq(i,k)*86400.*1000.,outqc(i,k)*86400.*1000.,outu(i,k)*86400.,outv(i,k)*86400.

                  enddo
                endif
                write(15,*) "=====end=========="
             endif;endif
           enddo
         ENDIF
         close(15)
         ENDIF
        ELSE

         IF(part == 1) THEN

            open(15, file=trim(rundata),status='old')
            i=1
            read(15,*) lixo
            read(15,*) lixo
            read(15,*) X_i,X_jcol,xlats(i),xlons(i),X_TIME
            read(15,*) lixo
            read(15,*) X_kte,z1(i),psur(i),tsur(i),X_xland(i)
            !-- check
            if(X_kte .ne. kte) stop " X_kte .ne. kte "
            read(15,*) lixo
            read(15,*) h_sfc_flux(i),le_sfc_flux(i),ztexec(i),zqexec(i)
            read(15,*) lixo
            read(15,*) X_stochastic_sig(i), dx(i),zws(i),kpbl(i)
              read(15,*) lixo
            do k = kts,kte
               read(15,100)  X_k,  zo(i,k), po(i,k) ,t  (i,k) ,tn  (i,k) ,q  (i,k) ,qo  (i,k), us  (i,k) ,vs(i,k) &
                                , qes(i,k) ,he(i,k) ,hes(i,k) ,qeso(i,k) ,heo(i,k) ,heso(i,k), dhdt(i,k) &
                                ,omeg(i,k,1:ens4)
            enddo
            close(15)
            !---settings
             tn    (i,:) = t  (i,:) + tn   (i,:) ! input is delta(T)
             qo    (i,:) = q  (i,:) + qo   (i,:) ! input is delta(Q)
             qeso  (i,:) = qes(i,:) + qeso (i,:) ! input is delta(Q)
             heo   (i,:) = he (i,:) + heo  (i,:) ! input is delta(H)
             heso  (i,:) = hes(i,:) + heso (i,:) ! input is delta(HO)
             xz    (i,:) = zo (i,:)
             z     (i,:) = zo (i,:)
             rho   (i,:) = 1.e2*po(i,:)/(rgas*t(i,:))

          ELSE

           do i=its,itf
             if(xlats(i)  > LATSND-DELTX .and. xlats(i) < LATSND+DELTX ) then
              if(xlons(i) > LONSND-DELTX .and. xlons(i) < LONSND+DELTX ) then

                print*, "====outputs======="
                print*, "A=",aa0(i),aa1(i),xaa0(i),sig(i)
                print*, "K=",k22(i),klcl(i),kpbl(i),kbcon(i),ktop(i)
                print*, "Z=",zo_cup(i,k22(i))-z1(i),zo_cup(i,klcl(i))-z1(i),zo_cup(i,kpbl(i))-z1(i)&
                            ,zo_cup(i,kbcon(i))-z1(i),zo_cup(i,ktop(i))-z1(i)
                print*, "H=",hkb(i)/cp,edto(i)
                print*, "T=",maxval(outt(i,1:ktop(i)))*86400.,maxval(outq(i,1:ktop(I)))*86400.*1000.,&
                             minval(outt(i,1:ktop(i)))*86400.,minval(outq(i,1:ktop(I)))*86400.*1000.
                print*, "P=",xmb(i)*1000.,'g/m2/s',3600*pre(i),'mm/h'
                if(xmb(i)>0.0) then
                  print*, "=> k zo po zuo,zdo,up_massentro,up_massdetro,outt, outq,outqc,outu,outv"
                  do k = kts,kte
                      write(*,101) k,zo(i,k), po(i,k) &
                                 ,zuo(i,k),zdo(i,k),up_massentro(i,k),up_massdetro(i,k),outt(i,k)*86400. &
                                   ,outq(i,k)*86400.*1000.,outqc(i,k)*86400.*1000.,outu(i,k)*86400.,outv(i,k)*86400.
                  enddo
                endif
             endif;endif
           enddo
         ENDIF
       ENDIF
    100 format(1x,i4,16E16.8)
    101 format(1x,i4,11E16.8)

 END SUBROUTINE SOUND
!------------------------------------------------------------------------------------
 SUBROUTINE cloud_dissipation(cumulus,itf,ktf, its,ite, kts,kte,ierr,kbcon,ktop,dtime,xmb,xland&
                             ,qo_cup,qeso_cup,po_cup,outt,outq,outqc,zuo,vvel2d,rho_hydr       &
                             ,qrco,sig,tempco,qco,tn_cup,heso_cup,zo)

 implicit none
 character*(*)                      ,intent(in)    :: cumulus
 integer                            ,intent(in)    :: itf,ktf, its,ite, kts,kte
 real                               ,intent(in)    :: dtime
 integer, dimension(its:ite)        ,intent(in)    :: ierr,kbcon,ktop
 real,    dimension(its:ite)        ,intent(in)    :: xmb,xland,sig
 real,    dimension(its:ite,kts:kte),intent(in)    :: po_cup,qo_cup,qeso_cup,zuo,vvel2d,rho_hydr&
                                                     ,tempco,qco,tn_cup,heso_cup,zo
 real,    dimension(its:ite,kts:kte),intent(inout) :: outt,outq,outqc,qrco

 !-- locals
 integer :: i,k
 real    ::  del_t,del_q,dp,frh
 real    :: qrc_diss,fractional_area,outqc_diss,outq_mix,outt_diss,outt_mix,tempx,qvx
 real, parameter :: cloud_lifetime= 1800.

 integer, parameter :: versionx = 2
 do i=its,itf

     if(ierr(i) /= 0) cycle

     do k=ktop(i),kbcon(i),-1

       !--- cloud liq/ice remained in the convection plume
       qrc_diss = max(0., qrco(i,k) - outqc(i,k) * dtime)

       !dp  = 100.*(po_cup(i,k)-po_cup(i,k+1))

       !--- get relative humidity
       frh = 0. !min(qo_cup(i,k)/qeso_cup(i,k),1.)

       !--- estimation of the fractional area
       fractional_area = (xmb(i)/sig(i)) * zuo(i,k) / (rho_hydr(i,k)*vvel2d(i,k))

       !--- source of enviroment moistening/cooling due to the 'remained' cloud dissipation into it.
       outqc_diss = ( qrc_diss * (1.-frh) ) / cloud_lifetime

       if(versionx==1 .or. COUPL_MPHYSICS .eqv. .false.) then

         outt_diss  = -outqc_diss*(xlv/cp) !--- cooling

         !--- source of enviroment moistening/warming due to the 'remained' in-cloud water vapor mixing into it.
         !  qvx   = qco   (i,k)
         !  tempx = tempco(i,k)
            qvx   = qeso_cup(i,k)
            tempx =(heso_cup(i,k)-g*zo(i,k)-xlv*qeso_cup(i,k))/cp

         outq_mix = ( qvx   - qo_cup(i,k) ) / cloud_lifetime

         outt_mix = ( tempx - tn_cup(i,k) ) / cloud_lifetime

         !-- feedback
         del_q = (outqc_diss + outq_mix) * use_cloud_dissipation * fractional_area ! units: kg[water]/kg[air}/sec
         del_t = (outt_diss  + outt_mix) * use_cloud_dissipation * fractional_area ! units: K/sec

         outq (i,k) = outq (i,k) + del_q
         outt (i,k) = outt (i,k) + del_t

        else

         outqc(i,k) = outqc(i,k) + outqc_diss* fractional_area * use_cloud_dissipation

        endif

        !print*,"diss2=",k,real(outqc_diss*86400.*1000),real(sqrt(1.-sig(i)),4),real( fractional_area*100.,4)

        qrco    (i,k) = max(0., qrco(i,k) - outqc_diss * use_cloud_dissipation * fractional_area *dtime)
       !if(qrco (i,k) <0.) print*,"qrc<0",trim(cumulus),qrco(i,k)

     enddo
 enddo

 END SUBROUTINE cloud_dissipation
!------------------------------------------------------------------------------------
 subroutine GF_convpar_init (mynum)
    implicit none
    integer,intent(in) :: mynum
    integer :: nlunit = 4
    character (len = 64) :: fn_nml = 'GF_ConvPar_nml'
    logical :: exists
    namelist /GF_NML/ icumulus_gf, closure_choice,use_scale_dep,dicycle                                  &
                 ,use_tracer_transp, use_tracer_scaven,use_flux_form, use_tracer_evap,downdraft,use_fct  &
                 ,use_rebcb, vert_discr, satur_calc, clev_grid,apply_sub_mp, alp1         &
                 ,sgs_w_timescale, lightning_diag,autoconv, bc_meth,overshoot,use_wetbulb &
                 ,c1,c0_deep, qrc_crit,lambau_deep,lambau_shdn,c0_mid                     &
                 ,cum_max_edt_land  ,cum_max_edt_ocean, cum_hei_down_land                 &
                 ,cum_hei_down_ocean,cum_hei_updf_land, cum_hei_updf_ocean                &
                 ,cum_entr_rate ,tau_deep,tau_mid                                         &
                 ,zero_diff ,use_momentum_transp ,moist_trigger,frac_modis                &
                 ,cum_use_excess,cum_ave_layer,adv_trigger, evap_fix      &
                 ,use_cloud_dissipation,use_smooth_tend,use_gustiness, use_random_num     &
                 ,dcape_threshold,beta_sh,c0_shal,use_linear_subcl_mf

    inquire (file = trim (fn_nml), exist = exists)
    if (.not. exists) then
        write (6, *) 'GF_convpar_nml :: namelist file: ', trim (fn_nml), ' does not exist'
        stop 31415
    else
        open (nlunit,file=fn_nml,status='old',form='formatted')
        read (nlunit,nml=GF_NML)
        close(nlunit)
    endif
    if(mynum==1) then
      !- print the namelist
      print*,"           "
      print*,"------------- GF ConvPar namelist -------------"
      print*, 'icumulus_gf        ' , icumulus_gf
      print*, 'closure_choice     ' , closure_choice
      print*, 'clev_grid          ' , clev_grid
      print*, 'use_rebcb          ' , use_rebcb
      print*, 'vert_discr         ' , vert_discr
      print*, 'satur_calc         ' , satur_calc
      print*, 'autoconv           ' , autoconv
      print*, 'bc_meth            ' , bc_meth
      print*, 'overshoot          ' , real(overshoot          ,4)
      print*, 'use_wetbulb        ' , use_wetbulb
      print*, 'hei_down_land      ' , real(cum_hei_down_land  ,4)
      print*, 'hei_down_ocean     ' , real(cum_hei_down_ocean ,4)
      print*, 'hei_updf_land      ' , real(cum_hei_updf_land  ,4)
      print*, 'hei_updf_ocean     ' , real(cum_hei_updf_ocean ,4)
      print*, 'max_edt_land       ' , real(cum_max_edt_land   ,4)
      print*, 'max_edt_ocean      ' , real(cum_max_edt_ocean  ,4)
      print*, 'cum_use_excess     ' , cum_use_excess
      print*, 'c0_deep            ' , real(c0_deep            ,4)
      print*, 'c0_mid             ' , real(c0_mid             ,4)
      print*, 'c0_shal            ' , real(c0_shal            ,4)
      print*, 'c1                 ' , real(c1                 ,4)
      print*, 'qrc_crit           ' , real(qrc_crit           ,4)
      print*, 'lambau_deep        ' , real(lambau_deep        ,4)
      print*, 'lambau_shdn        ' , real(lambau_shdn        ,4)
      print*, 'use_tracer_transp  ' , use_tracer_transp
      print*, 'use_tracer_scaven  ' , use_tracer_scaven
      print*, 'use_flux_form      ' , use_flux_form
      print*, 'use_fct            ' , use_fct
      print*, 'use_tracer_evap    ' , use_tracer_evap
      print*, 'use_momentum_trans ' , use_momentum_transp
      print*, 'downdraft          ' , downdraft
      print*, 'sgs_w_timescale    ' , sgs_w_timescale
      print*, 'apply_sub_mp       ' , apply_sub_mp
      print*, 'alp1               ' , real(alp1               ,4)
      print*, 'lightning_diag     ' , lightning_diag
      print*, 'use_scale_dep      ' , use_scale_dep
      print*, 'dicycle            ' , dicycle
      print*, 'zero_diff          ' , zero_diff
      print*, 'cum_entr           ' , real(cum_entr_rate      ,4)
      print*, 'frac_modis         ' , frac_modis
      print*, 'moist_trigger      ' , moist_trigger
      print*, 'adv_trigger        ' , adv_trigger
      print*, 'dcape_threshold    ' , real(dcape_threshold    ,4)
      print*, 'tau_deep,tau_mid   ' , real(tau_deep,4),real(tau_mid,4)
      print*, 'evap_fix           ' , evap_fix
      print*, 'use_cloud_dissipat ' , real(use_cloud_dissipation,4)
      print*, 'use_gustiness      ' , use_gustiness
      print*, 'use_random_num     ' , use_random_num
      print*, 'beta_sh            ' , real(beta_sh            , 4)
      print*, 'use_linear_subcl_mf' , use_linear_subcl_mf
      print*,"==============================================="
      call flush(6)
    endif
 end subroutine GF_convpar_init
!------------------------------------------------------------------------------------
 subroutine get_liq_ice_number_conc(itf,ktf,its,ite, kts,kte,ierr,ktop&
                                  ,cnvfrc,srftype,dtime,rho,outqc,tempco,outnliq,outnice)

    implicit none
    integer,   intent (in )  :: itf,ktf,its,ite,kts,kte
    real,      intent (in )  :: dtime

    real,    dimension (its:ite)         ,intent (in )  :: cnvfrc,srftype
    integer, dimension (its:ite)         ,intent (in )  :: ierr,ktop
    real,    dimension (its:ite,kts:kte) ,intent (in )  :: outqc,tempco,rho
    real,    dimension (its:ite,kts:kte) ,intent (out)  :: outnliq,outnice

    integer :: i,k
    real    :: fr,tqliq,tqice,dtinv

    real,   dimension (its:ite,kts:kte) :: nwfa   ! in the future set this as NCPL
    real,   dimension (its:ite,kts:kte) :: nifa   ! in the future set this as NCPI


    nwfa(:,:) =  99.e7  ! in the future set this as NCPL
    nifa(:,:) =  0.     ! in the future set this as NCPI
    dtinv    = 1./dtime
    do i=its,itf
         if(ierr(i) /= 0) cycle

         do k=kts,ktop(i)+1

            fr    = fract_liq_f(tempco(i,k),cnvfrc(i),srftype(i))
            tqliq = dtime * outqc(i,k)* rho(i,k) * fr
            tqice = dtime * outqc(i,k)* rho(i,k) * (1.-fr)

            outnice(i,k) = max(0.0,  make_IceNumber    (tqice, tempco(i,k))/rho(i,k))

            outnliq(i,k) = max(0.0,  make_DropletNumber(tqliq, nwfa  (i,k))/rho(i,k))

         enddo
         !-- convert in tendencies
         outnice = outnice * dtinv ! unit [1/s]
         outnliq = outnliq * dtinv ! unit [1/s]
         !--- for update
         ! nwfa =nwfa + outnliq*dtime
         ! nifa =nifa + outnice*dtime

    enddo

  end subroutine get_liq_ice_number_conc
!DSM {
  pure function intfuncgamma(x, y) result(z)
    real :: z
    real, intent(in) :: x, y

    z = x**(y-1.0) * exp(-x)
  end function intfuncgamma

  function gammaBrams(a) result(g)
    real :: g
    real, intent(in) :: a

    real, parameter :: small = 1.0e-4
    integer, parameter :: points = 100000

    real :: infty, dx, p, sp(2, points), x
    integer :: i
    logical :: correction

    x = a

    correction = .false.
    ! value with x<1 gives \infty, so we use
    ! \Gamma(x+1) = x\Gamma(x)
    ! to avoid the problem
    if ( x < 1.0 ) then
      correction = .true.
      x = x + 1
    end if

    ! find a "reasonable" infinity...
    ! we compute this integral indeed
    ! \int_0^M dt t^{x-1} e^{-t}
    ! where M is such that M^{x-1} e^{-M} ≤ \epsilon
    infty = 1.0e4
    do while ( intfuncgamma(infty, x) > small )
      infty = infty * 10.0
    end do

    ! using simpson
    dx = infty/real(points)
    sp = 0.0
    forall(i=1:points/2-1) sp(1, 2*i) = intfuncgamma(2.0*(i)*dx, x)
    forall(i=1:points/2) sp(2, 2*i - 1) = intfuncgamma((2.0*(i)-1.0)*dx, x)
    g = (intfuncgamma(0.0, x) + 2.0*sum(sp(1,:)) + 4.0*sum(sp(2,:)) + &
    intfuncgamma(infty, x))*dx/3.0

    if ( correction ) g = g/a

  end function gammaBrams
!DSM}
!----------------------------------------------------------------------
  subroutine gen_random(its,ite,use_random_num,random)
   implicit none
   integer, intent(in)  :: its,ite
   real,    intent(in)  :: use_random_num
   real,    intent(out) :: random(its:ite)

   !-local vars
   integer   :: i
   integer(8) :: iran, ranseed = 0

   call system_clock(ranseed)
   ranseed=mod(ranseed,2147483646)+1 !seed between 1 and 2^31-2
   iran = -ranseed

   !-- ran1 produces numbers between [ 0,1]
   !-- random        will be between [-1,1]
   !-- with use_random_num the interval will be [-use_random_num,+use_random_num]
   do i=its,ite
     random(i) = use_random_num * 2.0*(0.5-real(RAN1(IRAN),4))
     !print*,"ran=",i,random(i)
   enddo

   if(maxval(abs(random)) > use_random_num) stop "random > use_random_num"

  end subroutine gen_random
!----------------------------------------------------------------------

  real(8) function ran1(idum)

! This is contributed code standardized by Yong Wang
! Random number generator taken from Press et al.
!
! Returns numbers in the range 0-->1
!
! Their description...
! "Minimal" random number generator of Park and Miller with Bays-Durham
! shuffle and added safeguards. Returns a uniform deviate between 0.0 and 1.0
! (exclusive of the endpoint values). Call with idum a negative integer to
! initialize; thereafter, do not alter idum between successive calls in a
! sequence. RNMX should approximate the largest floating value that is less
! than 1.

   !use shr_kind_mod,         only: r8 => shr_kind_r8, i8 => shr_kind_i8
   implicit none
   integer(8), parameter:: ntab = 32,iq = 127773,ia = 16807,ir = 2836, &
       im = 2147483647,ndiv = 1+(im-1)/ntab

   real(8), parameter:: am = 1.0/im,eps = 1.2e-7,rnmx = 1.0-eps

   integer(8), intent(inout):: idum

   integer(8):: iy
   integer(8), dimension(ntab):: iv
   !save iv,iy
   data iv /ntab*0/, iy /0/
   integer(8):: j,k

   !
   if (idum.le.0.or.iy.eq.0) then
   ! initalize
       idum = max(-idum,1)
       do j = ntab+8,1,-1
          k = idum/iq
          idum = ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum = idum+im
          if (j.le.ntab) iv(j) = idum
       end do
       iy = iv(1)
   end if
   !
   k = idum/iq
   ! compute idum = mod(ia*idum,im) without overflows by schrage's method
   idum = ia*(idum-k*iq)-ir*k
   if (idum.lt.0) idum = idum+im
   ! j will be in the range 1-->ntab
   j = 1+iy/ndiv
   ! output previously stored value and refill the shuffle table
   iy = iv(j)
   iv(j) = idum
   ran1 = min(am*iy,rnmx)

 end function ran1
!----------------------------------------------------------------------

 SUBROUTINE get_delmix(cumulus,kts,kte,ktf,xland,subcl_level,po,ain,aout)
    implicit none
    character *(*)   ,intent (in) :: cumulus
    integer,intent(in)            :: kts,kte,ktf,subcl_level
    real   ,intent(in)            :: ain(kts:kte),po(kts:kte),xland
    real   ,intent(inout)         :: aout(kts:kte)

    !-- local var
    real :: x1,x2,dp,del,qc
    integer :: k

    !-
    qc = aout(kts)

    x2=0. ; x1=0.
    do k = kts,subcl_level
      dp = po(k+1)-po(k)
      x2 = x2 + dp
      x1 = x1 + dp*ain(k)
    enddo
    del = abs(qc-x1/(x2+1.e-12))
    aout(kts:subcl_level) =  ain(kts:subcl_level) + del
!----------------------------------------------------------------------

 end SUBROUTINE get_delmix


END  MODULE ConvPar_GF2020
