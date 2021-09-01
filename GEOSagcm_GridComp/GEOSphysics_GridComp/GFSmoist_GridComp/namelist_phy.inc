!--- BEGIN NAMELIST VARIABLES
    real(kind=kind_phys) :: fhzero         = 0.0             !< hours between clearing of diagnostic buckets
    logical              :: ldiag3d        = .true.          !< flag for 3d diagnostic fields
    logical              :: qdiag3d        = .true.          !< flag for 3d tracer diagnostic fields
    logical              :: lssav          = .true.         !< logical flag for storing diagnostics
    integer              :: naux2d         = 0               !< number of auxiliary 2d arrays to output (for debugging)
    integer              :: naux3d         = 0               !< number of auxiliary 3d arrays to output (for debugging)
    logical              :: aux2d_time_avg(1:20) = .false. !< flags for time averaging of auxiliary 2d arrays
    logical              :: aux3d_time_avg(1:20) = .false. !< flags for time averaging of auxiliary 3d arrays

    real(kind=kind_phys) :: fhcyc          = 0.              !< frequency for surface data cycling (hours)
    integer              :: thermodyn_id   =  1              !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id    =  1              !< valid for GFS only for get_prs/phi

    !--- coupling parameters
    logical              :: cplflx         = .false.         !< default no cplflx collection
    logical              :: cplwav         = .false.         !< default no cplwav collection
    logical              :: cplwav2atm     = .false.         !< default no cplwav2atm coupling
    logical              :: cplchm         = .false.         !< default no cplchm collection

!--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         = .false.

!--- radiation parameters
    real(kind=kind_phys) :: fhswr          = 3600.           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr          = 3600.           !< frequency for longwave radiation (secs)
    integer              :: nhfrad         = 0               !< number of timesteps for which to call radiation on physics timestep (coldstarts)
    integer              :: levr           = -99             !< number of vertical levels for radiation calculations
    integer              :: nfxr           = 39+6            !< second dimension of input/output array fluxr   
    logical              :: iaerclm        = .false.         !< flag for initializing aero data 
    integer              :: iccn           =  0              !< logical to use IN CCN forcing for MG2/3
    integer              :: iflip          =  1              !< iflip - is not the same as flipv
    integer              :: isol           =  0              !< use prescribed solar constant
    integer              :: ico2           =  0              !< prescribed global mean value (old opernl)
    integer              :: ialb           =  0              !< use climatology alb, based on sfc type
                                                             !< 1 => use modis based alb
    integer              :: iems           =  0              !< use fixed value of 1.0
    integer              :: iaer           =  1              !< default aerosol effect in sw only
    integer              :: icliq_sw       =  1              !< sw optical property for liquid clouds
    integer              :: iovr_sw        =  1              !< sw: max-random overlap clouds
    integer              :: iovr_lw        =  1              !< lw: max-random overlap clouds
    integer              :: ictm           =  1              !< ictm=0 => use data at initial cond time, if not
                                                             !<           available; use latest; no extrapolation.
                                                             !< ictm=1 => use data at the forecast time, if not
                                                             !<           available; use latest; do extrapolation.
                                                             !< ictm=yyyy0 => use yyyy data for the forecast time;
                                                             !<           no extrapolation.
                                                             !< ictm=yyyy1 = > use yyyy data for the fcst. If needed, 
                                                             !<           do extrapolation to match the fcst time.
                                                             !< ictm=-1 => use user provided external data for
                                                             !<           the fcst time; no extrapolation.
                                                             !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                                             !<           from climatology; no extrapolation.
    integer              :: isubc_sw          =  0           !< sw clouds without sub-grid approximation
    integer              :: isubc_lw          =  0           !< lw clouds without sub-grid approximation
                                                             !< =1 => sub-grid cloud with prescribed seeds
                                                             !< =2 => sub-grid cloud with randomly generated
                                                             !< seeds
    logical              :: crick_proof       = .false.      !< CRICK-Proof cloud water
    logical              :: ccnorm            = .false.      !< Cloud condensate normalized by cloud cover 
    logical              :: norad_precip      = .false.      !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr             = .true.       !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr             = .true.       !< flag to output sw heating rate (Radtend%swhc)
    ! RRTMGP
    logical              :: do_RRTMGP        = .false.       !< Use RRTMGP?
    character(len=128)   :: active_gases    = ''             !< Character list of active gases used in RRTMGP
    integer              :: nGases          = 0              !< Number of active gases
    character(len=128)   :: rrtmgp_root     = ''             !< Directory of rte+rrtmgp source code
    character(len=128)   :: lw_file_gas     = ''             !< RRTMGP K-distribution file, coefficients to compute optics for gaseous atmosphere
    character(len=128)   :: lw_file_clouds  = ''             !< RRTMGP file containing coefficients used to compute clouds optical properties
    integer              :: rrtmgp_nBandsLW = 16             !< Number of RRTMGP LW bands.
    integer              :: rrtmgp_nGptsLW  = 256            !< Number of RRTMGP LW spectral points.
    character(len=128)   :: sw_file_gas     = ''             !< RRTMGP K-distribution file, coefficients to compute optics for gaseous atmosphere
    character(len=128)   :: sw_file_clouds  = ''             !< RRTMGP file containing coefficients used to compute clouds optical properties
    integer              :: rrtmgp_nBandsSW = 14             !< Number of RRTMGP SW bands.
    integer              :: rrtmgp_nGptsSW  = 224            !< Number of RRTMGP SW spectral points.
    logical              :: doG_cldoptics       = .false.    !< Use legacy RRTMG cloud-optics?                                             
    logical              :: doGP_cldoptics_PADE = .false.    !< Use RRTMGP cloud-optics: PADE approximation?
    logical              :: doGP_cldoptics_LUT  = .false.    !< Use RRTMGP cloud-optics: LUTs?     
    integer              :: rrtmgp_nrghice = 0               !< Number of ice-roughness categories
    integer              :: rrtmgp_nGauss_ang=1              !< Number of angles used in Gaussian quadrature
    logical              :: do_GPsw_Glw    = .false.     
    logical              :: use_LW_jacobian = .false.        !< Use Jacobian of LW to update LW radiation tendencies. 

!--- Z-C microphysical parameters
    integer              :: ncld              =  1                 !< choice of cloud scheme
    integer              :: imp_physics       =  99                !< choice of cloud scheme
    real(kind=kind_phys) :: psautco(2)        = (/6.0d-4,3.0d-4/)  !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)        = (/1.0d-4,1.0d-4/)  !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco             = 2.0d-5             !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)         = (/1.0d-5,1.0d-5/)  !< [in] water and ice minimum threshold for Zhao
!---Max hourly
    real(kind=kind_phys) :: avg_max_length = 3600.                 !< reset value in seconds for max hourly
!--- Ferrier-Aligo microphysical parameters
    real(kind=kind_phys) :: rhgrd             = 0.98               !< fer_hires microphysics only; for 3-km domain
    logical              :: spec_adv          = .true.             !< Individual cloud species advected
    integer              :: icloud            = 0                  !< cloud effect to the optical depth in radiation; this also controls the cloud fraction options
                                                                   !<  3: with cloud effect from FA, and use cloud fraction option 3, based on Sundqvist et al. (1989)
!--- M-G microphysical parameters
    integer              :: fprcp             =  0                 !< no prognostic rain and snow (MG)
    integer              :: pdfflag           =  4                 !< pdf flag for MG macro physics
    real(kind=kind_phys) :: mg_dcs            = 200.0              !< Morrison-Gettelman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar          = 1.0
    real(kind=kind_phys) :: mg_ts_auto_ice(2) = (/180.0,180.0/)    !< ice auto conversion time scale
    real(kind=kind_phys) :: mg_rhmini         = 1.01               !< relative humidity threshold parameter for nucleating ice
    real(kind=kind_phys) :: mg_ncnst          = 100.e6             !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst          = 0.15e6             !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst          = 0.10e6             !< constant graupel/hail num concentration (m-3) = 0.1e6_r8
    real(kind=kind_phys) :: mg_alf            = 1.0                !< tuning factor for alphs in MG macrophysics
    real(kind=kind_phys) :: mg_qcmin(2)       = (/1.0d-9,1.0d-9/)  !< min liquid and ice mixing ratio in Mg macro clouds
    real(kind=kind_phys) :: mg_berg_eff_factor = 2.0               !< berg efficiency factor
    character(len=16)    :: mg_precip_frac_method = 'max_overlap'  !< type of precipitation fraction method
    real(kind=kind_phys) :: tf              = 258.16d0
    real(kind=kind_phys) :: tcr             = 273.16d0
!
    logical              :: effr_in         = .true.              !< flag to use effective radii of cloud species in radiation
    logical              :: microp_uniform  = .true.
    logical              :: do_cldliq       = .true.
    logical              :: do_cldice       = .true.
    logical              :: hetfrz_classnuc = .false.
    logical              :: mg_nccons       = .false.           !< set .true. to specify constant cloud droplet number
    logical              :: mg_nicons       = .false.           !< set .true. to specify constant cloud ice number
    logical              :: mg_ngcons       = .false.           !< set .true. to specify constant graupel/hail number
    logical              :: sed_supersat    = .true.
    logical              :: do_sb_physics   = .true.
    logical              :: mg_do_graupel   = .true.            !< set .true. to turn on prognostic grapuel (with fprcp=2)
    logical              :: mg_do_hail      = .false.           !< set .true. to turn on prognostic hail (with fprcp=2)
    logical              :: mg_do_ice_gmao  = .false.           !< set .true. to turn on gmao ice formulation
    logical              :: mg_do_liq_liu   = .true.            !< set .true. to turn on liu liquid treatment


    !--- Thompson microphysical parameters
    logical              :: ltaerosol      = .false.            !< flag for aerosol version
    logical              :: lradar         = .false.            !< flag for radar reflectivity 
    real(kind=kind_phys) :: nsradar_reset  = -999.0             !< seconds between resetting radar reflectivity calculation, set to <0 for every time step
    real(kind=kind_phys) :: ttendlim       = -999.0             !< temperature tendency limiter, set to <0 to deactivate

    !--- GFDL microphysical parameters
    logical              :: lgfdlmprad     = .false.            !< flag for GFDLMP radiation interaction 

    !--- Thompson,GFDL microphysical parameter
    logical              :: lrefres        = .false.            !< flag for radar reflectivity in restart file

    !--- land/surface model parameters
    integer              :: lsm            =  1              !< flag for land surface model to use =0  for osu lsm; =1  for noah lsm; =2  for noah mp lsm; =3  for RUC lsm
    integer              :: lsoil          =  4              !< number of soil layers
    integer              :: lsoil_lsm      =  -1             !< number of soil layers internal to land surface model; -1 use lsoil
    integer              :: lsnow_lsm      =  3              !< maximum number of snow layers internal to land surface model
    logical              :: rdlai          = .false.         !< read LAI from input file (for RUC LSM or NOAH LSM WRFv4)
    logical              :: ua_phys        = .false.         !< flag for using University of Arizona? extension to NOAH LSM WRFv4
    logical              :: usemonalb      = .true.          !< flag to read surface diffused shortwave albedo from input file for NOAH LSM WRFv4
    real(kind=kind_phys) :: aoasis         = 1.0             !< potential evaporation multiplication factor for NOAH LSM WRFv4
    integer              :: fasdas         = 0               !< flag to use "flux-adjusting surface data assimilation system"; 0 = OFF, 1 = ON
    integer              :: ivegsrc        =  2              !< ivegsrc = 0   => USGS,
                                                             !< ivegsrc = 1   => IGBP (20 category)
                                                             !< ivegsrc = 2   => UMD  (13 category)
    integer              :: isot           =  0              !< isot = 0   => Zobler soil type  ( 9 category)
                                                             !< isot = 1   => STATSGO soil type (19 category)
    integer              :: iopt_thcnd     =  1              !< option to treat thermal conductivity in Noah LSM (new in 3.8)
                                                             !< = 1, original (default)
                                                             !< = 2, McCumber and Pielke for silt loam and sandy loam
    ! -- to use Noah MP, lsm needs to be set to 2 and both ivegsrc and isot are set
    ! to 1 - MODIS IGBP and STATSGO - the defaults are the same as in the
    ! scripts;change from namelist

    integer              :: iopt_dveg      =  4  ! 4 -> off (use table lai; use maximum vegetation fraction)
    integer              :: iopt_crs       =  1  !canopy stomatal resistance (1-> ball-berry; 2->jarvis)
    integer              :: iopt_btr       =  1  !soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
    integer              :: iopt_run       =  3  !runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
    integer              :: iopt_sfc       =  1  !surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
    integer              :: iopt_frz       =  1  !supercooled liquid water (1-> ny06; 2->koren99)
    integer              :: iopt_inf       =  1  !frozen soil permeability (1-> ny06; 2->koren99)
    integer              :: iopt_rad       =  3  !radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
    integer              :: iopt_alb       =  2  !snow surface albedo (1->bats; 2->class)
    integer              :: iopt_snf       =  1  !rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
    integer              :: iopt_tbot      =  2  !lower boundary of soil temperature (1->zero-flux; 2->noah)
    integer              :: iopt_stc       =  1  !snow/soil temperature time scheme (only layer 1)

    logical              :: use_ufo        = .false.         !< flag for gcycle surface option

!--- flake model parameters
    integer              :: lkm            =  0  !< flag for flake model

    logical              :: lcurr_sf       = .false.                  !< flag for taking ocean currents into account in GFDL surface layer
    logical              :: pert_cd        = .false.                  !< flag for perturbing the surface drag coefficient for momentum in surface layer scheme
    integer              :: ntsflg         = 0                        !< flag for updating skin temperature in the GFDL surface layer scheme
    real(kind=kind_phys) :: sfenth         = 0.0                      !< enthalpy flux factor 0 zot via charnock ..>0 zot enhanced>15m/s
    
!--- tuning parameters for physical parameterizations
    logical              :: ras            = .false.                  !< flag for ras convection scheme
    logical              :: flipv          = .true.                   !< flag for vertical direction flip (ras)
                                                                      !< .true. implies surface at k=1
    logical              :: trans_trac     = .false.                  !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin      = .false.                  !< flag for diff monin schemes
    logical              :: cnvgwd         = .false.                  !< flag for conv gravity wave drag
    integer              :: gwd_opt        =  1                       !< flag for configuring gwd scheme
                                                                      !< gwd_opt = 3 : GSDdrag suite
                                                                      !< gwd_opt = 33: GSDdrag suite with extra output
!--- vay-2018
    logical              :: ldiag_ugwp     = .false.                  !< flag for UGWP diag fields
    logical              :: do_ugwp        = .false.                  !< flag do UGWP+RF
    logical              :: do_tofd        = .false.                  !< flag do Turb oro Form Drag

    logical              :: do_gwd         = .false.                  !< flag for running gravity wave drag
    logical              :: do_cnvgwd      = .false.                  !< flag for running conv gravity wave drag

    logical              :: mstrat         = .false.                  !< flag for moorthi approach for stratus
    logical              :: moist_adj      = .false.                  !< flag for moist convective adjustment
    logical              :: cscnv          = .false.                  !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre        = .false.                  !< flag controls precip type algorithm
    logical              :: do_aw          = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_awdd        = .false.                  !< AW scale-aware option in cs convection
    logical              :: flx_form       = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_shoc        = .false.                  !< flag for SHOC
    logical              :: shocaftcnv     = .false.                  !< flag for SHOC
    logical              :: shoc_cld       = .false.                  !< flag for SHOC in grrad
    logical              :: oz_phys        = .true.                   !< flag for old (2006) ozone physics
    logical              :: oz_phys_2015   = .false.                  !< flag for new (2015) ozone physics
    logical              :: h2o_phys       = .false.                  !< flag for stratosphere h2o
    logical              :: pdfcld         = .false.                  !< flag for pdfcld
    logical              :: shcnvcw        = .false.                  !< flag for shallow convective cloud
    logical              :: redrag         = .false.                  !< flag for reduced drag coeff. over sea
    logical              :: hybedmf        = .false.                  !< flag for hybrid edmf pbl scheme
    logical              :: satmedmf       = .false.                  !< flag for scale-aware TKE-based moist edmf
                                                                      !< vertical turbulent mixing scheme
    logical              :: shinhong       = .false.                  !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu         = .false.                  !< flag for YSU vertical turbulent mixing scheme
    logical              :: acm            = .false.                  !< flag for ACM vertical turbulent mixing scheme
    logical              :: dspheat        = .false.                  !< flag for tke dissipative heating
    logical              :: hurr_pbl       = .false.                  !< flag for hurricane-specific options in PBL scheme
    logical              :: lheatstrg      = .false.                  !< flag for canopy heat storage parameterization
    logical              :: cnvcld         = .false.
    logical              :: random_clds    = .false.                  !< flag controls whether clouds are random
    logical              :: shal_cnv       = .false.                  !< flag for calling shallow convection
    integer              :: imfshalcnv     =  1                       !< flag for mass-flux shallow convection scheme
                                                                      !<     1: July 2010 version of mass-flux shallow conv scheme
                                                                      !<         current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
                                                                      !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                                                      !<    -1: no shallow convection used
    integer              :: imfdeepcnv     =  1                       !< flag for mass-flux deep convection scheme
                                                                      !<     1: July 2010 version of SAS conv scheme
                                                                      !<           current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
    integer              :: isatmedmf      =  0                       !< flag for scale-aware TKE-based moist edmf scheme
                                                                      !<     0: initial version of satmedmf (Nov. 2018)
                                                                      !<     1: updated version of satmedmf (as of May 2019)
    logical              :: do_deep        = .true.                   !< whether to do deep convection
    logical              :: hwrf_samfdeep     = .false.               !< flag for HWRF SAMF deepcnv scheme 
    logical              :: hwrf_samfshal     = .false.               !< flag for HWRF SAMF shalcnv scheme 
    logical              :: do_mynnedmf       = .false.               !< flag for MYNN-EDMF
    logical              :: do_mynnsfclay     = .false.               !< flag for MYNN Surface Layer Scheme
    integer              :: grav_settling     = 0
    integer              :: bl_mynn_tkebudget = 0
    logical              :: bl_mynn_tkeadvect = .false.
    integer              :: bl_mynn_cloudpdf  = 2
    integer              :: bl_mynn_mixlength = 2
    integer              :: bl_mynn_edmf      = 0
    integer              :: bl_mynn_edmf_mom  = 1
    integer              :: bl_mynn_edmf_tke  = 0
    integer              :: bl_mynn_edmf_part = 0
    integer              :: bl_mynn_cloudmix  = 1
    integer              :: bl_mynn_mixqt     = 0
    integer              :: bl_mynn_output    = 0
    integer              :: icloud_bl         = 1
    real(kind=kind_phys) :: var_ric           = 1.0
    real(kind=kind_phys) :: coef_ric_l        = 0.16
    real(kind=kind_phys) :: coef_ric_s        = 0.25
    logical              :: do_myjsfc         = .false.               !< flag for MYJ surface layer scheme
    logical              :: do_myjpbl         = .false.               !< flag for MYJ PBL scheme

    integer              :: nmtvr          = 14                       !< number of topographic variables such as variance etc
                                                                      !< used in the GWD parameterization
    integer              :: jcap           =  1              !< number of spectral wave trancation used only by sascnv shalcnv
!   real(kind=kind_phys) :: cs_parm(10) = (/5.0,2.5,1.0e3,3.0e3,20.0,-999.,-999.,0.,0.,0./)
    real(kind=kind_phys) :: cs_parm(10) = (/8.0,4.0,1.0e3,3.5e3,20.0,1.0,-999.,1.,0.6,0./)
    real(kind=kind_phys) :: flgmin(2)      = (/0.180,0.220/)          !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)        = (/0.5d0,0.05d0/)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)        = (/1.0d0,1.0d0/)          !< multiplication factor for critical cloud
                                                                      !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(4)     = (/2.0d0,0.25d0,1.0d0,1.0d0/)   !< multiplication factors for cdmb, gwd, and NS gwd, tke based enhancement
    real(kind=kind_phys) :: sup            = 1.0                      !< supersaturation in pdf cloud (IMP_physics=98) when t is very low
                                                                      !< or ice super saturation in SHOC (when do_shoc=.true.)
    real(kind=kind_phys) :: ctei_rm(2)     = (/10.0d0,10.0d0/)        !< critical cloud top entrainment instability criteria 
                                                                      !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)       = (/0.90d0,0.90d0,0.90d0/) !< critical relative humidity at the surface
                                                                      !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)        = (/0.0d0,0.0d0/)          !< factor for cloud condensate detrainment 
                                                                      !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)     = (/1.0d-3,1.0d-3/)        !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)     = (/2.0d-3,2.0d-3/)        !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)     = (/1.0d-5,1.0d-5/)        !< [in] water and ice minimum threshold for ras
    integer              :: nrcmax         = 32                       !< number of random numbers used in RAS
    
    real(kind=kind_phys) :: rbcr           = 0.25                     !< Critical Richardson Number in PBL scheme
    real(kind=kind_phys) :: shoc_parm(5)   = (/7000.0,1.0,4.2857143,0.7,-999.0/)  !< some tunable parameters for shoc

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0        = 0.0d0           !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts         = 0.0d0           !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep      = 0.1             !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep       = 0.002           !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep        = 0.002           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evfact_deep    = 0.3             !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep   = 0.3             !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land 

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal      = 0.3             !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal       = 0.002           !< conversion parameter of detrainment from liquid water into convetive precipitaiton
    real(kind=kind_phys) :: c1_shal        = 5.e-4           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land 

!--- near surface sea temperature model
    logical              :: nst_anl        = .false.         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           = 0
    integer              :: nstf_name(5)   = (/0,0,1,0,5/)   !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                                             !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                                             !<                2 = NSSTM on and coupled
                                                             !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                                             !< nstf_name(3) : 1 = NSSTM analysis on, 0 = NSSTM analysis off
                                                             !< nstf_name(4) : zsea1 in mm
                                                             !< nstf_name(5) : zsea2 in mm
!--- fractional grid
    logical              :: frac_grid       = .false.         !< flag for fractional grid
    logical              :: ignore_lake     = .true.          !< flag for ignoring lakes
    real(kind=kind_phys) :: min_lakeice     = 0.15d0          !< minimum lake ice value
    real(kind=kind_phys) :: min_seaice      = 1.0d-11         !< minimum sea  ice value
    real(kind=kind_phys) :: min_lake_height = 250.0           !< minimum lake height value
    real(kind=kind_phys) :: rho_h2o         = 1000.d0 !rhowater        !< fresh water density

!--- surface layer z0 scheme
    integer              :: sfc_z0_type    = 0               !< surface roughness options over ocean
                                                             !< 0=no change
                                                             !< 6=areodynamical roughness over water with input 10-m wind
                                                             !< 7=slightly decrease Cd for higher wind speed compare to 6
                                                             !< negative when cplwav2atm=.true. - i.e. two way wave coupling

!--- vertical diffusion
    real(kind=kind_phys) :: xkzm_m         = 1.0d0           !< [in] bkgd_vdif_m  background vertical diffusion for momentum  
    real(kind=kind_phys) :: xkzm_h         = 1.0d0           !< [in] bkgd_vdif_h  background vertical diffusion for heat q  
    real(kind=kind_phys) :: xkzm_s         = 1.0d0           !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion  
    real(kind=kind_phys) :: xkzminv        = 0.3             !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac     = 1.0             !< turbulence diffusion coefficient factor
    real(kind=kind_phys) :: dspfac         = 1.0             !< tke dissipative heating factor
    real(kind=kind_phys) :: bl_upfr        = 0.13            !< updraft fraction in boundary layer mass flux scheme
    real(kind=kind_phys) :: bl_dnfr        = 0.1             !< downdraft fraction in boundary layer mass flux scheme

!--- parameters for canopy heat storage (CHS) parameterization
    real(kind=kind_phys) :: z0fac          = 0.3
    real(kind=kind_phys) :: e0fac          = 0.5

 
!---Cellular automaton options
    integer              :: nca            = 1
    integer              :: ncells         = 5
    integer              :: nlives         = 30
    integer              :: nca_g          = 1
    integer              :: ncells_g       = 1
    integer              :: nlives_g       = 100
    real(kind=kind_phys) :: nfracseed      = 0.5
    integer              :: nseed          = 100000
    integer              :: nseed_g        = 100
    integer              :: iseed_ca       = 0
    integer              :: nspinup        = 1
    logical              :: do_ca          = .false.
    logical              :: ca_sgs         = .false.
    logical              :: ca_global      = .false.
    logical              :: ca_smooth      = .false.
    real(kind=kind_phys) :: nthresh        = 0.0
    real                 :: ca_amplitude   = 500.
    integer              :: nsmooth        = 100 
    logical              :: ca_closure     = .false.
    logical              :: ca_entr        = .false.
    logical              :: ca_trigger     = .false.

!--- IAU options
    real(kind=kind_phys)  :: iau_delthrs      = 0           !< iau time interval (to scale increments)
    character(len=240)    :: iau_inc_files(7) = ''          !< list of increment files
    real(kind=kind_phys)  :: iaufhrs(7)       = -1          !< forecast hours associated with increment files
    logical  :: iau_filter_increments         = .false.     !< filter IAU increments
    logical  :: iau_drymassfixer              = .false.     !< IAU dry mass fixer

!--- debug flag
    logical              :: debug          = .false.
    logical              :: pre_rad        = .false.         !< flag for testing purpose

!  max and min lon and lat for critical relative humidity
    integer :: max_lon=5000, max_lat=2000, min_lon=192, min_lat=94
    real(kind=kind_phys) :: rhcmax = 0.9999999               !< max critical rel. hum.

!--- stochastic physics control parameters
    logical :: do_sppt      = .false.
    logical :: use_zmtnblck = .false.
    logical :: do_shum      = .false.
    logical :: do_skeb      = .false.
    integer :: skeb_npass = 11
    logical :: do_sfcperts = .false.   ! mg, sfc-perts
    integer :: nsfcpert    =  6        ! mg, sfc-perts
    real(kind=kind_phys) :: pertz0   = -999.
    real(kind=kind_phys) :: pertzt   = -999.
    real(kind=kind_phys) :: pertshc  = -999.
    real(kind=kind_phys) :: pertlai  = -999.
    real(kind=kind_phys) :: pertalb  = -999.
    real(kind=kind_phys) :: pertvegf = -999.

!--- aerosol scavenging factors
    character(len=20) :: fscav_aero(20) = 'default'

!--- END NAMELIST VARIABLES

    NAMELIST /gfs_physics_nml/                                                              &
                          !--- general parameters
                               fhzero, ldiag3d, qdiag3d, lssav, naux2d, naux3d,             &
                               aux2d_time_avg, aux3d_time_avg, fhcyc,                       &
                               thermodyn_id, sfcpress_id,                                   &
                          !--- coupling parameters
                               cplflx, cplwav, cplwav2atm, cplchm, lsidea,                  &
                          !--- radiation parameters
                               fhswr, fhlwr, levr, nfxr, iaerclm, iflip, isol, ico2, ialb,  &
                               isot, iems, iaer, icliq_sw, iovr_sw, iovr_lw, ictm, isubc_sw,&
                               isubc_lw, crick_proof, ccnorm, lwhtr, swhtr,                 &
                               nhfrad,                                                      &
                          ! --- RRTMGP
                               do_RRTMGP, active_gases, nGases, rrtmgp_root,                &
                               lw_file_gas, lw_file_clouds, rrtmgp_nBandsLW, rrtmgp_nGptsLW,&
                               sw_file_gas, sw_file_clouds, rrtmgp_nBandsSW, rrtmgp_nGptsSW,&
                               doG_cldoptics, doGP_cldoptics_PADE, doGP_cldoptics_LUT,      &
                               rrtmgp_nrghice, rrtmgp_nGauss_ang, do_GPsw_Glw,              &
                               use_LW_jacobian,                                             &
                          ! IN CCN forcing
                               iccn,                                                        &
                          !--- microphysical parameterizations
                               ncld, imp_physics, psautco, prautco, evpco, wminco,          &
                               fprcp, pdfflag, mg_dcs, mg_qcvar, mg_ts_auto_ice, mg_rhmini, &
                               effr_in, tf, tcr,                                            &
                               microp_uniform, do_cldice, hetfrz_classnuc,                  &
                               mg_do_graupel, mg_do_hail, mg_nccons, mg_nicons, mg_ngcons,  &
                               mg_ncnst, mg_ninst, mg_ngnst, sed_supersat, do_sb_physics,   &
                               mg_alf,   mg_qcmin, mg_do_ice_gmao, mg_do_liq_liu,           &
                               ltaerosol, lradar, nsradar_reset, lrefres, ttendlim,         &
                               lgfdlmprad,                                                  &
                          !--- max hourly
                               avg_max_length,                                              &
                          !--- land/surface model control
                               lsm, lsoil, lsoil_lsm, lsnow_lsm, rdlai,                     &
                               nmtvr, ivegsrc, use_ufo, iopt_thcnd, ua_phys, usemonalb,     &
                               aoasis, fasdas,                                              &
                          !    Noah MP options
                               iopt_dveg,iopt_crs,iopt_btr,iopt_run,iopt_sfc, iopt_frz,     &
                               iopt_inf, iopt_rad,iopt_alb,iopt_snf,iopt_tbot,iopt_stc,     &
                          !--- lake model control
                               lkm,                                                         &
                          !    GFDL surface layer options
                               lcurr_sf, pert_cd, ntsflg, sfenth,                           &
                          !--- physical parameterizations
                               ras, trans_trac, old_monin, cnvgwd, mstrat, moist_adj,       &
                               cscnv, cal_pre, do_aw, do_shoc, shocaftcnv, shoc_cld,        &
                               oz_phys, oz_phys_2015,                                       &
                               do_mynnedmf, do_mynnsfclay,                                  &
                               bl_mynn_cloudpdf, bl_mynn_edmf, bl_mynn_edmf_mom,            &
                               bl_mynn_edmf_tke, bl_mynn_edmf_part, bl_mynn_cloudmix,       &
                               bl_mynn_mixqt, bl_mynn_output, icloud_bl, bl_mynn_tkeadvect, &
                               gwd_opt,                                                     &
                               var_ric, coef_ric_l, coef_ric_s,                             &
                               do_myjsfc, do_myjpbl,                                        &
                               hwrf_samfdeep, hwrf_samfshal,                                &
                               h2o_phys, pdfcld, shcnvcw, redrag, hybedmf, satmedmf,        &
                               shinhong, do_ysu, acm, dspheat, hurr_pbl, lheatstrg, cnvcld, &
                               random_clds, shal_cnv, imfshalcnv, imfdeepcnv, isatmedmf,    &
                               do_deep, jcap,                                               &
                               cs_parm, flgmin, cgwf, ccwf, cdmbgwd, sup, ctei_rm, crtrh,   &
                               dlqf, rbcr, shoc_parm, psauras, prauras, wminras,            &
                               do_sppt, do_shum, do_skeb, do_sfcperts,                      &
                          !--- Rayleigh friction
                               prslrd0, ral_ts,  ldiag_ugwp, do_ugwp, do_tofd,              &
                          ! --- Ferrier-Aligo
                               spec_adv, rhgrd, icloud,                                     &
                          !--- mass flux deep convection
                               clam_deep, c0s_deep, c1_deep, betal_deep,                    &
                               betas_deep, evfact_deep, evfactl_deep, pgcon_deep,           &
                               asolfac_deep,                                                &
                          !--- mass flux shallow convection
                               clam_shal, c0s_shal, c1_shal, pgcon_shal, asolfac_shal,      &
                          !--- near surface sea temperature model
                               nst_anl, lsea, nstf_name,                                    &
                               frac_grid, min_lakeice, min_seaice, min_lake_height,         &
                               ignore_lake,                                                 &
                          !--- surface layer
                               sfc_z0_type,                                                 &
                          !    vertical diffusion
                               xkzm_m, xkzm_h, xkzm_s, xkzminv, moninq_fac, dspfac,         &
                               bl_upfr, bl_dnfr,                                            &
                          !--- canopy heat storage parameterization
                               z0fac, e0fac,                                                &
                          !--- cellular automata
                               nca, ncells, nlives, nca_g, ncells_g, nlives_g, nfracseed,   &
                               nseed, nseed_g, nthresh, do_ca,                              &
                               ca_sgs, ca_global,iseed_ca,ca_smooth,                        &
                               nspinup,ca_amplitude,nsmooth,ca_closure,ca_entr,ca_trigger,  &
                          !--- IAU
                               iau_delthrs,iaufhrs,iau_inc_files,iau_filter_increments,     &
                               iau_drymassfixer,                                            &
                          !--- debug options
                               debug, pre_rad,                                              &
                          !--- parameter range for critical relative humidity
                               max_lon, max_lat, min_lon, min_lat, rhcmax,                  &
!                               phys_version,                                                &
                          !--- aerosol scavenging factors ('name:value' string array)
                               fscav_aero
