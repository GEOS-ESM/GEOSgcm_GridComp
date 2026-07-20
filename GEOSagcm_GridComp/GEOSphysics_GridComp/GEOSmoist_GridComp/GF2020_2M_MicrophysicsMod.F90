!===============================================================================
! GF2020_2M_MicrophysicsMod.F90
!
! Diagnostic two-moment updraft microphysics for GF2020 convective plumes.
!
! - Follows the standard GF updraft moisture transport.
! - Activates droplets from a local cloud-base aerosol reservoir.
! - Transports liquid, ice, dust, soot, and sea-salt number as # kg-1.
! - Diagnoses immersion freezing, homogeneous freezing, and deposition ice growth.
! - Uses MG liquid/ice autoconversion for convective precipitation removal.
! - Outputs pice plus liquid/ice mass and number updraft profiles; parent applies the detrainment operator.
!===============================================================================


module GF2020_2M_MicrophysicsMod


  use module_gate
  use ConvPar_GF_SharedParams
  use GEOSmoist_Process_Library, only: AerPropsNew, ERFAPP, ICE_AUTO_TSC_CNV, DCS_CNV, &
       DEBUG_GF2M, ACC_ENH_CNV, ACC_ENH_ICE, FDROPDUST, FDROPSOOT, GF2M_USE_CORRECTOR, &
       FHETSOOT, FHETDUST, AUT_SCALE_CNV, GF2M_W_OPTION, BKG_INP_SC_CNV


  use micro_mg_utils, only: r8, MGHydrometeorProps, size_dist_param_liq, &
       size_dist_param_ice, liu2006_liq_autoconversion, ice_autoconversion, &
       ice_deposition_sublimation,                                      &
       accrete_cloud_water_rain, accrete_cloud_water_snow,              &
       accrete_rain_snow, accrete_cloud_ice_snow,                       &
       secondary_ice_production, self_collection_rain, snow_self_aggregation
  
  
  !use ConvPar_GF_SharedParams, only: cp, g, xlv, xlf
                                   
  
  implicit none
  private

  public :: cup_up_moisture_2M

  !---------------------------------------------------------------------------
  ! Diagnostic updraft-velocity option for full 2M microphysics.
  !   1 = use the internal 2M buoyancy/loading velocity estimate.
  !   2 = use the parent legacy cup_up_vvel profile passed into this routine.
  !   3 = use max(internal, parent legacy) for a conservative sensitivity test.
  ! The parent computes/passes the legacy profile only when needed.
  !---------------------------------------------------------------------------
  !integer :: GF2M_W_OPTION = 1

  ! Full 2M precipitation coupling mode.
  ! .true.  : predictor + fixed-precipitation corrector; repeats the upward plume pass.
  ! .false. : single upward pass for activation/freezing/autoconversion, followed by
  !           one downward precipitation/accretion pass.  This avoids replaying the
  !           activation/aerosol path and is the cleanest non-corrector diagnostic.
  !logical :: GF2M_USE_CORRECTOR = .true.

  !---------------------------------------------------------------------------
  ! Module constants.
  ! These are module-scope because both cup_up_moisture_2M and the helper
  ! procedures below use them.
  !---------------------------------------------------------------------------

  real, parameter :: qsmall_2m = 1.0e-12
  real, parameter :: nsmall_2m = 1.0e4        ! # m-3, numerical floor only
  real, parameter :: ncond_floor_cm3_liq = 5.0 ! # cm-3 floor if condensate exists
  real, parameter :: ncond_floor_cm3_ice = 0.01 ! # cm-3 floor if condensate exists

  real, parameter :: wmin_2m   = 0.5          ! m s-1
  real, parameter :: wmax_2m   = 20.0         ! m s-1
  real, parameter :: dtmax_2m  = 600.0        ! s
  real, parameter :: dpmin_2m  = 1.0e-6       ! hPa protection before Pa conversion

  real, parameter :: dust_frac_min_2m = 0.10
  real, parameter :: soot_frac_min_2m = 0.10
  real, parameter :: seasalt_kappa_min_2m = 1.00
  real, parameter :: kappa_wfa_min_2m = 0.20
  real, parameter :: eps_vap_2m       = 0.622

  ! INP source tuning parameters.  FDROP_* and FHET* come from
  ! GEOSmoist_Process_Library so the finite aerosol reservoirs can be
  ! adjusted without editing this diagnostic plume code.
  real :: fdrop_dust_2m
  real :: fdrop_soot_2m
  real, parameter :: fdrop_seasalt_2m = 1.0
  real, parameter :: acorr_dust_2m    = 2.7e7
  real, parameter :: acorr_soot_2m    = 8.0e7
  !real, parameter :: bkg_inp_scaling  = 0.001
  
  real(r8) :: bkg_inp_scaling 
  

  ! If enabled, droplet autoconversion removes the same fraction of the
  ! drop-mediated finite INP aerosol reservoirs that is removed from droplets.
  logical, parameter :: remove_inp_aerosol_with_liq_auto_2m = .true.

  real, parameter :: thom_2m  = 238.15        ! K homogeneous freezing T
  real, parameter :: dthom_2m = 0.05           ! K sigmoid width
  
  real(r8), parameter :: pi_r8       = 3.14159265358979323846_r8
  real(r8), parameter :: rho_liq_mg  = 1000._r8
  real(r8), parameter :: rho_ice_mg  = 917._r8
  real(r8), parameter :: rho_snow_2m = 600._r8

  ! Homogeneous-freezing ice-number closure.
  ! Homogeneous freezing can consume all selected liquid mass, but the
  ! resulting ice crystals are not assumed one-to-one with frozen droplets.
  ! Approximate the new ice number source using the mass ratio of a 12 um
  ! liquid drop to a 40 um ice particle.
  real(r8), parameter :: dliq_hom_number_ref_2m = 12.0e-6_r8  ! m diameter
  real(r8), parameter :: dice_hom_number_ref_2m = 40.0e-6_r8  ! m diameter
  real(r8), parameter :: hom_ice_number_ratio_2m =                  &
       (rho_liq_mg/rho_ice_mg) *                                    &
       (dliq_hom_number_ref_2m/dice_hom_number_ref_2m)**3

  ! Diagnostic convective precipitation distributions.
  ! Rain/snow use exponential Marshall-Palmer form, N(D)=N0 exp(-lambda D).
  real(r8), parameter :: lam_bnd_rain_2m(2) = 1._r8 / (/ 500.e-6_r8, 10.e-6_r8 /)
  real(r8), parameter :: lam_bnd_snow_2m(2) = 1._r8 / (/ 2000.e-6_r8, 20.e-6_r8 /)

  ! Terminal velocity V(D)=a D**b, SI units.
  real(r8), parameter :: arain_2m = 892._r8
  real(r8), parameter :: brain_2m = 0.80_r8
  real(r8), parameter :: asnow_2m = 11.8_r8
  real(r8), parameter :: bsnow_2m = 0.41_r8
  real(r8), parameter :: vt_min_2m = 0.10_r8

  ! Autoconversion precip-number source closure.
  ! Keep cloud-droplet/ice-number sinks separate from rain/snow number sources.
  ! Source rain/snow number is diagnosed from the autoconverted mass using an
  ! exponential Marshall-Palmer-like PSD with fixed intercepts.
  logical, parameter :: use_mp_auto_number_2m = .true.
  real(r8), parameter :: n0r_mp_auto_2m = 8.0e6_r8   ! m-4, rain source intercept
  real(r8), parameter :: n0s_mp_auto_2m = 2.0e7_r8   ! m-4, snow source intercept

  ! Liquid-only autoconversion efficiency.
  real(r8) :: auto_eff_liq_2m

  ! Existing fixed-precipitation accretion rates scaling
  logical, parameter :: use_accretion_eff_tuning_2m = .true.
  real(r8) :: accre_eff_rain_2m
  real(r8) :: accre_eff_snow_liq_2m
  real(r8) :: accre_eff_snow_ice_2m
  real(r8), parameter :: accre_eff_rain_snow_2m = 1.0_r8
  real(r8), parameter :: accre_eff_sip_2m       = 1.0_r8

  ! Aerosol activation source option.
  !   .true. : dNdrop_act = max(Ndrop_act - Ndrop, 0)
  !   .false.: additive behavior, dNdrop_act = Ndrop_act 
  logical, parameter :: use_activation_target_increment_2m = .true.

  ! Optional precip-number-only coalescence/aggregation.
  ! These conserve falling precip mass flux but reduce rain/snow number flux,
  logical, parameter :: use_mg_rain_self_collection_2m = .true.
  logical, parameter :: use_mg_snow_self_aggregation_2m = .true.
  real(r8), parameter :: rain_selfcollect_eff_2m = 1.00_r8
  real(r8), parameter :: snow_selfagg_eff_2m     = 1.00_r8
  real(r8), parameter :: rain_selfcollect_max_frac_2m = 0.95_r8
  real(r8), parameter :: snow_selfagg_max_frac_2m     = 0.95_r8
 
  
  

 
 contains

!-------------------------------------------------------------------------------------------

SUBROUTINE cup_up_moisture_2M(name,start_level,                                         &
                              ierr,ierrc,z_cup,qc,qrc,pw,pwav,hc,tempc,xland,                   &
                              AerP,jcol,flip,po,ktop,cd,dby,clw_all,                            &
                              t_cup,q,gamma_cup,zu,qes_cup,k22,qe_cup,                          &
                              zqexec,rho,                                                       &
                              up_massentr,up_massdetr,psum,psumh,x_add_buoy,                    &
                              zws,entr_rate,vvel2d_parent,vvel1d_parent,                        &
                              pice_2m,nliq_up_2m,nice_up_2m,qliq_up_2m,qice_up_2m,              &
                              nact_up_m3_out,                                                   &
                              itf,ktf,its,ite,kts,kte, use_linear_subcl_mf)

  implicit none

  !---------------------------------------------------------------------------
  ! Arguments
  !---------------------------------------------------------------------------

  integer, intent(in) :: itf, ktf, its, ite, kts, kte, use_linear_subcl_mf
  integer, intent(in) :: jcol
  integer, dimension(:), intent(in) :: flip

  character*(*)                    , intent(in) :: name
  integer, dimension(its:ite)      , intent(in) :: ktop, k22, start_level

  real, dimension(its:ite,kts:kte), intent(in) :: t_cup, rho, q, zu, gamma_cup,               &
                                                  qe_cup, hc, po, up_massentr, up_massdetr,  &
                                                  dby, qes_cup, z_cup, cd

  real, dimension(its:ite)        , intent(in) :: zqexec, xland, x_add_buoy
  real, dimension(its:ite)        , intent(in) :: zws
  real, dimension(its:ite,kts:kte), intent(in) :: entr_rate, vvel2d_parent
  real, dimension(its:ite)        , intent(in) :: vvel1d_parent

  type(AerPropsNew), dimension(:), intent(in) :: AerP

  integer, dimension(its:ite), intent(inout) :: ierr
  character*128, dimension(its:ite), intent(inout) :: ierrc

  real, dimension(its:ite,kts:kte), intent(out) :: qc, qrc, pw, clw_all, tempc
  real, dimension(its:ite)        , intent(out) :: pwav, psum, psumh

  real, dimension(its:ite,kts:kte), intent(out) :: pice_2m
  real, dimension(its:ite,kts:kte), intent(out) :: nliq_up_2m
  real, dimension(its:ite,kts:kte), intent(out) :: nice_up_2m
  real, dimension(its:ite,kts:kte), intent(out) :: qliq_up_2m
  real, dimension(its:ite,kts:kte), intent(out) :: qice_up_2m
  real, dimension(its:ite,kts:kte), intent(out) :: nact_up_m3_out

  !---------------------------------------------------------------------------
  ! Local variables
  !---------------------------------------------------------------------------

  integer :: i, k
  integer :: khost, nmodes
  integer :: kbase_2m
  integer :: m
  logical :: aer_reservoir_initialized
  logical :: debug_2m = .false.

  ! Local mode-sized aerosol reservoir used only for activation.
  ! AerP%num and aer_num_work_m3 are number concentrations (# m-3).
  ! This avoids copying the full AerP object while still preventing
  ! upper levels from reactivating aerosol already consumed below.
  real, dimension(size(AerP)) :: aer_num_work_m3   ! # m-3
  real, dimension(size(AerP)) :: aer_dpg_work_m
  real, dimension(size(AerP)) :: aer_sig_work
  real, dimension(size(AerP)) :: aer_kap_work

  real :: dz
  real :: qrch
  real :: qaver
  real :: denom
  real :: dp_pa

  real :: pice_mass
  real :: pliq_mass
  real :: fhom

  real :: w2_cup
  real :: w_cup
  real :: qrch_prev
  real :: qv_up_prev
  real :: qv_up_curr
  real :: qcond_prev
  real :: qcond_curr
  real :: nact_lev_mix
  real :: nact_target_mix
  real :: ncond_floor_mix
  real :: ndust_cb_mix
  real :: ddust_cb_m
  real :: sigdust_cb
  real :: ddust_mean_m
  real :: sigdust_mean
  real :: nsoot_cb_mix
  real :: dsoot_cb_m
  real :: sigsoot_cb
  real :: nseasalt_cb_mix
  real :: dseasalt_cb_m
  real :: sigseasalt_cb
  real :: nsoot_mix
  real :: dsoot_mean_m
  real :: sigsoot_mean
  real :: nseasalt_mix
  real :: dseasalt_mean_m
  real :: sigseasalt_mean
  real :: qliq_for_inp
  real :: nwfa_cb_m3
  real :: nwfa_now_m3
  real :: seasalt_imm_scale
  real :: nseasalt_for_inp_mix

  ! Number variables transported as mixing ratios: # kg-1.
  real, dimension(its:ite,kts:kte) :: nact_up_mix
  real, dimension(its:ite,kts:kte) :: nliq_up_mix
  real, dimension(its:ite,kts:kte) :: nice_up_mix
  real, dimension(its:ite,kts:kte) :: ndust_up_mix
  real, dimension(its:ite,kts:kte) :: ddust_up_m
  real, dimension(its:ite,kts:kte) :: sigdust_up
  real, dimension(its:ite,kts:kte) :: nsoot_up_mix
  real, dimension(its:ite,kts:kte) :: dsoot_up_m
  real, dimension(its:ite,kts:kte) :: sigsoot_up
  real, dimension(its:ite,kts:kte) :: nseasalt_up_mix
  real, dimension(its:ite,kts:kte) :: dseasalt_up_m
  real, dimension(its:ite,kts:kte) :: sigseasalt_up

  ! Explicit cloud-ice mass mixing ratio transported in the updraft: kg kg-1.
  real, dimension(its:ite,kts:kte) :: qice_up_mix

  real :: ncloud_mix0

  real :: qcloud_after_liq_auto
  real :: qcloud_after_all_auto

  real :: qice_mass
  real :: qliq_mass
  real :: qice_adv_mix
  real :: qfreeze_mass
  real :: qhom_mass
  real :: qdep_mass
  real :: qliq_avail

  real :: qc_liq_auto
  real :: qc_ice_auto

  real :: nliq_for_auto_mix
  real :: nliq_mix
  real :: nice_mix
  real :: nice_mix_before_ice_auto
  real :: ndust_mix
  real :: ninp_imm_mix
  real :: ninp_dust_mix
  real :: ninp_soot_mix
  real :: ninp_seasalt_mix
  real :: ninp_bkg_mix
  real :: dn_freeze_mix
  real :: dn_freeze_dust_mix
  real :: dn_freeze_soot_mix
  real :: dn_freeze_seasalt_mix
  real :: dn_freeze_bkg_mix
  real :: dn_hom_mix
  real :: dn_hom_ice_src_mix
  real :: nliq_before_freeze_mix
  real :: nliq_freeze_avail_mix

  real :: dq_auto_liq
  real :: dq_auto_ice
  ! Number bookkeeping for autoconversion:
  !   dn_auto_liq_mix      = cloud droplet sink (#/kg)
  !   dn_rain_auto_liq_mix = rain number source (#/kg)
  !   dn_auto_ice_mix      = cloud ice sink (#/kg)
  !   dn_snow_auto_ice_mix = snow number source (#/kg)
  real :: dn_auto_liq_mix
  real :: dn_rain_auto_liq_mix
  real :: dn_auto_ice_mix
  real :: dn_snow_auto_ice_mix
  real :: frac_liq_auto

  real :: prc_liu_kgkg_s
  real :: nprc1_liu

  real :: dt_plume

  real :: pres_pa_2m
  real :: qvl_sat_2m
  real :: qvi_sat_2m
  real :: esi_pa_2m
  real :: qliq_det_2m
  real :: qice_det_2m
  
  
  
  !diagnostic vars
  real :: qprecip_mix
  real :: qtotal_for_ratio
  real :: precip_to_total
  real :: nliq_cm3
  real :: nice_L
  real :: nact_lev_cm3
  real :: nwfa_cm3
  real :: ndust_cm3
  real :: ninp_imm_L

  ! Column storage for split precipitation/accretion pass.
  real, dimension(kts:kte) :: qliq_col
  real, dimension(kts:kte) :: qice_col
  real, dimension(kts:kte) :: nliq_col
  real, dimension(kts:kte) :: nice_col

  real, dimension(kts:kte) :: dq_auto_liq_col
  real, dimension(kts:kte) :: dq_auto_ice_col
  real, dimension(kts:kte) :: dn_auto_liq_col
  real, dimension(kts:kte) :: dn_auto_ice_col
  real, dimension(kts:kte) :: dn_rain_auto_liq_col
  real, dimension(kts:kte) :: dn_snow_auto_ice_col

  real, dimension(kts:kte) :: dt_col
  real, dimension(kts:kte) :: dz_col
  real, dimension(kts:kte) :: qrch_col
  real, dimension(kts:kte) :: w_col

  real, dimension(kts:kte) :: qrain_diag_col
  real, dimension(kts:kte) :: qsnow_diag_col
  real, dimension(kts:kte) :: nrain_diag_col
  real, dimension(kts:kte) :: nsnow_diag_col
  real, dimension(kts:kte) :: lamr_col
  real, dimension(kts:kte) :: n0r_col
  real, dimension(kts:kte) :: lams_col
  real, dimension(kts:kte) :: n0s_col
  real, dimension(kts:kte) :: umr_col
  real, dimension(kts:kte) :: ums_col
  real, dimension(kts:kte) :: unr_col
  real, dimension(kts:kte) :: uns_col
  real, dimension(kts:kte) :: nact_lev_mix_col
  real, dimension(kts:kte) :: nact_target_mix_col
  real, dimension(kts:kte) :: ndust_mix_col
  real, dimension(kts:kte) :: ninp_imm_mix_col
  real, dimension(kts:kte) :: nwfa_cm3_col
  real, dimension(kts:kte) :: seasalt_imm_scale_col

  real, dimension(kts:kte) :: dq_acc_liq_rain_col
  real, dimension(kts:kte) :: dq_acc_liq_snow_col
  real, dimension(kts:kte) :: dq_acc_ice_snow_col
  real, dimension(kts:kte) :: dq_sip_col
  real, dimension(kts:kte) :: dq_rain_to_snow_col

  real :: rho_safe
  real :: dt_eff
  real :: dz_eff

  real :: rain_flux_m
  real :: snow_flux_m
  real :: rain_flux_n
  real :: snow_flux_n
  real :: rain_flux_n_before_self
  real :: snow_flux_n_before_self
  real :: dn_rain_self_mix
  real :: dn_snow_self_mix

  real :: qrain_diag
  real :: qsnow_diag
  real :: nrain_diag
  real :: nsnow_diag

  real :: lamr_2m
  real :: n0r_2m
  real :: lams_2m
  real :: n0s_2m
  real :: umr_2m
  real :: ums_2m
  real :: unr_2m
  real :: uns_2m

  real :: dq_cw_rain
  real :: dn_cw_rain
  real :: dq_cw_snow_total
  real :: dq_cw_snow_to_snow
  real :: dn_cw_snow
  real :: dq_ci_snow
  real :: dn_ci_snow
  real :: dq_sip
  real :: dn_sip
  real :: dq_rs
  real :: dn_rs

  real :: liq_sink_mass
  real :: liq_sink_num
  real :: ice_sink_mass
  real :: ice_sink_num
  real :: scale_sink
  real :: dq_precip_level
  logical :: microp_uniform_acc
  
  !---------------------------------------------------------------------------
  ! MG one-column local interface
  !---------------------------------------------------------------------------

  integer, parameter :: mgncol_1 = 1
  logical :: microp_uniform_liu

  type(MGHydrometeorProps) :: liq_props_mg
  type(MGHydrometeorProps) :: ice_props_mg
  real(r8), dimension(2) :: liq_lambda_bounds_mg
  real(r8), dimension(2) :: ice_lambda_bounds_mg


  real :: ac_time_mg
  real(r8) :: dcs_mg

  ! Liquid 2M variables
  real(r8) :: qcic_mg
  real(r8) :: ncic_mg
  real(r8) :: rho_mg
  real(r8) :: pgam_mg
  real(r8) :: lamc_mg
  real(r8) :: relvar_mg

  real(r8), dimension(mgncol_1) :: pgam_arr
  real(r8), dimension(mgncol_1) :: qcic_arr
  real(r8), dimension(mgncol_1) :: ncic_arr
  real(r8), dimension(mgncol_1) :: rho_arr
  real(r8), dimension(mgncol_1) :: relvar_arr
  real(r8), dimension(mgncol_1) :: prc_arr
  real(r8), dimension(mgncol_1) :: nprc_arr
  real(r8), dimension(mgncol_1) :: nprc1_arr

  ! Ice 2M variables
  real(r8) :: qiic_mg
  real(r8) :: niic_mg
  real(r8) :: lami_mg
  real(r8) :: tice_mg

  real(r8), dimension(mgncol_1) :: tice_arr
  real(r8), dimension(mgncol_1) :: qiic_arr
  real(r8), dimension(mgncol_1) :: niic_arr
  real(r8), dimension(mgncol_1) :: lami_arr
  real(r8), dimension(mgncol_1) :: prci_arr
  real(r8), dimension(mgncol_1) :: nprci_arr

  ! MG ice deposition/sublimation arrays
  real(r8), dimension(mgncol_1) :: t_dep_arr
  real(r8), dimension(mgncol_1) :: qv_dep_arr
  real(r8), dimension(mgncol_1) :: qi_dep_arr
  real(r8), dimension(mgncol_1) :: ni_dep_arr
  real(r8), dimension(mgncol_1) :: icldm_dep_arr
  real(r8), dimension(mgncol_1) :: rho_dep_arr
  real(r8), dimension(mgncol_1) :: dv_dep_arr
  real(r8), dimension(mgncol_1) :: qvl_dep_arr
  real(r8), dimension(mgncol_1) :: qvi_dep_arr
  real(r8), dimension(mgncol_1) :: berg_dep_arr
  real(r8), dimension(mgncol_1) :: vap_dep_arr
  real(r8), dimension(mgncol_1) :: ice_sublim_arr

  ! MG accretion/collection arrays
  real(r8), dimension(mgncol_1) :: t_acc_arr
  real(r8), dimension(mgncol_1) :: mu_acc_arr
  real(r8), dimension(mgncol_1) :: asn_arr
  real(r8), dimension(mgncol_1) :: accre_enhan_arr

  real(r8), dimension(mgncol_1) :: qric_arr
  real(r8), dimension(mgncol_1) :: qsic_arr
  real(r8), dimension(mgncol_1) :: qliq_acc_arr
  real(r8), dimension(mgncol_1) :: nliq_acc_arr
  real(r8), dimension(mgncol_1) :: qice_acc_arr
  real(r8), dimension(mgncol_1) :: nice_acc_arr

  real(r8), dimension(mgncol_1) :: lamc_arr
  real(r8), dimension(mgncol_1) :: lamr_arr
  real(r8), dimension(mgncol_1) :: n0r_arr
  real(r8), dimension(mgncol_1) :: lams_arr
  real(r8), dimension(mgncol_1) :: n0s_arr

  real(r8), dimension(mgncol_1) :: umr_arr
  real(r8), dimension(mgncol_1) :: ums_arr
  real(r8), dimension(mgncol_1) :: unr_arr
  real(r8), dimension(mgncol_1) :: uns_arr

  real(r8), dimension(mgncol_1) :: pra_arr
  real(r8), dimension(mgncol_1) :: npra_arr
  real(r8), dimension(mgncol_1) :: psacws_arr
  real(r8), dimension(mgncol_1) :: npsacws_arr
  real(r8), dimension(mgncol_1) :: pracs_arr
  real(r8), dimension(mgncol_1) :: npracs_arr
  real(r8), dimension(mgncol_1) :: prai_arr
  real(r8), dimension(mgncol_1) :: nprai_arr
  real(r8), dimension(mgncol_1) :: msacwi_arr
  real(r8), dimension(mgncol_1) :: nsacwi_arr


     !---------------------------------------------------------------------------
    ! MG hydrometeor properties.
    !
    ! Use the MG constructor so shape_coef is initialized exactly as in MG:
    !   shape_coef = rho*pi*gamma(eff_dim+1)/6
    !
    ! For eff_dim = 3, this gives shape_coef = rho*pi.
    !---------------------------------------------------------------------------

    dcs_mg = real(DCS_CNV, r8)
    auto_eff_liq_2m       = real(AUT_SCALE_CNV, r8)
    accre_eff_rain_2m     = real(ACC_ENH_CNV, r8)
    accre_eff_snow_liq_2m = real(ACC_ENH_CNV, r8)
    accre_eff_snow_ice_2m = real(ACC_ENH_ICE, r8)
    fdrop_dust_2m         = FDROPDUST*FHETDUST
    fdrop_soot_2m         = FDROPSOOT*FHETSOOT
    bkg_inp_scaling       = real(BKG_INP_SC_CNV, r8) 
    
    liq_lambda_bounds_mg = 1._r8 / (/ 50.e-6_r8, 2.e-6_r8 /)

    ! mean ice diameter cannot grow bigger than twice the snow-autoconversion
    ! threshold.
    ice_lambda_bounds_mg = 1._r8 / (/ 2._r8*dcs_mg, 1.e-6_r8 /)

    liq_props_mg = MGHydrometeorProps( &
         rho_liq_mg, 3._r8, liq_lambda_bounds_mg, 0._r8 )

    ice_props_mg = MGHydrometeorProps( &
         rho_ice_mg, 3._r8, ice_lambda_bounds_mg, 0._r8 )

      microp_uniform_liu = .false.
      microp_uniform_acc = .true.
      ac_time_mg         = real(ICE_AUTO_TSC_CNV, r8)   
      debug_2m = DEBUG_GF2M

    ! The full 2M physics is controlled by GEOSmoist_Process_Library tunings
    ! imported above.  The call signature intentionally excludes legacy Kessler
    ! controls such as C0, AUTOCONV, and qrc_crit so they cannot accidentally
    ! override the 2M mass pathway.

  !---------------------------------------------------------------------------
  ! Initialization
  !---------------------------------------------------------------------------

  do i = its, itf
     pwav(i)  = 0.0
     psum(i)  = 0.0
     psumh(i) = 0.0
  enddo

  do k = kts, ktf
     do i = its, itf
        pw(i,k)             = 0.0
        clw_all(i,k)        = 0.0
        tempc(i,k)          = t_cup(i,k)
        qrc(i,k)            = 0.0
        qc(i,k)             = qe_cup(i,k)

        nact_up_mix(i,k)    = 0.0
        nliq_up_mix(i,k)    = 0.0
        nice_up_mix(i,k)    = 0.0
        ndust_up_mix(i,k)   = 0.0
        ddust_up_m(i,k)     = 0.0
        sigdust_up(i,k)     = 0.0
        nsoot_up_mix(i,k)   = 0.0
        dsoot_up_m(i,k)     = 0.0
        sigsoot_up(i,k)     = 0.0
        nseasalt_up_mix(i,k)= 0.0
        dseasalt_up_m(i,k)  = 0.0
        sigseasalt_up(i,k)  = 0.0
        qice_up_mix(i,k)    = 0.0
        nact_up_m3_out(i,k) = 0.0

        pice_2m(i,k)        = 0.0
        nliq_up_2m(i,k)     = 0.0
        nice_up_2m(i,k)     = 0.0
        qliq_up_2m(i,k)     = 0.0
        qice_up_2m(i,k)     = 0.0
     enddo
  enddo

  !---------------------------------------------------------------------------
  ! Boundary condition for total water in the updraft.
  !---------------------------------------------------------------------------

  do i = its, itf
     if(ierr(i) /= 0) cycle

     call get_cloud_bc(name, kts, kte, ktf, xland(i), po(i,kts:kte), &
                       qe_cup(i,kts:kte), qaver, k22(i))

     qc(i,kts:start_level(i))  = qaver + zqexec(i) + 0.5*x_add_buoy(i)/xlv
     qrc(i,kts:start_level(i)) = 0.0
  enddo

  !---------------------------------------------------------------------------
  ! Optional linear subcloud mixing for shallow convection.
  !---------------------------------------------------------------------------

  if(trim(name) == 'shallow' .and. use_linear_subcl_mf == 1) then
     do i = its, itf
        if(ierr(i) /= 0) cycle

        call get_delmix_2m(name, kts, kte, ktf, xland(i), start_level(i), &
                        po(i,kts:kte), qe_cup(i,kts:kte), qc(i,kts:kte))
     enddo
  endif

  !---------------------------------------------------------------------------
  ! Main updraft moisture / number loop.
  !---------------------------------------------------------------------------

  do i = its, itf
     if(ierr(i) /= 0) cycle

         !-----------------------------------------------------------------------
         ! Effective 2M cloud base.
         !
         ! cup_up_moisture_2M only visits levels from start_level(i)+1 upward.
         ! Do not initialize aerosol/number state at kbcon(i), because kbcon(i)
         ! can be below the first level visited here.  For this diagnostic 2M
         ! microphysics, the first loop level is the local cloud base.
         !-----------------------------------------------------------------------

     kbase_2m = start_level(i) + 1
     if(kbase_2m > ktop(i) + 1) cycle

     qliq_col(:) = 0.0
     qice_col(:) = 0.0
     nliq_col(:) = 0.0
     nice_col(:) = 0.0

     dq_auto_liq_col(:) = 0.0
     dq_auto_ice_col(:) = 0.0
     dn_auto_liq_col(:) = 0.0
     dn_auto_ice_col(:) = 0.0
     dn_rain_auto_liq_col(:) = 0.0
     dn_snow_auto_ice_col(:) = 0.0

     dt_col(:)   = 0.0
     dz_col(:)   = 0.0
     qrch_col(:) = 0.0
     w_col(:)    = 0.0

     qrain_diag_col(:) = 0.0
     qsnow_diag_col(:) = 0.0
     nrain_diag_col(:) = 0.0
     nsnow_diag_col(:) = 0.0
     lamr_col(:) = 0.0
     n0r_col(:) = 0.0
     lams_col(:) = 0.0
     n0s_col(:) = 0.0
     umr_col(:) = 0.0
     ums_col(:) = 0.0
     unr_col(:) = 0.0
     uns_col(:) = 0.0
     nact_lev_mix_col(:) = 0.0
     nact_target_mix_col(:) = 0.0
     ndust_mix_col(:) = 0.0
     ninp_imm_mix_col(:) = 0.0
     nwfa_cm3_col(:) = 0.0
     seasalt_imm_scale_col(:) = 1.0

     dq_acc_liq_rain_col(:) = 0.0
     dq_acc_liq_snow_col(:) = 0.0
     dq_acc_ice_snow_col(:) = 0.0
     dq_sip_col(:)          = 0.0
     dq_rain_to_snow_col(:) = 0.0

     nact_lev_mix    = 0.0
     nact_target_mix = 0.0
     nwfa_cb_m3      = 0.0

     if(size(AerP) > 0) then
        nmodes = min(size(AerP), max(0, AerP(1)%nmods))
     else
        nmodes = 0
     endif

     aer_reservoir_initialized = .false.
     aer_num_work_m3(:) = 0.0
     aer_dpg_work_m(:)  = 0.0
     aer_sig_work(:)    = 0.0
     aer_kap_work(:)    = 0.0

     ! Initialize local updraft kinetic energy.  The first 2M loop level is
     ! treated as cloud base.
     w2_cup = max(wmin_2m*wmin_2m, zws(i)*zws(i))
     if((GF2M_W_OPTION == 2 .or. GF2M_W_OPTION == 3) .and. vvel1d_parent(i) > 0.0) then
        w2_cup = max(w2_cup, vvel1d_parent(i)*vvel1d_parent(i))
     endif
     w_cup  = cup_2m_w_from_w2(w2_cup)

     do k = kbase_2m, ktop(i) + 1

        dz = z_cup(i,k) - z_cup(i,k-1)

        !--------------------------------------------------------------------
        ! Saturation vapor value in the updraft.
        !--------------------------------------------------------------------

        qrch = qes_cup(i,k) + (1.0/xlv) * &
               (gamma_cup(i,k)/(1.0 + gamma_cup(i,k))) * dby(i,k)

        !--------------------------------------------------------------------
        ! Steady plume equation for total water and cloud condensate.
        !--------------------------------------------------------------------

        denom = zu(i,k-1) - 0.5*up_massdetr(i,k-1) + up_massentr(i,k-1)

        if(denom > 0.0) then

           qc(i,k) = ( qc(i,k-1)*zu(i,k-1)                              &
                     - 0.5*up_massdetr(i,k-1)*qc(i,k-1)                  &
                     + up_massentr(i,k-1)*q(i,k-1) ) / denom

           if(k == kbase_2m) then
              qc(i,k) = qc(i,k) + zqexec(i)*up_massentr(i,k-1)/denom
           endif

           qrc(i,k) = ( qrc(i,k-1)*zu(i,k-1)                             &
                      - 0.5*up_massdetr(i,k-1)*qrc(i,k-1) ) / denom

           qice_up_mix(i,k) = ( qice_up_mix(i,k-1)*zu(i,k-1)              &
                              - 0.5*up_massdetr(i,k-1)*qice_up_mix(i,k-1) &
                              ) / denom

        else

           qc(i,k)          = qc(i,k-1)
           qrc(i,k)         = qrc(i,k-1)
           qice_up_mix(i,k) = qice_up_mix(i,k-1)

        endif

        ! Updraft temperature.
        tempc(i,k) = (1.0/cp) * (hc(i,k) - g*z_cup(i,k) - xlv*qrch)

        !--------------------------------------------------------------------
        ! Local diagnostic updraft velocity.
        !
        ! w2_cup carries w^2 upward.  At cloud base the speed is the
        ! initialized subcloud value.  Above cloud base, advance w^2 from
        ! k-1 to k using buoyancy, condensate loading, and
        ! entrainment/detrainment damping.
        !--------------------------------------------------------------------

        if(k <= kbase_2m) then
           w_cup = cup_2m_w_from_w2(w2_cup)
        else
           qrch_prev = qes_cup(i,k-1) + (1.0/xlv) *                       &
                       (gamma_cup(i,k-1)/(1.0 + gamma_cup(i,k-1))) *       &
                       dby(i,k-1)

           qv_up_prev  = min(qc(i,k-1), qrch_prev)
           qv_up_curr  = min(qc(i,k  ), qrch)
           qcond_prev  = max(0.0, qc(i,k-1) - qrch_prev)
           qcond_curr  = max(0.0, qc(i,k  ) - qrch)

           w2_cup = cup_2m_next_w2( w2_cup,                                &
                                    z_cup(i,k-1), z_cup(i,k),              &
                                    t_cup(i,k-1), t_cup(i,k),              &
                                    q(i,k-1),     q(i,k),                  &
                                    tempc(i,k-1), tempc(i,k),              &
                                    qv_up_prev,   qv_up_curr,              &
                                    qcond_prev,   qcond_curr,              &
                                    entr_rate(i,k-1), cd(i,k-1) )

           w_cup = cup_2m_w_from_w2(w2_cup)
        endif

        ! Optional use of the parent legacy updraft velocity profile.
        ! This affects activation/residence-time-dependent 2M microphysics only;
        ! the parent still computes the authoritative mass tendencies.
        if((GF2M_W_OPTION == 2 .or. GF2M_W_OPTION == 3) .and. vvel2d_parent(i,k) > 0.0) then
           if(GF2M_W_OPTION == 2) then
              w_cup = max(wmin_2m, min(wmax_2m, vvel2d_parent(i,k)))
           else
              w_cup = max(w_cup, max(wmin_2m, min(wmax_2m, vvel2d_parent(i,k))))
           endif
        endif

        !--------------------------------------------------------------------
        ! Vertical transport of cloud/dust number plus in-plume activation.
        !
        ! Transported number variables are # kg-1:
        !   nliq_up_mix
        !   nice_up_mix
        !   ndust_up_mix
        !
        ! Activation is a point source of liquid droplets at each level.
        ! Only the droplet reservoir is transported upward. Newly activated
        ! droplets are added locally after transport:
        !
        !   nliq(k) = transported nliq(k-1) + activated_from_remaining_aerosol(k)
        !
        ! AerP itself remains intent(in).  The only mutable aerosol state is
        ! aer_num_work_m3(:), a small local mode-sized reservoir initialized at
        ! cloud base and depleted by the activation helper.
        !
        ! Aerosol, not INP, is transported. INP is diagnosed locally from the
        ! transported dust number, transported dust size, and local temperature.
        !--------------------------------------------------------------------

        nact_lev_mix = 0.0

        !--------------------------------------------------------------------
        ! Number transport plus local aerosol activation.
        !
        ! kbase_2m = start_level(i)+1 is the effective cloud base for this
        ! routine.  This version
        ! intentionally does not entrain fresh aerosol aloft.
        !
        ! Transported hydrometeor/dust numbers are # kg-1.
        ! aer_num_work_m3 is # m-3 and is only the local activation reservoir.
        !--------------------------------------------------------------------

        if(k == kbase_2m) then

           khost = flip(kbase_2m)

           if(nmodes > 0) then
              do m = 1, nmodes
                 aer_num_work_m3(m) = max(0.0,    AerP(m)%num(i,jcol,khost))
                 aer_dpg_work_m(m)  = max(1.0e-9, AerP(m)%dpg(i,jcol,khost))
                 aer_sig_work(m)    = max(1.0e-6, AerP(m)%sig(i,jcol,khost))
                 aer_kap_work(m)    = max(0.0,    AerP(m)%kap(i,jcol,khost))
              enddo

              if(nmodes < size(AerP)) then
                 aer_num_work_m3(nmodes+1:size(AerP)) = 0.0
                 aer_dpg_work_m (nmodes+1:size(AerP)) = 0.0
                 aer_sig_work   (nmodes+1:size(AerP)) = 0.0
                 aer_kap_work   (nmodes+1:size(AerP)) = 0.0
              endif

              aer_reservoir_initialized = .true.
              ! Cloud-base water-friendly aerosol reservoir before activation.
              ! Used below as a proxy for sea-salt removal by droplet activation.
              nwfa_cb_m3 = cup_2m_nwfa_m3(nmodes, aer_num_work_m3, aer_kap_work)
           else
              aer_reservoir_initialized = .false.
              nwfa_cb_m3 = 0.0
           endif

           ! At the effective cloud base, start with no transported cloud
           ! number.  Activation below adds droplets locally.
           nliq_up_mix(i,k)  = 0.0
           nice_up_mix(i,k)  = 0.0
           qice_up_mix(i,k)  = 0.0

           ! Finite INP-active aerosol reservoirs are initialized at the same
           ! effective cloud base and then transported upward.  INP is diagnosed
           ! locally from transported dust, soot, and sea-salt reservoirs.
           if(nmodes > 0) then
              call aerosol_ice_species_from_aerp( &
                   AerP, nmodes, i, jcol, khost, rho(i,k), &
                   ndust_cb_mix, ddust_cb_m, sigdust_cb, &
                   nsoot_cb_mix, dsoot_cb_m, sigsoot_cb, &
                   nseasalt_cb_mix, dseasalt_cb_m, sigseasalt_cb )
           else
              ndust_cb_mix    = 0.0
              ddust_cb_m      = 0.0
              sigdust_cb      = 0.0
              nsoot_cb_mix    = 0.0
              dsoot_cb_m      = 0.0
              sigsoot_cb      = 0.0
              nseasalt_cb_mix = 0.0
              dseasalt_cb_m   = 0.0
              sigseasalt_cb   = 0.0
           endif

           ndust_up_mix(i,k)    = max(0.0, ndust_cb_mix)
           ddust_up_m(i,k)      = max(0.0, ddust_cb_m)
           sigdust_up(i,k)      = max(0.0, sigdust_cb)
           nsoot_up_mix(i,k)    = max(0.0, nsoot_cb_mix)
           dsoot_up_m(i,k)      = max(0.0, dsoot_cb_m)
           sigsoot_up(i,k)      = max(0.0, sigsoot_cb)
           nseasalt_up_mix(i,k) = max(0.0, nseasalt_cb_mix)
           dseasalt_up_m(i,k)   = max(0.0, dseasalt_cb_m)
           sigseasalt_up(i,k)   = max(0.0, sigseasalt_cb)

        elseif(k > kbase_2m) then

           if(denom > 0.0) then

              nliq_up_mix(i,k) = ( nliq_up_mix(i,k-1)*zu(i,k-1)              &
                                 - 0.5*up_massdetr(i,k-1)*nliq_up_mix(i,k-1) &
                                 ) / denom

              nice_up_mix(i,k) = ( nice_up_mix(i,k-1)*zu(i,k-1)              &
                                 - 0.5*up_massdetr(i,k-1)*nice_up_mix(i,k-1) &
                                 ) / denom

              ndust_up_mix(i,k) = ( ndust_up_mix(i,k-1)*zu(i,k-1)              &
                                  - 0.5*up_massdetr(i,k-1)*ndust_up_mix(i,k-1) &
                                  ) / denom

              nsoot_up_mix(i,k) = ( nsoot_up_mix(i,k-1)*zu(i,k-1)              &
                                  - 0.5*up_massdetr(i,k-1)*nsoot_up_mix(i,k-1) &
                                  ) / denom

              nseasalt_up_mix(i,k) = ( nseasalt_up_mix(i,k-1)*zu(i,k-1)              &
                                     - 0.5*up_massdetr(i,k-1)*nseasalt_up_mix(i,k-1) &
                                     ) / denom

           else

              nliq_up_mix(i,k)  = nliq_up_mix(i,k-1)
              nice_up_mix(i,k)  = nice_up_mix(i,k-1)
              ndust_up_mix(i,k)    = ndust_up_mix(i,k-1)
              nsoot_up_mix(i,k)    = nsoot_up_mix(i,k-1)
              nseasalt_up_mix(i,k) = nseasalt_up_mix(i,k-1)

           endif

           nliq_up_mix(i,k)  = max(0.0, nliq_up_mix(i,k))
           nice_up_mix(i,k)  = max(0.0, nice_up_mix(i,k))
           ndust_up_mix(i,k)    = max(0.0, ndust_up_mix(i,k))
           nsoot_up_mix(i,k)    = max(0.0, nsoot_up_mix(i,k))
           nseasalt_up_mix(i,k) = max(0.0, nseasalt_up_mix(i,k))

           ddust_up_m(i,k)      = ddust_up_m(i,k-1)
           sigdust_up(i,k)      = sigdust_up(i,k-1)
           dsoot_up_m(i,k)      = dsoot_up_m(i,k-1)
           sigsoot_up(i,k)      = sigsoot_up(i,k-1)
           dseasalt_up_m(i,k)   = dseasalt_up_m(i,k-1)
           sigseasalt_up(i,k)   = sigseasalt_up(i,k-1)

        else

           nliq_up_mix(i,k)  = 0.0
           nice_up_mix(i,k)  = 0.0
           ndust_up_mix(i,k)    = 0.0
           ddust_up_m(i,k)      = 0.0
           sigdust_up(i,k)      = 0.0
           nsoot_up_mix(i,k)    = 0.0
           dsoot_up_m(i,k)      = 0.0
           sigsoot_up(i,k)      = 0.0
           nseasalt_up_mix(i,k) = 0.0
           dseasalt_up_m(i,k)   = 0.0
           sigseasalt_up(i,k)   = 0.0
           qice_up_mix(i,k)     = 0.0

        endif

        ! Point-source activation from the remaining cloud-base aerosol
        ! reservoir.  The helper depletes aer_num_work_m3(:), not AerP.
        ! Activation behavior is controlled by use_activation_target_increment_2m:
        !   .true. : dNdrop_act = max(Ndrop_act - Ndrop, 0)
        !   .false.: dNdrop_act = Ndrop_act  (legacy additive source)
        nact_target_mix = 0.0
        if(k >= kbase_2m .and. aer_reservoir_initialized .and. nmodes > 0) then

           call activate_from_aerp_arg2002_local( &
                nmodes, aer_num_work_m3, aer_dpg_work_m, aer_sig_work, aer_kap_work, &
                w_cup, tempc(i,k), rho(i,k), nliq_up_mix(i,k), &
                nact_target_mix, nact_lev_mix )

           nliq_up_mix(i,k) = max(0.0, nliq_up_mix(i,k) + nact_lev_mix)

        else

           nact_lev_mix = 0.0

        endif

        ! Scale the finite sea-salt reservoir available for immersion freezing by
        ! the remaining soluble/water-friendly aerosol reservoir.  The transported
        ! sea-salt reservoir itself is depleted only by actual freezing below.
        if(aer_reservoir_initialized .and. nmodes > 0 .and. nwfa_cb_m3 > 0.0) then
           nwfa_now_m3 = cup_2m_nwfa_m3(nmodes, aer_num_work_m3, aer_kap_work)
           seasalt_imm_scale = max(0.0, min(1.0, nwfa_now_m3/nwfa_cb_m3))
        else
           nwfa_now_m3 = 0.0
           seasalt_imm_scale = 1.0
        endif
        seasalt_imm_scale_col(k) = seasalt_imm_scale

        !--------------------------------------------------------------------
        ! Total condensate before autoconversion.
        !--------------------------------------------------------------------

        clw_all(i,k) = max(0.0, qc(i,k) - qrch)
        qrc(i,k)     = min(clw_all(i,k), qrc(i,k))

        qice_up_mix(i,k) = max(0.0, min(qice_up_mix(i,k), qrc(i,k), clw_all(i,k)))

        !--------------------------------------------------------------------
        ! Safeguard: condensate cannot exist with zero cloud number.
        !
        ! If activation failed or the cloud-base aerosol reservoir was zero,
        ! impose a small droplet-number floor so the plume never carries
        ! qcond > 0 with Nl+Ni = 0.  
        !--------------------------------------------------------------------

        if(clw_all(i,k) > qsmall_2m) then
           ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / max(rho(i,k), 1.0e-12)

           if(nliq_up_mix(i,k) + nice_up_mix(i,k) < ncond_floor_mix) then
              nliq_up_mix(i,k) = max(nliq_up_mix(i,k), ncond_floor_mix)
           endif
           
           
        endif

        !--------------------------------------------------------------------
        ! Residence time for this layer.
        ! Used by deposition growth and autoconversion.
        !--------------------------------------------------------------------

        dt_plume = dz / max(w_cup, wmin_2m)
        dt_plume = max(0.0, min(dtmax_2m, dt_plume))

        !--------------------------------------------------------------------
        ! Number mixing ratios before autoconversion and immersion freezing.
        !--------------------------------------------------------------------

        nliq_mix       = nliq_up_mix(i,k)
        nice_mix       = nice_up_mix(i,k)
        ndust_mix      = ndust_up_mix(i,k)
        ddust_mean_m   = ddust_up_m(i,k)
        sigdust_mean   = sigdust_up(i,k)
        nsoot_mix      = nsoot_up_mix(i,k)
        dsoot_mean_m   = dsoot_up_m(i,k)
        sigsoot_mean   = sigsoot_up(i,k)
        nseasalt_mix   = nseasalt_up_mix(i,k)
        dseasalt_mean_m= dseasalt_up_m(i,k)
        sigseasalt_mean= sigseasalt_up(i,k)

        nliq_before_freeze_mix = nliq_mix

        !--------------------------------------------------------------------
        ! Heterogeneous freezing number diagnosed in place.
        !
        ! INP is diagnosed from transported dust, soot, and sea-salt aerosol
        ! reservoirs plus an infinite background source.  Dust/soot include
        ! immersion plus a pseudo-contact source.  Sea salt is a proxy for biological INP.
        ! Finite aerosol reservoirs are depleted by the number nucleated.
        !--------------------------------------------------------------------

        qliq_for_inp = max(0.0, clw_all(i,k) - qice_up_mix(i,k))
        nseasalt_for_inp_mix = max(0.0, nseasalt_mix * seasalt_imm_scale)

        ! Do not let the artificial condensate-number floor trigger freezing.
        ! INP freezing is evaluated only against droplets above the numerical
        ! liquid-number floor, and only if liquid mass is actually present.
        ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / max(rho(i,k), 1.0e-12)
        nliq_freeze_avail_mix = max(0.0, nliq_mix - ncond_floor_mix)

        if(qliq_for_inp > qsmall_2m .and. nliq_freeze_avail_mix > 1.0) then
           call cup_2m_inp_freezing_sources(                               &
                tempc(i,k), po(i,k), rho(i,k), w_cup, dt_plume,            &
                qliq_for_inp, nliq_freeze_avail_mix,                       &
                ndust_mix, ddust_mean_m, sigdust_mean,                     &
                nsoot_mix, dsoot_mean_m, sigsoot_mean,                     &
                nseasalt_for_inp_mix, dseasalt_mean_m, sigseasalt_mean,    &
                ninp_dust_mix, ninp_soot_mix, ninp_seasalt_mix,            &
                ninp_bkg_mix, ninp_imm_mix )
        else
           ninp_dust_mix    = 0.0
           ninp_soot_mix    = 0.0
           ninp_seasalt_mix = 0.0
           ninp_bkg_mix     = 0.0
           ninp_imm_mix     = 0.0
        endif

        dn_freeze_dust_mix    = ninp_dust_mix
        dn_freeze_soot_mix    = ninp_soot_mix
        dn_freeze_seasalt_mix = ninp_seasalt_mix
        dn_freeze_bkg_mix     = ninp_bkg_mix
        dn_freeze_mix         = ninp_imm_mix

        if(dn_freeze_mix > nliq_mix .and. dn_freeze_mix > 1.0e-30) then
           scale_sink = max(0.0, min(1.0, nliq_mix/dn_freeze_mix))
           dn_freeze_dust_mix    = dn_freeze_dust_mix    * scale_sink
           dn_freeze_soot_mix    = dn_freeze_soot_mix    * scale_sink
           dn_freeze_seasalt_mix = dn_freeze_seasalt_mix * scale_sink
           dn_freeze_bkg_mix     = dn_freeze_bkg_mix     * scale_sink
           dn_freeze_mix         = nliq_mix
        endif

        nliq_mix       = max(0.0, nliq_mix  - dn_freeze_mix)
        nice_mix       = max(0.0, nice_mix  + dn_freeze_mix)
        ndust_mix      = max(0.0, ndust_mix - dn_freeze_dust_mix)
        nsoot_mix      = max(0.0, nsoot_mix - dn_freeze_soot_mix)
        nseasalt_mix   = max(0.0, nseasalt_mix - dn_freeze_seasalt_mix)

        nliq_up_mix(i,k)     = nliq_mix
        nice_up_mix(i,k)     = nice_mix
        ndust_up_mix(i,k)    = ndust_mix
        nsoot_up_mix(i,k)    = nsoot_mix
        nseasalt_up_mix(i,k) = nseasalt_mix

        ncloud_mix0 = max(0.0, nliq_mix + nice_mix)
        nact_up_mix(i,k) = ncloud_mix0
        nact_up_m3_out(i,k) = nact_up_mix(i,k) * max(rho(i,k), 1.0e-12)

        !--------------------------------------------------------------------
        ! Explicit ice mass fraction.
        !
        ! Ice mass is carried as qice_up_mix and grown by deposition.
        ! This replaces the T-based mass partition.
        !
        ! Initial ice mass from immersion freezing:
        !   freeze the same fraction of remaining liquid condensate mass as
        !   the frozen fraction of liquid droplet number.
        !
        ! Homogeneous freezing:
        !   complete at and below thom_2m, sigmoid transition above.
        !
        ! Deposition growth:
        !   uses MG ice_deposition_sublimation assuming water saturation.
        !   qv  = qvl
        !   qvl = water saturation
        !   qvi = ice saturation
        !   sublimation is clipped to zero.
        !--------------------------------------------------------------------

        qice_adv_mix = max(0.0, min(qice_up_mix(i,k), clw_all(i,k)))
        qliq_mass    = max(0.0, clw_all(i,k) - qice_adv_mix)

        qfreeze_mass = 0.0
        if(nliq_before_freeze_mix > 1.0e-30 .and. qliq_mass > qsmall_2m) then
           qfreeze_mass = qliq_mass * dn_freeze_mix / nliq_before_freeze_mix
        endif
        qfreeze_mass = max(0.0, min(qfreeze_mass, qliq_mass))

        qice_mass = qice_adv_mix + qfreeze_mass
        qice_mass = max(0.0, min(qice_mass, clw_all(i,k)))
        qliq_mass = max(0.0, clw_all(i,k) - qice_mass)

        if(tempc(i,k) <= thom_2m) then
           fhom = 1.0
        else
           fhom = 1.0 / (1.0 + exp((tempc(i,k) - thom_2m)/dthom_2m))
        endif
        fhom = max(0.0, min(1.0, fhom))

        dn_hom_mix = fhom * nliq_mix
        dn_hom_mix = max(0.0, min(dn_hom_mix, nliq_mix))

        ! Homogeneous freezing consumes the selected liquid mass, but the
        ! resulting ice crystal number is not one-to-one with the frozen
        ! droplet sink. Use the 12 um liquid / 40 um ice mass-ratio closure.
        dn_hom_ice_src_mix = real(hom_ice_number_ratio_2m) * dn_hom_mix
        dn_hom_ice_src_mix = max(0.0, min(dn_hom_ice_src_mix, dn_hom_mix))

        qhom_mass = fhom * qliq_mass
        qhom_mass = max(0.0, min(qhom_mass, qliq_mass))

        nliq_mix = max(0.0, nliq_mix - dn_hom_mix)
        nice_mix = max(0.0, nice_mix + dn_hom_ice_src_mix)

        qice_mass = max(0.0, min(qice_mass + qhom_mass, clw_all(i,k)))
        qliq_mass = max(0.0, clw_all(i,k) - qice_mass)

        !--------------------------------------------------------------------
        ! Deposition-only ice growth.
        !--------------------------------------------------------------------

        qdep_mass = 0.0

        if(tempc(i,k) < 273.15 .and. qice_mass > qsmall_2m .and. nice_mix > 1.0) then

           pres_pa_2m = max(100.0*po(i,k), 1.0)

           ! Use the same in-plume water saturation value used by the closure.
           qvl_sat_2m = max(qrch, 0.0)

           ! Ice saturation vapor pressure over ice [Pa], Murphy-Koop style.
           esi_pa_2m = exp(9.550426 - 5723.265/tempc(i,k)                  &
                         + 3.53068*log(tempc(i,k))                         &
                         - 0.00728332*tempc(i,k))

           esi_pa_2m = max(0.0, min(esi_pa_2m, 0.99*pres_pa_2m))

           qvi_sat_2m = eps_vap_2m*esi_pa_2m / max(pres_pa_2m - esi_pa_2m, 1.0)
           qvi_sat_2m = max(0.0, min(qvi_sat_2m, qvl_sat_2m))

           t_dep_arr(1)     = real(tempc(i,k), r8)
           qv_dep_arr(1)    = real(qvl_sat_2m, r8)
           qi_dep_arr(1)    = real(qice_mass, r8)
           ni_dep_arr(1)    = real(nice_mix, r8)
           icldm_dep_arr(1) = 1._r8
           rho_dep_arr(1)   = real(max(rho(i,k), 1.0e-12), r8)

           ! vapor diffusivity.
           ! p is pressure in Pa.
           dv_dep_arr(1)    = 8.794e-5_r8 * t_dep_arr(1)**1.81_r8          &
                            / real(pres_pa_2m, r8)

           qvl_dep_arr(1)   = real(qvl_sat_2m, r8)
           qvi_dep_arr(1)   = real(qvi_sat_2m, r8)

           call ice_deposition_sublimation(                                &
                t_dep_arr, qv_dep_arr, qi_dep_arr, ni_dep_arr,              &
                icldm_dep_arr, rho_dep_arr, dv_dep_arr,                    &
                qvl_dep_arr, qvi_dep_arr,                                  &
                berg_dep_arr, vap_dep_arr, ice_sublim_arr, mgncol_1 )

           ! Deposition only. Do not allow sublimation in this plume step.
           ice_sublim_arr(1) = 0._r8
           vap_dep_arr(1)    = max(vap_dep_arr(1), 0._r8)

           qdep_mass = real(vap_dep_arr(1)) * dt_plume
           qdep_mass = max(0.0, qdep_mass)

           ! Treat deposition growth as Bergeron transfer from the remaining
           ! cloud-liquid reservoir to cloud ice, conserving total condensate.
           qliq_avail = max(0.0, clw_all(i,k) - qice_mass)
           qdep_mass  = min(qdep_mass, qliq_avail)

           qice_mass = max(0.0, min(qice_mass + qdep_mass, clw_all(i,k)))
           qliq_mass = max(0.0, clw_all(i,k) - qice_mass)

        endif

        !--------------------------------------------------------------------
        ! Phase-aware number safeguard before autoconversion.
        !
        ! If retained liquid mass exists, liquid droplet number must exist.
        ! If retained ice mass exists, ice crystal number must exist.
        !--------------------------------------------------------------------

        if(clw_all(i,k) > qsmall_2m) then
           ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / max(rho(i,k), 1.0e-12)

           if(qliq_mass > qsmall_2m .and. nliq_mix < ncond_floor_mix) then
              nliq_mix = ncond_floor_mix
           endif

            ncond_floor_mix = ncond_floor_cm3_ice * 1.0e6 / max(rho(i,k), 1.0e-12)
            
           if(qice_mass > qsmall_2m .and. nice_mix < ncond_floor_mix) then
              nice_mix = ncond_floor_mix
           endif
        endif

        !--------------------------------------------------------------------
        ! Final number and mass fractions before autoconversion.
        !--------------------------------------------------------------------

        if(clw_all(i,k) > qsmall_2m) then
           pice_mass = qice_mass / clw_all(i,k)
        else
           pice_mass = 0.0
        endif

        pice_mass = max(0.0, min(1.0, pice_mass))
        pliq_mass = 1.0 - pice_mass

        qice_up_mix(i,k) = qice_mass

        nliq_up_mix(i,k) = nliq_mix
        nice_up_mix(i,k) = nice_mix
        nact_up_mix(i,k) = max(0.0, nliq_mix + nice_mix)

        nact_up_m3_out(i,k) = nact_up_mix(i,k) * max(rho(i,k), 1.0e-12)

        pice_2m(i,k) = pice_mass
        qice_up_2m(i,k) = max(0.0, min(qice_up_mix(i,k), clw_all(i,k)))
        qliq_up_2m(i,k) = max(0.0, clw_all(i,k) - qice_up_2m(i,k))
        nliq_up_2m(i,k) = max(0.0, nliq_up_mix(i,k))
        nice_up_2m(i,k) = max(0.0, nice_up_mix(i,k))

        !--------------------------------------------------------------------
        ! Initialize autoconversion sinks.
        !--------------------------------------------------------------------

        dq_auto_liq          = 0.0
        dq_auto_ice          = 0.0
        dn_auto_liq_mix      = 0.0
        dn_rain_auto_liq_mix = 0.0
        dn_auto_ice_mix      = 0.0
        dn_snow_auto_ice_mix = 0.0
        nliq_for_auto_mix    = 0.0

        !--------------------------------------------------------------------
        ! 1. Liquid autoconversion: cloud liquid -> warm precipitation.
        !
        !   qcic = kg kg-1
        !   ncic = # kg-1
        !   rho  = kg m-3
        !--------------------------------------------------------------------

        if(clw_all(i,k) > qsmall_2m) then

           qc_liq_auto       = pliq_mass * clw_all(i,k)
           nliq_for_auto_mix = nliq_mix

           if(qc_liq_auto > qsmall_2m .and. nliq_for_auto_mix > 0.0) then

              qcic_mg = real(qc_liq_auto, r8)
              rho_mg  = real(max(rho(i,k), 1.0e-12), r8)

              ! nsmall_2m is # m-3. Convert inline to # kg-1.
              ncic_mg = real(max(nsmall_2m/max(rho(i,k), 1.0e-12), &
                                 nliq_for_auto_mix), r8)

              call size_dist_param_liq( &
                   liq_props_mg, qcic_mg, ncic_mg, rho_mg, &
                   pgam_mg, lamc_mg )

              if(pgam_mg > 0._r8 .and. qcic_mg > 0._r8 .and. ncic_mg > 0._r8) then

                 pgam_arr(1) = pgam_mg
                 qcic_arr(1) = qcic_mg
                 ncic_arr(1) = ncic_mg
                 rho_arr(1)  = rho_mg

                 relvar_mg     = 1._r8 / (pgam_mg + 1._r8)
                 relvar_arr(1) = relvar_mg

                 call liu2006_liq_autoconversion( &
                      pgam_arr, microp_uniform_liu, qcic_arr, ncic_arr, &
                      rho_arr, relvar_arr, prc_arr, nprc_arr, nprc1_arr, &
                      mgncol_1 )

                 prc_liu_kgkg_s = real(prc_arr(1))
                 nprc1_liu      = real(nprc1_arr(1))

              else
                 prc_liu_kgkg_s = 0.0
                 nprc1_liu      = 0.0

              endif

              dq_auto_liq = prc_liu_kgkg_s * dt_plume
              dq_auto_liq = max(0.0, min(qc_liq_auto, dq_auto_liq))

              dn_auto_liq_mix = nprc1_liu * dt_plume
              dn_auto_liq_mix = max(0.0, min(nliq_for_auto_mix, dn_auto_liq_mix))

              ! Keep liquid droplet-number loss consistent with liquid mass loss.
              ! This is a cloud-droplet sink, not a rain-number source.
              if(qc_liq_auto > qsmall_2m .and. nliq_for_auto_mix > 0.0) then
                 dn_auto_liq_mix = min(dn_auto_liq_mix,                    &
                                        nliq_for_auto_mix *                 &
                                        dq_auto_liq / qc_liq_auto)
              endif

              dq_auto_liq     = real(auto_eff_liq_2m) * dq_auto_liq
              dn_auto_liq_mix = real(auto_eff_liq_2m) * dn_auto_liq_mix

              if(use_mp_auto_number_2m) then
                 call cup_2m_precip_number_from_auto_mass(                 &
                      dq_auto_liq, dn_auto_liq_mix, max(rho(i,k),1.0e-12), &
                      rho_liq_mg, n0r_mp_auto_2m, lam_bnd_rain_2m,          &
                      dn_rain_auto_liq_mix )
              else
                 dn_rain_auto_liq_mix = max(0.0, real(auto_eff_liq_2m) * &
                                                  real(nprc_arr(1)) * dt_plume)
                 if(dn_auto_liq_mix > 0.0) then
                    dn_rain_auto_liq_mix = min(dn_rain_auto_liq_mix, dn_auto_liq_mix)
                 endif
              endif

           endif
           
           

        endif


        !--------------------------------------------------------------------
        ! State after liquid autoconversion.
        !--------------------------------------------------------------------

        qcloud_after_liq_auto = max(0.0, clw_all(i,k) - dq_auto_liq)

        if(remove_inp_aerosol_with_liq_auto_2m) then
           ! Optional drop-mediated INP aerosol sink tied to droplet autoconversion.
           ! If autoconversion removes a fraction of droplets, remove the same
           ! fraction from finite aerosol reservoirs assumed to be inside droplets.
           ! Background INP is not depleted.
           frac_liq_auto = 0.0
           if(nliq_for_auto_mix > 1.0e-30) then
              frac_liq_auto = max(0.0, min(1.0, dn_auto_liq_mix/nliq_for_auto_mix))
           endif

           if(frac_liq_auto > 0.0) then
              ndust_mix    = max(0.0, ndust_mix   *(1.0 - frac_liq_auto))
              nsoot_mix    = max(0.0, nsoot_mix   *(1.0 - frac_liq_auto))
              nseasalt_mix = max(0.0, nseasalt_mix*(1.0 - frac_liq_auto))
           endif

           ndust_up_mix(i,k)    = ndust_mix
           nsoot_up_mix(i,k)    = nsoot_mix
           nseasalt_up_mix(i,k) = nseasalt_mix
        endif

        nliq_mix = max(0.0, nliq_mix - dn_auto_liq_mix)
       
        
        nice_mix_before_ice_auto = nice_mix

        !--------------------------------------------------------------------
        ! 2. Ice autoconversion: cloud ice -> cold precipitation/snow.
        !--------------------------------------------------------------------

        qc_ice_auto = min(qice_mass, qcloud_after_liq_auto)

        if(qc_ice_auto > qsmall_2m .and. nice_mix_before_ice_auto > 0.0) then

           tice_mg = real(tempc(i,k), r8)
           qiic_mg = real(qc_ice_auto, r8)
           niic_mg = real(nice_mix_before_ice_auto, r8)
           rho_mg  = real(max(rho(i,k), 1.0e-12), r8)

           call size_dist_param_ice( &
                ice_props_mg, qiic_mg, niic_mg, lami_mg )

           if(lami_mg > 0._r8 .and. qiic_mg > 0._r8 .and. niic_mg > 0._r8) then

              tice_arr(1) = tice_mg
              qiic_arr(1) = qiic_mg
              lami_arr(1) = lami_mg
              niic_arr(1) = niic_mg

              call ice_autoconversion( &
                   tice_arr, qiic_arr, lami_arr, niic_arr, &
                   dcs_mg, prci_arr, nprci_arr, mgncol_1, ac_time_mg )

              dq_auto_ice = real(prci_arr(1)) * dt_plume
              dq_auto_ice = max(0.0, min(qc_ice_auto, dq_auto_ice))

              dn_auto_ice_mix = real(nprci_arr(1)) * dt_plume
              dn_auto_ice_mix = max(0.0, min(nice_mix_before_ice_auto, dn_auto_ice_mix))

              if(use_mp_auto_number_2m) then
                 call cup_2m_precip_number_from_auto_mass(                 &
                      dq_auto_ice, dn_auto_ice_mix, max(rho(i,k),1.0e-12), &
                      rho_snow_2m, n0s_mp_auto_2m, lam_bnd_snow_2m,         &
                      dn_snow_auto_ice_mix )
              else
                 dn_snow_auto_ice_mix = dn_auto_ice_mix
              endif

           endif

        endif

        !--------------------------------------------------------------------
        ! Final retained cloud condensate and number after warm + cold precip.
        !--------------------------------------------------------------------

        qcloud_after_all_auto = max(0.0, qcloud_after_liq_auto - dq_auto_ice)

        nliq_mix = max(0.0, nliq_mix)
        nice_mix = max(0.0, nice_mix_before_ice_auto - dn_auto_ice_mix)

        qice_up_mix(i,k) = max(0.0, qice_mass - dq_auto_ice)
        qice_up_mix(i,k) = min(qice_up_mix(i,k), qcloud_after_all_auto)

        !--------------------------------------------------------------------
        ! Final post-autoconversion safeguard.
        !
        ! Autoconversion can remove number more aggressively than mass in
        ! edge cases.  Do not detrain retained condensate with zero number.
        !--------------------------------------------------------------------

        if(qcloud_after_all_auto > qsmall_2m) then
           ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / max(rho(i,k), 1.0e-12)
           qliq_mass = max(0.0, qcloud_after_all_auto - qice_mass)

           if(qliq_mass > qsmall_2m .and. nliq_mix < ncond_floor_mix) then
              nliq_mix = ncond_floor_mix
           endif

           ncond_floor_mix = ncond_floor_cm3_ice * 1.0e6 / max(rho(i,k), 1.0e-12)
           qice_mass = max(0.0, min(qice_up_mix(i,k), qcloud_after_all_auto))
         
           if(qice_mass > qsmall_2m .and. nice_mix < ncond_floor_mix) then
              nice_mix = ncond_floor_mix
           endif
        endif

        nliq_up_mix(i,k) = nliq_mix
        nice_up_mix(i,k) = nice_mix
        nact_up_mix(i,k) = max(0.0, nliq_mix + nice_mix)

        nact_up_m3_out(i,k) = nact_up_mix(i,k) * max(rho(i,k), 1.0e-12)

        !--------------------------------------------------------------------
        ! Store provisional post-autoconversion cloud state.
        !
        ! The upward transport uses this rainout-adjusted state.  The
        ! precipitation-dependent collection processes are applied afterward in
        ! a separate top-down pass because local rain/snow depend on production
        ! from levels aloft.
        !--------------------------------------------------------------------

        qice_col(k) = max(0.0, min(qice_up_mix(i,k), qcloud_after_all_auto))
        qliq_col(k) = max(0.0, qcloud_after_all_auto - qice_col(k))

        nliq_col(k) = max(0.0, nliq_mix)
        nice_col(k) = max(0.0, nice_mix)

        dq_auto_liq_col(k) = max(0.0, dq_auto_liq)
        dq_auto_ice_col(k) = max(0.0, dq_auto_ice)
        dn_auto_liq_col(k)      = max(0.0, dn_auto_liq_mix)
        dn_auto_ice_col(k)      = max(0.0, dn_auto_ice_mix)
        dn_rain_auto_liq_col(k) = max(0.0, dn_rain_auto_liq_mix)
        dn_snow_auto_ice_col(k) = max(0.0, dn_snow_auto_ice_mix)

        dt_col(k)   = max(0.0, dt_plume)
        dz_col(k)   = max(0.0, dz)
        qrch_col(k) = qrch
        w_col(k)    = w_cup

        nact_lev_mix_col(k)    = nact_lev_mix
        nact_target_mix_col(k) = nact_target_mix
        ndust_mix_col(k)       = ndust_mix
        ninp_imm_mix_col(k)    = ninp_imm_mix
        if(aer_reservoir_initialized .and. nmodes > 0) then
           nwfa_cm3_col(k) = nwfa_now_m3 * 1.0e-6
        else
           nwfa_cm3_col(k) = 0.0
        endif

        ! Provisional state needed by the next upward level.
        qrc(i,k)         = qcloud_after_all_auto
        qice_up_mix(i,k) = qice_col(k)
        qc(i,k)          = qrc(i,k) + min(qc(i,k), qrch)

        nliq_up_mix(i,k) = nliq_col(k)
        nice_up_mix(i,k) = nice_col(k)
        nact_up_mix(i,k) = max(0.0, nliq_col(k) + nice_col(k))
        nact_up_m3_out(i,k) = nact_up_mix(i,k) * max(rho(i,k), 1.0e-12)

        if(qcloud_after_all_auto > qsmall_2m) then
           pice_2m(i,k) = max(0.0, min(1.0, qice_col(k)/qcloud_after_all_auto))
        else
           pice_2m(i,k) = 0.0
        endif

     enddo

         !-----------------------------------------------------------------------
         ! Predictor-corrector precipitation treatment.
         !
         ! The upward loop above is  the predictor.  It supplies provisional
         ! liquid/ice autoconversion sources, which are used here to diagnose fixed
         ! rain/snow profiles.  The final pass then replays the plume upward with
         ! those fixed precipitation profiles.  Thus rain/snow accretion competes
         ! with freezing before the final ice-autoconversion calculation.
         !-----------------------------------------------------------------------

     ! Full 2M precipitation processes are always evaluated here.
     ! Individual autoconversion/accretion rates are controlled by
     ! GEOSmoist_Process_Library tunings, not by legacy c0 constants.

        !--------------------------------------------------------------------
        ! 1. Diagnose fixed rain/snow profiles from provisional autoconversion
        !    and provisional accretion.
        !
        ! This downward predictor lets falling precipitation grow by
        ! collecting the provisional cloud field.  The resulting qrain/qsnow
        ! profiles are then held fixed during the final upward corrector.
        !
        ! For the current layer, qrain_diag_col/qsnow_diag_col represent the
        ! precipitating water available to collect cloud at that level before
        ! same-layer accretion growth.  The accreted mass is added to the flux
        ! leaving the bottom of the layer and therefore affects lower levels.
        !--------------------------------------------------------------------

        rain_flux_m = 0.0
        snow_flux_m = 0.0
        rain_flux_n = 0.0
        snow_flux_n = 0.0

        do k = ktop(i) + 1, kbase_2m, -1

           rho_safe = max(rho(i,k), 1.0e-12)
           dz_eff   = max(dz_col(k), 0.0)
           dt_eff   = max(dt_col(k), 1.0e-6)

           ! Add provisional autoconversion production in this layer.
           rain_flux_m = rain_flux_m + rho_safe * dq_auto_liq_col(k) / dt_eff * dz_eff
           snow_flux_m = snow_flux_m + rho_safe * dq_auto_ice_col(k) / dt_eff * dz_eff

           ! Use precip-number sources, not cloud-number sinks, for falling flux PSD.
           rain_flux_n = rain_flux_n + rho_safe * dn_rain_auto_liq_col(k) / dt_eff * dz_eff
           snow_flux_n = snow_flux_n + rho_safe * dn_snow_auto_ice_col(k) / dt_eff * dz_eff

           ! Diagnose rain/snow seen by the provisional cloud at this level.
           call cup_2m_mp_from_flux(                                      &
                real(rain_flux_m, r8), real(rain_flux_n, r8),             &
                real(rho_safe, r8), rho_liq_mg,                           &
                arain_2m, brain_2m, lam_bnd_rain_2m,                      &
                qrain_diag, nrain_diag, lamr_2m, n0r_2m, umr_2m, unr_2m )

           ! Optional number-only rain self-collection.  Mass flux is conserved;
           ! rain number flux is reduced and the MP rain PSD is rebuilt.
           rain_flux_n_before_self = rain_flux_n
           dn_rain_self_mix = 0.0

           call cup_2m_apply_mg_rain_self_collection(                     &
                qrain_diag, nrain_diag, rho_safe, dt_eff,                 &
                rain_flux_n, dn_rain_self_mix )

           if(rain_flux_n < rain_flux_n_before_self) then
              call cup_2m_mp_from_flux(                                   &
                   real(rain_flux_m, r8), real(rain_flux_n, r8),          &
                   real(rho_safe, r8), rho_liq_mg,                        &
                   arain_2m, brain_2m, lam_bnd_rain_2m,                   &
                   qrain_diag, nrain_diag, lamr_2m, n0r_2m, umr_2m, unr_2m )
           endif

           call cup_2m_mp_from_flux(                                      &
                real(snow_flux_m, r8), real(snow_flux_n, r8),             &
                real(rho_safe, r8), rho_snow_2m,                          &
                asnow_2m, bsnow_2m, lam_bnd_snow_2m,                      &
                qsnow_diag, nsnow_diag, lams_2m, n0s_2m, ums_2m, uns_2m )

           ! Optional number-only snow self-aggregation.  Mass flux is conserved;
           ! snow number flux is reduced and the MP snow PSD is rebuilt.
           snow_flux_n_before_self = snow_flux_n
           dn_snow_self_mix = 0.0

           call cup_2m_apply_mg_snow_self_aggregation(                    &
                tempc(i,k), qsnow_diag, nsnow_diag, rho_safe, dt_eff,      &
                snow_flux_n, dn_snow_self_mix )

           if(snow_flux_n < snow_flux_n_before_self) then
              call cup_2m_mp_from_flux(                                   &
                   real(snow_flux_m, r8), real(snow_flux_n, r8),          &
                   real(rho_safe, r8), rho_snow_2m,                       &
                   asnow_2m, bsnow_2m, lam_bnd_snow_2m,                   &
                   qsnow_diag, nsnow_diag, lams_2m, n0s_2m, ums_2m, uns_2m )
           endif

           qrain_diag_col(k) = qrain_diag
           qsnow_diag_col(k) = qsnow_diag
           nrain_diag_col(k) = nrain_diag
           nsnow_diag_col(k) = nsnow_diag

           lamr_col(k) = lamr_2m
           n0r_col(k)  = n0r_2m
           lams_col(k) = lams_2m
           n0s_col(k)  = n0s_2m

           umr_col(k) = umr_2m
           ums_col(k) = ums_2m
           unr_col(k) = unr_2m
           uns_col(k) = uns_2m

           ! Provisional accretion predictor.  This does not finalize cloud
           ! tendencies; it only grows the downward rain/snow flux profiles.
           qliq_mass = max(0.0, qliq_col(k))
           qice_mass = max(0.0, qice_col(k))
           nliq_mix  = max(0.0, nliq_col(k))
           nice_mix  = max(0.0, nice_col(k))

           call cup_2m_apply_fixed_precip_accretion(                       &
                tempc(i,k), rho_safe, dt_eff,                              &
                qrain_diag, nrain_diag,                                    &
                qsnow_diag, nsnow_diag,                                    &
                lamr_2m, n0r_2m, lams_2m, n0s_2m,                          &
                umr_2m, ums_2m, unr_2m, uns_2m,                            &
                liq_props_mg,                                              &
                qliq_mass, nliq_mix, qice_mass, nice_mix,                  &
                dq_cw_rain, dq_cw_snow_to_snow, dq_ci_snow,                &
                dq_sip, dq_rs )

           if(.not. GF2M_USE_CORRECTOR) then
              ! No-corrector mode: this downward pass is the only precipitation-
              ! dependent cloud update.  Store the post-accretion cloud state;
              ! do not replay activation/freezing/autoconversion upward.
              dq_acc_liq_rain_col(k) = max(0.0, dq_cw_rain)
              dq_acc_liq_snow_col(k) = max(0.0, dq_cw_snow_to_snow)
              dq_acc_ice_snow_col(k) = max(0.0, dq_ci_snow)
              dq_sip_col(k)          = max(0.0, dq_sip)
              dq_rain_to_snow_col(k) = max(0.0, dq_rs)

              qliq_col(k) = max(0.0, qliq_mass)
              qice_col(k) = max(0.0, qice_mass)
              nliq_col(k) = max(0.0, nliq_mix)
              nice_col(k) = max(0.0, nice_mix)
           endif

           ! Estimate rain number converted to snow using the diagnosed rain
           ! mean mass.  The helper caps dq_rs against available rain mass and
           ! number internally, but currently returns only the mass transfer.
           if(qrain_diag > 1.0e-30 .and. nrain_diag > 0.0) then
              dn_rs = min(nrain_diag, nrain_diag * dq_rs / qrain_diag)
           else
              dn_rs = 0.0
           endif
           dn_rs = max(0.0, dn_rs)

           ! Update fluxes leaving the bottom of this layer.  Cloud-water
           ! accretion by rain grows rain mass but not rain number.  Snow
           ! riming and cloud-ice collection grow snow mass but not snow
           ! number.  Rain-snow collection transfers rain mass/number to snow
           ! mass without creating new snow particles.
           rain_flux_m = rain_flux_m                                      &
                       + rho_safe * dq_cw_rain / dt_eff * dz_eff          &
                       - rho_safe * dq_rs      / dt_eff * dz_eff

           snow_flux_m = snow_flux_m                                      &
                       + rho_safe * (dq_cw_snow_to_snow + dq_ci_snow)     &
                                   / dt_eff * dz_eff                      &
                       + rho_safe * dq_rs / dt_eff * dz_eff

           rain_flux_n = rain_flux_n                                      &
                       - rho_safe * dn_rs / dt_eff * dz_eff

           rain_flux_m = max(0.0, rain_flux_m)
           snow_flux_m = max(0.0, snow_flux_m)
           rain_flux_n = max(0.0, rain_flux_n)
           snow_flux_n = max(0.0, snow_flux_n)

        enddo

        if(.not. GF2M_USE_CORRECTOR) then

           !-----------------------------------------------------------------
           ! No-corrector finalization.
           !
           ! The upward pass supplied activation/freezing/deposition plus
           ! cloud liquid/ice autoconversion.  The downward pass above supplied
           ! precipitation-dependent accretion/collection.  Finalize the plume
           ! state directly from those two passes; do not repeat the upward
           ! transport or mutate the aerosol activation reservoir a second time.
           !-----------------------------------------------------------------

           do k = kbase_2m, ktop(i) + 1

              dz_eff   = max(dz_col(k), 0.0)
              rho_safe = max(rho(i,k), 1.0e-12)
              qrch     = qrch_col(k)

              qliq_mass = max(0.0, qliq_col(k))
              qice_mass = max(0.0, qice_col(k))
              nliq_mix  = max(0.0, nliq_col(k))
              nice_mix  = max(0.0, nice_col(k))

              qcloud_after_all_auto = max(0.0, qliq_mass + qice_mass)

              qrc(i,k)         = qcloud_after_all_auto
              qice_up_mix(i,k) = max(0.0, min(qice_mass, qcloud_after_all_auto))

              nliq_up_mix(i,k) = nliq_mix
              nice_up_mix(i,k) = nice_mix
              nact_up_mix(i,k) = max(0.0, nliq_mix + nice_mix)
              nact_up_m3_out(i,k) = nact_up_mix(i,k) * rho_safe

              if(qcloud_after_all_auto > qsmall_2m) then
                 pice_2m(i,k) = max(0.0, min(1.0, qice_up_mix(i,k)/qcloud_after_all_auto))
              else
                 pice_2m(i,k) = 0.0
              endif

              qice_up_2m(i,k) = max(0.0, min(qice_up_mix(i,k), qcloud_after_all_auto))
              qliq_up_2m(i,k) = max(0.0, qcloud_after_all_auto - qice_up_2m(i,k))
              nliq_up_2m(i,k) = max(0.0, nliq_up_mix(i,k))
              nice_up_2m(i,k) = max(0.0, nice_up_mix(i,k))

              ! Return the full-2M phase-resolved retained condensate as the
              ! unique cloud condensate state used by the parent detrainment path.
              qrc(i,k) = qliq_up_2m(i,k) + qice_up_2m(i,k)
              qice_up_mix(i,k) = qice_up_2m(i,k)
              if(qrc(i,k) > qsmall_2m) then
                 pice_2m(i,k) = max(0.0, min(1.0, qice_up_2m(i,k)/qrc(i,k)))
              else
                 pice_2m(i,k) = 0.0
              endif
              qc(i,k) = qrc(i,k) + min(qc(i,k), qrch)

              dq_precip_level = dq_auto_liq_col(k) + dq_auto_ice_col(k)       &
                              + dq_acc_liq_rain_col(k)                       &
                              + dq_acc_liq_snow_col(k)                       &
                              + dq_acc_ice_snow_col(k)

              pw(i,k) = max(0.0, dq_precip_level) * zu(i,k)
              qc(i,k) = qrc(i,k) + min(qc(i,k), qrch)

              pwav(i) = pwav(i) + pw(i,k)
              psum(i) = psum(i) + clw_all(i,k)*zu(i,k)*dz_eff

           enddo

           if(pwav(i) < 0.0) then
              ierr(i)  = 66
              ierrc(i) = "pwav negative"
           endif

           cycle

        endif

        !--------------------------------------------------------------------
        ! 2. Final combined plume pass with fixed rain/snow.
        !
        ! This pass repeats the upward transport and final microphysics.  The
        ! fixed precipitation profiles diagnosed above are collection targets
        ! only; they are not updated by this final pass.
        !--------------------------------------------------------------------

        do k = kbase_2m, ktop(i) + 1

           dz_eff   = max(dz_col(k), 0.0)
           dt_eff   = max(dt_col(k), 1.0e-6)
           rho_safe = max(rho(i,k), 1.0e-12)
           qrch     = qrch_col(k)

           denom = zu(i,k-1) - 0.5*up_massdetr(i,k-1) + up_massentr(i,k-1)

           if(denom > 0.0) then

              qc(i,k) = ( qc(i,k-1)*zu(i,k-1)                              &
                        - 0.5*up_massdetr(i,k-1)*qc(i,k-1)                  &
                        + up_massentr(i,k-1)*q(i,k-1) ) / denom

              if(k == kbase_2m) then
                 qc(i,k) = qc(i,k) + zqexec(i)*up_massentr(i,k-1)/denom
              endif

              qrc(i,k) = ( qrc(i,k-1)*zu(i,k-1)                             &
                         - 0.5*up_massdetr(i,k-1)*qrc(i,k-1) ) / denom

              qice_up_mix(i,k) = ( qice_up_mix(i,k-1)*zu(i,k-1)              &
                                 - 0.5*up_massdetr(i,k-1)*qice_up_mix(i,k-1) &
                                 ) / denom

           else

              qc(i,k)          = qc(i,k-1)
              qrc(i,k)         = qrc(i,k-1)
              qice_up_mix(i,k) = qice_up_mix(i,k-1)

           endif

           ! Reconstruct transported number state for the final pass.
           if(k == kbase_2m) then

              nliq_up_mix(i,k)  = 0.0
              nice_up_mix(i,k)  = 0.0
              qice_up_mix(i,k)  = 0.0

              khost = flip(kbase_2m)
              if(nmodes > 0) then
                 call aerosol_ice_species_from_aerp(                       &
                      AerP, nmodes, i, jcol, khost, rho(i,k),              &
                      ndust_cb_mix, ddust_cb_m, sigdust_cb,                &
                      nsoot_cb_mix, dsoot_cb_m, sigsoot_cb,                &
                      nseasalt_cb_mix, dseasalt_cb_m, sigseasalt_cb )
              else
                 ndust_cb_mix    = 0.0
                 ddust_cb_m      = 0.0
                 sigdust_cb      = 0.0
                 nsoot_cb_mix    = 0.0
                 dsoot_cb_m      = 0.0
                 sigsoot_cb      = 0.0
                 nseasalt_cb_mix = 0.0
                 dseasalt_cb_m   = 0.0
                 sigseasalt_cb   = 0.0
              endif

              ndust_up_mix(i,k)    = max(0.0, ndust_cb_mix)
              ddust_up_m(i,k)      = max(0.0, ddust_cb_m)
              sigdust_up(i,k)      = max(0.0, sigdust_cb)
              nsoot_up_mix(i,k)    = max(0.0, nsoot_cb_mix)
              dsoot_up_m(i,k)      = max(0.0, dsoot_cb_m)
              sigsoot_up(i,k)      = max(0.0, sigsoot_cb)
              nseasalt_up_mix(i,k) = max(0.0, nseasalt_cb_mix)
              dseasalt_up_m(i,k)   = max(0.0, dseasalt_cb_m)
              sigseasalt_up(i,k)   = max(0.0, sigseasalt_cb)

           else

              if(denom > 0.0) then

                 nliq_up_mix(i,k) = ( nliq_up_mix(i,k-1)*zu(i,k-1)              &
                                    - 0.5*up_massdetr(i,k-1)*nliq_up_mix(i,k-1) &
                                    ) / denom

                 nice_up_mix(i,k) = ( nice_up_mix(i,k-1)*zu(i,k-1)              &
                                    - 0.5*up_massdetr(i,k-1)*nice_up_mix(i,k-1) &
                                    ) / denom

                 ndust_up_mix(i,k) = ( ndust_up_mix(i,k-1)*zu(i,k-1)              &
                                     - 0.5*up_massdetr(i,k-1)*ndust_up_mix(i,k-1) &
                                     ) / denom

                 nsoot_up_mix(i,k) = ( nsoot_up_mix(i,k-1)*zu(i,k-1)              &
                                     - 0.5*up_massdetr(i,k-1)*nsoot_up_mix(i,k-1) &
                                     ) / denom

                 nseasalt_up_mix(i,k) = ( nseasalt_up_mix(i,k-1)*zu(i,k-1)              &
                                        - 0.5*up_massdetr(i,k-1)*nseasalt_up_mix(i,k-1) &
                                        ) / denom

              else

                 nliq_up_mix(i,k)  = nliq_up_mix(i,k-1)
                 nice_up_mix(i,k)  = nice_up_mix(i,k-1)
                 ndust_up_mix(i,k)    = ndust_up_mix(i,k-1)
                 nsoot_up_mix(i,k)    = nsoot_up_mix(i,k-1)
                 nseasalt_up_mix(i,k) = nseasalt_up_mix(i,k-1)

              endif

              nliq_up_mix(i,k)  = max(0.0, nliq_up_mix(i,k))
              nice_up_mix(i,k)  = max(0.0, nice_up_mix(i,k))
              ndust_up_mix(i,k)    = max(0.0, ndust_up_mix(i,k))
              nsoot_up_mix(i,k)    = max(0.0, nsoot_up_mix(i,k))
              nseasalt_up_mix(i,k) = max(0.0, nseasalt_up_mix(i,k))

              ddust_up_m(i,k)      = ddust_up_m(i,k-1)
              sigdust_up(i,k)      = sigdust_up(i,k-1)
              dsoot_up_m(i,k)      = dsoot_up_m(i,k-1)
              sigsoot_up(i,k)      = sigsoot_up(i,k-1)
              dseasalt_up_m(i,k)   = dseasalt_up_m(i,k-1)
              sigseasalt_up(i,k)   = sigseasalt_up(i,k-1)

           endif

           ! Reuse the predictor activation result.  This avoids a second
           ! mutation of the local aerosol reservoir.  With target-increment
           ! activation, recompute the positive increment against the corrected
           ! transported droplet number.  With legacy additive activation, reuse
           ! the predictor-pass source directly.
           nact_target_mix = max(0.0, nact_target_mix_col(k))
           if(use_activation_target_increment_2m) then
              nact_lev_mix = max(0.0, nact_target_mix - nliq_up_mix(i,k))
           else
              nact_lev_mix = max(0.0, nact_lev_mix_col(k))
           endif
           nliq_up_mix(i,k) = max(0.0, nliq_up_mix(i,k) + nact_lev_mix)
           nact_lev_mix_col(k) = nact_lev_mix

           clw_all(i,k) = max(0.0, qc(i,k) - qrch)
           qrc(i,k)     = min(clw_all(i,k), qrc(i,k))

           qice_adv_mix = max(0.0, min(qice_up_mix(i,k), qrc(i,k), clw_all(i,k)))
           qice_mass    = qice_adv_mix
           qliq_mass    = max(0.0, clw_all(i,k) - qice_mass)

           nliq_mix       = max(0.0, nliq_up_mix(i,k))
           nice_mix       = max(0.0, nice_up_mix(i,k))
           ndust_mix      = max(0.0, ndust_up_mix(i,k))
           ddust_mean_m   = max(0.0, ddust_up_m(i,k))
           sigdust_mean   = max(0.0, sigdust_up(i,k))
           nsoot_mix      = max(0.0, nsoot_up_mix(i,k))
           dsoot_mean_m   = max(0.0, dsoot_up_m(i,k))
           sigsoot_mean   = max(0.0, sigsoot_up(i,k))
           nseasalt_mix   = max(0.0, nseasalt_up_mix(i,k))
           dseasalt_mean_m= max(0.0, dseasalt_up_m(i,k))
           sigseasalt_mean= max(0.0, sigseasalt_up(i,k))

           ! Initial condensate-number safeguard before collection/freezing.
           if(clw_all(i,k) > qsmall_2m) then
              ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / rho_safe
              if(qliq_mass > qsmall_2m .and. nliq_mix < ncond_floor_mix) then
                 nliq_mix = ncond_floor_mix
              endif

              ncond_floor_mix = ncond_floor_cm3_ice * 1.0e6 / rho_safe
              if(qice_mass > qsmall_2m .and. nice_mix < ncond_floor_mix) then
                 nice_mix = ncond_floor_mix
              endif
           endif

           !---------------------------------------------------------------
           ! Collection/accretion against fixed rain/snow, before freezing.
           !---------------------------------------------------------------

           call cup_2m_apply_fixed_precip_accretion(                       &
                tempc(i,k), rho_safe, dt_eff,                              &
                qrain_diag_col(k), nrain_diag_col(k),                      &
                qsnow_diag_col(k), nsnow_diag_col(k),                      &
                lamr_col(k), n0r_col(k), lams_col(k), n0s_col(k),          &
                umr_col(k), ums_col(k), unr_col(k), uns_col(k),            &
                liq_props_mg,                                              &
                qliq_mass, nliq_mix, qice_mass, nice_mix,                  &
                dq_cw_rain, dq_cw_snow_to_snow, dq_ci_snow,                &
                dq_sip, dq_rs )

           dq_acc_liq_rain_col(k) = dq_cw_rain
           dq_acc_liq_snow_col(k) = dq_cw_snow_to_snow
           dq_acc_ice_snow_col(k) = dq_ci_snow
           dq_sip_col(k)          = dq_sip
           dq_rain_to_snow_col(k) = dq_rs

           !---------------------------------------------------------------
           ! Immersion freezing after accretion.
           !---------------------------------------------------------------

           nliq_before_freeze_mix = nliq_mix

           seasalt_imm_scale = max(0.0, min(1.0, seasalt_imm_scale_col(k)))
           nseasalt_for_inp_mix = max(0.0, nseasalt_mix * seasalt_imm_scale)

           ! Do not let the artificial condensate-number floor trigger freezing.
           ! INP freezing is evaluated only against droplets above the numerical
           ! liquid-number floor, and only if liquid mass is actually present.
           ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / rho_safe
           nliq_freeze_avail_mix = max(0.0, nliq_mix - ncond_floor_mix)

           if(qliq_mass > qsmall_2m .and. nliq_freeze_avail_mix > 1.0) then
              call cup_2m_inp_freezing_sources(                               &
                   tempc(i,k), po(i,k), rho_safe, w_col(k), dt_eff,           &
                   qliq_mass, nliq_freeze_avail_mix,                          &
                   ndust_mix, ddust_mean_m, sigdust_mean,                     &
                   nsoot_mix, dsoot_mean_m, sigsoot_mean,                     &
                   nseasalt_for_inp_mix, dseasalt_mean_m, sigseasalt_mean,    &
                   ninp_dust_mix, ninp_soot_mix, ninp_seasalt_mix,            &
                   ninp_bkg_mix, ninp_imm_mix )
           else
              ninp_dust_mix    = 0.0
              ninp_soot_mix    = 0.0
              ninp_seasalt_mix = 0.0
              ninp_bkg_mix     = 0.0
              ninp_imm_mix     = 0.0
           endif

           dn_freeze_dust_mix    = ninp_dust_mix
           dn_freeze_soot_mix    = ninp_soot_mix
           dn_freeze_seasalt_mix = ninp_seasalt_mix
           dn_freeze_bkg_mix     = ninp_bkg_mix
           dn_freeze_mix         = ninp_imm_mix

           if(dn_freeze_mix > nliq_mix .and. dn_freeze_mix > 1.0e-30) then
              scale_sink = max(0.0, min(1.0, nliq_mix/dn_freeze_mix))
              dn_freeze_dust_mix    = dn_freeze_dust_mix    * scale_sink
              dn_freeze_soot_mix    = dn_freeze_soot_mix    * scale_sink
              dn_freeze_seasalt_mix = dn_freeze_seasalt_mix * scale_sink
              dn_freeze_bkg_mix     = dn_freeze_bkg_mix     * scale_sink
              dn_freeze_mix         = nliq_mix
           endif

           qfreeze_mass = 0.0
           if(nliq_before_freeze_mix > 1.0e-30 .and. qliq_mass > qsmall_2m) then
              qfreeze_mass = qliq_mass * dn_freeze_mix / nliq_before_freeze_mix
           endif
           qfreeze_mass = max(0.0, min(qfreeze_mass, qliq_mass))

           nliq_mix       = max(0.0, nliq_mix  - dn_freeze_mix)
           nice_mix       = max(0.0, nice_mix  + dn_freeze_mix)
           ndust_mix      = max(0.0, ndust_mix - dn_freeze_dust_mix)
           nsoot_mix      = max(0.0, nsoot_mix - dn_freeze_soot_mix)
           nseasalt_mix   = max(0.0, nseasalt_mix - dn_freeze_seasalt_mix)

           qliq_mass = max(0.0, qliq_mass - qfreeze_mass)
           qice_mass = max(0.0, qice_mass + qfreeze_mass)

           !---------------------------------------------------------------
           ! Homogeneous freezing of remaining droplets.
           !---------------------------------------------------------------

           if(tempc(i,k) <= thom_2m) then
              fhom = 1.0
           else
              fhom = 1.0 / (1.0 + exp((tempc(i,k) - thom_2m)/dthom_2m))
           endif
           fhom = max(0.0, min(1.0, fhom))

           dn_hom_mix = fhom * nliq_mix
           dn_hom_mix = max(0.0, min(dn_hom_mix, nliq_mix))

           ! Homogeneous freezing consumes the selected liquid mass, but the
           ! resulting ice crystal number is not one-to-one with the frozen
           ! droplet sink. Use the 12 um liquid / 40 um ice mass-ratio closure.
           dn_hom_ice_src_mix = real(hom_ice_number_ratio_2m) * dn_hom_mix
           dn_hom_ice_src_mix = max(0.0, min(dn_hom_ice_src_mix, dn_hom_mix))

           qhom_mass = fhom * qliq_mass
           qhom_mass = max(0.0, min(qhom_mass, qliq_mass))

           nliq_mix = max(0.0, nliq_mix - dn_hom_mix)
           nice_mix = max(0.0, nice_mix + dn_hom_ice_src_mix)

           qliq_mass = max(0.0, qliq_mass - qhom_mass)
           qice_mass = max(0.0, qice_mass + qhom_mass)

           !---------------------------------------------------------------
           ! Deposition-only ice growth, treated as liquid-to-ice transfer.
           !---------------------------------------------------------------

           qdep_mass = 0.0

           if(tempc(i,k) < 273.15 .and. qice_mass > qsmall_2m .and. nice_mix > 1.0) then

              pres_pa_2m = max(100.0*po(i,k), 1.0)
              qvl_sat_2m = max(qrch, 0.0)

              esi_pa_2m = exp(9.550426 - 5723.265/tempc(i,k)                 &
                            + 3.53068*log(tempc(i,k))                        &
                            - 0.00728332*tempc(i,k))

              esi_pa_2m = max(0.0, min(esi_pa_2m, 0.99*pres_pa_2m))

              qvi_sat_2m = eps_vap_2m*esi_pa_2m / max(pres_pa_2m - esi_pa_2m, 1.0)
              qvi_sat_2m = max(0.0, min(qvi_sat_2m, qvl_sat_2m))

              t_dep_arr(1)     = real(tempc(i,k), r8)
              qv_dep_arr(1)    = real(qvl_sat_2m, r8)
              qi_dep_arr(1)    = real(qice_mass, r8)
              ni_dep_arr(1)    = real(nice_mix, r8)
              icldm_dep_arr(1) = 1._r8
              rho_dep_arr(1)   = real(rho_safe, r8)

              dv_dep_arr(1)    = 8.794e-5_r8 * t_dep_arr(1)**1.81_r8          &
                               / real(pres_pa_2m, r8)

              qvl_dep_arr(1)   = real(qvl_sat_2m, r8)
              qvi_dep_arr(1)   = real(qvi_sat_2m, r8)

              call ice_deposition_sublimation(                               &
                   t_dep_arr, qv_dep_arr, qi_dep_arr, ni_dep_arr,             &
                   icldm_dep_arr, rho_dep_arr, dv_dep_arr,                    &
                   qvl_dep_arr, qvi_dep_arr,                                  &
                   berg_dep_arr, vap_dep_arr, ice_sublim_arr, mgncol_1 )

              ice_sublim_arr(1) = 0._r8
              vap_dep_arr(1)    = max(vap_dep_arr(1), 0._r8)

              qdep_mass = real(vap_dep_arr(1)) * dt_eff
              qdep_mass = max(0.0, qdep_mass)
              qdep_mass = min(qdep_mass, qliq_mass)

              qliq_mass = max(0.0, qliq_mass - qdep_mass)
              qice_mass = max(0.0, qice_mass + qdep_mass)

           endif

           !---------------------------------------------------------------
           ! Safeguard before final autoconversion.
           !---------------------------------------------------------------

           if(qliq_mass + qice_mass > qsmall_2m) then
              ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / rho_safe
              if(qliq_mass > qsmall_2m .and. nliq_mix < ncond_floor_mix) then
                 nliq_mix = ncond_floor_mix
              endif

              ncond_floor_mix = ncond_floor_cm3_ice * 1.0e6 / rho_safe
              if(qice_mass > qsmall_2m .and. nice_mix < ncond_floor_mix) then
                 nice_mix = ncond_floor_mix
              endif
           endif

           !---------------------------------------------------------------
           ! Final liquid autoconversion.
           !---------------------------------------------------------------

           dq_auto_liq          = 0.0
           dq_auto_ice          = 0.0
           dn_auto_liq_mix      = 0.0
           dn_rain_auto_liq_mix = 0.0
           dn_auto_ice_mix      = 0.0
           dn_snow_auto_ice_mix = 0.0

           qc_liq_auto       = qliq_mass
           nliq_for_auto_mix = nliq_mix

           if(qc_liq_auto > qsmall_2m .and. nliq_for_auto_mix > 0.0) then

              qcic_mg = real(qc_liq_auto, r8)
              rho_mg  = real(rho_safe, r8)
              ncic_mg = real(max(nsmall_2m/rho_safe, nliq_for_auto_mix), r8)

              call size_dist_param_liq(                                    &
                   liq_props_mg, qcic_mg, ncic_mg, rho_mg,                 &
                   pgam_mg, lamc_mg )

              if(pgam_mg > 0._r8 .and. qcic_mg > 0._r8 .and. ncic_mg > 0._r8) then

                 pgam_arr(1) = pgam_mg
                 qcic_arr(1) = qcic_mg
                 ncic_arr(1) = ncic_mg
                 rho_arr(1)  = rho_mg

                 relvar_mg     = 1._r8 / (pgam_mg + 1._r8)
                 relvar_arr(1) = relvar_mg

                 call liu2006_liq_autoconversion(                         &
                      pgam_arr, microp_uniform_liu, qcic_arr, ncic_arr,    &
                      rho_arr, relvar_arr, prc_arr, nprc_arr, nprc1_arr,   &
                      mgncol_1 )

                 dq_auto_liq     = real(prc_arr(1))  * dt_eff
                 dn_auto_liq_mix = real(nprc1_arr(1)) * dt_eff

                 dq_auto_liq     = max(0.0, min(qc_liq_auto, dq_auto_liq))
                 dn_auto_liq_mix = max(0.0, min(nliq_for_auto_mix, dn_auto_liq_mix))

                 if(qc_liq_auto > qsmall_2m .and. nliq_for_auto_mix > 0.0) then
                    dn_auto_liq_mix = min(dn_auto_liq_mix,                 &
                                           nliq_for_auto_mix *              &
                                           dq_auto_liq / qc_liq_auto)
                 endif

                 dq_auto_liq     = real(auto_eff_liq_2m) * dq_auto_liq
                 dn_auto_liq_mix = real(auto_eff_liq_2m) * dn_auto_liq_mix

                 if(use_mp_auto_number_2m) then
                    call cup_2m_precip_number_from_auto_mass(              &
                         dq_auto_liq, dn_auto_liq_mix, rho_safe,          &
                         rho_liq_mg, n0r_mp_auto_2m, lam_bnd_rain_2m,     &
                         dn_rain_auto_liq_mix )
                 else
                    dn_rain_auto_liq_mix = max(0.0, real(auto_eff_liq_2m) * &
                                                     real(nprc_arr(1)) * dt_eff)
                    if(dn_auto_liq_mix > 0.0) then
                       dn_rain_auto_liq_mix = min(dn_rain_auto_liq_mix, dn_auto_liq_mix)
                    endif
                 endif

              endif

           endif

           qliq_mass = max(0.0, qliq_mass - dq_auto_liq)
           nliq_mix  = max(0.0, nliq_mix  - dn_auto_liq_mix)

           !---------------------------------------------------------------
           ! Final ice autoconversion.
           !---------------------------------------------------------------

           qc_ice_auto = qice_mass
           nice_mix_before_ice_auto = nice_mix

           if(qc_ice_auto > qsmall_2m .and. nice_mix_before_ice_auto > 0.0) then

              tice_mg = real(tempc(i,k), r8)
              qiic_mg = real(qc_ice_auto, r8)
              niic_mg = real(nice_mix_before_ice_auto, r8)

              call size_dist_param_ice(                                    &
                   ice_props_mg, qiic_mg, niic_mg, lami_mg )

              if(lami_mg > 0._r8 .and. qiic_mg > 0._r8 .and. niic_mg > 0._r8) then

                 tice_arr(1) = tice_mg
                 qiic_arr(1) = qiic_mg
                 lami_arr(1) = lami_mg
                 niic_arr(1) = niic_mg

                 call ice_autoconversion(                                  &
                      tice_arr, qiic_arr, lami_arr, niic_arr,              &
                      dcs_mg, prci_arr, nprci_arr, mgncol_1, ac_time_mg )

                 dq_auto_ice     = real(prci_arr(1))  * dt_eff
                 dn_auto_ice_mix = real(nprci_arr(1)) * dt_eff

                 dq_auto_ice     = max(0.0, min(qc_ice_auto, dq_auto_ice))
                 dn_auto_ice_mix = max(0.0, min(nice_mix_before_ice_auto, dn_auto_ice_mix))

                 if(use_mp_auto_number_2m) then
                    call cup_2m_precip_number_from_auto_mass(              &
                         dq_auto_ice, dn_auto_ice_mix, rho_safe,          &
                         rho_snow_2m, n0s_mp_auto_2m, lam_bnd_snow_2m,    &
                         dn_snow_auto_ice_mix )
                 else
                    dn_snow_auto_ice_mix = dn_auto_ice_mix
                 endif

              endif

           endif

           qice_mass = max(0.0, qice_mass - dq_auto_ice)
           nice_mix  = max(0.0, nice_mix  - dn_auto_ice_mix)

           qcloud_after_all_auto = max(0.0, qliq_mass + qice_mass)

           ! Final post-autoconversion safeguards.
           if(qcloud_after_all_auto > qsmall_2m) then
              ncond_floor_mix = ncond_floor_cm3_liq * 1.0e6 / rho_safe
              if(qliq_mass > qsmall_2m .and. nliq_mix < ncond_floor_mix) then
                 nliq_mix = ncond_floor_mix
              endif

              ncond_floor_mix = ncond_floor_cm3_ice * 1.0e6 / rho_safe
              if(qice_mass > qsmall_2m .and. nice_mix < ncond_floor_mix) then
                 nice_mix = ncond_floor_mix
              endif
           endif

           qrc(i,k)         = qcloud_after_all_auto
           qice_up_mix(i,k) = max(0.0, min(qice_mass, qcloud_after_all_auto))

           nliq_up_mix(i,k)  = max(0.0, nliq_mix)
           nice_up_mix(i,k)  = max(0.0, nice_mix)
           ndust_up_mix(i,k)    = max(0.0, ndust_mix)
           nsoot_up_mix(i,k)    = max(0.0, nsoot_mix)
           nseasalt_up_mix(i,k) = max(0.0, nseasalt_mix)
           nact_up_mix(i,k)     = max(0.0, nliq_up_mix(i,k) + nice_up_mix(i,k))

           nact_up_m3_out(i,k) = nact_up_mix(i,k) * rho_safe

           if(qcloud_after_all_auto > qsmall_2m) then
              pice_2m(i,k) = max(0.0, min(1.0, qice_up_mix(i,k)/qcloud_after_all_auto))
           else
              pice_2m(i,k) = 0.0
           endif

           qice_up_2m(i,k) = max(0.0, min(qice_up_mix(i,k), qcloud_after_all_auto))
           qliq_up_2m(i,k) = max(0.0, qcloud_after_all_auto - qice_up_2m(i,k))
           nliq_up_2m(i,k) = max(0.0, nliq_up_mix(i,k))
           nice_up_2m(i,k) = max(0.0, nice_up_mix(i,k))

           ! Return the full-2M phase-resolved retained condensate as the
           ! unique cloud condensate state used by the parent detrainment path.
           qrc(i,k) = qliq_up_2m(i,k) + qice_up_2m(i,k)
           qice_up_mix(i,k) = qice_up_2m(i,k)
           if(qrc(i,k) > qsmall_2m) then
              pice_2m(i,k) = max(0.0, min(1.0, qice_up_2m(i,k)/qrc(i,k)))
           else
              pice_2m(i,k) = 0.0
           endif
           qc(i,k) = qrc(i,k) + min(qc(i,k), qrch)

           dq_auto_liq_col(k) = max(0.0, dq_auto_liq)
           dq_auto_ice_col(k) = max(0.0, dq_auto_ice)
           dn_auto_liq_col(k)      = max(0.0, dn_auto_liq_mix)
           dn_auto_ice_col(k)      = max(0.0, dn_auto_ice_mix)
           dn_rain_auto_liq_col(k) = max(0.0, dn_rain_auto_liq_mix)
           dn_snow_auto_ice_col(k) = max(0.0, dn_snow_auto_ice_mix)

           dq_precip_level = dq_auto_liq + dq_auto_ice                    &
                           + dq_cw_rain + dq_cw_snow_to_snow + dq_ci_snow

           pw(i,k) = dq_precip_level * zu(i,k)

           qc(i,k) = qrc(i,k) + min(qc(i,k), qrch)

           pwav(i) = pwav(i) + pw(i,k)
           psum(i) = psum(i) + clw_all(i,k)*zu(i,k)*dz_eff

           if(debug_2m) then

              qprecip_mix = max(0.0, dq_precip_level)
              qtotal_for_ratio = max(qc(i,k) + qprecip_mix, 1.0e-30)
              precip_to_total  = qprecip_mix / qtotal_for_ratio

              nact_lev_cm3 = nact_lev_mix_col(k) * rho_safe * 1.0e-6
              nwfa_cm3     = nwfa_cm3_col(k)
              ndust_cm3    = ndust_mix        * rho_safe * 1.0e-6
              nliq_cm3     = nliq_up_mix(i,k) * rho_safe * 1.0e-6
              nice_L       = nice_up_mix(i,k) * rho_safe * 1.0e-3
              ninp_imm_L   = ninp_imm_mix     * rho_safe * 1.0e-3

              if(k == kbase_2m) then
                 write(*,'(a)') ' '
                 write(*,'(a,i6,a,i6,a,a)')                              &
                      'CUP_2M_DIAG: i=', i, ' j=', jcol, ' plume=', trim(name)
                 write(*,'(a)')                                           &
                      '   k    T[K]   W[m/s]  pice Nact[#/cm3] NWFA[#/cm3]' // &
                      ' Dust[#/cm3] INP[#/L]  Ni[#/L]  Nl[#/cm3]' //      &
                      ' qcond[kg/kg] qprecip/qtot qrain qsnow'
              endif

              write(*,'(i5,1x,f7.2,1x,f7.2,1x,f5.2,1x,es10.1,1x,es10.1,' // &
                       '1x,es10.1,1x,es10.1,1x,es10.1,1x,es10.1,1x,es10.1,' // &
                       '1x,es10.1,1x,es10.1,1x,es10.1)')                 &
                   k, tempc(i,k), w_col(k), pice_2m(i,k), nact_lev_cm3,  &
                   nwfa_cm3, ndust_cm3, ninp_imm_L, nice_L, nliq_cm3,    &
                   qrc(i,k), precip_to_total, qrain_diag_col(k), qsnow_diag_col(k)

           endif

           ! Do not compute mass or number detrainment tendencies here.
           ! The parent applies the exact legacy detrainment operator to
           ! qliq_up_2m/qice_up_2m/nliq_up_2m/nice_up_2m.

        enddo

     if(pwav(i) < 0.0) then
        ierr(i)  = 66
        ierrc(i) = "pwav negative"
     endif

  enddo

  !---------------------------------------------------------------------------
  ! Convert qc back from total water to water vapor in cloud.
  !---------------------------------------------------------------------------

  do i = its, itf
     if(ierr(i) /= 0) cycle

     do k = kts, ktop(i) + 1
        qc(i,k) = qc(i,k) - qrc(i,k)
     enddo
  enddo


END SUBROUTINE cup_up_moisture_2M


!============================ helpers ========================= 

 SUBROUTINE get_delmix_2m(cumulus,kts,kte,ktf,xland,subcl_level,po,ain,aout)
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

 end SUBROUTINE get_delmix_2m
  
!----------------------------------------------------------------------

  real FUNCTION cup_2m_w_from_w2(w2_in)

    implicit none

    real, intent(in) :: w2_in

    cup_2m_w_from_w2 = sqrt(max(wmin_2m*wmin_2m, w2_in))
    cup_2m_w_from_w2 = max(wmin_2m, cup_2m_w_from_w2)
    cup_2m_w_from_w2 = min(wmax_2m, cup_2m_w_from_w2)

  END FUNCTION cup_2m_w_from_w2


  real FUNCTION cup_2m_next_w2( w2_in,                         &
                                zc_km1, zc_k,                  &
                                tenv_km1, tenv_k,              &
                                qenv_km1, qenv_k,              &
                                tcld_km1, tcld_k,              &
                                qcld_km1, qcld_k,              &
                                qcond_km1, qcond_k,            &
                                entr_km1, cd_km1 )

    implicit none

    real, intent(in) :: w2_in
    real, intent(in) :: zc_km1, zc_k
    real, intent(in) :: tenv_km1, tenv_k
    real, intent(in) :: qenv_km1, qenv_k
    real, intent(in) :: tcld_km1, tcld_k
    real, intent(in) :: qcld_km1, qcld_k
    real, intent(in) :: qcond_km1, qcond_k
    real, intent(in) :: entr_km1, cd_km1

    real, parameter :: f_w       = 2.0
    real, parameter :: gam_w     = 0.5
    real, parameter :: ftun_load = 0.5

    real :: dz_w
    real :: tve
    real :: tvc
    real :: buoy
    real :: dw_buoy
    real :: kx
    real :: w2_raw

    dz_w = zc_k - zc_km1

    if(dz_w <= 0.0) then
       cup_2m_next_w2 = w2_in
       return
    endif

    ! Use the standard virtual-temperature form:
    !   Tv = T * (1 + q/eps) / (1 + q)
    tve = 0.5 * ( cup_2m_virtual_temperature(tenv_km1, qenv_km1) +     &
                  cup_2m_virtual_temperature(tenv_k,   qenv_k  ) )

    tvc = 0.5 * ( cup_2m_virtual_temperature(tcld_km1, qcld_km1) +     &
                  cup_2m_virtual_temperature(tcld_k,   qcld_k  ) )

    ! Buoyancy minus condensate loading.
    buoy = g * ( (tvc - tve) / max(tve, 1.0e-12) -                    &
                 ftun_load * 0.5 * (qcond_k + qcond_km1) )

    dw_buoy = 2.0 / (f_w * (1.0 + gam_w)) * buoy * dz_w

    ! Entrainment/detrainment damping.
    kx = max(0.0, max(entr_km1, cd_km1) * dz_w)

    w2_raw = (dw_buoy + w2_in - 2.0*kx*w2_in) / (1.0 + kx)

    if(w2_raw < 0.0) then
       w2_raw = 0.5*w2_in
    endif

    cup_2m_next_w2 = max(0.0, w2_raw)

  END FUNCTION cup_2m_next_w2


  real FUNCTION cup_2m_virtual_temperature(temp_k, qv_mix)

    implicit none

    real, intent(in) :: temp_k
    real, intent(in) :: qv_mix

    real :: qsafe

    qsafe = max(0.0, qv_mix)

    ! Standard virtual temperature using water-vapor mixing ratio qv:
    !   Tv = T * (1 + qv/eps) / (1 + qv)
    cup_2m_virtual_temperature = temp_k *                              &
                                 (1.0 + qsafe/eps_vap_2m) /            &
                                 (1.0 + qsafe)

  END FUNCTION cup_2m_virtual_temperature


  real FUNCTION cup_2m_nwfa_m3(nmodes_loc, aer_num_m3, aer_kap)

    implicit none

    integer, intent(in) :: nmodes_loc
    real, dimension(:), intent(in) :: aer_num_m3
    real, dimension(:), intent(in) :: aer_kap

    integer :: imode

    cup_2m_nwfa_m3 = 0.0

    do imode = 1, nmodes_loc
       if(aer_kap(imode) > kappa_wfa_min_2m) then
          cup_2m_nwfa_m3 = cup_2m_nwfa_m3 + max(0.0, aer_num_m3(imode))
       endif
    enddo

  END FUNCTION cup_2m_nwfa_m3



  SUBROUTINE aerosol_ice_species_from_aerp(AerP_loc, nmodes_loc, ii, jj, kk, rho_air, &
                                           ndust_mix, ddust_mean_m, sigdust_mean,     &
                                           nsoot_mix, dsoot_mean_m, sigsoot_mean,     &
                                           nseasalt_mix, dseasalt_mean_m, sigseasalt_mean)

    ! Diagnose finite aerosol reservoirs used by heterogeneous ice nucleation.
    ! All returned number reservoirs are mixing ratios [# kg-1].  Diameter and
    ! sigma are number-weighted over valid modes, using the same averaging for
    ! sigma as for diameter.

    implicit none

    type(AerPropsNew), dimension(:), intent(in) :: AerP_loc
    integer, intent(in) :: nmodes_loc
    integer, intent(in) :: ii, jj, kk
    real, intent(in)  :: rho_air

    real, intent(out) :: ndust_mix, ddust_mean_m, sigdust_mean
    real, intent(out) :: nsoot_mix, dsoot_mean_m, sigsoot_mean
    real, intent(out) :: nseasalt_mix, dseasalt_mean_m, sigseasalt_mean

    integer :: im
    real :: rho_safe
    real :: nmode_m3
    real :: dmode_m
    real :: sigmode
    real :: fdust_mode
    real :: fsoot_mode
    real :: kap_mode
    real :: ndust_m3, nsoot_m3, nseasalt_m3
    real :: dsum_dust, dsum_soot, dsum_seasalt
    real :: ssum_dust, ssum_soot, ssum_seasalt

    rho_safe = max(rho_air, 1.0e-12)

    ndust_m3 = 0.0
    nsoot_m3 = 0.0
    nseasalt_m3 = 0.0

    dsum_dust = 0.0
    dsum_soot = 0.0
    dsum_seasalt = 0.0

    ssum_dust = 0.0
    ssum_soot = 0.0
    ssum_seasalt = 0.0

    do im = 1, nmodes_loc

       nmode_m3    = max(0.0,    AerP_loc(im)%num(ii,jj,kk))
       dmode_m     = max(1.0e-9, AerP_loc(im)%dpg(ii,jj,kk))
       sigmode     = max(1.0e-6, AerP_loc(im)%sig(ii,jj,kk))
       fdust_mode  = max(0.0,    AerP_loc(im)%fdust(ii,jj,kk))
       fsoot_mode  = max(0.0,    AerP_loc(im)%fsoot(ii,jj,kk))
       kap_mode    = max(0.0,    AerP_loc(im)%kap(ii,jj,kk))

       if(nmode_m3 > 0.0 .and. fdust_mode > dust_frac_min_2m) then
          ndust_m3  = ndust_m3  + nmode_m3
          dsum_dust = dsum_dust + nmode_m3*dmode_m
          ssum_dust = ssum_dust + nmode_m3*sigmode
       endif

       if(nmode_m3 > 0.0 .and. fsoot_mode > soot_frac_min_2m) then
          nsoot_m3  = nsoot_m3  + nmode_m3
          dsum_soot = dsum_soot + nmode_m3*dmode_m
          ssum_soot = ssum_soot + nmode_m3*sigmode
       endif

       if(nmode_m3 > 0.0 .and. kap_mode > seasalt_kappa_min_2m) then
          nseasalt_m3  = nseasalt_m3  + nmode_m3
          dsum_seasalt = dsum_seasalt + nmode_m3*dmode_m
          ssum_seasalt = ssum_seasalt + nmode_m3*sigmode
       endif

    enddo

    if(ndust_m3 > 0.0) then
       ndust_mix     = ndust_m3/rho_safe
       ddust_mean_m  = dsum_dust/ndust_m3
       sigdust_mean  = ssum_dust/ndust_m3
    else
       ndust_mix     = 0.0
       ddust_mean_m  = 0.0
       sigdust_mean  = 0.0
    endif

    if(nsoot_m3 > 0.0) then
       nsoot_mix     = nsoot_m3/rho_safe
       dsoot_mean_m  = dsum_soot/nsoot_m3
       sigsoot_mean  = ssum_soot/nsoot_m3
    else
       nsoot_mix     = 0.0
       dsoot_mean_m  = 0.0
       sigsoot_mean  = 0.0
    endif

    if(nseasalt_m3 > 0.0) then
       nseasalt_mix     = nseasalt_m3/rho_safe
       dseasalt_mean_m  = dsum_seasalt/nseasalt_m3
       sigseasalt_mean  = ssum_seasalt/nseasalt_m3
    else
       nseasalt_mix     = 0.0
       dseasalt_mean_m  = 0.0
       sigseasalt_mean  = 0.0
    endif

  END SUBROUTINE aerosol_ice_species_from_aerp


  real FUNCTION cup_2m_aerosol_area(dmean_m, sig_ln, acorr)

    ! Effective heterogeneous freezing surface area per particle [m2].
    ! Uses the same lognormal moment structure as the dust-only helper, but
    ! with species-specific correction factors and the mode-averaged sigma.

    implicit none

    real, intent(in) :: dmean_m
    real, intent(in) :: sig_ln
    real, intent(in) :: acorr
    real :: sig_safe

    if(dmean_m <= 0.0) then
       cup_2m_aerosol_area = 0.0
       return
    endif

    sig_safe = sig_ln
    if(sig_safe <= 0.0) sig_safe = log(2.0)

    cup_2m_aerosol_area = 0.52*dmean_m**3 * acorr * exp(min(50.0, 4.5*sig_safe*sig_safe))
    cup_2m_aerosol_area = max(0.0, cup_2m_aerosol_area)

  END FUNCTION cup_2m_aerosol_area


  SUBROUTINE cup_2m_inp_freezing_sources(temp_k, pres_hpa, rho_air, w_up, dt_step, &
                                         qliq, nliq,                               &
                                         ndust, ddust_m, sigdust,                  &
                                         nsoot, dsoot_m, sigsoot,                  &
                                         nseasalt, dseasalt_m, sigseasalt,         &
                                         ninp_dust, ninp_soot, ninp_seasalt,       &
                                         ninp_bkg, ninp_total)

    ! Heterogeneous freezing source for the convective 2M plume [# kg-1 over dt].
    !
    ! immersion and contact modes
    ! following Barahona et al. GMD (2014).  
    !
    ! Components:
    !   - Dust immersion: active-site ns(T) form from Niemand et al. (2012).
    !   - Soot immersion: active-site ns(T) form from Murray et al. review (2012).
    !   - Dust/soot contact: pseudo-contact source, using immersion ns evaluated
    !     at T-3 K,  plus Brownian diffusional collection.
    !   - Sea-salt immersion: number-based DeMott et al. PNAS (2015) expression.
    !   - Background INP: infinite reservoir Phillips
    !     et al. (2007)-style rate, coupled only to the presence of liquid.
    !
    ! Supersaturation over ice at water saturation follows the Murphy and Koop
    ! (2005)-based linear approximation used in the original helper.

    implicit none

    real, intent(in)  :: temp_k
    real, intent(in)  :: pres_hpa
    real, intent(in)  :: rho_air
    real, intent(in)  :: w_up
    real, intent(in)  :: dt_step
    real, intent(in)  :: qliq
    real, intent(in)  :: nliq
    real, intent(in)  :: ndust, ddust_m, sigdust
    real, intent(in)  :: nsoot, dsoot_m, sigsoot
    real, intent(in)  :: nseasalt, dseasalt_m, sigseasalt
    real, intent(out) :: ninp_dust, ninp_soot, ninp_seasalt
    real, intent(out) :: ninp_bkg, ninp_total

    real :: tx, taux
    real :: nsdust, nssoot, nsss
    real :: dnsd, dnss, dnsss
    real :: nsdust_c, nssoot_c
    real :: min_ns_dust, min_ns_soot, min_ns_seasalt
    real :: coolr, wsub
    real :: ahet
    real :: viscosity, mfp, nslip, ndfaer, lam_r
    real :: pres_pa, rho_safe
    real :: dust_imm, dust_con
    real :: soot_imm, soot_con
    real :: ss_imm
    real :: si, bkg_rate_L_s
    real :: scale_local

    ninp_dust    = 0.0
    ninp_soot    = 0.0
    ninp_seasalt = 0.0
    ninp_bkg     = 0.0
    ninp_total   = 0.0

    if(qliq <= 1.0e-10 .or. nliq <= 1.0 .or. dt_step <= 0.0) return
    if(temp_k >= 272.0 .or. temp_k <= thom_2m) return

    rho_safe = max(rho_air, 1.0e-12)
    pres_pa  = max(100.0*pres_hpa, 1.0)

    ! Minimum active-site concentrations:
    !   dust    threshold chosen to limit nucleation to about -12 C,
    !   soot    threshold chosen to limit nucleation to about -18 C,
    !   sea-salt threshold chosen to limit nucleation to about -5 C.
    min_ns_dust    = 3.75e6
    min_ns_soot    = 3.75e9
    min_ns_seasalt = 4.0e2

    ! Saturated cooling-rate approximation.
    wsub  = max(min(w_up, 10.0), 0.8)
    coolr = 5.0e-3*wsub

    tx = max(temp_k - 273.16, -38.0)

    ! Droplet distribution slope used by the contact-freezing collection term.
    lam_r = (pi_r8*950.0*nliq/(rho_safe*max(qliq,1.0e-30)))**(1.0/3.0)
    lam_r = min(max(lam_r, 1.0e-12), 1.0e8)

    ! Air viscosity and molecular mean free path used for Brownian contact
    ! freezing collection, as in the original contact-nucleation block.
    viscosity = 1.8e-5*(temp_k/298.0)**0.85
    mfp = 2.0*viscosity/(pres_pa*sqrt(8.0*28.96e-3/(pi_r8*8.314409*temp_k)))

    ! Immersion active-site densities:
    !   dust: Niemand et al. (2012)
    !   soot: Murray et al. review (2012)
    nsdust = max(exp(min(80.0, -0.517*tx + 9.3585)) - min_ns_dust, 0.0)
    nssoot = max(1.0e4*exp(min(80.0, -0.0101*tx*tx - 0.8525*tx + 0.7667)) - min_ns_soot, 0.0)
    dnsd   = 0.517*nsdust
    dnss   = max(-(-2.0*0.0101*tx - 0.8525)*nssoot, 0.0)

    ! Contact nucleation ns is
    ! equivalent to immersion ns evaluated at T-3 K.
    taux = max(tx - 3.0, -35.0)
    nsdust_c = max(exp(min(80.0, -0.517*taux + 9.3585)) - min_ns_dust, 0.0)
    nssoot_c = max(1.0e4*exp(min(80.0, -0.0101*taux*taux - 0.8525*taux + 0.7667)) - min_ns_soot, 0.0)

    dust_imm = 0.0
    dust_con = 0.0
    if(ndust > 0.0 .and. ddust_m > 0.0) then
       ! Dust effective area uses the Murray et al. (2011) nonsphericity /
       ! aggregation correction, equivalent to about 10 m2 g-1 
       ! The lognormal-width factor is evaluated with the transported
       ! species-mean sigma.
       ahet = cup_2m_aerosol_area(ddust_m, sigdust, acorr_dust_2m)
       if(ahet > 0.0) then
          dust_imm = ndust*exp(-min(80.0,nsdust*ahet))*dnsd*coolr*ahet*fdrop_dust_2m
          ! Slip-corrected Brownian diffusivity and contact collection.
          nslip = 1.0 + (2.0*mfp/ddust_m)*(1.257 + 0.4*exp(-(1.1*ddust_m*0.5/mfp)))
          ndfaer = 1.381e-23*temp_k*nslip*(1.0-exp(-min(80.0,nsdust_c*ahet))) / &
                   (12.0*pi_r8*viscosity*ddust_m)
          dust_con = 2.0*pi_r8*ndfaer*ndust*nliq/lam_r
       endif
    endif
    ninp_dust = max(0.0, min(ndust, (dust_imm + dust_con)*dt_step))

    soot_imm = 0.0
    soot_con = 0.0
    if(nsoot > 0.0 .and. dsoot_m > 0.0) then
       ! Soot effective area uses the Popovicheva (2003)
       ! equivalent to about 50 m2 g-1.
       ahet = cup_2m_aerosol_area(dsoot_m, sigsoot, acorr_soot_2m)
       if(ahet > 0.0) then
          soot_imm = nsoot*exp(-min(80.0,nssoot*ahet))*dnss*coolr*ahet*fdrop_soot_2m
          ! Slip-corrected Brownian diffusivity and contact collection.
          nslip = 1.0 + (2.0*mfp/dsoot_m)*(1.257 + 0.4*exp(-(1.1*dsoot_m*0.5/mfp)))
          ndfaer = 1.381e-23*temp_k*nslip*(1.0-exp(-min(80.0,nssoot_c*ahet))) / &
                   (12.0*pi_r8*viscosity*dsoot_m)
          soot_con = 2.0*pi_r8*ndfaer*nsoot*nliq/lam_r
       endif
    endif
    ninp_soot = max(0.0, min(nsoot, (soot_imm + soot_con)*dt_step))

    ss_imm = 0.0
    if(nseasalt > 0.0) then
       ! Sea-salt INP follows DeMott et al. PNAS (2015)
       ! tendency:
       !   INsea_rate =  dnsss * cooling_rate
       ! It intentionally does not use aerosol surface area nor number, only limited by available sea salt.
       nsss  = max(exp(min(80.0, -0.459*temp_k + 128.6235)) - min_ns_seasalt, 0.0)
       dnsss = max(0.459*nsss, 0.0)
       ss_imm = dnsss*coolr*fdrop_seasalt_2m
    endif
    ninp_seasalt = max(0.0, min(nseasalt, ss_imm*dt_step))

    ! Infinite background INP reservoir.  Only requires liquid to be present;
    ! it is not depleted.  Ice supersaturation at water saturation uses the
    ! Murphy and Koop (2005)-based approximation. Phillips et al. (2007)-style
    ! background source, converted from # L-1 s-1 to # kg-1 s-1.
    if(temp_k < 260.0) then
       si = max(0.0, -1.2379e-2*temp_k + 3.3595)
       bkg_rate_L_s = coolr*42.8*exp(min(50.0,3.88*si))*bkg_inp_scaling
       ninp_bkg = max(0.0, bkg_rate_L_s*1000.0/rho_safe*dt_step)
    endif

    ninp_total = ninp_dust + ninp_soot + ninp_seasalt + ninp_bkg

    !print *, temp_k, ninp_total , ninp_dust , ninp_soot , ninp_seasalt , ninp_bkg, nliq
    
    if(ninp_total > nliq .and. ninp_total > 1.0e-30) then
       scale_local = max(0.0, min(1.0, nliq/ninp_total))
       ninp_dust    = ninp_dust*scale_local
       ninp_soot    = ninp_soot*scale_local
       ninp_seasalt = ninp_seasalt*scale_local
       ninp_bkg     = ninp_bkg*scale_local
       ninp_total   = nliq
    endif

  END SUBROUTINE cup_2m_inp_freezing_sources


  SUBROUTINE activate_from_aerp_arg2002_local( &
       NMODES, AER_NUM_M3, AER_DPG_M, AER_SIG, AER_KAP, &
       W, TEMP, RHO_AIR, NDROP_CURR_MIX, NCPL_ACT_TARGET_MIX, NCPL_ACT_MIX )

    implicit none

    integer, intent(in) :: NMODES
    real, dimension(:), intent(inout) :: AER_NUM_M3  ! remaining aerosol number, # m-3
    real, dimension(:), intent(in)    :: AER_DPG_M   ! modal geometric mean diameter, m
    real, dimension(:), intent(in)    :: AER_SIG     ! ln(sigma_g)
    real, dimension(:), intent(in)    :: AER_KAP     ! hygroscopicity kappa
    real, intent(in)  :: W              ! local updraft velocity, m s-1
    real, intent(in)  :: TEMP           ! K
    real, intent(in)  :: RHO_AIR        ! kg m-3
    real, intent(in)  :: NDROP_CURR_MIX ! transported droplet population before activation, # kg-1
    real, intent(out) :: NCPL_ACT_TARGET_MIX ! ARG2002 activated target, # kg-1
    real, intent(out) :: NCPL_ACT_MIX        ! newly activated tendency, # kg-1

    integer :: INDEX
    real :: auxx, aux1, T
    real :: kappa, fi, gi, nui, citai, ui, smax
    real :: WPOS
    real :: Akoh, alfa_arg, beta_arg, growth_arg
    real :: rho_safe
    real :: ncurr_m3
    real :: ncpl_act_target_sum_m3
    real :: ncpl_act_delta_sum_m3
    real :: deplete_scale
    real :: aux1_m3

    real, dimension(NMODES) :: TPI
    real, dimension(NMODES) :: DPGI
    real, dimension(NMODES) :: SIGI
    real, dimension(NMODES) :: KAPI
    real, dimension(NMODES) :: SMI
    real, dimension(NMODES) :: ACT_POT_M3

    !---------------------------------------------------------------
    ! AER_NUM_M3 is the local remaining aerosol concentration
    ! reservoir for this plume column.  ARG2002 uses number concentration
    ! (# m-3), so no unit conversion is applied before calculating smax.
    ! The external AerP argument is not modified.
    !
    ! Activation behavior is controlled by use_activation_target_increment_2m:
    !   .true. : dNdrop_act = max(Ndrop_act - Ndrop, 0)
    !   .false.: dNdrop_act = Ndrop_act
    ! Only the returned positive increment is depleted from the local aerosol
    ! reservoir.  The .false. option restores the older additive behavior.
    !---------------------------------------------------------------

    rho_safe = max(RHO_AIR, 1.0e-12)
    ncurr_m3 = max(0.0, NDROP_CURR_MIX) * rho_safe

    do INDEX = 1, NMODES
       TPI (INDEX) = max(0.0,    AER_NUM_M3(INDEX))
       DPGI(INDEX) = max(1.0e-9, AER_DPG_M(INDEX))
       SIGI(INDEX) = max(1.0e-6, AER_SIG(INDEX))
       KAPI(INDEX) = max(0.0,    AER_KAP(INDEX))
       SMI (INDEX) = 0.0
       ACT_POT_M3(INDEX) = 0.0
    enddo

    WPOS = max(W, 1.0e-3)

    !---------------------------------------------------------------
    ! Calculate constants. These fits were obtained from detailed
    ! correlations of physical properties. growth_arg is actually 1/G
    ! in the ARG2002-style formulation.
    !---------------------------------------------------------------

    T = min(max(TEMP, 243.0), 323.0)

    alfa_arg   = 2.8915E-08*(T**2) - 2.1328E-05*T + 4.2523E-03
    beta_arg   = exp(3.49996E-04*T**2 - 2.27938E-01*T + 4.20901E+01)
    growth_arg = exp(-2.94362E-06*T**3 + 2.77941E-03*T**2 - 8.92889E-01*T + 1.18787E+02)
    Akoh       = 0.66e-6/T

    !===============================================================
    ! Calculate maximum supersaturation following ARG2002-like form.
    !===============================================================

    auxx = 0.0

    DO INDEX = 1, NMODES

       kappa = max(KAPI(INDEX), 0.001)

       SMI(INDEX) = ((0.667*Akoh/DPGI(INDEX))**1.5) / sqrt(2.0*kappa)
       SMI(INDEX) = max(SMI(INDEX), 1.0e-5)

       if ((TPI(INDEX) .gt. 1.0e4) .and. (kappa .gt. 0.1)) then

          fi = 0.5 * exp(2.5*SIGI(INDEX))
          gi = 1.0 + 0.25*SIGI(INDEX)

          nui = ((alfa_arg*WPOS*growth_arg)**1.5) / &
                (2.0*pi_r8*980.0*beta_arg*TPI(INDEX))

          citai = 0.667*Akoh*sqrt(alfa_arg*WPOS*growth_arg)

          aux1 = fi*((citai/nui)**1.5) + &
                 gi*(SMI(INDEX)*SMI(INDEX)/(nui+(3.0*citai)))**0.75

          aux1 = aux1/(SMI(INDEX)*SMI(INDEX))

          auxx = auxx + aux1

       endif

    ENDDO

    !===============================================================
    ! Calculate the activated target, then only apply/deplete the
    ! positive increment above the transported droplet population.
    !===============================================================

    ncpl_act_target_sum_m3 = 0.0
    ncpl_act_delta_sum_m3  = 0.0

    if (auxx .gt. 0.0) then

       smax = 1.0/sqrt(auxx)

       DO INDEX = 1, NMODES

          if ((TPI(INDEX) .gt. 1.0e4) .and. (KAPI(INDEX) .gt. 0.1)) then

             ui = sqrt(2.0)*log(SMI(INDEX)/smax)/3.0

             aux1_m3 = 0.5*TPI(INDEX)*(1.0 - erfapp(ui))
             aux1_m3 = max(0.0, min(aux1_m3, AER_NUM_M3(INDEX)))

             ACT_POT_M3(INDEX) = aux1_m3
             ncpl_act_target_sum_m3 = ncpl_act_target_sum_m3 + aux1_m3

          endif

       ENDDO

       if(use_activation_target_increment_2m) then
          ncpl_act_delta_sum_m3 = max(0.0, ncpl_act_target_sum_m3 - ncurr_m3)
          ncpl_act_delta_sum_m3 = min(ncpl_act_delta_sum_m3, ncpl_act_target_sum_m3)
       else
          ! Legacy behavior: the diagnosed activated number is the actual
          ! additive source, independent of the currently transported Ndrop.
          ncpl_act_delta_sum_m3 = ncpl_act_target_sum_m3
       endif

       if(ncpl_act_target_sum_m3 > 1.0e-30) then
          deplete_scale = max(0.0, min(1.0, ncpl_act_delta_sum_m3/ncpl_act_target_sum_m3))
       else
          deplete_scale = 0.0
       endif

       DO INDEX = 1, NMODES
          if(ACT_POT_M3(INDEX) > 0.0) then
             ! Deplete only the local activation reservoir and only by the
             ! newly activated increment.  Do not touch dust/INP reservoirs.
             AER_NUM_M3(INDEX) = max(0.0, AER_NUM_M3(INDEX) - deplete_scale*ACT_POT_M3(INDEX))
          endif
       ENDDO

    endif

    NCPL_ACT_TARGET_MIX = max(0.0, ncpl_act_target_sum_m3 / rho_safe)
    NCPL_ACT_MIX        = max(0.0, ncpl_act_delta_sum_m3  / rho_safe)

  END SUBROUTINE activate_from_aerp_arg2002_local



  SUBROUTINE cup_2m_apply_fixed_precip_accretion(                         &
       temp_k, rho_air, dt_step,                                           &
       qrain, nrain, qsnow, nsnow,                                         &
       lamr, n0r, lams, n0s, umr, ums, unr, uns,                           &
       liq_props,                                                          &
       qliq, nliq, qice, nice,                                             &
       dq_cw_rain, dq_cw_snow_to_snow, dq_ci_snow, dq_sip, dq_rs )

    !-----------------------------------------------------------------------
    ! Apply cloud-precipitation collection against fixed diagnostic rain/snow.
    !
    ! The rain/snow fields are held fixed in this final corrector pass.  This
    ! helper updates cloud liquid/ice mass and number only.  Rain-snow
    ! collection is returned diagnostically but does not change the fixed
    ! precipitation profile.
    !-----------------------------------------------------------------------

    implicit none

    real, intent(in) :: temp_k
    real, intent(in) :: rho_air
    real, intent(in) :: dt_step
    real, intent(in) :: qrain
    real, intent(in) :: nrain
    real, intent(in) :: qsnow
    real, intent(in) :: nsnow
    real, intent(in) :: lamr
    real, intent(in) :: n0r
    real, intent(in) :: lams
    real, intent(in) :: n0s
    real, intent(in) :: umr
    real, intent(in) :: ums
    real, intent(in) :: unr
    real, intent(in) :: uns

    type(MGHydrometeorProps), intent(in) :: liq_props

    real, intent(inout) :: qliq
    real, intent(inout) :: nliq
    real, intent(inout) :: qice
    real, intent(inout) :: nice

    real, intent(out) :: dq_cw_rain
    real, intent(out) :: dq_cw_snow_to_snow
    real, intent(out) :: dq_ci_snow
    real, intent(out) :: dq_sip
    real, intent(out) :: dq_rs

    integer, parameter :: mgncol_loc = 1
    logical :: microp_uniform_acc_loc

    real :: rho_safe
    real :: dt_safe
    real :: dq_cw_snow_total
    real :: dn_cw_rain
    real :: dn_cw_snow
    real :: dn_ci_snow
    real :: dn_sip
    real :: dn_rs
    real :: liq_sink_mass
    real :: liq_sink_num
    real :: ice_sink_mass
    real :: ice_sink_num
    real :: scale_sink

    real(r8) :: qcic_mg_loc
    real(r8) :: ncic_mg_loc
    real(r8) :: rho_mg_loc
    real(r8) :: pgam_mg_loc
    real(r8) :: lamc_mg_loc

    real(r8), dimension(mgncol_loc) :: t_arr
    real(r8), dimension(mgncol_loc) :: rho_arr_l
    real(r8), dimension(mgncol_loc) :: mu_arr
    real(r8), dimension(mgncol_loc) :: asn_arr_l
    real(r8), dimension(mgncol_loc) :: accre_enhan_arr_l
    real(r8), dimension(mgncol_loc) :: relvar_arr_l

    real(r8), dimension(mgncol_loc) :: qric_arr_l
    real(r8), dimension(mgncol_loc) :: qsic_arr_l
    real(r8), dimension(mgncol_loc) :: qliq_arr_l
    real(r8), dimension(mgncol_loc) :: nliq_arr_l
    real(r8), dimension(mgncol_loc) :: qice_arr_l
    real(r8), dimension(mgncol_loc) :: nice_arr_l

    real(r8), dimension(mgncol_loc) :: pgam_arr_l
    real(r8), dimension(mgncol_loc) :: lamc_arr_l
    real(r8), dimension(mgncol_loc) :: lamr_arr_l
    real(r8), dimension(mgncol_loc) :: n0r_arr_l
    real(r8), dimension(mgncol_loc) :: lams_arr_l
    real(r8), dimension(mgncol_loc) :: n0s_arr_l

    real(r8), dimension(mgncol_loc) :: umr_arr_l
    real(r8), dimension(mgncol_loc) :: ums_arr_l
    real(r8), dimension(mgncol_loc) :: unr_arr_l
    real(r8), dimension(mgncol_loc) :: uns_arr_l

    real(r8), dimension(mgncol_loc) :: pra_arr_l
    real(r8), dimension(mgncol_loc) :: npra_arr_l
    real(r8), dimension(mgncol_loc) :: psacws_arr_l
    real(r8), dimension(mgncol_loc) :: npsacws_arr_l
    real(r8), dimension(mgncol_loc) :: pracs_arr_l
    real(r8), dimension(mgncol_loc) :: npracs_arr_l
    real(r8), dimension(mgncol_loc) :: prai_arr_l
    real(r8), dimension(mgncol_loc) :: nprai_arr_l
    real(r8), dimension(mgncol_loc) :: msacwi_arr_l
    real(r8), dimension(mgncol_loc) :: nsacwi_arr_l

    dq_cw_rain         = 0.0
    dq_cw_snow_to_snow = 0.0
    dq_ci_snow         = 0.0
    dq_sip             = 0.0
    dq_rs              = 0.0

    dn_cw_rain         = 0.0
    dn_cw_snow         = 0.0
    dn_ci_snow         = 0.0
    dn_sip             = 0.0
    dn_rs              = 0.0
    dq_cw_snow_total   = 0.0

    rho_safe = max(rho_air, 1.0e-12)
    dt_safe  = max(dt_step, 1.0e-6)

    if((qrain <= qsmall_2m .and. qsnow <= qsmall_2m) .or.                 &
       (qliq <= qsmall_2m .and. qice <= qsmall_2m)) then
       return
    endif

    microp_uniform_acc_loc = .true.

    pgam_mg_loc = 0._r8
    lamc_mg_loc = 0._r8

    if(qliq > qsmall_2m .and. nliq > 0.0) then
       qcic_mg_loc = real(qliq, r8)
       rho_mg_loc  = real(rho_safe, r8)
       ncic_mg_loc = real(max(nsmall_2m/rho_safe, nliq), r8)

       call size_dist_param_liq(                                          &
            liq_props, qcic_mg_loc, ncic_mg_loc, rho_mg_loc,              &
            pgam_mg_loc, lamc_mg_loc )
    endif

    t_arr(1)             = real(temp_k, r8)
    rho_arr_l(1)         = real(rho_safe, r8)
    mu_arr(1)            = cup_2m_air_viscosity(real(temp_k, r8))
    asn_arr_l(1)         = asnow_2m
    accre_enhan_arr_l(1) = 1._r8
    relvar_arr_l(1)      = 0._r8

    qric_arr_l(1)        = real(max(0.0, qrain), r8)
    qsic_arr_l(1)        = real(max(0.0, qsnow), r8)

    qliq_arr_l(1)        = real(max(0.0, qliq), r8)
    nliq_arr_l(1)        = real(max(0.0, nliq), r8)
    qice_arr_l(1)        = real(max(0.0, qice), r8)
    nice_arr_l(1)        = real(max(0.0, nice), r8)

    pgam_arr_l(1)        = pgam_mg_loc
    lamc_arr_l(1)        = max(lamc_mg_loc, 1.e-12_r8)

    lamr_arr_l(1)        = real(max(0.0, lamr), r8)
    n0r_arr_l(1)         = real(max(0.0, n0r), r8)
    lams_arr_l(1)        = real(max(0.0, lams), r8)
    n0s_arr_l(1)         = real(max(0.0, n0s), r8)

    umr_arr_l(1)         = real(max(0.0, umr), r8)
    ums_arr_l(1)         = real(max(0.0, ums), r8)
    unr_arr_l(1)         = real(max(0.0, unr), r8)
    uns_arr_l(1)         = real(max(0.0, uns), r8)

    pra_arr_l(1)         = 0._r8
    npra_arr_l(1)        = 0._r8
    psacws_arr_l(1)      = 0._r8
    npsacws_arr_l(1)     = 0._r8
    pracs_arr_l(1)       = 0._r8
    npracs_arr_l(1)      = 0._r8
    prai_arr_l(1)        = 0._r8
    nprai_arr_l(1)       = 0._r8
    msacwi_arr_l(1)      = 0._r8
    nsacwi_arr_l(1)      = 0._r8

    if(qrain > qsmall_2m .and. qliq > qsmall_2m) then
       call accrete_cloud_water_rain(                                      &
            microp_uniform_acc_loc, qric_arr_l, qliq_arr_l,               &
            nliq_arr_l, relvar_arr_l, accre_enhan_arr_l,                  &
            pra_arr_l, npra_arr_l, mgncol_loc )
    endif

    if(qsnow > qsmall_2m .and. qliq > qsmall_2m .and. temp_k <= 273.15 .and. &
       pgam_mg_loc > 0._r8) then

       call accrete_cloud_water_snow(                                      &
            t_arr, rho_arr_l, asn_arr_l, uns_arr_l, mu_arr,               &
            qliq_arr_l, nliq_arr_l, qsic_arr_l,                           &
            pgam_arr_l, lamc_arr_l, lams_arr_l, n0s_arr_l,                &
            psacws_arr_l, npsacws_arr_l, mgncol_loc )

       dq_cw_snow_total = real(psacws_arr_l(1)) * dt_safe

       call secondary_ice_production(                                      &
            t_arr, psacws_arr_l, msacwi_arr_l, nsacwi_arr_l, mgncol_loc )

    endif

    if(qrain > qsmall_2m .and. qsnow > qsmall_2m .and. temp_k <= 273.15) then
       call accrete_rain_snow(                                             &
            t_arr, rho_arr_l,                                             &
            umr_arr_l, ums_arr_l, unr_arr_l, uns_arr_l,                   &
            qric_arr_l, qsic_arr_l,                                       &
            lamr_arr_l, n0r_arr_l, lams_arr_l, n0s_arr_l,                 &
            pracs_arr_l, npracs_arr_l, mgncol_loc )
    endif

    if(qsnow > qsmall_2m .and. qice > qsmall_2m .and. temp_k <= 273.15) then
       call accrete_cloud_ice_snow(                                        &
            t_arr, rho_arr_l, asn_arr_l,                                  &
            qice_arr_l, nice_arr_l, qsic_arr_l,                           &
            lams_arr_l, n0s_arr_l, prai_arr_l, nprai_arr_l,               &
            mgncol_loc, accre_enhan_arr_l )
    endif

    if(use_accretion_eff_tuning_2m) then
       pra_arr_l(1)    = pra_arr_l(1)    * max(0._r8, accre_eff_rain_2m)
       npra_arr_l(1)   = npra_arr_l(1)   * max(0._r8, accre_eff_rain_2m)
       psacws_arr_l(1) = psacws_arr_l(1) * max(0._r8, accre_eff_snow_liq_2m)
       npsacws_arr_l(1)= npsacws_arr_l(1)* max(0._r8, accre_eff_snow_liq_2m)
       msacwi_arr_l(1) = msacwi_arr_l(1) * max(0._r8, accre_eff_sip_2m)
       nsacwi_arr_l(1) = nsacwi_arr_l(1) * max(0._r8, accre_eff_sip_2m)
       pracs_arr_l(1)  = pracs_arr_l(1)  * max(0._r8, accre_eff_rain_snow_2m)
       npracs_arr_l(1) = npracs_arr_l(1) * max(0._r8, accre_eff_rain_snow_2m)
       prai_arr_l(1)   = prai_arr_l(1)   * max(0._r8, accre_eff_snow_ice_2m)
       nprai_arr_l(1)  = nprai_arr_l(1)  * max(0._r8, accre_eff_snow_ice_2m)
    endif

    dq_cw_rain = max(0.0, real(pra_arr_l(1))      * dt_safe)
    dn_cw_rain = max(0.0, real(npra_arr_l(1))     * dt_safe)

    dq_cw_snow_to_snow = max(0.0, real(psacws_arr_l(1))  * dt_safe)
    dn_cw_snow         = max(0.0, real(npsacws_arr_l(1)) * dt_safe)

    dq_sip = max(0.0, real(msacwi_arr_l(1)) * dt_safe)
    dn_sip = max(0.0, real(nsacwi_arr_l(1)) * dt_safe)

    dq_cw_snow_total = max(dq_cw_snow_total, dq_cw_snow_to_snow + dq_sip)

    dq_ci_snow = max(0.0, real(prai_arr_l(1))  * dt_safe)
    dn_ci_snow = max(0.0, real(nprai_arr_l(1)) * dt_safe)

    dq_rs = max(0.0, real(pracs_arr_l(1))  * dt_safe)
    dn_rs = max(0.0, real(npracs_arr_l(1)) * dt_safe)

    ! Cap liquid sinks consistently in mass and number.
    liq_sink_mass = dq_cw_rain + dq_cw_snow_total
    liq_sink_num  = dn_cw_rain + dn_cw_snow

    scale_sink = 1.0
    if(liq_sink_mass > 1.0e-30) scale_sink = min(scale_sink, qliq/liq_sink_mass)
    if(liq_sink_num  > 1.0e-30) scale_sink = min(scale_sink, nliq/liq_sink_num)
    scale_sink = max(0.0, min(1.0, scale_sink))

    dq_cw_rain         = dq_cw_rain         * scale_sink
    dn_cw_rain         = dn_cw_rain         * scale_sink
    dq_cw_snow_total   = dq_cw_snow_total   * scale_sink
    dq_cw_snow_to_snow = dq_cw_snow_to_snow * scale_sink
    dn_cw_snow         = dn_cw_snow         * scale_sink
    dq_sip             = dq_sip             * scale_sink
    dn_sip             = dn_sip             * scale_sink

    ! Cap cloud-ice sink.
    ice_sink_mass = dq_ci_snow
    ice_sink_num  = dn_ci_snow

    scale_sink = 1.0
    if(ice_sink_mass > 1.0e-30) scale_sink = min(scale_sink, qice/ice_sink_mass)
    if(ice_sink_num  > 1.0e-30) scale_sink = min(scale_sink, nice/ice_sink_num)
    scale_sink = max(0.0, min(1.0, scale_sink))

    dq_ci_snow = dq_ci_snow * scale_sink
    dn_ci_snow = dn_ci_snow * scale_sink

    ! Rain-to-snow is diagnostic only in this fixed-precip corrector.
    scale_sink = 1.0
    if(dq_rs > 1.0e-30) scale_sink = min(scale_sink, max(qrain,0.0)/dq_rs)
    if(dn_rs > 1.0e-30) scale_sink = min(scale_sink, max(nrain,0.0)/dn_rs)
    scale_sink = max(0.0, min(1.0, scale_sink))

    dq_rs = dq_rs * scale_sink
    dn_rs = dn_rs * scale_sink

    qliq = max(0.0, qliq - dq_cw_rain - dq_cw_snow_total)
    nliq = max(0.0, nliq - dn_cw_rain - dn_cw_snow)

    qice = max(0.0, qice - dq_ci_snow + dq_sip)
    nice = max(0.0, nice - dn_ci_snow + dn_sip)

    if(qliq > qsmall_2m) then
       nliq = max(nliq, ncond_floor_cm3_liq * 1.0e6 / rho_safe)
    endif

    if(qice > qsmall_2m) then
       nice = max(nice, ncond_floor_cm3_ice * 1.0e6 / rho_safe)
    endif

  END SUBROUTINE cup_2m_apply_fixed_precip_accretion


  SUBROUTINE cup_2m_precip_number_from_auto_mass( dq_auto, dn_cloud_sink, &
                                                           rho_air, rho_hyd, n0_src, &
                                                           lam_bnd, dn_precip_src )

    implicit none

    real, intent(in) :: dq_auto          ! kg kg-1 generated over the layer step
    real, intent(in) :: dn_cloud_sink    ! # kg-1 removed from cloud category
    real, intent(in) :: rho_air          ! kg m-3
    real(r8), intent(in) :: rho_hyd      ! kg m-3
    real(r8), intent(in) :: n0_src       ! m-4
    real(r8), dimension(2), intent(in) :: lam_bnd
    real, intent(out) :: dn_precip_src   ! # kg-1 precip number source

    real(r8) :: dq_r8
    real(r8) :: dn_sink_r8
    real(r8) :: rho_r8
    real(r8) :: lwc_src
    real(r8) :: lam_src
    real(r8) :: nvol_src
    real(r8) :: dn_src_mix

    dn_precip_src = 0.0

    dq_r8      = max(0._r8, real(dq_auto, r8))
    dn_sink_r8 = max(0._r8, real(dn_cloud_sink, r8))
    rho_r8     = max(1.e-12_r8, real(rho_air, r8))

    if(dq_r8 <= 0._r8 .or. n0_src <= 0._r8) return

    ! Mass concentration generated by autoconversion over this layer step.
    lwc_src = rho_r8 * dq_r8
    if(lwc_src <= 0._r8) return

    ! Exponential PSD source closure consistent with cup_2m_mp_from_flux:
    ! q*rho = pi*rho_hyd*N0/lambda^4, N = N0/lambda.
    lam_src = (pi_r8 * rho_hyd * n0_src / lwc_src)**0.25_r8
    lam_src = max(lam_bnd(1), min(lam_bnd(2), lam_src))

    nvol_src   = n0_src / max(lam_src, 1.e-30_r8)
    dn_src_mix = nvol_src / rho_r8

    ! A precip source cannot exceed the cloud particles removed.
    ! This is mainly protective; for source PSDs the cap is usually inactive.
    if(dn_sink_r8 > 0._r8) dn_src_mix = min(dn_src_mix, dn_sink_r8)

    dn_precip_src = real(max(0._r8, dn_src_mix))

  END SUBROUTINE cup_2m_precip_number_from_auto_mass


  SUBROUTINE cup_2m_apply_mg_rain_self_collection( qrain, nrain, rho_air,  &
                                                               dt_step, rain_flux_n,   &
                                                               dn_self_mix )

    !-----------------------------------------------------------------------
    ! Optional MG Beheng-style rain self-collection for the diagnostic falling
    ! rain flux.  This is a number-only sink:
    !
    !   rain + rain -> larger rain
    !
    ! It conserves rain mass flux, reduces rain number flux, and lets the
    ! subsequent MP reconstruction produce larger/faster rain collectors.
    !-----------------------------------------------------------------------

    implicit none

    real, intent(in)    :: qrain
    real, intent(in)    :: nrain
    real, intent(in)    :: rho_air
    real, intent(in)    :: dt_step
    real, intent(inout) :: rain_flux_n
    real, intent(out)   :: dn_self_mix

    integer, parameter :: mgncol_loc = 1
    real(r8), dimension(mgncol_loc) :: rho_arr_l
    real(r8), dimension(mgncol_loc) :: qrain_arr_l
    real(r8), dimension(mgncol_loc) :: nrain_arr_l
    real(r8), dimension(mgncol_loc) :: nragg_arr_l

    real(r8) :: n_old
    real(r8) :: n_new
    real(r8) :: dt_safe_l
    real(r8) :: rate_l
    real(r8) :: dn_l
    real(r8) :: dn_max_l
    real(r8) :: flux_scale_l

    dn_self_mix = 0.0

    if(.not. use_mg_rain_self_collection_2m) return
    if(qrain <= qsmall_2m) return
    if(nrain <= 1.0) return
    if(rain_flux_n <= 0.0) return
    if(rho_air <= 1.0e-12) return
    if(dt_step <= 0.0) return

    rho_arr_l(1)   = real(max(rho_air, 1.0e-12), r8)
    qrain_arr_l(1) = real(max(qrain, 0.0), r8)
    nrain_arr_l(1) = real(max(nrain, 0.0), r8)
    nragg_arr_l(1) = 0._r8

    call self_collection_rain( rho_arr_l, qrain_arr_l, nrain_arr_l,        &
                               nragg_arr_l, mgncol_loc )

    ! MG returns a negative tendency [# kg-1 s-1].  Use an implicit update
    ! derived from dN/dt = -k N, where k = -nragg/N, to avoid zeroing rain
    ! number when qrain is very large.
    n_old     = nrain_arr_l(1)
    dt_safe_l = real(max(dt_step, 0.0), r8)

    if(n_old <= 0._r8 .or. nragg_arr_l(1) >= 0._r8) return

    rate_l = max(0._r8, -nragg_arr_l(1)/max(n_old, 1.0e-30_r8))
    rate_l = rate_l * max(0._r8, rain_selfcollect_eff_2m)

    n_new = n_old / (1._r8 + rate_l*dt_safe_l)
    n_new = max(0._r8, min(n_old, n_new))

    dn_l     = n_old - n_new
    dn_max_l = max(0._r8, rain_selfcollect_max_frac_2m) * n_old
    dn_l     = max(0._r8, min(dn_l, dn_max_l))
    n_new    = max(0._r8, n_old - dn_l)

    if(n_old > 0._r8) then
       flux_scale_l = n_new/n_old
    else
       flux_scale_l = 1._r8
    endif

    rain_flux_n = max(0.0, rain_flux_n * real(flux_scale_l))
    dn_self_mix = real(dn_l)

  END SUBROUTINE cup_2m_apply_mg_rain_self_collection


  SUBROUTINE cup_2m_apply_mg_snow_self_aggregation( temp_k, qsnow, nsnow,   &
                                                    rho_air, dt_step,       &
                                                    snow_flux_n,            &
                                                    dn_self_mix )

    !-----------------------------------------------------------------------
    ! Optional MG snow self-aggregation.  This is also a number-only sink:
    !
    !   snow + snow -> larger snow
    !
    ! It conserves snow mass flux, reduces snow number flux, and lets the
    ! subsequent MP reconstruction produce larger/faster snow collectors.
    !-----------------------------------------------------------------------

    implicit none

    real, intent(in)    :: temp_k
    real, intent(in)    :: qsnow
    real, intent(in)    :: nsnow
    real, intent(in)    :: rho_air
    real, intent(in)    :: dt_step
    real, intent(inout) :: snow_flux_n
    real, intent(out)   :: dn_self_mix

    integer, parameter :: mgncol_loc = 1
    real(r8), dimension(mgncol_loc) :: t_arr_l
    real(r8), dimension(mgncol_loc) :: rho_arr_l
    real(r8), dimension(mgncol_loc) :: asn_arr_l
    real(r8), dimension(mgncol_loc) :: qsnow_arr_l
    real(r8), dimension(mgncol_loc) :: nsnow_arr_l
    real(r8), dimension(mgncol_loc) :: nsagg_arr_l

    real(r8) :: n_old
    real(r8) :: n_new
    real(r8) :: dt_safe_l
    real(r8) :: rate_l
    real(r8) :: dn_l
    real(r8) :: dn_max_l
    real(r8) :: flux_scale_l

    dn_self_mix = 0.0

    if(.not. use_mg_snow_self_aggregation_2m) return
    if(temp_k > 273.15) return
    if(qsnow <= qsmall_2m) return
    if(nsnow <= 1.0) return
    if(snow_flux_n <= 0.0) return
    if(rho_air <= 1.0e-12) return
    if(dt_step <= 0.0) return

    t_arr_l(1)     = real(temp_k, r8)
    rho_arr_l(1)   = real(max(rho_air, 1.0e-12), r8)
    asn_arr_l(1)   = asnow_2m
    qsnow_arr_l(1) = real(max(qsnow, 0.0), r8)
    nsnow_arr_l(1) = real(max(nsnow, 0.0), r8)
    nsagg_arr_l(1) = 0._r8

    call snow_self_aggregation( t_arr_l, rho_arr_l, asn_arr_l,             &
                                rho_snow_2m, qsnow_arr_l, nsnow_arr_l,     &
                                nsagg_arr_l, mgncol_loc )

    n_old     = nsnow_arr_l(1)
    dt_safe_l = real(max(dt_step, 0.0), r8)

    if(n_old <= 0._r8 .or. nsagg_arr_l(1) >= 0._r8) return

    rate_l = max(0._r8, -nsagg_arr_l(1)/max(n_old, 1.0e-30_r8))
    rate_l = rate_l * max(0._r8, snow_selfagg_eff_2m)

    n_new = n_old / (1._r8 + rate_l*dt_safe_l)
    n_new = max(0._r8, min(n_old, n_new))

    dn_l     = n_old - n_new
    dn_max_l = max(0._r8, snow_selfagg_max_frac_2m) * n_old
    dn_l     = max(0._r8, min(dn_l, dn_max_l))
    n_new    = max(0._r8, n_old - dn_l)

    if(n_old > 0._r8) then
       flux_scale_l = n_new/n_old
    else
       flux_scale_l = 1._r8
    endif

    snow_flux_n = max(0.0, snow_flux_n * real(flux_scale_l))
    dn_self_mix = real(dn_l)

  END SUBROUTINE cup_2m_apply_mg_snow_self_aggregation


  SUBROUTINE cup_2m_mp_from_flux( mflux, nflux, rho_air, rho_hyd,          &
                                  a_vt, b_vt, lam_bnd,                    &
                                  qmix, nmix, lam, n0, um, un )

    !-----------------------------------------------------------------------
    ! Diagnose an exponential Marshall-Palmer precipitation distribution from
    ! downward mass and number fluxes.
    !
    !   N(D) = N0 exp(-lambda D)
    !   V(D) = a_vt D**b_vt
    !
    ! Inputs:
    !   mflux   kg m-2 s-1
    !   nflux   #  m-2 s-1
    !   rho_air kg m-3
    !   rho_hyd kg m-3
    !
    ! Outputs:
    !   qmix    kg kg-1
    !   nmix    # kg-1
    !   lam     m-1
    !   n0      # m-4
    !   um      mass-weighted terminal velocity, m s-1
    !   un      number-weighted terminal velocity, m s-1
    !-----------------------------------------------------------------------

    implicit none

    real(r8), intent(in)  :: mflux
    real(r8), intent(in)  :: nflux
    real(r8), intent(in)  :: rho_air
    real(r8), intent(in)  :: rho_hyd
    real(r8), intent(in)  :: a_vt
    real(r8), intent(in)  :: b_vt
    real(r8), intent(in)  :: lam_bnd(2)

    real, intent(out) :: qmix
    real, intent(out) :: nmix
    real, intent(out) :: lam
    real, intent(out) :: n0
    real, intent(out) :: um
    real, intent(out) :: un

    real(r8) :: lam_r8
    real(r8) :: n0_r8
    real(r8) :: um_r8
    real(r8) :: un_r8
    real(r8) :: vel_ratio
    real(r8) :: qmix_r8
    real(r8) :: nmix_r8
    real(r8) :: rho_safe
    real(r8) :: lam_min
    real(r8) :: lam_max

    qmix = 0.0
    nmix = 0.0
    lam  = real(lam_bnd(2))
    n0   = 0.0
    um   = 0.0
    un   = 0.0

    if(mflux <= 0._r8 .or. rho_air <= 0._r8) return

    rho_safe = max(rho_air, 1.e-12_r8)

    lam_min = min(lam_bnd(1), lam_bnd(2))
    lam_max = max(lam_bnd(1), lam_bnd(2))

    ! For an exponential distribution:
    !   un = a Gamma(b+1) / lambda**b
    !   um = a Gamma(b+4) / (Gamma(4) lambda**b)
    ! and um/un is independent of lambda.
    vel_ratio = gamma(b_vt + 4._r8) / (6._r8 * gamma(b_vt + 1._r8))

    if(nflux > 0._r8) then
       lam_r8 = (pi_r8 * rho_hyd * vel_ratio * nflux / max(mflux, 1.e-30_r8))**(1._r8/3._r8)
    else
       ! Fallback if the source has mass but no reliable precip number.
       lam_r8 = sqrt(lam_min * lam_max)
    endif

    lam_r8 = max(lam_min, min(lam_max, lam_r8))

    un_r8 = a_vt * gamma(b_vt + 1._r8) / max(lam_r8**b_vt, 1.e-30_r8)
    um_r8 = a_vt * gamma(b_vt + 4._r8) / (6._r8 * max(lam_r8**b_vt, 1.e-30_r8))

    un_r8 = max(vt_min_2m, un_r8)
    um_r8 = max(vt_min_2m, um_r8)

    qmix_r8 = mflux / (rho_safe * um_r8)

    ! Keep N0, lambda, q, and n internally consistent after bounding lambda.
    n0_r8   = qmix_r8 * rho_safe * lam_r8**4 / (pi_r8 * rho_hyd)
    nmix_r8 = n0_r8 / (rho_safe * lam_r8)

    qmix = real(max(0._r8, qmix_r8))
    nmix = real(max(0._r8, nmix_r8))
    lam  = real(max(0._r8, lam_r8))
    n0   = real(max(0._r8, n0_r8))
    um   = real(max(0._r8, um_r8))
    un   = real(max(0._r8, un_r8))

  END SUBROUTINE cup_2m_mp_from_flux


  real(r8) FUNCTION cup_2m_air_viscosity(t)

    ! Dynamic viscosity of air [kg m-1 s-1].
    ! Sutherland-style fit, adequate for the snow-riming collection efficiency.

    implicit none

    real(r8), intent(in) :: t
    real(r8) :: tsafe

    tsafe = max(180._r8, min(330._r8, t))

    cup_2m_air_viscosity = 1.496e-6_r8 * tsafe**1.5_r8 / (tsafe + 120._r8)

  END FUNCTION cup_2m_air_viscosity




end module GF2020_2M_MicrophysicsMod
