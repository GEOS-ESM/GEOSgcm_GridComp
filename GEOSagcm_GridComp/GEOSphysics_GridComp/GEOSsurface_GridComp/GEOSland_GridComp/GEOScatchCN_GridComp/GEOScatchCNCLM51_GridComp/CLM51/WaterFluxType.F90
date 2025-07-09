module WaterFluxType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use clm_varpar       , only : nlevsno, nlevsoi
  use clm_varcon       , only : spval
  use LandunitType     , only : lun
  use ColumnType       , only : col
  use netcdf
  use MAPL_ExceptionHandling
  use decompMod        , only : bounds_type
  use AnnualFluxDribbler, only : annual_flux_dribbler_type

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

 !
  type, public :: waterflux_type

     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_through_snow_patch  (:)   ! patch canopy throughfall of snow (mm H2O/s)
     real(r8), pointer :: qflx_through_liq_patch  (:)    ! patch canopy throughfal of liquid (rain+irrigation) (mm H2O/s)
     real(r8), pointer :: qflx_intercepted_snow_patch(:) ! patch canopy interception of snow (mm H2O/s)
     real(r8), pointer :: qflx_intercepted_liq_patch(:)  ! patch canopy interception of liquid (rain+irrigation) (mm H2O/s)
     real(r8), pointer :: qflx_snocanfall_patch(:)       ! patch rate of excess canopy snow falling off canopy (mm H2O/s)
     real(r8), pointer :: qflx_liqcanfall_patch(:)       ! patch rate of excess canopy liquid falling off canopy (mm H2O/s)
     real(r8), pointer :: qflx_snow_unload_patch(:)      ! patch rate of canopy snow unloading (mm H2O/s)
     real(r8), pointer :: qflx_liq_grnd_col        (:)   ! col liquid (rain+irrigation) on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_snow_grnd_col       (:)   ! col snow on ground after interception (mm H2O/s) [+]
     real(r8), pointer :: qflx_rain_plus_snomelt_col(:)  ! col rain plus snow melt falling on the soil (mm/s)
     real(r8), pointer :: qflx_solidevap_from_top_layer_patch(:) ! patch rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]
     real(r8), pointer :: qflx_solidevap_from_top_layer_col(:)   ! col rate of ice evaporated from top soil or snow layer (sublimation) (mm H2O /s) [+]
     real(r8), pointer :: qflx_evap_soi_patch      (:)   ! patch soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_soi_col        (:)   ! col soil evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_patch      (:)   ! patch vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_veg_col        (:)   ! col vegetation evaporation (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_patch      (:)   ! patch evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_can_col        (:)   ! col evaporation from leaves and stems (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_evap_tot_patch      (:)   ! patch pft_qflx_evap_soi + pft_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_evap_tot_col        (:)   ! col col_qflx_evap_soi + col_qflx_evap_veg + qflx_tran_veg
     real(r8), pointer :: qflx_liqevap_from_top_layer_patch(:) ! patch rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]
     real(r8), pointer :: qflx_liqevap_from_top_layer_col(:)   ! col rate of liquid water evaporated from top soil or snow layer (mm H2O/s) [+]

     ! In the snow capping parametrization excess mass above h2osno_max is removed.  A breakdown of mass into liquid 
     ! and solid fluxes is done, these are represented by qflx_snwcp_liq_col and qflx_snwcp_ice_col. 
     real(r8), pointer :: qflx_snwcp_liq_col       (:)   ! col excess liquid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_ice_col       (:)   ! col excess solid h2o due to snow capping (outgoing) (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_liq_col(:) ! col excess liquid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_snwcp_discarded_ice_col(:) ! col excess solid h2o due to snow capping, which we simply discard in order to reset the snow pack (mm H2O /s)
     real(r8), pointer :: qflx_glcice_col(:)              ! col net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC; only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_frz_col (:)         ! col ice growth (positive definite) (mm H2O/s); only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_melt_col(:)         ! col ice melt (positive definite) (mm H2O/s); only valid inside the do_smb_c filter
     real(r8), pointer :: qflx_glcice_dyn_water_flux_col(:) ! col water flux needed for balance check due to glc_dyn_runoff_routing (mm H2O/s) (positive means addition of water to the system); valid for all columns

     real(r8), pointer :: qflx_tran_veg_patch      (:)   ! patch vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_tran_veg_col        (:)   ! col vegetation transpiration (mm H2O/s) (+ = to atm)
     real(r8), pointer :: qflx_soliddew_to_top_layer_patch(:) ! patch rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+]
     real(r8), pointer :: qflx_soliddew_to_top_layer_col(:)   ! col rate of solid water deposited on top soil or snow layer (frost) (mm H2O /s) [+] (+ = to atm); usually eflx_bot >= 0)
     real(r8), pointer :: qflx_liqdew_to_top_layer_patch(:)   ! patch rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]
     real(r8), pointer :: qflx_liqdew_to_top_layer_col(:)     ! col rate of liquid water deposited on top soil or snow layer (dew) (mm H2O /s) [+]

     real(r8), pointer :: qflx_infl_col            (:)   ! col infiltration (mm H2O /s)
     real(r8), pointer :: qflx_surf_col            (:)   ! col total surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_col           (:)   ! col sub-surface runoff (mm H2O /s)
     real(r8), pointer :: qflx_drain_perched_col   (:)   ! col sub-surface runoff from perched wt (mm H2O /s)                                                                                                      
     real(r8), pointer :: qflx_top_soil_col        (:)   ! col net water input into soil from top (mm/s)
     real(r8), pointer :: qflx_floodc_col          (:)   ! col flood water flux at column level
     real(r8), pointer :: qflx_sl_top_soil_col     (:)   ! col liquid water + ice from layer above soil to top soil layer or sent to qflx_qrgwl (mm H2O/s)
     real(r8), pointer :: qflx_snomelt_col         (:)   ! col snow melt (mm H2O /s)
     real(r8), pointer :: qflx_qrgwl_col           (:)   ! col qflx_surf at glaciers, wetlands, lakes
     real(r8), pointer :: qflx_runoff_col          (:)   ! col total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_r_col        (:)   ! col Rural total runoff (qflx_drain+qflx_surf+qflx_qrgwl) (mm H2O /s)
     real(r8), pointer :: qflx_runoff_u_col        (:)   ! col urban total runoff (qflx_drain+qflx_surf) (mm H2O /s) 
     real(r8), pointer :: qflx_rsub_sat_col        (:)   ! col soil saturation excess [mm/s]
     real(r8), pointer :: qflx_snofrz_lyr_col      (:,:) ! col snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
     real(r8), pointer :: qflx_snofrz_col          (:)   ! col column-integrated snow freezing rate (positive definite) (col) [kg m-2 s-1]
     real(r8), pointer :: qflx_snow_drain_col      (:)   ! col drainage from snow pack
     real(r8), pointer :: qflx_ice_runoff_snwcp_col(:)   ! col solid runoff from snow capping (mm H2O /s)
     real(r8), pointer :: qflx_ice_runoff_xs_col   (:)   ! col solid runoff from excess ice in soil (mm H2O /s)

     real(r8), pointer :: qflx_h2osfc_to_ice_col   (:)   ! col conversion of h2osfc to ice
     real(r8), pointer :: qflx_snow_h2osfc_col     (:)   ! col snow falling on surface water
     real(r8), pointer :: qflx_too_small_h2osfc_to_soil_col(:) ! col h2osfc transferred to soil if h2osfc is below some threshold (mm H2O /s)
     real(r8), pointer :: qflx_snow_percolation_col(:,:) ! col liquid percolation out of the bottom of snow layer j (mm H2O /s)

     ! Dynamic land cover change
     real(r8), pointer :: qflx_liq_dynbal_grc      (:)   ! grc liq dynamic land cover change conversion runoff flux
     real(r8), pointer :: qflx_ice_dynbal_grc      (:)   ! grc ice dynamic land cover change conversion runoff flux

     real(r8), pointer :: qflx_sfc_irrig_col        (:)   ! col surface irrigation flux (mm H2O/s) [+]             
     real(r8), pointer :: qflx_gw_uncon_irrig_col   (:)   ! col unconfined groundwater irrigation flux (mm H2O/s)
     real(r8), pointer :: qflx_gw_uncon_irrig_lyr_col(:,:) ! col unconfined groundwater irrigation flux, separated by layer (mm H2O/s)
     real(r8), pointer :: qflx_gw_con_irrig_col     (:)   ! col confined groundwater irrigation flux (mm H2O/s)
     real(r8), pointer :: qflx_irrig_drip_patch     (:)   ! patch drip irrigation
     real(r8), pointer :: qflx_irrig_sprinkler_patch(:)   ! patch sprinkler irrigation

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: qflx_liq_dynbal_dribbler
     type(annual_flux_dribbler_type) :: qflx_ice_dynbal_dribbler

   contains

     procedure, public :: Init

  end type waterflux_type
  !type(waterflux_type), public, target, save :: waterflux_inst

contains

!---------------------------------------------
  subroutine Init(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM type for water flux variables that just apply to bulk water and are needed for calling CTSM routines                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !

  ! !USES:
    use landunit_varcon, only : istsoil, istcrop                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(waterflux_type), intent(inout)         :: this

    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: c, l 
    !--------------------
    
    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
    begg = bounds%begg ; endg = bounds%endg

    allocate(this%qflx_through_liq_patch(begp:endp))
    allocate(this%qflx_through_snow_patch(begp:endp))
    allocate(this%qflx_liqcanfall_patch(begp:endp))
    allocate(this%qflx_snocanfall_patch(begp:endp))
    allocate(this%qflx_snow_unload_patch(begp:endp))
    allocate(this%qflx_top_soil_col(begc:endc))
    allocate(this%qflx_infl_col(begc:endc))
    allocate(this%qflx_surf_col(begc:endc))
    allocate(this%qflx_qrgwl_col(begc:endc))
    allocate(this%qflx_drain_col(begc:endc))
    allocate(this%qflx_drain_perched_col(begc:endc))
    allocate(this%qflx_liq_dynbal_grc(begg:endg))
    allocate(this%qflx_ice_dynbal_grc(begg:endg))
    allocate(this%qflx_runoff_col(begc:endc))
    allocate(this%qflx_runoff_u_col(begc:endc))
    allocate(this%qflx_runoff_r_col(begc:endc))
    allocate(this%qflx_snomelt_col(begc:endc))
    allocate(this%qflx_snofrz_col(begc:endc))
    allocate(this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0))
    allocate(this%qflx_snow_drain_col(begc:endc))
    allocate(this%qflx_evap_soi_patch(begp:endp))
    allocate(this%qflx_evap_can_patch(begp:endp))
    allocate(this%qflx_tran_veg_patch(begp:endp))
    allocate(this%qflx_snwcp_liq_col(begc:endc))
    allocate(this%qflx_snwcp_ice_col(begc:endc))
    allocate(this%qflx_glcice_col(begc:endc))
    allocate(this%qflx_glcice_frz_col(begc:endc))
    allocate(this%qflx_glcice_melt_col(begc:endc))
    allocate(this%qflx_liq_grnd_col(begc:endc))
    allocate(this%qflx_snow_grnd_col(begc:endc))
    allocate(this%qflx_liqevap_from_top_layer_patch(begp:endp))
    allocate(this%qflx_evap_veg_patch(begp:endp))
    allocate(this%qflx_evap_tot_patch(begp:endp))
    allocate(this%qflx_liqdew_to_top_layer_patch(begp:endp))
    allocate(this%qflx_solidevap_from_top_layer_patch(begp:endp))
    allocate(this%qflx_soliddew_to_top_layer_patch(begp:endp))
    allocate(this%qflx_rsub_sat_col(begc:endc))
    allocate(this%qflx_h2osfc_to_ice_col(begc:endc))
    allocate(this%qflx_sfc_irrig_col(begc:endc))
    allocate(this%qflx_gw_uncon_irrig_col(begc:endc))
    allocate(this%qflx_gw_con_irrig_col(begc:endc))
    allocate(this%qflx_irrig_drip_patch(begp:endp))
    allocate(this%qflx_irrig_sprinkler_patch(begp:endp))

    allocate(this%qflx_liqevap_from_top_layer_col(begc:endc))
    allocate(this%qflx_liqdew_to_top_layer_col(begc:endc))
    allocate(this%qflx_soliddew_to_top_layer_col(begc:endc))
    allocate(this%qflx_ice_runoff_xs_col(begc:endc))
    allocate(this%qflx_glcice_dyn_water_flux_col(begc:endc))
    allocate(this%qflx_gw_uncon_irrig_lyr_col(begc:endc,1:nlevsoi))

    this%qflx_through_liq_patch(begp:endp) = spval
    this%qflx_through_snow_patch(begp:endp) = spval
    this%qflx_liqcanfall_patch(begp:endp) = spval
    this%qflx_snocanfall_patch(begp:endp) = spval
    this%qflx_snow_unload_patch(begp:endp) = spval
    this%qflx_top_soil_col(begc:endc) = spval
    this%qflx_infl_col(begc:endc) = spval
    this%qflx_surf_col(begc:endc) = spval
    this%qflx_qrgwl_col(begc:endc) = spval
    this%qflx_drain_col(begc:endc) = spval
    this%qflx_drain_perched_col(begc:endc) = spval
    this%qflx_liq_dynbal_grc(begg:endg) = spval
    this%qflx_ice_dynbal_grc(begg:endg) = spval
    this%qflx_runoff_col(begc:endc) = spval
    this%qflx_runoff_u_col(begc:endc) = spval
    this%qflx_runoff_r_col(begc:endc) = spval
    this%qflx_snomelt_col(begc:endc) = spval
    this%qflx_snofrz_col(begc:endc) = spval
    this%qflx_snofrz_lyr_col(begc:endc,-nlevsno+1:0) = spval
    this%qflx_snow_drain_col(begc:endc) = spval
    this%qflx_evap_soi_patch(begp:endp) = spval
    this%qflx_evap_can_patch(begp:endp) = spval
    this%qflx_tran_veg_patch(begp:endp) = spval
    this%qflx_snwcp_liq_col(begc:endc) = spval
    this%qflx_snwcp_ice_col(begc:endc) = spval
    this%qflx_glcice_col(begc:endc) = spval
    this%qflx_glcice_frz_col(begc:endc) = spval
    this%qflx_glcice_melt_col(begc:endc) = spval
    this%qflx_liq_grnd_col(begc:endc) = spval
    this%qflx_snow_grnd_col(begc:endc) = spval
    this%qflx_liqevap_from_top_layer_patch(begp:endp) = spval
    this%qflx_evap_veg_patch(begp:endp) = spval
    this%qflx_evap_tot_patch(begp:endp) = spval
    this%qflx_liqdew_to_top_layer_patch(begp:endp) = spval
    this%qflx_solidevap_from_top_layer_patch(begp:endp) = spval
    this%qflx_soliddew_to_top_layer_patch(begp:endp) = spval
    this%qflx_rsub_sat_col(begc:endc) = spval
    this%qflx_h2osfc_to_ice_col(begc:endc) = spval
    this%qflx_sfc_irrig_col(begc:endc) = spval
    this%qflx_gw_uncon_irrig_col(begc:endc) = spval
    this%qflx_gw_con_irrig_col(begc:endc) = spval
    this%qflx_irrig_drip_patch(begp:endp) = spval
    this%qflx_irrig_sprinkler_patch(begp:endp) = spval

    ! assign cold start values for variables where it is needed

    this%qflx_snocanfall_patch(bounds%begp:bounds%endp)       = 0.0_r8
    this%qflx_liqcanfall_patch(bounds%begp:bounds%endp)       = 0.0_r8
    this%qflx_snow_unload_patch(bounds%begp:bounds%endp)      = 0.0_r8

    this%qflx_liqevap_from_top_layer_patch(bounds%begp:bounds%endp) = 0.0_r8
    this%qflx_liqdew_to_top_layer_patch(bounds%begp:bounds%endp)    = 0.0_r8
    this%qflx_soliddew_to_top_layer_patch (bounds%begp:bounds%endp) = 0.0_r8

    this%qflx_sfc_irrig_col (bounds%begc:bounds%endc)         = 0.0_r8
    this%qflx_gw_uncon_irrig_col (bounds%begc:bounds%endc)    = 0.0_r8
    this%qflx_gw_uncon_irrig_lyr_col(bounds%begc:bounds%endc,:) = 0.0_r8
    this%qflx_gw_con_irrig_col (bounds%begc:bounds%endc)      = 0.0_r8
    this%qflx_irrig_drip_patch (bounds%begp:bounds%endp)      = 0.0_r8
    this%qflx_irrig_sprinkler_patch (bounds%begp:bounds%endp) = 0.0_r8

    this%qflx_liqevap_from_top_layer_col(bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_liqdew_to_top_layer_col(bounds%begc:bounds%endc)    = 0.0_r8
    this%qflx_soliddew_to_top_layer_col (bounds%begc:bounds%endc) = 0.0_r8
    this%qflx_snow_drain_col(bounds%begc:bounds%endc)  = 0._r8
    this%qflx_ice_runoff_xs_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_glcice_dyn_water_flux_col(bounds%begc:bounds%endc) = 0._r8
    this%qflx_tran_veg_patch(bounds%begp:bounds%endp) = 0._r8
    this%qflx_evap_veg_patch(bounds%begp:bounds%endp) = 0._r8

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%qflx_drain_col(c) = 0._r8
          this%qflx_surf_col(c)  = 0._r8
       end if
    end do

  end subroutine Init

end module WaterFluxType
