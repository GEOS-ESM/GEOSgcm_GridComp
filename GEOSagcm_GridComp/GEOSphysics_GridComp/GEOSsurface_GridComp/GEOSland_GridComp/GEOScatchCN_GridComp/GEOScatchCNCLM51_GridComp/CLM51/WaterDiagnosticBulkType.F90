module WaterDiagnosticBulkType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use clm_varpar       , only : nlevgrnd, nlevsno, nlevcan
  use clm_varcon       , only : spval
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use WaterDiagnosticType, only : waterdiagnostic_type

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  !
  type, extends(waterdiagnostic_type), public :: waterdiagnosticbulk_type

     real(r8), pointer :: h2osno_total_col       (:)   ! col total snow water (mm H2O)
     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
     real(r8), pointer :: snow_5day_col          (:)   ! col snow height 5 day avg
     real(r8), pointer :: snowdp_col             (:)   ! col area-averaged snow height (m)
     real(r8), pointer :: snow_layer_unity_col   (:,:) ! value 1 for each snow layer, used for history diagnostics
     real(r8), pointer :: bw_col                 (:,:) ! col partial density of water in the snow pack (ice + liquid) [kg/m3] 

     real(r8), pointer :: h2osoi_liq_tot_col     (:)   ! vertically summed col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_tot_col     (:)   ! vertically summed col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: air_vol_col            (:,:) ! col air filled porosity
     real(r8), pointer :: h2osoi_liqvol_col      (:,:) ! col volumetric liquid water content (v/v)
     real(r8), pointer :: swe_old_col            (:,:) ! col initial snow water

     real(r8), pointer :: snw_rds_col            (:,:) ! col snow grain radius (col,lyr)    [m^-6, microns]
     real(r8), pointer :: snw_rds_top_col        (:)   ! col snow grain radius (top layer)  [m^-6, microns]
     real(r8), pointer :: h2osno_top_col         (:)   ! col top-layer mass of snow  [kg]
     real(r8), pointer :: sno_liq_top_col        (:)   ! col snow liquid water fraction (mass), top layer  [fraction]

     real(r8), pointer :: iwue_ln_patch          (:)   ! patch intrinsic water use efficiency near local noon (umolCO2/molH2O)
     real(r8), pointer :: vpd_ref2m_patch        (:)   ! patch 2 m height surface vapor pressure deficit (Pa)
     real(r8), pointer :: rh_ref2m_patch         (:)   ! patch 2 m height surface relative humidity (%)
     real(r8), pointer :: rh_ref2m_r_patch       (:)   ! patch 2 m height surface relative humidity - rural (%)
     real(r8), pointer :: rh_ref2m_u_patch       (:)   ! patch 2 m height surface relative humidity - urban (%)
     real(r8), pointer :: rh_af_patch            (:)   ! patch fractional humidity of canopy air (dimensionless) ! private
     real(r8), pointer :: rh10_af_patch          (:)   ! 10-day mean patch fractional humidity of canopy air (dimensionless)
     real(r8), pointer :: dqgdT_col              (:)   ! col d(qg)/dT

     ! Fractions
     real(r8), pointer :: frac_sno_col           (:)   ! col fraction of ground covered by snow (0 to 1)
     real(r8), pointer :: frac_sno_eff_col       (:)   ! col fraction of ground covered by snow (0 to 1) (note: this can be 1 even if there is no snow, but should be ignored in the no-snow case)
     real(r8), pointer :: frac_iceold_col        (:,:) ! col fraction of ice relative to the tot water (new) (-nlevsno+1:nlevgrnd) 
     real(r8), pointer :: frac_h2osfc_col        (:)   ! col fractional area with surface water greater than zero
     real(r8), pointer :: frac_h2osfc_nosnow_col (:)   ! col fractional area with surface water greater than zero (if no snow present)
     real(r8), pointer :: wf_col                 (:)   ! col soil water as frac. of whc for top 0.05 m (0-1) 
     real(r8), pointer :: wf2_col                (:)   ! col soil water as frac. of whc for top 0.17 m (0-1) 
     real(r8), pointer :: fwet_patch             (:)   ! patch canopy fraction that is wet (0 to 1)
     real(r8), pointer :: fcansno_patch          (:)   ! patch canopy fraction that is snow covered (0 to 1)
     real(r8), pointer :: fdry_patch             (:)   ! patch canopy fraction of foliage that is green and dry [-] (new)

     ! Summed fluxes
     real(r8), pointer :: qflx_prec_intr_patch   (:)   ! patch interception of precipitation (mm H2O/s)
     real(r8), pointer :: qflx_prec_grnd_col     (:)   ! col water onto ground including canopy runoff (mm H2O/s)

   contains

     procedure, public :: InitBulk

end type waterdiagnosticbulk_type
!type(waterdiagnosticbulk_type), public, target, save :: waterdiagnosticbulk_inst

contains

!-----------------------------------------------
  subroutine InitBulk(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM type for water diagnostic variables that just apply to bulk water and are needed for calling CTSM routines                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT
    type(bounds_type), intent(in) :: bounds
    class(waterdiagnosticbulk_type), intent(inout) :: this
   
    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !----------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
    begl = bounds%begl ; endl = bounds%endl
    begg = bounds%begg ; endg = bounds%endg

    call this%Init(bounds)    

    allocate(this%h2osno_total_col       (begc:endc))                     ; this%h2osno_total_col       (:)   = nan
    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%snow_5day_col          (begc:endc))                     ; this%snow_5day_col          (:)   = nan
    allocate(this%snowdp_col             (begc:endc))                     ; this%snowdp_col             (:)   = nan
    allocate(this%snow_layer_unity_col   (begc:endc,-nlevsno+1:0))        ; this%snow_layer_unity_col   (:,:) = nan
    allocate(this%bw_col                 (begc:endc,-nlevsno+1:0))        ; this%bw_col                 (:,:) = nan
    allocate(this%air_vol_col            (begc:endc, 1:nlevgrnd))         ; this%air_vol_col            (:,:) = nan
    allocate(this%h2osoi_liqvol_col      (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liqvol_col      (:,:) = nan
    allocate(this%h2osoi_ice_tot_col     (begc:endc))                     ; this%h2osoi_ice_tot_col     (:)   = nan
    allocate(this%h2osoi_liq_tot_col     (begc:endc))                     ; this%h2osoi_liq_tot_col     (:)   = nan
    allocate(this%swe_old_col            (begc:endc,-nlevsno+1:0))        ; this%swe_old_col            (:,:) = nan

    allocate(this%snw_rds_col            (begc:endc,-nlevsno+1:0))        ; this%snw_rds_col            (:,:) = nan
    allocate(this%snw_rds_top_col        (begc:endc))                     ; this%snw_rds_top_col        (:)   = nan
    allocate(this%h2osno_top_col         (begc:endc))                     ; this%h2osno_top_col         (:)   = nan
    allocate(this%sno_liq_top_col        (begc:endc))                     ; this%sno_liq_top_col        (:)   = nan

    allocate(this%dqgdT_col              (begc:endc))                     ; this%dqgdT_col              (:)   = nan
    allocate(this%iwue_ln_patch          (begp:endp))                     ; this%iwue_ln_patch          (:)   = nan
    allocate(this%vpd_ref2m_patch        (begp:endp))                     ; this%vpd_ref2m_patch        (:)   = nan
    allocate(this%rh_ref2m_patch         (begp:endp))                     ; this%rh_ref2m_patch         (:)   = nan
    allocate(this%rh_ref2m_u_patch       (begp:endp))                     ; this%rh_ref2m_u_patch       (:)   = nan
    allocate(this%rh_ref2m_r_patch       (begp:endp))                     ; this%rh_ref2m_r_patch       (:)   = nan
    allocate(this%rh_af_patch            (begp:endp))                     ; this%rh_af_patch            (:)   = nan
    allocate(this%rh10_af_patch          (begp:endp))                     ; this%rh10_af_patch          (:)   = spval

    allocate(this%frac_sno_col           (begc:endc))                     ; this%frac_sno_col           (:)   = nan
    allocate(this%frac_sno_eff_col       (begc:endc))                     ; this%frac_sno_eff_col       (:)   = nan
    allocate(this%frac_iceold_col        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%frac_iceold_col        (:,:) = nan
    allocate(this%frac_h2osfc_col        (begc:endc))                     ; this%frac_h2osfc_col        (:)   = nan
    allocate(this%frac_h2osfc_nosnow_col (begc:endc))                     ; this%frac_h2osfc_nosnow_col        (:)   = nan
    allocate(this%wf_col                 (begc:endc))                     ; this%wf_col                 (:)   = nan
    allocate(this%wf2_col                (begc:endc))                     ; this%wf2_col                (:)   = nan
    allocate(this%fwet_patch             (begp:endp))                     ; this%fwet_patch             (:)   = nan
    allocate(this%fcansno_patch          (begp:endp))                     ; this%fcansno_patch          (:)   = nan
    allocate(this%fdry_patch             (begp:endp))                     ; this%fdry_patch             (:)   = nan
    allocate(this%qflx_prec_intr_patch   (begp:endp))                     ; this%qflx_prec_intr_patch   (:)   = nan
    allocate(this%qflx_prec_grnd_col     (begc:endc))                     ; this%qflx_prec_grnd_col     (:)   = nan

 end subroutine InitBulk

end module WaterDiagnosticBulkType
