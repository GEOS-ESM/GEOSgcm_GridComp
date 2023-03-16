module WaterFluxBulkType

  use MAPL_ConstantsMod    , ONLY : r8 => MAPL_R8
  use nanMod               , only : nan
  use clm_varpar           , only : nlevsno, nlevsoi
  use clm_varcon           , only : spval
  use MAPL_ExceptionHandling
  use WaterFluxType  , only : waterflux_type
  use decompMod      , only : bounds_type

  implicit none

  ! !PUBLIC TYPES:
  type, extends(waterflux_type), public :: waterfluxbulk_type
     ! water fluxes are in units or mm/s

     real(r8), pointer :: qflx_phs_neg_col         (:)   ! col sum of negative hydraulic redistribution fluxes (mm H2O/s) [+]

     real(r8), pointer :: qflx_snowindunload_patch (:)   ! patch canopy snow wind unloading (mm H2O /s)
     real(r8), pointer :: qflx_snotempunload_patch (:)   ! patch canopy snow temp unloading (mm H2O /s) 

     real(r8), pointer :: qflx_ev_snow_patch       (:)   ! patch evaporation heat flux from snow       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_snow_col         (:)   ! col evaporation heat flux from snow         (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_patch       (:)   ! patch evaporation heat flux from soil       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_soil_col         (:)   ! col evaporation heat flux from soil         (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_patch     (:)   ! patch evaporation heat flux from soil       (mm H2O/s) [+ to atm]
     real(r8), pointer :: qflx_ev_h2osfc_col       (:)   ! col evaporation heat flux from soil         (mm H2O/s) [+ to atm]

     real(r8), pointer :: qflx_adv_col             (:,:) ! col advective flux across different soil layer interfaces [mm H2O/s] [+ downward]
     real(r8), pointer :: qflx_rootsoi_col         (:,:) ! col root and soil water exchange [mm H2O/s] [+ into root]
     real(r8), pointer :: qflx_hydr_redist_patch   (:)   ! patch hydraulic redistribution [mm H2O/s]
     real(r8), pointer :: qflx_sat_excess_surf_col (:)   ! col surface runoff due to saturated surface (mm H2O /s)
     real(r8), pointer :: qflx_infl_excess_col     (:)   ! col infiltration excess runoff (mm H2O /s)
     real(r8), pointer :: qflx_infl_excess_surf_col(:)   ! col surface runoff due to infiltration excess (mm H2O /s)
     real(r8), pointer :: qflx_h2osfc_surf_col     (:)   ! col surface water runoff (mm H2O /s)
     real(r8), pointer :: qflx_in_soil_col         (:)   ! col surface input to soil (mm/s)
     real(r8), pointer :: qflx_in_soil_limited_col (:)   ! col surface input to soil, limited by max infiltration rate (mm/s)
     real(r8), pointer :: qflx_h2osfc_drain_col    (:)   ! col bottom drainage from h2osfc (mm/s)
     real(r8), pointer :: qflx_top_soil_to_h2osfc_col(:) ! col portion of qflx_top_soil going to h2osfc, minus evaporation (mm/s)
     real(r8), pointer :: qflx_in_h2osfc_col(:)          ! col total surface input to h2osfc
     real(r8), pointer :: qflx_deficit_col         (:)   ! col water deficit to keep non-negative liquid water content (mm H2O)   
     real(r8), pointer :: qflx_snomelt_lyr_col     (:,:) ! col snow melt in each layer (mm H2O /s)
     real(r8), pointer :: qflx_drain_vr_col        (:,:) ! col liquid water losted as drainage (m /time step)

     ! ET accumulation
     real(r8), pointer :: AnnEt                    (:)   ! Annual average ET flux mmH20/s      

   contains

     procedure , public :: InitBulk

  end type waterfluxbulk_type
!  type(waterfluxbulk_type), public, target, save :: waterfluxbulk_inst

contains

!---------------------------------------------
  subroutine InitBulk(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM type for water flux bulk variables that just apply to bulk water and are needed for calling CTSM routines                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(waterfluxbulk_type), intent(inout)     :: this

    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    call this%Init(bounds)

    allocate(this%qflx_snowindunload_patch (begp:endp))              ; this%qflx_snowindunload_patch (:)   = nan
    allocate(this%qflx_snotempunload_patch (begp:endp))              ; this%qflx_snotempunload_patch (:)   = nan
    allocate(this%qflx_hydr_redist_patch   (begp:endp))              ; this%qflx_hydr_redist_patch   (:)   = nan
    allocate(this%qflx_phs_neg_col         (begc:endc))              ; this%qflx_phs_neg_col         (:)   = nan

    allocate( this%qflx_ev_snow_patch      (begp:endp))              ; this%qflx_ev_snow_patch       (:)   = nan
    allocate( this%qflx_ev_snow_col        (begc:endc))              ; this%qflx_ev_snow_col         (:)   = nan
    allocate( this%qflx_ev_soil_patch      (begp:endp))              ; this%qflx_ev_soil_patch       (:)   = nan
    allocate( this%qflx_ev_soil_col        (begc:endc))              ; this%qflx_ev_soil_col         (:)   = nan
    allocate( this%qflx_ev_h2osfc_patch    (begp:endp))              ; this%qflx_ev_h2osfc_patch     (:)   = nan
    allocate( this%qflx_ev_h2osfc_col      (begc:endc))              ; this%qflx_ev_h2osfc_col       (:)   = nan

    allocate(this%qflx_drain_vr_col      (begc:endc,1:nlevsoi))      ; this%qflx_drain_vr_col        (:,:) = nan
    allocate(this%qflx_adv_col             (begc:endc,0:nlevsoi))    ; this%qflx_adv_col             (:,:) = nan
    allocate(this%qflx_rootsoi_col         (begc:endc,1:nlevsoi))    ; this%qflx_rootsoi_col         (:,:) = nan
    allocate(this%qflx_sat_excess_surf_col (begc:endc))              ; this%qflx_sat_excess_surf_col (:)   = nan
    allocate(this%qflx_infl_excess_col     (begc:endc))              ; this%qflx_infl_excess_col     (:)   = nan
    allocate(this%qflx_in_soil_col         (begc:endc))              ; this%qflx_in_soil_col         (:)   = nan
    allocate(this%qflx_in_soil_limited_col (begc:endc))              ; this%qflx_in_soil_limited_col (:)   = nan
    allocate(this%qflx_h2osfc_drain_col    (begc:endc))              ; this%qflx_h2osfc_drain_col    (:)   = nan
    allocate(this%qflx_top_soil_to_h2osfc_col(begc:endc))            ; this%qflx_top_soil_to_h2osfc_col(:) = nan
    allocate(this%qflx_in_h2osfc_col       (begc:endc))              ; this%qflx_in_h2osfc_col(:)          = nan
    allocate(this%qflx_infl_excess_surf_col(begc:endc))              ; this%qflx_infl_excess_surf_col(:)   = nan
    allocate(this%qflx_h2osfc_surf_col     (begc:endc))              ; this%qflx_h2osfc_surf_col     (:)   = nan
    allocate(this%qflx_snomelt_lyr_col     (begc:endc,-nlevsno+1:0)) ; this%qflx_snomelt_lyr_col     (:,:) = nan
    allocate(this%qflx_deficit_col         (begc:endc))              ; this%qflx_deficit_col         (:)   = nan
    allocate(this%AnnET                    (begc:endc))              ; this%AnnET                    (:)   = nan

  end subroutine InitBulk
end module WaterFluxBulkType
