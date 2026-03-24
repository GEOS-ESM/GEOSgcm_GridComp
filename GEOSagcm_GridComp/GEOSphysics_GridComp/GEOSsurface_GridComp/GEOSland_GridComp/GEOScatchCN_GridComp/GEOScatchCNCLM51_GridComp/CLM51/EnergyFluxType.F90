module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! Energy flux data structure
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use nanMod         , only : nan
  use clm_varcon     , only : spval
  use clm_varctl     , only : use_biomass_heat_storage
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use PatchType      , only : patch
  use clm_varpar   , only: nlevgrnd

 !
  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:

  !
  type, public :: energyflux_type

     ! Fluxes
     real(r8), pointer :: eflx_sh_stem_patch      (:)   ! patch sensible heat flux from stem (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_h2osfc_to_snow_col (:)   ! col snow melt to h2osfc heat flux (W/m**2)
     real(r8), pointer :: eflx_sh_grnd_patch      (:)   ! patch sensible heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_veg_patch       (:)   ! patch sensible heat flux from leaves (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_snow_patch      (:)   ! patch sensible heat flux from snow (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_soil_patch      (:)   ! patch sensible heat flux from soil  (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_h2osfc_patch    (:)   ! patch sensible heat flux from surface water (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_patch       (:)   ! patch total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_u_patch     (:)   ! patch urban total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_tot_r_patch     (:)   ! patch rural total sensible heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_sh_precip_conversion_col(:) ! col sensible heat flux from precipitation conversion (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_patch       (:)   ! patch total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_u_patch     (:)   ! patch urban total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_tot_r_patch     (:)   ! patch rural total latent heat flux (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vegt_patch      (:)   ! patch transpiration heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_vege_patch      (:)   ! patch evaporation heat flux from veg (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_lh_grnd_patch      (:)   ! patch evaporation heat flux from ground (W/m**2) [+ to atm]
     real(r8), pointer :: eflx_soil_grnd_patch    (:)   ! patch soil heat flux (W/m**2) [+ = into soil] 
     real(r8), pointer :: eflx_soil_grnd_u_patch  (:)   ! patch urban soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_soil_grnd_r_patch  (:)   ! patch rural soil heat flux (W/m**2) [+ = into soil]
     real(r8), pointer :: eflx_lwrad_net_patch    (:)   ! patch net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_net_r_patch  (:)   ! patch rural net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_net_u_patch  (:)   ! patch urban net infrared (longwave) rad (W/m**2) [+ = to atm]
     real(r8), pointer :: eflx_lwrad_out_patch    (:)   ! patch emitted infrared (longwave) radiation (W/m**2)
     real(r8), pointer :: eflx_lwrad_out_r_patch  (:)   ! patch rural emitted infrared (longwave) rad (W/m**2)
     real(r8), pointer :: eflx_lwrad_out_u_patch  (:)   ! patch urban emitted infrared (longwave) rad (W/m**2)
     real(r8), pointer :: eflx_snomelt_col        (:)   ! col snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_snomelt_r_col      (:)   ! col rural snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_snomelt_u_col      (:)   ! col urban snow melt heat flux (W/m**2)
     real(r8), pointer :: eflx_gnet_patch         (:)   ! patch net heat flux into ground  (W/m**2)
     real(r8), pointer :: eflx_grnd_lake_patch    (:)   ! patch net heat flux into lake / snow surface, excluding light transmission (W/m**2)
     real(r8), pointer :: eflx_dynbal_grc         (:)   ! grc dynamic land cover change conversion energy flux (W/m**2)
     real(r8), pointer :: eflx_bot_col            (:)   ! col heat flux from beneath the soil or ice column (W/m**2)
     real(r8), pointer :: eflx_fgr12_col          (:)   ! col ground heat flux between soil layers 1 and 2 (W/m**2)
     real(r8), pointer :: eflx_fgr_col            (:,:) ! col (rural) soil downward heat flux (W/m2) (1:nlevgrnd)  (pos upward; usually eflx_bot >= 0)
     real(r8), pointer :: eflx_building_heat_errsoi_col(:) ! col heat flux to interior surface of walls and roof for errsoi check (W m-2)
     real(r8), pointer :: eflx_urban_ac_col       (:)   ! col urban air conditioning flux (W/m**2)
     real(r8), pointer :: eflx_urban_heat_col     (:)   ! col urban heating flux (W/m**2)
     real(r8), pointer :: eflx_anthro_patch       (:)   ! patch total anthropogenic heat flux (W/m**2)
     real(r8), pointer :: eflx_traffic_patch      (:)   ! patch traffic sensible heat flux (W/m**2)
     real(r8), pointer :: eflx_wasteheat_patch    (:)   ! patch sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
     real(r8), pointer :: eflx_heat_from_ac_patch (:)   ! patch sensible heat flux put back into canyon due to removal by AC (W/m**2)
     real(r8), pointer :: eflx_traffic_lun        (:)   ! lun traffic sensible heat flux (W/m**2)
     real(r8), pointer :: eflx_wasteheat_lun      (:)   ! lun sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
     real(r8), pointer :: eflx_heat_from_ac_lun   (:)   ! lun sensible heat flux to be put back into canyon due to removal by AC (W/m**2)
     real(r8), pointer :: eflx_building_lun       (:)   ! lun building heat flux from change in interior building air temperature (W/m**2)
     real(r8), pointer :: eflx_urban_ac_lun       (:)   ! lun urban air conditioning flux (W/m**2)
     real(r8), pointer :: eflx_urban_heat_lun     (:)   ! lun urban heating flux (W/m**2)

     ! Derivatives of energy fluxes
     real(r8), pointer :: dgnetdT_patch           (:)   ! patch derivative of net ground heat flux wrt soil temp  (W/m**2 K)
     real(r8), pointer :: netrad_patch            (:)   ! col net radiation (W/m**2) [+ = to sfc]
     real(r8), pointer :: cgrnd_patch             (:)   ! col deriv. of soil energy flux wrt to soil temp [W/m2/k]
     real(r8), pointer :: cgrndl_patch            (:)   ! col deriv. of soil latent heat flux wrt soil temp  [W/m**2/k]
     real(r8), pointer :: cgrnds_patch            (:)   ! col deriv. of soil sensible heat flux wrt soil temp [W/m2/k]

     ! Canopy radiation
     real(r8), pointer :: dlrad_patch             (:)   ! col downward longwave radiation below the canopy [W/m2]
     real(r8), pointer :: ulrad_patch             (:)   ! col upward longwave radiation above the canopy [W/m2]

     ! Wind Stress
     real(r8), pointer :: taux_patch              (:)   ! patch wind (shear) stress: e-w (kg/m/s**2)
     real(r8), pointer :: tauy_patch              (:)   ! patch wind (shear) stress: n-s (kg/m/s**2)

     ! Conductance
     real(r8), pointer :: canopy_cond_patch       (:)   ! patch tracer conductance for canopy [m/s] 

     ! Transpiration
     real(r8), pointer :: btran_patch             (:)   ! patch transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_patch         (:)   ! patch daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: btran_min_inst_patch    (:)   ! patch instantaneous daily minimum transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsun_patch              (:)   ! patch sunlit canopy transpiration wetness factor (0 to 1)
     real(r8), pointer :: bsha_patch              (:)   ! patch shaded canopy transpiration wetness factor (0 to 1)

     ! Roots
     real(r8), pointer :: rresis_patch            (:,:) ! patch root resistance by layer (0-1)  (nlevgrnd)

     ! Latent heat
     real(r8), pointer :: htvp_col                (:)   ! latent heat of vapor of water (or sublimation) [j/kg]

     ! Canopy heat
     real(r8), pointer :: dhsdt_canopy_patch      (:)   ! patch change in heat content of canopy (leaf+stem) (W/m**2) [+ to atm]

     ! Balance Checks
     real(r8), pointer :: errsoi_patch            (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errsoi_col              (:)   ! soil/lake energy conservation error   (W/m**2)
     real(r8), pointer :: errseb_patch            (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errseb_col              (:)   ! surface energy conservation error     (W/m**2)
     real(r8), pointer :: errsol_patch            (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errsol_col              (:)   ! solar radiation conservation error    (W/m**2)
     real(r8), pointer :: errlon_patch            (:)   ! longwave radiation conservation error (W/m**2)
     real(r8), pointer :: errlon_col              (:)   ! longwave radiation conservation error (W/m**2)


   contains

    procedure , public :: Init

  end type energyflux_type
  type(energyflux_type), public, target, save :: energyflux_inst

contains

!---------------------------------------------
  subroutine Init(this, bounds)
    
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(energyflux_type)        :: this

    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !--------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
    begl = bounds%begl ; endl= bounds%endl
    begg = bounds%begg ; endg = bounds%endg

    allocate( this%eflx_h2osfc_to_snow_col (begc:endc))             ; this%eflx_h2osfc_to_snow_col (:)   = nan
    allocate( this%eflx_sh_snow_patch      (begp:endp))             ; this%eflx_sh_snow_patch      (:)   = nan
    allocate( this%eflx_sh_soil_patch      (begp:endp))             ; this%eflx_sh_soil_patch      (:)   = nan
    allocate( this%eflx_sh_h2osfc_patch    (begp:endp))             ; this%eflx_sh_h2osfc_patch    (:)   = nan
    allocate( this%eflx_sh_tot_patch       (begp:endp))             ; this%eflx_sh_tot_patch       (:)   = nan
    allocate( this%eflx_sh_tot_u_patch     (begp:endp))             ; this%eflx_sh_tot_u_patch     (:)   = nan
    allocate( this%eflx_sh_tot_r_patch     (begp:endp))             ; this%eflx_sh_tot_r_patch     (:)   = nan
    allocate( this%eflx_sh_grnd_patch      (begp:endp))             ; this%eflx_sh_grnd_patch      (:)   = nan
    allocate( this%eflx_sh_stem_patch      (begp:endp))             ; this%eflx_sh_stem_patch      (:)   = nan
    allocate( this%eflx_sh_veg_patch       (begp:endp))             ; this%eflx_sh_veg_patch       (:)   = nan
    allocate( this%eflx_sh_precip_conversion_col(begc:endc))        ; this%eflx_sh_precip_conversion_col(:) = nan
    allocate( this%eflx_lh_tot_u_patch     (begp:endp))             ; this%eflx_lh_tot_u_patch     (:)   = nan
    allocate( this%eflx_lh_tot_patch       (begp:endp))             ; this%eflx_lh_tot_patch       (:)   = nan
    allocate( this%eflx_lh_tot_r_patch     (begp:endp))             ; this%eflx_lh_tot_r_patch     (:)   = nan
    allocate( this%eflx_lh_grnd_patch      (begp:endp))             ; this%eflx_lh_grnd_patch      (:)   = nan
    allocate( this%eflx_lh_vege_patch      (begp:endp))             ; this%eflx_lh_vege_patch      (:)   = nan
    allocate( this%eflx_lh_vegt_patch      (begp:endp))             ; this%eflx_lh_vegt_patch      (:)   = nan
    allocate( this%eflx_soil_grnd_patch    (begp:endp))             ; this%eflx_soil_grnd_patch    (:)   = nan
    allocate( this%eflx_soil_grnd_u_patch  (begp:endp))             ; this%eflx_soil_grnd_u_patch  (:)   = nan
    allocate( this%eflx_soil_grnd_r_patch  (begp:endp))             ; this%eflx_soil_grnd_r_patch  (:)   = nan
    allocate( this%eflx_lwrad_net_patch    (begp:endp))             ; this%eflx_lwrad_net_patch    (:)   = nan
    allocate( this%eflx_lwrad_net_u_patch  (begp:endp))             ; this%eflx_lwrad_net_u_patch  (:)   = nan
    allocate( this%eflx_lwrad_net_r_patch  (begp:endp))             ; this%eflx_lwrad_net_r_patch  (:)   = nan
    allocate( this%eflx_lwrad_out_patch    (begp:endp))             ; this%eflx_lwrad_out_patch    (:)   = nan
    allocate( this%eflx_lwrad_out_u_patch  (begp:endp))             ; this%eflx_lwrad_out_u_patch  (:)   = nan
    allocate( this%eflx_lwrad_out_r_patch  (begp:endp))             ; this%eflx_lwrad_out_r_patch  (:)   = nan
    allocate( this%eflx_gnet_patch         (begp:endp))             ; this%eflx_gnet_patch         (:)   = nan
    allocate( this%eflx_grnd_lake_patch    (begp:endp))             ; this%eflx_grnd_lake_patch    (:)   = nan
    allocate( this%eflx_dynbal_grc         (begg:endg))             ; this%eflx_dynbal_grc         (:)   = nan
    allocate( this%eflx_bot_col            (begc:endc))             ; this%eflx_bot_col            (:)   = nan
    allocate( this%eflx_snomelt_col        (begc:endc))             ; this%eflx_snomelt_col        (:)   = nan
    allocate( this%eflx_snomelt_r_col      (begc:endc))             ; this%eflx_snomelt_r_col      (:)   = nan
    allocate( this%eflx_snomelt_u_col      (begc:endc))             ; this%eflx_snomelt_u_col      (:)   = nan
    allocate( this%eflx_fgr12_col          (begc:endc))             ; this%eflx_fgr12_col          (:)   = nan
    allocate( this%eflx_fgr_col            (begc:endc, 1:nlevgrnd)) ; this%eflx_fgr_col            (:,:) = nan
    allocate( this%eflx_building_heat_errsoi_col  (begc:endc))      ; this%eflx_building_heat_errsoi_col(:)= nan
    allocate( this%eflx_urban_ac_col       (begc:endc))             ; this%eflx_urban_ac_col       (:)   = nan
    allocate( this%eflx_urban_heat_col     (begc:endc))             ; this%eflx_urban_heat_col     (:)   = nan
    allocate( this%eflx_wasteheat_patch    (begp:endp))             ; this%eflx_wasteheat_patch    (:)   = nan
    allocate( this%eflx_traffic_patch      (begp:endp))             ; this%eflx_traffic_patch      (:)   = nan
    allocate( this%eflx_heat_from_ac_patch (begp:endp))             ; this%eflx_heat_from_ac_patch (:)   = nan
    allocate( this%eflx_heat_from_ac_lun   (begl:endl))             ; this%eflx_heat_from_ac_lun   (:)   = nan
    allocate( this%eflx_building_lun       (begl:endl))             ; this%eflx_building_lun       (:)   = nan
    allocate( this%eflx_urban_ac_lun       (begl:endl))             ; this%eflx_urban_ac_lun       (:)   = nan
    allocate( this%eflx_urban_heat_lun     (begl:endl))             ; this%eflx_urban_heat_lun     (:)   = nan
    allocate( this%eflx_traffic_lun        (begl:endl))             ; this%eflx_traffic_lun        (:)   = nan
    allocate( this%eflx_wasteheat_lun      (begl:endl))             ; this%eflx_wasteheat_lun      (:)   = nan
    allocate( this%eflx_anthro_patch       (begp:endp))             ; this%eflx_anthro_patch       (:)   = nan

    allocate( this%dgnetdT_patch           (begp:endp))             ; this%dgnetdT_patch           (:)   = nan
    allocate( this%cgrnd_patch             (begp:endp))             ; this%cgrnd_patch             (:)   = nan
    allocate( this%cgrndl_patch            (begp:endp))             ; this%cgrndl_patch            (:)   = nan
    allocate( this%cgrnds_patch            (begp:endp))             ; this%cgrnds_patch            (:)   = nan
    allocate( this%dlrad_patch             (begp:endp))             ; this%dlrad_patch             (:)   = nan
    allocate( this%ulrad_patch             (begp:endp))             ; this%ulrad_patch             (:)   = nan
    allocate( this%netrad_patch            (begp:endp))             ; this%netrad_patch            (:)   = nan

    allocate( this%taux_patch              (begp:endp))             ; this%taux_patch              (:)   = nan
    allocate( this%tauy_patch              (begp:endp))             ; this%tauy_patch              (:)   = nan

    allocate( this%canopy_cond_patch       (begp:endp))             ; this%canopy_cond_patch       (:)   = nan

    allocate( this%htvp_col                (begc:endc))             ; this%htvp_col                (:)   = nan

    allocate( this%dhsdt_canopy_patch      (begp:endp))             ; this%dhsdt_canopy_patch      (:)   = nan

    allocate(this%rresis_patch             (begp:endp,1:nlevgrnd))  ; this%rresis_patch            (:,:) = nan
    allocate(this%btran_patch              (begp:endp))             ; this%btran_patch             (:)   = nan
    allocate(this%btran_min_patch          (begp:endp))             ; this%btran_min_patch         (:)   = nan
    allocate(this%btran_min_inst_patch     (begp:endp))             ; this%btran_min_inst_patch    (:)   = nan
    allocate( this%bsun_patch              (begp:endp))             ; this%bsun_patch              (:)   = nan
    allocate( this%bsha_patch              (begp:endp))             ; this%bsha_patch              (:)   = nan
    allocate( this%errsoi_patch            (begp:endp))             ; this%errsoi_patch            (:)   = nan
    allocate( this%errsoi_col              (begc:endc))             ; this%errsoi_col              (:)   = nan
    allocate( this%errseb_patch            (begp:endp))             ; this%errseb_patch            (:)   = nan
    allocate( this%errseb_col              (begc:endc))             ; this%errseb_col              (:)   = nan
    allocate( this%errsol_patch            (begp:endp))             ; this%errsol_patch            (:)   = nan
    allocate( this%errsol_col              (begc:endc))             ; this%errsol_col              (:)   = nan
    allocate( this%errlon_patch            (begp:endp))             ; this%errlon_patch            (:)   = nan
    allocate( this%errlon_col              (begc:endc))             ; this%errlon_col              (:)   = nan


  end subroutine Init

end module EnergyFluxType

