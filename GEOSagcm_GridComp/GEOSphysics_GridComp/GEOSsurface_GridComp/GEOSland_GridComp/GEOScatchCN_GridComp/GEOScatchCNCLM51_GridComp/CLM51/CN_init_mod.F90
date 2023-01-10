#include "MAPL_Generic.h" 

module CN_initMod

  use ESMF

  use clm_varcon        , only : clm_varcon_init
  use clm_varpar        , only : VAR_COL, VAR_PFT, clm_varpar_init
  use clm_varctl        , only : use_century_decomp
  use decompMod
  use filterMod
  use CNVegNitrogenStateType
  use CNVegCarbonStateType
  use atm2lndType
  use TemperatureType
  use SoilStateType
  use WaterDiagnosticBulkType
  use CanopyStateType
  use SolarAbsorbedType
  use SurfaceAlbedoType
  use OzoneBaseMod
  use pftconMod
  use WaterFluxType
  use SoilBiogeochemCarbonStateType
  use SoilBiogeochemNitrogenStateType
  use CNProductsMod 
  use SoilBiogeochemStateType
  use CNVegStateType
  use CNVegCarbonFluxType
  use CNVegNitrogenFluxType
  use GridcellType
  use WaterFluxBulkType
  use SoilBiogeochemCarbonFluxType
  use SoilBiogeochemNitrogenFluxType
  use PatchType
  use ColumnType
  use ch4Mod
  use SoilBiogeochemDecompCascadeConType
  use ActiveLayerMod
  use CropType
  use CNDVType
  use LandunitType
  use RootBiophysMod
  use CNMRespMod         , only : readCNMRespParams => readParams
  use CNSharedParamsMod  , only : CNParamsReadShared
  use spmdMod
  use Wateratm2lndBulkType
  use WaterDiagnosticBulkType
  use Wateratm2lndType
  use EnergyFluxType
  use SaturatedExcessRunoffMod
  use WaterStateBulkType
  use WaterStateType
  use FrictionVelocityMod
  use PhotosynthesisMod
  use CNVegetationFacade, only : cn_vegetation_type

  use SoilBiogeochemDecompCascadeBGCMod  , only : init_decompcascade_bgc
  use SoilBiogeochemDecompCascadeCNMod   , only : init_decompcascade_cn
  use SoilBiogeochemDecompCascadeCNMod   , only : readSoilBiogeochemDecompCnParams       => readParams
  use NutrientCompetitionFactoryMod      , only : create_nutrient_competition_method
  use NutrientCompetitionMethodMod       , only : nutrient_competition_method_type
  use SoilBiogeochemDecompMod            , only : readSoilBiogeochemDecompParams         => readParams
  use CNPhenologyMod                     , only : readCNPhenolParams                     => readParams
  use SoilBiogeochemLittVertTranspMod    , only : readSoilBiogeochemLittVertTranspParams => readParams
  use CNPhenologyMod                     , only : CNPhenologyReadNML
  use dynSubgridControlMod               , only : dynSubgridControl_init
  use CNFireFactoryMod                   , only : CNFireReadNML, create_cnfire_method
  use FireMethodType                     , only : fire_method_type
  use SoilBiogeochemNLeachingMod         , only : readSoilBiogeochemNLeachingParams      => readParams
  use SoilBiogeochemCompetitionMod       , only : readSoilBiogeochemCompetitionParams    => readParams
  use SoilBiogeochemPotentialMod         , only : readSoilBiogeochemPotentialParams      => readParams

 
  use clm_varpar       , only : numpft, num_zon, num_veg, var_pft, var_col, &
                                nlevgrnd, nlevsoi

  use MAPL             , only : NetCDF4_FileFormatter

 implicit none 
 private 

 type(photosyns_type), public            :: photosyns_inst
 class(nutrient_competition_method_type), public,  allocatable :: nutrient_competition_method
 class(fire_method_type), allocatable :: cnfire_method

 contains

!------------------------------------------------------
 subroutine CN_init(nch,ityp,fveg,cncol,cnpft,lats,lons,cn5_cold_start)

  !ARGUMENTS
  implicit none
  !INPUT/OUTPUT
  integer,                                         intent(in) :: nch   ! number of tiles
  integer, dimension(nch,num_veg,num_zon),         intent(in) :: ityp  ! PFT index
  real,    dimension(nch,num_veg,num_zon),         intent(in) :: fveg  ! PFT fraction
  real,    dimension(nch,num_zon,var_col),         intent(in) :: cncol ! column-level CN restart variables
  real,    dimension(nch,num_zon,num_veg,var_pft), intent(in) :: cnpft ! patch/pft-level CN restart variables
  real,    dimension(nch),                         intent(in) :: lats  ! Catchment tile latitudes [rad]
  real,    dimension(nch),                         intent(in) :: lons  ! Catchment tile longitudes [rad]
  logical, optional,                               intent(in) :: cn5_cold_start  ! cold start for the CLM variables that are new in Catchment-CN5.0

                                                                                                        
  !LOCAL

  type(bounds_type)                       :: bounds
  type(patch_type)                        :: patch
  type(column_type)                       :: col
  type(landunit_type)                     :: lun
  type(cnveg_nitrogenstate_type)          :: cnveg_nitrogenstate_inst
  type(cnveg_carbonstate_type)            :: cnveg_carbonstate_inst
  type(atm2lnd_type)                      :: atm2lnd_inst
  type(temperature_type)                  :: temperature_inst
  type(soilstate_type)                    :: soilstate_inst
  type(waterdiagnosticbulk_type)          :: waterdiagnosticbulk_inst
  type(wateratm2lndbulk_type)             :: wateratm2lndbulk_inst
  type(wateratm2lnd_type)                 :: wateratm2lnd_inst
  type(canopystate_type)                  :: canopystate_inst
  type(solarabs_type)                     :: solarabs_inst
  type(surfalb_type)                      :: surfalb_inst
  type(ozone_base_type)                   :: ozone_inst
  type(pftcon_type)                       :: pftcon
  type(waterflux_type)                    :: waterflux_inst
  type(soilbiogeochem_carbonstate_type)   :: soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_nitrogenstate_type) :: soilbiogeochem_nitrogenstate_inst
  type(cn_products_type)                  :: c_products_inst
  type(cn_products_type)                  :: n_products_inst
  type(soilbiogeochem_state_type)         :: soilbiogeochem_state_inst
  type(cnveg_state_type)                  :: cnveg_state_inst
  type(cnveg_carbonflux_type)             :: cnveg_carbonflux_inst
  type(cnveg_nitrogenflux_type)           :: cnveg_nitrogenflux_inst
  type(gridcell_type)                     :: grc
  type(soilbiogeochem_carbonflux_type)    :: soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_nitrogenflux_type)  :: soilbiogeochem_nitrogenflux_inst
  type(ch4_type)                          :: ch4_inst
  type(crop_type)                         :: crop_inst
  type(dgvs_type)                         :: dgvs_inst
  type(saturated_excess_runoff_type)      :: saturated_excess_runoff_inst
  type(energyflux_type)                   :: energyflux_inst
  type(waterstatebulk_type)               :: waterstatebulk_inst
  type(waterstate_type)                   :: waterstate_inst
  type(frictionvel_type)                  :: frictionvel_inst
  type(cn_vegetation_type)               :: bgc_vegetation_inst

  character(300)     :: paramfile
  character(300)     :: NLFilename
  type(Netcdf4_fileformatter) :: ncid
  integer            :: rc, status

  integer, parameter :: zeng_2001_root    = 0 !the zeng 2001 root profile function
  !-----------------------------------------

! initialize CN model
! -------------------

    call spmd_init()

    call clm_varpar_init()

    call clm_varcon_init()

    call init_clm_varctl()

    call init_bounds                    (nch, bounds)

    ! initialize subrgid types

    call init_patch_type                (bounds, nch, ityp, fveg, patch)

    call init_column_type               (bounds, nch, col)

    call init_landunit_type             (bounds, lun)

    call init_gridcell_type             (bounds, nch, cnpft, lats, lons, grc)

    ! create subgrid structure

    call clm_ptrs_compdown              (bounds)

    ! initialize filters

    call allocFilters                  (bounds, nch, ityp, fveg)

    ! read parameters and configurations from namelist file

    call CNPhenologyReadNML       ( NLFilename )
    call dynSubgridControl_init   ( )
    call CNFireReadNML            ( NLFilename )

    ! initialize states and fluxes

    call init_cnveg_nitrogenstate_type  (bounds, nch, ityp, fveg, cncol, cnpft, cnveg_nitrogenstate_inst) 

    call init_cnveg_carbonstate_type    (bounds, nch, ityp, fveg, cncol, cnpft, cnveg_carbonstate_inst)

    call init_atm2lnd_type              (bounds, atm2lnd_inst)

    call init_temperature_type          (bounds, temperature_inst)

    call init_soilstate_type            (bounds, soilstate_inst)

    call init_waterdiagnosticbulk_type  (bounds, waterdiagnosticbulk_inst)

    call init_wateratm2lndbulk_type     (bounds, wateratm2lndbulk_inst)

    call init_wateratm2lnd_type         (bounds, wateratm2lnd_inst)

    call init_canopystate_type          (bounds, nch, ityp, fveg, cncol, cnpft, canopystate_inst, cn5_cold_start)

    call init_solarabs_type             (bounds, solarabs_inst)

    call init_surfalb_type              (bounds, nch, cncol, cnpft, surfalb_inst)

    call init_ozone_base_type           (bounds, ozone_inst)

    call photosyns_inst%Init            (bounds, nch, ityp, fveg, cncol, cnpft, cn5_cold_start)

    call init_pftcon_type               (pftcon)

    call init_waterflux_type            (bounds, waterflux_inst)

    call init_soilbiogeochem_carbonstate_type(bounds, nch, cncol,  soilbiogeochem_carbonstate_inst)

    call init_soilbiogeochem_nitrogenstate_type(bounds, nch, cncol,  soilbiogeochem_nitrogenstate_inst)

    call init_cn_products_type          (bounds, nch, cncol, 'C',  c_products_inst)

    call init_cn_products_type          (bounds, nch, cncol, 'N',  n_products_inst)

    call init_soilbiogeochem_state_type (bounds, nch, cncol,  soilbiogeochem_state_inst)

    call init_cnveg_state_type          (bounds, nch, ityp, fveg, cncol, cnpft, cnveg_state_inst)

    call init_cnveg_carbonflux_type     (bounds, nch, ityp, fveg, cncol, cnpft, cnveg_carbonflux_inst, cn5_cold_start)
 
    call init_cnveg_nitrogenflux_type   (bounds, nch, ityp, fveg, cncol, cnpft, cnveg_nitrogenflux_inst)

    call init_waterfluxbulk_type        (bounds, waterfluxbulk_inst)

    call init_soilbiogeochem_carbonflux_type(bounds,soilbiogeochem_carbonflux_inst)

    call init_soilbiogeochem_nitrogenflux_type(bounds,soilbiogeochem_nitrogenflux_inst)

    call init_ch4_type                  (bounds, ch4_inst)

    call init_decomp_cascade_constants  (use_century_decomp)
  
    call init_active_layer_type         (bounds, active_layer_inst)

    call init_crop_type                 (bounds, crop_inst)

    call init_dgvs_type                 (bounds, dgvs_inst)

    call init_saturated_excess_runoff_type(bounds, saturated_excess_runoff_inst)

    call init_energyflux_type           (bounds, energyflux_inst)  

    call init_waterstatebulk_type       (bounds, waterstatebulk_inst)

    call init_waterstate_type           (bounds, waterstate_inst)

    call init_frictionvel_type          (bounds, frictionvel_inst)

    call CNPhenologyInit                (bounds)

    call bgc_vegetation_inst%cn_balance_inst%Init      (bounds)
    call create_cnfire_method( bgc_vegetation_inst%cnfire_method)
 

    ! calls to original CTSM initialization routines

    ! initialize rooting profile with default values
    rooting_profile_method_water    = zeng_2001_root
    rooting_profile_method_carbon   = zeng_2001_root
    rooting_profile_varindex_water  = 1 
    rooting_profile_varindex_carbon = 2 


    ! initialize root fractions
    call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%rootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd),'water')
    call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%crootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd),'carbon')

   ! allocate CLM arrays that are not allocated in their modules

    allocate(nutrient_competition_method, &
         source=create_nutrient_competition_method(bounds)) ! jkolassa: this allocates and initializes the nutrient_competition_method_type

   ! initialize CLM parameters from parameter file

   paramfile = '/discover/nobackup/jkolassa/CLM/parameter_files/ctsm51_params.c210923.nc'

   call ncid%open(trim(paramfile),pFIO_READ, __RC__)

   call readCNMRespParams(ncid)
   call CNParamsReadShared(ncid)  ! this is called CN params but really is for the soil biogeochem parameters
   call readSoilBiogeochemDecompCnParams(ncid) 
   call nutrient_competition_method%readParams(ncid)
   call readSoilBiogeochemDecompParams(ncid)
   call readCNPhenolParams(ncid)
   call readSoilBiogeochemLittVertTranspParams(ncid)
   call photosyns_inst%ReadParams( ncid )
   call readSoilBiogeochemNLeachingParams(ncid)
   call readSoilBiogeochemCompetitionParams(ncid)
   call readSoilBiogeochemPotentialParams(ncid)
   call bgc_vegetation_inst%cnfire_method%CNFireReadParams( ncid )

   call ncid%close(rc=status)

   call FireMethodInit(bounds,paramfile)

   if (use_century_decomp) then
      call init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, &
                                      soilstate_inst )
   else
      call init_decompcascade_cn(bounds, soilbiogeochem_state_inst)
   end if

 end  subroutine CN_init
end module CN_initMod


