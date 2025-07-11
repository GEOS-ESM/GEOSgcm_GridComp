#include "MAPL_Generic.h" 

module CNCLM_initMod

  use ESMF

  use clm_varcon        , only : clm_varcon_init
  use clm_varpar        , only : VAR_COL, VAR_PFT, clm_varpar_init
  use clm_varctl        , only : use_century_decomp, init_clm_varctl
  use clm_time_manager  , only : get_step_size, update_rad_dtime
  use decompMod         , only : bounds
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
  use SoilBiogeochemDecompCascadeConType, only : init_decomp_cascade_constants
  use ActiveLayerMod
  use CropType
  use CNDVType
  use LandunitType          , only : lun
  use RootBiophysMod
  use CNMRespMod         , only : readCNMRespParams => readParams, CNMRespReadNML
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
  use CNVegetationFacade
  use initSubgridMod
  use CN2CLMType
  use WaterType         , only : water_type
  use CNBalanceCheckMod

  use SoilBiogeochemDecompCascadeBGCMod  , only : init_decompcascade_bgc, DecompCascadeBGCreadNML
  use SoilBiogeochemDecompCascadeBGCMod  , only : readSoilBiogeochemDecompBgcParams      => readParams
  use SoilBiogeochemDecompCascadeCNMod   , only : init_decompcascade_cn
  use SoilBiogeochemDecompCascadeCNMod   , only : readSoilBiogeochemDecompCnParams       => readParams
  use NutrientCompetitionFactoryMod      , only : create_nutrient_competition_method
  use NutrientCompetitionMethodMod       , only : nutrient_competition_method_type
  use SoilBiogeochemDecompMod            , only : readSoilBiogeochemDecompParams         => readParams
  use CNPhenologyMod                     , only : readCNPhenolParams                     => readParams
  use SoilBiogeochemLittVertTranspMod    , only : readSoilBiogeochemLittVertTranspParams => readParams
  use CNPhenologyMod                     , only : CNPhenologyReadNML, CNPhenologyInit
  use dynSubgridControlMod               , only : dynSubgridControl_init
  use CNFireFactoryMod                   , only : CNFireReadNML, create_cnfire_method
  use FireMethodType                     , only : fire_method_type
  use SoilBiogeochemNLeachingMod         , only : readSoilBiogeochemNLeachingParams      => readParams
  use SoilBiogeochemCompetitionMod       , only : readSoilBiogeochemCompetitionParams    => readParams
  use SoilBiogeochemCompetitionMod       , only : SoilBiogeochemCompetitionInit
  use SoilBiogeochemPotentialMod         , only : readSoilBiogeochemPotentialParams      => readParams
  use CNGapMortalityMod                  , only : readCNGapMortalityParams               => readParams
  use CNFUNMod                           , only : readCNFUNParams                        => readParams, &
                                                  CNFUNInit
  use CNNDynamicsMod                     , only : CNNDynamicsReadNML
  use SurfaceAlbedoMod                   , only: SurfaceAlbedo_readnl
  use SoilBiogeochemPrecisionControlMod  , only: SoilBiogeochemPrecisionControlInit
  use SoilBiogeochemNitrifDenitrifMod    , only : readSoilBiogeochemNitrifDenitrifParams => readParams 
  use SoilStateInitTimeConstMod          , only : readParams_SoilStateInitTimeConst      => readParams, &
                                                  SoilStateInitTimeConst
  use clm_varpar       , only : numpft, num_zon, num_veg, var_pft, var_col, &
                                nlevgrnd, nlevsoi

  use MAPL             , only : NetCDF4_FileFormatter, pFIO_READ
  use ESMF

 implicit none 
 private 

 class(nutrient_competition_method_type), public,  allocatable :: nutrient_competition_method
 class(fire_method_type),                          allocatable :: cnfire_method
 type(saturated_excess_runoff_type),      public               :: saturated_excess_runoff_inst
 type(waterflux_type),                    public               :: waterflux_inst
 type(waterfluxbulk_type),                public               :: waterfluxbulk_inst



! !PUBLIC MEMBER FUNCTIONS:
 public :: CN_init

 contains

!------------------------------------------------------
 subroutine CN_init(nch,ityp,fveg,cncol,cnpft,lats,lons,dtcn,paramfile,water_inst,bgc_vegetation_inst,cn5_cold_start)

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
  real,                                            intent(in) :: dtcn  ! Catchment-CN step size
  character(len=ESMF_MAXSTR),                      intent(in) :: paramfile
  logical, optional,                               intent(in) :: cn5_cold_start  ! cold start for the CLM variables that are new in Catchment-CN5.0
  type(water_type), intent(out)                       :: water_inst
  type(cn_vegetation_type), intent(out)               :: bgc_vegetation_inst                                                                                                        
  !LOCAL

  character(300)     :: NLFilename
  type(Netcdf4_fileformatter) :: ncid
  integer            :: rc, status, ndt

  !-----------------------------------------

! initialize CN step size

  ndt = get_step_size( nint(dtcn) )

! initialize CN model
! -------------------

    call spmd_init()

    call clm_varpar_init()

    call clm_varcon_init()

    call init_clm_varctl()

    call bounds%Init                    (nch)

    ! initialize subrgid types

    call patch%Init                     (bounds, nch, ityp, fveg)

    call col%Init                       (bounds, nch)

    call lun%Init                       (bounds, nch)

    call grc%Init                       (bounds, nch, cnpft, lats, lons)

    ! create subgrid structure

    call clm_ptrs_compdown              (bounds)

    ! initialize filters

    call allocFilters                  (bounds, nch, ityp, fveg)

    ! read parameters and configurations from namelist file

    NLFilename = trim('./CN_CLM51.nml')
    
    call CNPhenologyReadNML       ( NLFilename )
    call dynSubgridControl_init   ( )
    call CNFireReadNML            ( NLFilename )
    call CNNDynamicsReadNML       ( NLFilename )
    call photosyns_inst%ReadNML   ( NLFilename )
    call canopystate_inst%ReadNML ( NLFilename )
    call DecompCascadeBGCreadNML  ( NLFilename )
    call CNMRespReadNML           ( NLFilename )
    call SurfaceAlbedo_readnl     ( NLFilename )

    ! initialize states and fluxes

    call pftcon%init_pftcon_type        (paramfile)

    call bgc_vegetation_inst%Init(bounds, NLFilename, nch, ityp, fveg, cncol, cnpft, paramfile, cn5_cold_start)

    call atm2lnd_inst%Init              (bounds)

    call temperature_inst%Init          (bounds)

    call soilstate_inst%Init            (bounds)

    call SoilStateInitTimeConst         (bounds, soilstate_inst, NLFilename) ! sets hydraulic and thermal soil properties

    call water_inst%Init                (bounds)

    call canopystate_inst%Init          (bounds, nch, ityp, fveg, cncol, cnpft, cn5_cold_start)

    call solarabs_inst%Init             (bounds)

    call surfalb_inst%Init              (bounds, nch, cncol, cnpft)

    call ozone_inst%Init           (bounds)

    call photosyns_inst%Init            (bounds, nch, ityp, fveg, cncol, cnpft, cn5_cold_start)

    call soilbiogeochem_carbonstate_inst%Init(bounds, nch, cncol)

    call soilbiogeochem_nitrogenstate_inst%Init(bounds, nch, cncol)

    call soilbiogeochem_state_inst%Init (bounds, nch, cncol, cnpft, ityp, fveg)

    call soilbiogeochem_carbonflux_inst%Init (bounds)

    call soilbiogeochem_nitrogenflux_inst%Init(bounds)

    call ch4_inst%Init                  (bounds)

    call init_decomp_cascade_constants  (use_century_decomp)
  
    call active_layer_inst%Init         (bounds)

    call crop_inst%Init                 (bounds)

    call dgvs_inst%Init                 (bounds)

    call saturated_excess_runoff_inst%Init(bounds)

    call energyflux_inst%Init           (bounds)  

    call frictionvel_inst%Init          (bounds)

    call cn_balance_inst%Init           (bounds)


    ! initialize rooting profile parameters from namelist

    call init_rootprof(NLFilename)

    ! initialize root fractions
    call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%rootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd),'water')
    call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
             soilstate_inst%crootfr_patch(bounds%begp:bounds%endp,1:nlevgrnd),'carbon')

   ! allocate CLM arrays that are not allocated in their modules

    allocate(nutrient_competition_method, &
         source=create_nutrient_competition_method(bounds)) ! jkolassa: this allocates and initializes the nutrient_competition_method_type

   ! initialize CLM parameters from parameter file

   call ncid%open(trim(paramfile),pFIO_READ, RC=status)

   call readCNMRespParams(ncid)
   call CNParamsReadShared(ncid, NLFilename)  ! this is called CN params but really is for the soil biogeochem parameters
   call readSoilBiogeochemDecompCnParams(ncid) 
   call readSoilBiogeochemDecompBgcParams(ncid)
   call nutrient_competition_method%readParams(ncid)
   call readSoilBiogeochemDecompParams(ncid)
   call readCNPhenolParams(ncid)
   call readSoilBiogeochemLittVertTranspParams(ncid)
   call photosyns_inst%ReadParams( ncid )
   call readSoilBiogeochemNLeachingParams(ncid)
   call readSoilBiogeochemCompetitionParams(ncid)
   call readSoilBiogeochemPotentialParams(ncid)
   call readCNGapMortalityParams(ncid)
   call readCNFUNParams(ncid)
   call readSoilBiogeochemNitrifDenitrifParams(ncid)
   !call readParams_SoilStateInitTimeConst(ncid)   

   call ncid%close(rc=status)

   ! initialize types that depend on parameters

   call CNPhenologyInit                (bounds)
   call SoilBiogeochemCompetitionInit  (bounds)
   call CNFUNInit(bounds,bgc_vegetation_inst%cnveg_state_inst,bgc_vegetation_inst%cnveg_carbonstate_inst,bgc_vegetation_inst%cnveg_nitrogenstate_inst)

   ! Initialize precision control for soil biogeochemistry (use soilbiogeochem_carbonstate three times, since we   do not currently use isotopes)
    call SoilBiogeochemPrecisionControlInit( soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonstate_inst, &         
                                             soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst)


  ! call FireMethodInit(bounds,paramfile)

   if (use_century_decomp) then
      call init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, &
                                      soilstate_inst )
   else
      call init_decompcascade_cn(bounds, soilbiogeochem_state_inst)
   end if

  ! initialize custom type used to pass Catchment information to nested CLM fire types

  call cn2clm_inst%Init         (bounds)

  ! initialize radiation time

  call update_rad_dtime(.true.)

 end subroutine CN_init
end module CNCLM_initMod


