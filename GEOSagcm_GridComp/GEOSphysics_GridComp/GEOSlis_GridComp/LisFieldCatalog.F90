module LisFieldCatalogMod

   !DESCRIPTION:
   ! Self-contained copy of the stdName/stateName/units metadata and role
   ! classification (adImport/adExport/directConn) that LisFieldsMod
   ! needs from LIS_NUOPC_Gluecode's LIS_FieldList/LIS_HookupInit -- kept here
   ! instead so the Plug does not have to compile/link LISF's NUOPC cap
   ! (LIS_NUOPC_Cap/LIS_NUOPC_Gluecode/LIS_NUOPC_DataCopy/LIS_ESMF_Extensions/
   ! LIS_NUOPC_Flags, and transitively LISWRFGridCompMod and NUOPC itself) or
   ! LIS_FORC_AttributesMod just to get this static data.
   !
   ! Roles mirror LIS_HookupInit's select case in
   ! @LISF/lis/runmodes/nuopc_cpl_mode/LIS_NUOPC_Gluecode.F90 for the default
   ! build only (no WRF_HYDRO/PARFLOW/GSM_EXTLND); keep in sync if that
   ! select case changes upstream.

   implicit none
   private

   public :: LIS_Field
   public :: LIS_FieldList

   type :: LIS_Field
      character(len=:), allocatable :: stdName
      character(len=:), allocatable :: stateName
      character(len=:), allocatable :: units
      logical :: adImport = .false.
      logical :: adExport = .false.
      logical :: directConn = .false.
   end type LIS_Field

   type :: LIS_FieldList
      type(LIS_Field), allocatable :: fields(:)
   end type LIS_FieldList

   interface LIS_FieldList
      module procedure new_LIS_FieldList
   end interface LIS_FieldList

contains

   ! Constructs the catalog. Roles (adImport/adExport/directConn) are baked
   ! in directly, mirroring LIS_HookupInit's select case (default build only,
   ! no WRF_HYDRO/PARFLOW/GSM_EXTLND -- those fields keep all-.FALSE. defaults).
   function new_LIS_FieldList() result(this)
      type(LIS_FieldList) :: this

      this%fields = [ &
           LIS_Field( &
                stdName='2m_air_temperature', &
                stateName='t2_f', &
                units='K', &
                adImport=.true.), &
           LIS_Field( &
                stdName='2m_heat_exchange_coefficient', &
                stateName='chs2_f', &
                units='m s-1', &
                adImport=.true., &
                adExport=.true.), &
           LIS_Field( &
                stdName='2m_moisture_exchange_coefficient', &
                stateName='cqs2_f', &
                units='m s-1', &
                adImport=.true., &
                adExport=.true.), &
           LIS_Field( &
                stdName='2m_potential_temperature', &
                stateName='th2_f', &
                units='K', &
                adImport=.true.), &
           LIS_Field( &
                stdName='2m_specific_humidity', &
                stateName='q2_f', &
                units='kg kg-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='air_temperature', &
                stateName='tair_f', &
                units='K', &
                adImport=.true.), &
           LIS_Field( &
                stdName='albedo', &
                stateName='albedo_f', &
                units='1', &
                adImport=.true., &
                adExport=.true.), &
           LIS_Field( &
                stdName='albedo_w_snow_effect', &
                stateName='albedo_snwff', &
                units='1', &
                adExport=.true.), &
           LIS_Field( &
                stdName='atmospheric_density', &
                stateName='density_f', &
                units='kg m-3', &
                adImport=.true.), &
           LIS_Field( &
                stdName='canopy_moisture', &
                stateName='canopmoist', &
                units='m', &
                adExport=.true.), &
           LIS_Field( &
                stdName='convective_available_potential_energy', &
                stateName='cape_f', &
                units='J kg-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='convective_rainfall_flux', &
                stateName='crainf_f', &
                units='kg m-2 s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='cosine_solar_zenith_angle', &
                stateName='coszenith_f', &
                units='-', &
                adImport=.true.), &
           LIS_Field( &
                stdName='downward_heat_flux_in_soil', &
                stateName='qg', &
                units='W m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='eastward_wind', &
                stateName='ewind_f', &
                units='m s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='effective_mixing_ratio', &
                stateName='effmixratio', &
                units='kg kg-1', &
                adExport=.true.), &
           LIS_Field( &
                stdName='emissivity', &
                stateName='emiss_f', &
                units='-', &
                adImport=.true., &
                adExport=.true.), &
           LIS_Field( &
                stdName='forcing_height', &
                stateName='fheight_f', &
                units='m', &
                adImport=.true.), &
           LIS_Field( &
                stdName='green_vegetation_fraction', &
                stateName='greenness_f', &
                units='1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='heat_exchange_coefficient_in_air', &
                stateName='ch_f', &
                units='-', &
                adImport=.true.), &
           LIS_Field( &
                stdName='latent_heat_flux_kinematic', &
                stateName='qlekinematic', &
                units='kg m-2 s-1', &
                adExport=.true.), &
           LIS_Field( &
                stdName='level_pressure', &
                stateName='lpressure_f', &
                units='Pa', &
                adImport=.true.), &
           LIS_Field( &
                stdName='liquid_fraction_of_soil_moisture_layer_1', &
                stateName='smliqfracl1', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='liquid_fraction_of_soil_moisture_layer_2', &
                stateName='smliqfracl2', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='liquid_fraction_of_soil_moisture_layer_3', &
                stateName='smliqfracl3', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='liquid_fraction_of_soil_moisture_layer_4', &
                stateName='smliqfracl4', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='liquid_water_content_of_surface_snow', &
                stateName='swe', &
                units='kg m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='mixing_ratio', &
                stateName='mixratio_f', &
                units='kg kg-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='momentum_exchange_coefficient_in_air', &
                stateName='cm_f', &
                units='m s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='northward_wind', &
                stateName='nwind_f', &
                units='m s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='ozone_concentration', &
                stateName='o3_f', &
                units='kg kg-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='porosity', &
                stateName='porosity', &
                units='-'), &
           LIS_Field( &
                stdName='potential_evaporation', &
                stateName='pet_f', &
                units='kg m-2 s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='rainfall_flux', &
                stateName='rainf_f', &
                units='kg m-2 s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='reference_et', &
                stateName='refet_f', &
                units='kg m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='relative_soil_moisture', &
                stateName='relsmc', &
                units='m3 m-3', &
                adExport=.true.), &
           LIS_Field( &
                stdName='root_zone_soil_moisture', &
                stateName='rootmoist', &
                units='m3 m-3', &
                adExport=.true.), &
           LIS_Field( &
                stdName='seaicemask', &
                stateName='xice_f', &
                units='-', &
                adImport=.true., &
                adExport=.true.), &
           LIS_Field( &
                stdName='snow_depth', &
                stateName='snowdepth', &
                units='m', &
                adExport=.true.), &
           LIS_Field( &
                stdName='snowfall_flux', &
                stateName='snowf_f', &
                units='kg m-2 s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='snowflag', &
                stateName='snowflag_f', &
                units='-', &
                adImport=.true.), &
           LIS_Field( &
                stdName='snowmelt', &
                stateName='qsm', &
                units='kg m-2 s-1', &
                adExport=.true.), &
           LIS_Field( &
                stdName='soil_moisture_content', &
                stateName='soilmoist', &
                units='kg m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='soil_moisture_fraction_layer_1', &
                stateName='smfracl1', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_moisture_fraction_layer_2', &
                stateName='smfracl2', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_moisture_fraction_layer_3', &
                stateName='smfracl3', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_moisture_fraction_layer_4', &
                stateName='smfracl4', &
                units='m3 m-3', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_temperature_layer_1', &
                stateName='soiltempl1', &
                units='K', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_temperature_layer_2', &
                stateName='soiltempl2', &
                units='K', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_temperature_layer_3', &
                stateName='soiltempl3', &
                units='K', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_temperature_layer_4', &
                stateName='soiltempl4', &
                units='K', &
                adImport=.true., &
                adExport=.true., &
                directConn=.true.), &
           LIS_Field( &
                stdName='soil_temperature_lower_boundary', &
                stateName='tmn_f', &
                units='K', &
                adImport=.true.), &
           LIS_Field( &
                stdName='specific_humidity', &
                stateName='qair_f', &
                units='kg kg-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='subsurface_runoff_amount', &
                stateName='qsb', &
                units='kg m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='surface_air_pressure', &
                stateName='psurf_f', &
                units='Pa', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_diffuse_downwelling_shortwave_flux_in_air', &
                stateName='diffusesw_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_direct_downwelling_shortwave_flux_in_air', &
                stateName='directsw_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_downward_par_diffuse', &
                stateName='pardf_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_downward_par_direct', &
                stateName='pardr_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_downwelling_longwave_flux_in_air', &
                stateName='lwdown_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_downwelling_shortwave_flux_in_air', &
                stateName='swdown_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_net_downward_shortwave_flux', &
                stateName='swnet_f', &
                units='W m-2', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_roughness_length', &
                stateName='roughness_f', &
                units='m', &
                adImport=.true., &
                adExport=.true.), &
           LIS_Field( &
                stdName='surface_runoff_amount', &
                stateName='qs', &
                units='kg m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='surface_snow_area_fraction', &
                stateName='snowcover', &
                units='-', &
                adExport=.true.), &
           LIS_Field( &
                stdName='surface_specific_humidity', &
                stateName='qsfc_f', &
                units='kg kg-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_temperature', &
                stateName='avgsurft', &
                units='K', &
                adExport=.true.), &
           LIS_Field( &
                stdName='surface_upward_latent_heat_flux', &
                stateName='qle', &
                units='W m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='surface_upward_sensible_heat_flux', &
                stateName='qh', &
                units='W m-2', &
                adExport=.true.), &
           LIS_Field( &
                stdName='vapor_pressure', &
                stateName='vaporpress_f', &
                units='-', &
                adImport=.true.), &
           LIS_Field( &
                stdName='vapor_pressure_deficit', &
                stateName='vaporpressdef_f', &
                units='-', &
                adImport=.true.), &
           LIS_Field( &
                stdName='vegetationtype', &
                stateName='vegtype', &
                units='-', &
                adExport=.true.), &
           LIS_Field( &
                stdName='wind_speed', &
                stateName='wind_f', &
                units='m s-1', &
                adImport=.true.), &
           LIS_Field( &
                stdName='surface_water_depth', &
                stateName='sfcheadrt_f', &
                units='mm'), &
           LIS_Field( &
                stdName='time_step_infiltration_excess', &
                stateName='infxsrt', &
                units='mm'), &
           LIS_Field( &
                stdName='soil_column_drainage', &
                stateName='soldrain', &
                units='mm'), &
           LIS_Field( &
                stdName='final_potential_evaporation', &
                stateName='etp', &
                units='W m-2'), &
           LIS_Field( &
                stdName='accum_plant_transpiration', &
                stateName='ett', &
                units='W m-2'), &
           LIS_Field( &
                stdName='total_water_flux_layer_1', &
                stateName='wtrflx1', &
                units='kg m-2 s-1'), &
           LIS_Field( &
                stdName='total_water_flux_layer_2', &
                stateName='wtrflx2', &
                units='kg m-2 s-1'), &
           LIS_Field( &
                stdName='total_water_flux_layer_3', &
                stateName='wtrflx3', &
                units='kg m-2 s-1'), &
           LIS_Field( &
                stdName='total_water_flux_layer_4', &
                stateName='wtrflx4', &
                units='kg m-2 s-1'), &
           LIS_Field( &
                stdName='ground_water_storage', &
                stateName='wa', &
                units='mm') ]
   end function new_LIS_FieldList

end module LisFieldCatalogMod
