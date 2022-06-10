module CN_DriverMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use nanMod           , only : nan
  use CNVegetationFacade

contains

!---------------------------------
 subroutine CN_Driver(nch,ndep)

 use CNCLM_decompMod, only : bounds
 use CNCLM_filterMod, only : filter
 use CNCLM_SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
 use CNCLM_SoilBiogeochemNitrogenFluxType, only : soilbiogeochem_nitrogenflux_type

 !ARGUMENTS
 implicit none
 
 !INPUT
 integer,              intent(in) :: nch     ! number of tiles 
 real, dimension(nch), intent(in) :: ndep    ! nitrogen deposition


 !LOCAL

 ! jkolassa: not sure the below type declarations are necessary or whether use statements 
 ! above are enough

 type(bounds_type) :: bounds
 type(clumpfilter_type) :: filter
 type(soilbiogeochem_carbonflux_type) :: soilbiogeochem_carbonflux_inst
 type(soilbiogeochem_nitrogenflux_type) :: soilbiogeochem_nitrogenflux_inst

 logical, save :: doalb = .true.         ! assume surface albedo calculation time step; jkolassa: following setting from previous CNCLM versions


 call bgc_vegetation_inst%EcosystemDynamicsPreDrainage(bounds,            &
      filter%num_soilc, filter%soilc,                       &
      filter%num_soilp, filter%soilp,                       &
      filter%num_actfirec, filter%actfirec,                 &
      filter%num_actfirep, filter%actfirep,                 &
      filter%num_pcropp, filter%pcropp, &
      filter%num_exposedvegp, filter%exposedvegp, &
      filter%num_noexposedvegp, filter%noexposedvegp, &
      doalb,              &
      soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,         &
      c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
      c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
      soilbiogeochem_state_inst,                                               &
      soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,     &
      active_layer_inst, &
      atm2lnd_inst, water_inst%waterstatebulk_inst, &
      water_inst%waterdiagnosticbulk_inst, water_inst%waterfluxbulk_inst,      &
      water_inst%wateratm2lndbulk_inst, canopystate_inst, soilstate_inst, temperature_inst, &
      soil_water_retention_curve, crop_inst, ch4_inst, &
      photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,          &
      nutrient_competition_method, fireemis_inst)

 end subroutine CN_Driver

end module CN_DriverMod
