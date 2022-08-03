module CN_DriverMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use nanMod           , only : nan
  use CNVegetationFacade
  use clm_varpar       , only : nlevsno, nlevmaxurbgrnd, num_zon
  use clm_varcon       , only : grav, denh2o


contains

!---------------------------------
 subroutine CN_Driver(nch,ndep,tp1,tairm,rzm,psis,bee,dayl)

 use CNCLM_decompMod, only : bounds
 use CNCLM_filterMod, only : filter
 use CNCLM_SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
 use CNCLM_SoilBiogeochemNitrogenFluxType, only : soilbiogeochem_nitrogenflux_type
 use CNCLM_ActiveLayerMod
 use CNCLM_GridcellType

 !ARGUMENTS
 implicit none
 
 !INPUT
 integer,              intent(in) :: nch     ! number of tiles 
 real, dimension(nch), intent(in) :: ndep    ! nitrogen deposition [g m^-2 s^-1]
 real, dimension(nch), intent(in) :: tp1     ! soil temperatures [K]
 real, dimension(nch), intent(in) :: tairm   ! surface air temperature [K] averaged over CN interval
 real, dimension(nch,nzone), intent(in) :: rzm ! weighted root-zone moisture content as frac of WHC
 real, dimension(nch), intent(in) :: bee     ! Clapp-Hornberger 'b' [-]
 real, dimension(nch), intent(in) :: psis    ! saturated matric potential [m]
 real, dimension(nch), intent(in) :: dayl    ! daylength [seconds]

 !LOCAL

 ! jkolassa: not sure the below type declarations are necessary or whether use statements 
 ! above are enough

 type(bounds_type)                      :: bounds
 type(clumpfilter_type)                 :: filter
 type(soilbiogeochem_carbonflux_type)   :: soilbiogeochem_carbonflux_inst
 type(soilbiogeochem_nitrogenflux_type) :: soilbiogeochem_nitrogenflux_inst
 type(gridcell_type)                    :: grc

 logical, save :: doalb = .true.         ! assume surface albedo calculation time step; jkolassa: following setting from previous CNCLM versions
 integer  :: n, p, nc, nz, np

 !-------------------------------

  ! update CLM types with current states

  n = 0
  p = 0
  do nc = 1,nch        ! catchment tile loop

     grc%dayl(nc) = dayl(nc)

     do nz = 1,num_zon    ! CN zone loop
        n = n + 1
        temperature_inst%t_soisno_col(n,-nlevsno+1:nlevmaxurbgrnd) = tp1(nc) ! jkolassa: only one soil and no snow column at this point (may change in future)
        soilstate_inst%soilpsi_col(n,nlevgrnd) = 1.e-6*psis(nc)*grav*denh2o*rzm(nc,nz)**(-bee(nc)) ! jkolassa: only one soil layer at this point
        do np = 0,numpft  ! PFT index loop
           p = p + 1
           temperature_inst%t_ref2m_patch(p) = tairm(nc)
        end do ! np
     end do ! nz
  end do ! nc

  ! call CLM routines that are needed prior to Ecosystem Dynamics call

  call active_layer_inst%alt_calc(num_soilc, filter_soilc, &
       temperature_inst)


  ! Ecosystem Dynamics calculations
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



      grc%prev_dayl = grc%dayl ! set previous day length for following time steps (dayl itself is computed in GridComp)
 end subroutine CN_Driver

end module CN_DriverMod
