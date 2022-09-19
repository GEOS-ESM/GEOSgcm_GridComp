module CN_DriverMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use nanMod           , only : nan
  use CNVegetationFacade
  use clm_varpar       , only : nlevsno, nlevmaxurbgrnd, num_zon, CN_zone_weight
  use clm_varcon       , only : grav, denh2o


contains

!---------------------------------
 subroutine CN_Driver(nch,ndep,tp1,tairm,psis,bee,dayl,btran_fire,car1m,&
                      rzm,sfm,tm,rhm,windm,rainfm,snowfm,prec10d,prec60d,gdp,&
                      abm,peatf,hdm,lnfm,poros,rh30)

 use CNCLM_decompMod, only : bounds
 use CNCLM_filterMod, only : filter
 use CNCLM_SoilBiogeochemCarbonFluxType, only : soilbiogeochem_carbonflux_type
 use CNCLM_SoilBiogeochemNitrogenFluxType, only : soilbiogeochem_nitrogenflux_type
 use CNCLM_ActiveLayerMod
 use CNCLM_GridcellType
 use FireMethodType              , only : fire_method_type
 use SaturatedExcessRunoffMod    , only : saturated_excess_runoff_type
 use CNCLM_WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
 use CNCLM_atm2lndType           , only : atm2lnd_type
 use Wateratm2lndBulkType        , only : wateratm2lndbulk_type
 use CNCLM_CNVegStateType        , only : cnveg_state_type

 !ARGUMENTS
 implicit none
 
 !INPUT
 integer,              intent(in) :: nch     ! number of tiles 
 real, dimension(nch), intent(in) :: ndep    ! nitrogen deposition [g m^-2 s^-1]
 real, dimension(nch), intent(in) :: tp1     ! soil temperatures [K]
 real, dimension(nch), intent(in) :: tairm   ! surface air temperature [K] averaged over CN interval
 real, dimension(nch), intent(in) :: bee     ! Clapp-Hornberger 'b' [-]
 real, dimension(nch), intent(in) :: psis    ! saturated matric potential [m]
 real, dimension(nch), intent(in) :: dayl    ! daylength [seconds]
 real, dimension(nch,num_zon), intent(in) :: btran_fire
 real, dimension(nch), intent(in) :: car1m     ! fraction of tile that is saturated area
 real, dimension(nch,num_zon), intent(in) :: rzm ! weighted root-zone moisture content as frac of WHC
 real, dimension(nch,num_zon), intent(in) :: sfm ! weighted surface moisture content as frac of WHC
 real, dimension(nch), intent(in) :: tm        ! air temperature (K)
 real, dimension(nch), intent(in) :: rhm       ! relative humidity (%)
 real, dimension(nch), intent(in) :: windm     ! wind speed (m/s)
 real, dimension(nch), intent(in) :: rainfm    ! rainfall (convective + largescale) (kg/m2/s)
 real, dimension(nch), intent(in) :: snowfm    ! snowfall (kg/m2/s)  
 real, dimension(nch), intent(in) :: prec10d   ! 10-day running mean of total precipitation (mm H2O/s)
 real, dimension(nch), intent(in) :: prec60d   ! 60-day running mean of total precipitation (mm H2O/s)
 real, dimension(nch), intent(in) :: gdp       ! Real GDP (K 1995US$/capita)
 real, dimension(nch), intent(in) :: abm       ! Peak month for agricultural fire, unitless
 real, dimension(nch), intent(in) :: peatf     ! Fraction of peatland, unitless (0-1)
 real, dimension(nch), intent(in) :: hdm       ! Human population density in 2010 (individual/km2)
 real, dimension(nch), intent(in) :: lnfm      ! Lightning frequency [Flashes/km^2/day]
 real, dimension(nch), intent(in) :: poros     ! porosity
 real, dimension(nch), intent(in) :: rh30      ! 30-day running mean of relative humidity

 !LOCAL

 ! jkolassa: not sure the below type declarations are necessary or whether use statements 
 ! above are enough

 type(bounds_type)                      :: bounds
 type(clumpfilter_type)                 :: filter
 type(soilbiogeochem_carbonflux_type)   :: soilbiogeochem_carbonflux_inst
 type(soilbiogeochem_nitrogenflux_type) :: soilbiogeochem_nitrogenflux_inst
 type(gridcell_type)                    :: grc
 type(cn_vegetation_type), public       :: bgc_vegetation_inst
 type(fire_method_type)                 :: cnfire_method
 type(saturated_excess_runoff_type)     :: saturated_excess_runoff_inst

 logical, save :: doalb = .true.         ! assume surface albedo calculation time step; jkolassa: following setting from previous CNCLM versions
 integer  :: n, p, nc, nz, np

 !-------------------------------

  ! update CLM types with current states

  n = 0
  p = 0
  do nc = 1,nch        ! catchment tile loop

     grc%dayl(nc)                          = dayl(nc)
     wateratm2lndbulk_inst%forc_rh_grc(nc) = rhm(nc)
     atm2lnd_inst%forc_wind_grc(nc)        = windm(nc)
     cnfire_method%forc_hdm(nc)            = hdm(nc)
     cnfire_method%forc_lnfm(nc)           = lnfm(nc)

     do nz = 1,num_zon    ! CN zone loop
        n = n + 1

        temperature_inst%t_soisno_col(n,-nlevsno+1:nlevmaxurbgrnd) = tp1(nc) ! jkolassa: only one soil and no snow column at this point (may change in future)
        temperature_inst%t_grnd_col(n) = temperature_inst%t_soisno_col(n)
        temperature_inst%t_soi17cm_col(n) = temperature_inst%t_grnd_col(n)
        soilstate_inst%soilpsi_col(n,1:nlevgrnd) = 1.e-6*psis(nc)*grav*denh2o*rzm(nc,nz)**(-bee(nc)) ! jkolassa: only one soil layer at this point
        soilstate_inst%watsat_col(n,1:nlevmaxurbgrnd)   = poros(nc) 
        atm2lnd_inst%forc_t_downscaled_col(n)  = tm(nc)
        wateratm2lndbulk_inst%forc_rain_downscaled_col(n) = rainfm(nc)
        wateratm2lndbulk_inst%forc_snow_downscaled_col(n) = snowfm(nc)
        waterdiagnosticbulk_inst%wf_col(n)  = sfm(nc,nz)
        waterdiagnosticbulk_inst%wf2_col(n) = rzm(nc,nz)
        cnveg_state_inst%gdp_lf_col(n) = gdp(nc)
        cnveg_state_inst%abm_lf_col(n) = abm(nc)
        cnveg_state_inst%peatf_lf_col(n) = peatf(nc)

        ! compute column-level saturated area fraction (water table at surface)
        if(nz==1) then
          saturated_excess_runoff_inst%fsat_col(n) = min(max(0.,car1m(nc)/CN_zone_weight(nz)),1.)
        elseif(nz==2) then
          saturated_excess_runoff_inst%fsat_col(n) = min(max(0.,(car1m(nc)-CN_zone_weight(1))/CN_zone_weight(nz)),1.)
        elseif(nz==3)
          saturated_excess_runoff_inst%fsat_col(n) = min(max(0.,(car1m(nc)-CN_zone_weight(1)-CN_zone_weight(2))/CN_zone_weight(nz)),1.)
        endif

        do np = 0,numpft  ! PFT index loop
           p = p + 1
           temperature_inst%t_ref2m_patch(p) = tairm(nc)
           cnfire_method%btran2_patch(p)     = btran_fire(nc,nz)  ! only needed if fire method is Li 2016
           wateratm2lndbulk_inst%prec60_patch(p) = prec60d(nc)
           wateratm2lndbulk_inst%prec10_patch(p) = prec10d(nc)
           wateratm2lndbulk_inst%rh30_patch(p) = rh30(nc)
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
