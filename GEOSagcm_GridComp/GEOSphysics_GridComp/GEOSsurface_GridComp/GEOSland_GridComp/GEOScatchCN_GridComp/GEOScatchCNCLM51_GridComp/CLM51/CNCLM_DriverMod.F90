module CNCLM_DriverMod

 use, intrinsic :: iso_fortran_env, only: INT64
 use nanMod           , only : nan
 use CNVegetationFacade
 use clm_varpar       , only : nlevsno, nlevmaxurbgrnd, num_veg, num_zon, CN_zone_weight,&
                                var_col, var_pft, nlevgrnd, numpft, ndecomp_pools
 use clm_varcon       , only : grav, denh2o
 use clm_time_manager , only : is_first_step, get_nstep
 use decompMod
 use filterMod
 use SoilBiogeochemCarbonFluxType  
 use SoilBiogeochemNitrogenFluxType
 use SoilBiogeochemCarbonStateType  
 use SoilBiogeochemNitrogenStateType
 use ActiveLayerMod               
 use GridcellType               
 use FireMethodType            
 use SaturatedExcessRunoffMod    
 use WaterDiagnosticBulkType    
 use atm2lndType                 
 use Wateratm2lndBulkType       
 use CNVegStateType             
 use WaterStateBulkType         
 use SoilStateType             
 use TemperatureType            
 use WaterDiagnosticBulkType    
 use WaterStateBulkType        
 use WaterFluxBulkType          
 use FrictionVelocityMod       
 use ActiveLayerMod              
 use SoilBiogeochemStateType    
 use CanopyStateType
 use CropType                    
 use ch4Mod                     
 use PhotosynthesisMod           
 use EnergyFluxType            
 use CNFireEmissionsMod        
 use CN_initMod                  
 use CNVegCarbonFluxType        
 use CNVegCarbonStateType        
 use CNVegNitrogenFluxType      
 use CNVegNitrogenStateType     
 use CNProductsMod            
 use CNFireFactoryMod           
 use FireDataBaseType                 
 use CNFireLi2014Mod           
 use CNFireLi2016Mod             
 use CNFireLi2021Mod            
 use CNFireBaseMod 
 use CN2CLMType
 use WaterType

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: CN_Driver
  public :: CN_exit
  public :: get_CN_LAI

contains

!---------------------------------
 subroutine CN_Driver(istep,nch,ityp,fveg,ndep,tp1,tairm,psis,bee,dayl,btran_fire,car1m,&
                      rzm,sfm,rhm,windm,rainfm,snowfm,prec10d,prec60d,et365d,gdp,&
                      abm,peatf,hdm,lnfm,poros,rh30,totwat,bflow,runsrf,sndzn,&
                      fsnow,tg10d,t2m5d,sndzn5d,water_inst,first, &
                      psnsunm, psnsham, lmrsunm, lmrsham, laisunm, laisham, wpwet,   &
                      zlai,zsai,ztai,colc,nppg,gppg,srg,arg,hrg,neeg,burn,closs,nfire,&
                      som_closs,root,vegc,xsmr,ndeployg,denitg,sminn_leachedg,sminng,&
                      col_fire_nlossg,leafng,leafcg,gross_nming,net_nming,&
                      nfix_to_sminng,actual_immobg,fpgg,fpig,sminn_to_plantg,&
                      sminn_to_npoolg,ndep_to_sminng,totvegng,totlitng,totsomng,&
                      retransng,retransn_to_npoolg,fuelcg,totlitcg,cwdcg,rootcg)


 !ARGUMENTS
 implicit none
 
 !INPUT
 integer(INT64),       intent(in) :: istep ! number of CN time steps run
 integer,              intent(in) :: nch     ! number of tiles 
 integer, dimension(nch,num_veg,num_zon),      intent(in) :: ityp ! PFT index
 real, dimension(nch,num_veg,num_zon),         intent(in) :: fveg    ! PFT fraction
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
 real, dimension(nch), intent(in) :: rhm       ! relative humidity (%)
 real, dimension(nch), intent(in) :: windm     ! wind speed (m/s)
 real, dimension(nch), intent(in) :: rainfm    ! rainfall (convective + largescale) (kg/m2/s)
 real, dimension(nch), intent(in) :: snowfm    ! snowfall (kg/m2/s)  
 real, dimension(nch), intent(in) :: prec10d   ! 10-day running mean of total precipitation (mm H2O/s)
 real, dimension(nch), intent(in) :: prec60d   ! 60-day running mean of total precipitation (mm H2O/s)
 real, dimension(nch), intent(in) :: et365d    ! 365-day running mean of total ET (EVPSOI + EVPINT + EVPVEG) (W m-2)
 real, dimension(nch), intent(in) :: gdp       ! Real GDP (K 1995US$/capita)
 real, dimension(nch), intent(in) :: abm       ! Peak month for agricultural fire, unitless
 real, dimension(nch), intent(in) :: peatf     ! Fraction of peatland, unitless (0-1)
 real, dimension(nch), intent(in) :: hdm       ! Human population density in 2010 (individual/km2)
 real, dimension(nch), intent(in) :: lnfm      ! Lightning frequency [Flashes/km^2/day]
 real, dimension(nch), intent(in) :: poros     ! porosity
 real, dimension(nch), intent(in) :: rh30      ! 30-day running mean of relative humidity
 real, dimension(nch), intent(in) :: totwat    ! soil liquid water content [kg m^-2]
 real, dimension(nch), intent(in) :: bflow     ! baseflow
 real, dimension(nch), intent(in) :: runsrf    ! surface runoff [kg m^-2 s^-1]
 real, dimension(nch), intent(in) :: sndzn     ! snow height of snow covered area (m)
 real, dimension(nch), intent(in) :: fsnow     ! snow cover fraction [0-1]
 real, dimension(nch), intent(in) :: tg10d     ! 10-day running mean of ground temperature [K]
 real, dimension(nch), intent(in) :: t2m5d     ! 5-day running mean of daily minimum 2m temperature [K]
 real, dimension(nch), intent(in) :: sndzn5d   ! 5-day running mean of total snow depth
 type(water_type),     intent(in) :: water_inst
 logical,              intent(in) :: first
 real, dimension(nch,num_veg,num_zon), intent(in) :: psnsunm
 real, dimension(nch,num_veg,num_zon), intent(in) :: psnsham
 real, dimension(nch,num_veg,num_zon), intent(in) :: lmrsunm
 real, dimension(nch,num_veg,num_zon), intent(in) :: lmrsham
 real, dimension(nch,num_veg,num_zon), intent(in) :: laisunm
 real, dimension(nch,num_veg,num_zon), intent(in) :: laisham
 real, dimension(nch),                 intent(in) :: wpwet    ! wetness at wilting point


 ! OUTPUT

  real, dimension(nch,num_veg,num_zon), intent(out) :: zlai ! leaf-area index for tile (subject to burying by snow)
  real, dimension(nch,num_veg,num_zon), intent(out) :: zsai ! stem-area index for tile
  real, dimension(nch,num_veg,num_zon), intent(out) :: ztai ! leaf-area index for tile (not buried by snow)

  real, dimension(nch,num_zon),         intent(out) :: colc ! column total carbon
  real, dimension(nch),                 intent(out) :: nppg ! (gC/m2/s) net primary production [PFT]
  real, dimension(nch),                 intent(out) :: gppg ! (gC/m2/s) gross primary production [PFT]
 
  real, dimension(nch),                 intent(out) :: srg  ! (gC/m2/s) total soil respiration (HR + root resp) [column]
  real, dimension(nch),                 intent(out) :: arg  ! (gC/m2/s) autotrophic respiration  [column]
  real, dimension(nch),                 intent(out) :: hrg  ! (gC/m2/s) heterotrophic respiration [column]
  real, dimension(nch),                 intent(out) :: neeg ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source [column]

  real, dimension(nch),                 intent(out) :: fuelcg          ! fuel avalability for non-crop areas outside tropical closed broadleaf evergreen closed forests (gC/m2)
  real, dimension(nch),                 intent(out) :: totlitcg        ! (gC/m2) total litter carbon
  real, dimension(nch),                 intent(out) :: cwdcg            ! (gC/m2) coarse woody debris C
  real, dimension(nch),                 intent(out) :: rootcg          ! (gC/m2) total root carbon
 
  real, dimension(nch),                 intent(out) :: burn  ! burn rate / fractional area burned (/sec)
  real, dimension(nch),                 intent(out) :: closs ! (gC/m2/s) total fire C loss
  real, dimension(nch),                 intent(out) :: nfire ! fire counts (count/km2/s)
  real, dimension(nch),                 intent(out) :: som_closs ! (gC/m2/s) carbon emissions due to peat burning

  real, dimension(nch), intent(out) :: root ! fine root carbon [gC/m2]
  real, dimension(nch), intent(out) :: vegc ! (gC/m2) total vegetation carbon, excluding cpool
  real, dimension(nch), intent(out) :: xsmr ! (gC/m2) abstract C pool to meet excess maintenance respiration (MR) demand

  real, dimension(nch), intent(out) :: ndeployg        ! total N deployed to growth and storage (gN/m2/s)
  real, dimension(nch), intent(out) :: denitg          ! total rate of denitrification (gN/m2/s)
  real, dimension(nch), intent(out) :: sminn_leachedg  ! soil mineral N pool loss to leaching (gN/m2/s)
  real, dimension(nch), intent(out) :: sminng          ! (gN/m2) soil mineral N
  real, dimension(nch), intent(out) :: col_fire_nlossg ! (gN/m2/s) total column-level fire N loss
  real, dimension(nch), intent(out) :: leafng          ! (gN/m2) leaf N
  real, dimension(nch), intent(out) :: leafcg          ! (gC/m2) leaf C
  real, dimension(nch), intent(out) :: gross_nming     ! gross rate of N mineralization (gN/m2/s)
  real, dimension(nch), intent(out) :: net_nming       ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
  real, dimension(nch), intent(out) :: nfix_to_sminng  ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
  real, dimension(nch), intent(out) :: actual_immobg   ! vert-int (diagnostic) actual N immobilization (gN/m2/s) 
  real, dimension(nch), intent(out) :: fpgg            ! fraction of potential gpp (no units)
  real, dimension(nch), intent(out) :: fpig            ! fraction of potential immobilization (no units)
  real, dimension(nch), intent(out) :: sminn_to_plantg ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
  real, dimension(nch), intent(out) :: sminn_to_npoolg ! deployment of soil mineral N uptake (gN/m2/s)
  real, dimension(nch), intent(out) :: ndep_to_sminng  ! atmospheric N deposition to soil mineral N (gN/m2/s)
  real, dimension(nch), intent(out) :: totvegng              ! (gN/m2) total vegetation nitrogen
  real, dimension(nch), intent(out) :: totlitng              ! (gN/m2) total litter nitrogen
  real, dimension(nch), intent(out) :: totsomng              ! (gN/m2) total soil organic matter nitrogen
  real, dimension(nch), intent(out) :: retransng             ! (gN/m2) plant pool of retranslocated N
  real, dimension(nch), intent(out) :: retransn_to_npoolg    ! deployment of retranslocated N (gN/m2/s)

 !LOCAL

 ! jkolassa: not sure the below type declarations are necessary or whether use statements 
 ! above are enough

! type(bounds_type)                      :: bounds
! type(clumpfilter)                      :: filter
! type(soilbiogeochem_carbonflux_type)   :: soilbiogeochem_carbonflux_inst
 type(soilbiogeochem_carbonflux_type)   :: c13_soilbiogeochem_carbonflux_inst
 type(soilbiogeochem_carbonflux_type)   :: c14_soilbiogeochem_carbonflux_inst
! type(soilbiogeochem_nitrogenflux_type) :: soilbiogeochem_nitrogenflux_inst
! type(gridcell_type)                    :: grc
! type(cn_vegetation_type)               :: bgc_vegetation_inst
! type(saturated_excess_runoff_type)     :: saturated_excess_runoff_inst
  type(wateratm2lndbulk_type)            :: wateratm2lndbulk_inst
! type(soilstate_type)                   :: soilstate_inst 
! type(atm2lnd_type)                     :: atm2lnd_inst
! type(temperature_type)                 :: temperature_inst
 type(waterdiagnosticbulk_type)         :: waterdiagnosticbulk_inst
! type(cnveg_state_type)                 :: cnveg_state_inst
 type(waterstatebulk_type)              :: waterstatebulk_inst
! type(waterfluxbulk_type)               :: waterfluxbulk_inst
! type(frictionvel_type)                 :: frictionvel_inst
! type(active_layer_type)                :: active_layer_inst
! type(soilbiogeochem_carbonstate_type)  :: soilbiogeochem_carbonstate_inst
 type(soilbiogeochem_carbonstate_type)  :: c13_soilbiogeochem_carbonstate_inst
 type(soilbiogeochem_carbonstate_type)  :: c14_soilbiogeochem_carbonstate_inst
! type(soilbiogeochem_nitrogenstate_type):: soilbiogeochem_nitrogenstate_inst
! type(soilbiogeochem_state_type)        :: soilbiogeochem_state_inst
! type(crop_type)                        :: crop_inst
! type(ch4_type)                         :: ch4_inst
! type(photosyns_type)                   :: photosyns_inst
! type(energyflux_type)                  :: energyflux_inst
 type(fireemis_type)                    :: fireemis_inst
! type(cnveg_carbonflux_type)            :: cnveg_carbonflux_inst
! type(cnveg_carbonstate_type)           :: cnveg_carbonstate_inst
! type(cnveg_nitrogenflux_type)          :: cnveg_nitrogenflux_inst
! type(cnveg_nitrogenstate_type)         :: cnveg_nitrogenstate_inst
 !type(cnfire_li2014_type)               :: cnfire_li2014_inst
! type(cnfire_li2016_type)               :: cnfire_li2016_inst
! type(cnfire_li2021_type)               :: cnfire_li2021_inst

 real :: pwtgcell
 logical, save :: doalb = .true.         ! assume surface albedo calculation time step; jkolassa: following setting from previous CNCLM versions
 integer  :: n, p, nc, nz, np, nv
 integer(INT64) :: nstep_cn   ! number of CN model steps run

 !-------------------------------

  ! update time step 
  nstep_cn = get_nstep(istep)

  ! update CLM types with current states

  n = 0
  p = 0
  do nc = 1,nch        ! catchment tile loop

     grc%dayl(nc)                          = dayl(nc)
     water_inst%wateratm2lndbulk_inst%forc_rh_grc(nc) = rhm(nc)
     atm2lnd_inst%forc_wind_grc(nc)        = windm(nc)
     atm2lnd_inst%forc_ndep_grc(nc)        = ndep(nc)

     cn2clm_inst%forc_hdm_cn2clm(nc)  = hdm(nc)
     cn2clm_inst%forc_lnfm_cn2clm(nc) = lnfm(nc) 
    ! cnfire_li2016_inst%forc_hdm(nc)  = hdm(nc)
    ! cnfire_li2016_inst%forc_lnfm(nc) = lnfm(nc)
    ! cnfire_li2021_inst%forc_hdm(nc)  = hdm(nc)
    ! cnfire_li2021_inst%forc_lnfm(nc) = lnfm(nc)


     do nz = 1,num_zon    ! CN zone loop
        n = n + 1

        temperature_inst%t_soisno_col(n,-nlevsno+1:nlevmaxurbgrnd) = tp1(nc) ! jkolassa: only one soil and no snow column at this point (may change in future)
        temperature_inst%t_grnd_col(n) = temperature_inst%t_soisno_col(n,1)
        temperature_inst%t_soi17cm_col(n) = temperature_inst%t_grnd_col(n)
        soilstate_inst%soilpsi_col(n,1:nlevgrnd) = 1.e-6*psis(nc)*grav*denh2o*rzm(nc,nz)**(-bee(nc)) ! jkolassa: only one soil layer at this point
        soilstate_inst%psiwilt_col(n,1:nlevgrnd) = 1.e-6*psis(nc)*grav*denh2o*wpwet(nc)**(-bee(nc)) ! jkolassa: soil water potential at wilting point (not a CLM variable, but added to use instead of constant threshold to determine water stress)
        soilstate_inst%watsat_col(n,1:nlevmaxurbgrnd)   = poros(nc) 
        soilstate_inst%bd_col(n,1:nlevmaxurbgrnd)       = (1. - soilstate_inst%watsat_col(n,1:nlevmaxurbgrnd))*2700.
        atm2lnd_inst%forc_t_downscaled_col(n)  = tairm(nc)
        water_inst%wateratm2lndbulk_inst%forc_rain_downscaled_col(n) = rainfm(nc)
        water_inst%wateratm2lndbulk_inst%forc_snow_downscaled_col(n) = snowfm(nc)
        water_inst%waterdiagnosticbulk_inst%wf_col(n)  = sfm(nc,nz)
        water_inst%waterdiagnosticbulk_inst%wf2_col(n) = rzm(nc,nz)
        water_inst%waterdiagnosticbulk_inst%frac_sno_col(n) = fsnow(nc)
        water_inst%waterdiagnosticbulk_inst%snow_depth_col(n) = sndzn(nc)
        water_inst%waterdiagnosticbulk_inst%snow_5day_col(n) = sndzn5d(nc)  
        bgc_vegetation_inst%cnveg_state_inst%gdp_lf_col(n) = gdp(nc)
        bgc_vegetation_inst%cnveg_state_inst%abm_lf_col(n) = abm(nc)
        bgc_vegetation_inst%cnveg_state_inst%peatf_lf_col(n) = peatf(nc)
        water_inst%waterstatebulk_inst%h2osoi_liq_col(n,-nlevsno+1:nlevgrnd) = totwat(nc)
        water_inst%waterfluxbulk_inst%qflx_drain_col(n) = bflow(nc)
        water_inst%waterfluxbulk_inst%qflx_surf_col(n) = runsrf(nc)
        water_inst%waterfluxbulk_inst%AnnET(n) = et365d(nc)*(0.0864*0.408/3600)     ! convert from W m-2 to mm/s

        ! compute column-level saturated area fraction (water table at surface)
        if(nz==1) then
          saturated_excess_runoff_inst%fsat_col(n) = min(max(0.,car1m(nc)/CN_zone_weight(nz)),1.)
        elseif(nz==2) then
          saturated_excess_runoff_inst%fsat_col(n) = min(max(0.,(car1m(nc)-CN_zone_weight(1))/CN_zone_weight(nz)),1.)
        elseif(nz==3) then
          saturated_excess_runoff_inst%fsat_col(n) = min(max(0.,(car1m(nc)-CN_zone_weight(1)-CN_zone_weight(2))/CN_zone_weight(nz)),1.)
        endif

        do np = 0,numpft  ! PFT index loop
           p = p + 1
           temperature_inst%t_ref2m_patch(p) = tairm(nc)
           temperature_inst%soila10_patch(p) = tg10d(nc)
           temperature_inst%t_a5min_patch(p) = t2m5d(nc)
           cn2clm_inst%btran2_patch_cn2clm(p)     = btran_fire(nc,nz)  
          ! cnfire_li2016_inst%cnfire_base_type%btran2_patch(p)     = btran_fire(nc,nz)
          ! cnfire_li2021_inst%cnfire_base_type%btran2_patch(p)     = btran_fire(nc,nz)
           water_inst%wateratm2lndbulk_inst%prec60_patch(p) = prec60d(nc)
           water_inst%wateratm2lndbulk_inst%prec10_patch(p) = prec10d(nc)
           water_inst%wateratm2lndbulk_inst%rh30_patch(p) = rh30(nc)
           frictionvel_inst%forc_hgt_u_patch(p) = 30. ! following CNCLM45 implementation, but this should be available from the GridComp

           do nv = 1,num_veg ! defined veg loop
              if(ityp(nc,nv,nz)==np .and. fveg(nc,nv,nz)>1.e-4) then
                 photosyns_inst%psnsun_patch(p) = psnsunm(nc,nv,nz)
                 photosyns_inst%psnsha_patch(p) = psnsham(nc,nv,nz)
                 photosyns_inst%lmrsun_patch(p) = lmrsunm(nc,nv,nz)
                 photosyns_inst%lmrsha_patch(p) = lmrsham(nc,nv,nz)
                 canopystate_inst%laisun_patch(p) = laisunm(nc,nv,nz)
                 canopystate_inst%laisha_patch(p) = laisham(nc,nv,nz)
              end if 
           end do ! nv
        end do ! np
     end do ! nz
  end do ! nc



  ! call CLM routines that are needed prior to Ecosystem Dynamics call

  call active_layer_inst%alt_calc(filter(1)%num_soilc, filter(1)%soilc, &
       temperature_inst)

  call bgc_vegetation_inst%InitEachTimeStep(bounds, filter(1)%num_soilc, filter(1)%soilc)

  call bgc_vegetation_inst%InitGridcellBalance(bounds, &
       filter(1)%num_allc, filter(1)%allc, &
       filter(1)%num_soilc, filter(1)%soilc, &
       filter(1)%num_soilp, filter(1)%soilp, &
       soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst)

  call bgc_vegetation_inst%InitColumnBalance(bounds, &
       filter(1)%num_allc, filter(1)%allc, &
       filter(1)%num_soilc, filter(1)%soilc, &
       filter(1)%num_soilp, filter(1)%soilp, &
       soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenstate_inst)

  ! Ecosystem Dynamics calculations
  ! jkolassa: This call contains most of the CLM ecosystem dynamics
  ! calculations, including soil biogeochemistry, carbon/nitrogen state and 
  ! flux updates, fire, etc.
 call bgc_vegetation_inst%EcosystemDynamicsPreDrainage(bounds,            &
      filter(1)%num_soilc, filter(1)%soilc,                       &
      filter(1)%num_soilp, filter(1)%soilp,                       &
      filter(1)%num_actfirec, filter(1)%actfirec,                 &
      filter(1)%num_actfirep, filter(1)%actfirep,                 &
      filter(1)%num_pcropp, filter(1)%pcropp, &
      filter(1)%num_exposedvegp, filter(1)%exposedvegp, &
      filter(1)%num_noexposedvegp, filter(1)%noexposedvegp, &
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
      crop_inst, ch4_inst, &
      photosyns_inst, saturated_excess_runoff_inst, energyflux_inst,          &
      nutrient_competition_method, fireemis_inst)


 ! jkolassa: This call is mostly to compute the nitrogen leaching, summary states and fluxes
 ! and the vegetation structural updates
 call bgc_vegetation_inst%EcosystemDynamicsPostDrainage(bounds, &
      filter(1)%num_allc, filter(1)%allc, &
      filter(1)%num_soilc, filter(1)%soilc, &
      filter(1)%num_soilp, filter(1)%soilp, &
      filter(1)%num_actfirec, filter(1)%actfirec,                 &
      filter(1)%num_actfirep, filter(1)%actfirep,                 &
      doalb, crop_inst, &
      soilstate_inst, soilbiogeochem_state_inst, &
      water_inst%waterstatebulk_inst, water_inst%waterdiagnosticbulk_inst, &
      water_inst%waterfluxbulk_inst, frictionvel_inst, canopystate_inst, &
      soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
      c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
      c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
      soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)


! check carbon and nitrogen balances except on first time step
  if(.not.first) then
    call bgc_vegetation_inst%BalanceCheck( &
         bounds, filter(1)%num_soilc, filter(1)%soilc, &
         soilbiogeochem_carbonflux_inst, &
         soilbiogeochem_nitrogenflux_inst, atm2lnd_inst )
  else
    !first = .false.
  end if

      grc%prev_dayl = grc%dayl ! set previous day length for following time steps (dayl itself is computed in GridComp)


! map CLM outputs to Catchment space

  n = 0
  p = 0
  do nc = 1,nch        ! catchment tile loop

     nppg(nc) = 0.
     gppg(nc) = 0.
     srg(nc) = 0. 
     arg(nc) = 0.
     hrg(nc) = 0.
     burn(nc) = 0.
     closs(nc) = 0.
     som_closs(nc) = 0.
     nfire(nc) = 0.
     root(nc) = 0.
     vegc(nc) = 0.
     ndeployg(nc) = 0.
     leafng(nc) = 0.
     leafcg(nc) = 0.
     sminn_to_npoolg(nc) = 0.
     totvegng(nc) = 0.
     retransng(nc) = 0.
     retransn_to_npoolg(nc) = 0.
     rootcg(nc) = 0.
     denitg(nc) = 0.
     sminn_leachedg(nc) = 0.
     sminng(nc) = 0.
     col_fire_nlossg(nc) = 0.
     gross_nming(nc) = 0.
     net_nming(nc) = 0.
     nfix_to_sminng(nc) = 0.
     actual_immobg(nc) = 0.
     fpgg(nc) = 0.
     fpig(nc) = 0.
     sminn_to_plantg(nc) = 0.
     ndep_to_sminng(nc) = 0.
     totlitng(nc) = 0.
     totsomng(nc) = 0.
     fuelcg(nc) = 0.
     totlitcg(nc) = 0.
     cwdcg(nc) = 0.
     xsmr(nc) = 0.

     neeg(nc) = bgc_vegetation_inst%cnveg_carbonflux_inst%nee_grc(nc)

     do nz = 1,num_zon    ! CN zone loop
        n = n + 1

        colc(nc,nz) = bgc_vegetation_inst%cnveg_carbonstate_inst%totc_col(n)
        srg(nc) = srg(nc) + bgc_vegetation_inst%cnveg_carbonflux_inst%sr_col(n)*CN_zone_weight(nz)         
        arg(nc) = arg(nc) + bgc_vegetation_inst%cnveg_carbonflux_inst%ar_col(n)*CN_zone_weight(nz)
        hrg(nc) = hrg(nc) + soilbiogeochem_carbonflux_inst%hr_col(n)*CN_zone_weight(nz)
        burn(nc) = burn(nc) + bgc_vegetation_inst%cnveg_state_inst%farea_burned_col(n)*CN_zone_weight(nz)
        closs(nc) = closs(nc) + bgc_vegetation_inst%cnveg_carbonflux_inst%fire_closs_col(n)*CN_zone_weight(nz)
        som_closs(nc) = som_closs(nc) + soilbiogeochem_carbonflux_inst%somc_fire_col(n)*CN_zone_weight(nz)
        nfire(nc) = nfire(nc) + bgc_vegetation_inst%cnveg_state_inst%nfire_col(n)*CN_zone_weight(nz)
        denitg(nc) = denitg(nc) + soilbiogeochem_nitrogenflux_inst%denit_col(n)*CN_zone_weight(nz)
        sminn_leachedg(nc) = sminn_leachedg(nc) + soilbiogeochem_nitrogenflux_inst%sminn_leached_col(n)*CN_zone_weight(nz)
        sminng(nc) = sminng(nc) + soilbiogeochem_nitrogenstate_inst%sminn_col(n)*CN_zone_weight(nz)
        col_fire_nlossg(nc) = col_fire_nlossg(nc) + bgc_vegetation_inst%cnveg_nitrogenflux_inst%fire_nloss_col(n)*CN_zone_weight(nz)
        gross_nming(nc) = gross_nming(nc) + soilbiogeochem_nitrogenflux_inst%gross_nmin_col(n)*CN_zone_weight(nz)
        net_nming(nc) = net_nming(nc) + soilbiogeochem_nitrogenflux_inst%net_nmin_col(n)*CN_zone_weight(nz)
        nfix_to_sminng(nc) = nfix_to_sminng(nc) + soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col(n)*CN_zone_weight(nz)
        actual_immobg(nc) = actual_immobg(nc) + soilbiogeochem_nitrogenflux_inst%actual_immob_col(n)*CN_zone_weight(nz)
        fpgg(nc) = fpgg(nc) + soilbiogeochem_state_inst%fpg_col(n)*CN_zone_weight(nz)
        fpig(nc) = fpig(nc) + soilbiogeochem_state_inst%fpi_col(n)*CN_zone_weight(nz)
        sminn_to_plantg(nc) = sminn_to_plantg(nc) + soilbiogeochem_nitrogenflux_inst%sminn_to_plant_col(n)*CN_zone_weight(nz)
        ndep_to_sminng(nc) = ndep_to_sminng(nc) + soilbiogeochem_nitrogenflux_inst%ndep_to_sminn_col(n)*CN_zone_weight(nz)
        totlitng(nc) = totlitng(nc) + soilbiogeochem_nitrogenstate_inst%totlitn_col(n)*CN_zone_weight(nz)
        totsomng(nc) = totsomng(nc) + soilbiogeochem_nitrogenstate_inst%totsomn_col(n)*CN_zone_weight(nz)
        fuelcg(nc) = fuelcg(nc) + bgc_vegetation_inst%cnveg_carbonstate_inst%fuelc_col(n)*CN_zone_weight(nz)
        totlitcg(nc) = totlitcg(nc) + soilbiogeochem_carbonstate_inst%totlitc_col(n)*CN_zone_weight(nz)
        cwdcg(nc) = cwdcg(nc) + soilbiogeochem_carbonstate_inst%cwdc_col(n)*CN_zone_weight(nz)

        do np = 0,numpft  ! PFT index loop
           p = p + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==np .and. fveg(nc,nv,nz)>1.e-4) then

                  zlai(nc,nv,nz) = canopystate_inst%elai_patch(p)
                  zsai(nc,nv,nz) = canopystate_inst%esai_patch(p)
                  ztai(nc,nv,nz) = canopystate_inst%tlai_patch(p)

                  pwtgcell = fveg(nc,nv,nz)*CN_zone_weight(nz) ! PFT weight in catchment tile
                  nppg(nc) = nppg(nc) + bgc_vegetation_inst%cnveg_carbonflux_inst%npp_patch(p)*pwtgcell
                  gppg(nc) = gppg(nc) + bgc_vegetation_inst%cnveg_carbonflux_inst%gpp_patch(p)*pwtgcell
                  root(nc) = root(nc) + (bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_patch(p) &
                                       + bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_storage_patch(p) &
                                       + bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_xfer_patch(p) &
                                        )*pwtgcell
                  vegc(nc) = vegc(nc) + bgc_vegetation_inst%cnveg_carbonstate_inst%totvegc_patch(p)*pwtgcell
                  ndeployg(nc) = ndeployg(nc) + bgc_vegetation_inst%cnveg_nitrogenflux_inst%ndeploy_patch(p)*pwtgcell
                  leafng(nc) = leafng(nc) + bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch(p)*pwtgcell
                  leafcg(nc) = leafcg(nc) + bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch(p)*pwtgcell
                  sminn_to_npoolg(nc) = sminn_to_npoolg(nc) + bgc_vegetation_inst%cnveg_nitrogenflux_inst%sminn_to_npool_patch(p)*pwtgcell
                  totvegng(nc) = totvegng(nc) + bgc_vegetation_inst%cnveg_nitrogenstate_inst%totvegn_patch(p)*pwtgcell
                  retransng(nc) = retransng(nc) + bgc_vegetation_inst%cnveg_nitrogenstate_inst%retransn_patch(p)*pwtgcell
                  retransn_to_npoolg(nc) = retransn_to_npoolg(nc) + bgc_vegetation_inst%cnveg_nitrogenflux_inst%retransn_to_npool_patch(p)*pwtgcell
                  rootcg(nc) = rootcg(nc) + (bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_storage_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_xfer_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_storage_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_xfer_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%deadcrootc_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%deadcrootc_storage_patch(p) &
                                          + bgc_vegetation_inst%cnveg_carbonstate_inst%deadcrootc_xfer_patch(p) &
                                          )*pwtgcell
            
                  xsmr(nc) = xsmr(nc) + bgc_vegetation_inst%cnveg_carbonstate_inst%xsmrpool_patch(p)*pwtgcell
                end if
             end do ! nv
        end do !np
     end do ! nz
  end do ! nc

 end subroutine CN_Driver

!------------------------------------------------
 subroutine CN_exit(nch,ityp,fveg,cncol,cnpft)

  ! INPUT
  integer,                                 intent(in) :: nch  ! number of tiles
  integer, dimension(nch,num_veg,num_zon), intent(in) :: ityp ! PFT index
  real, dimension(nch,num_veg,num_zon),    intent(in) :: fveg ! PFT fraction

  ! OUTPUT
  real, dimension(nch,num_zon,var_col),         intent(out) :: cncol ! column-level restart variables 
  real, dimension(nch,num_zon,num_veg,var_pft), intent(out) :: cnpft ! PFT-level restart variables

  ! LOCAL

! type(bounds_type)                      :: bounds
! type(soilbiogeochem_carbonflux_type)   :: soilbiogeochem_carbonflux_inst
! type(soilbiogeochem_nitrogenflux_type) :: soilbiogeochem_nitrogenflux_inst
! type(gridcell_type)                    :: grc
! type(cn_vegetation_type)               :: bgc_vegetation_inst
!! type(cnveg_state_type)                 :: cnveg_state_inst
! type(soilbiogeochem_carbonstate_type)  :: soilbiogeochem_carbonstate_inst
! type(soilbiogeochem_nitrogenstate_type):: soilbiogeochem_nitrogenstate_inst
! type(soilbiogeochem_state_type)        :: soilbiogeochem_state_inst
! type(cnveg_carbonflux_type)            :: cnveg_carbonflux_inst
! type(cnveg_carbonstate_type)           :: cnveg_carbonstate_inst
! type(cnveg_nitrogenflux_type)          :: cnveg_nitrogenflux_inst
! type(cnveg_nitrogenstate_type)          :: cnveg_nitrogenstate_inst
! type(cn_products_type)                  :: c_products_inst
! type(cn_products_type)                  :: n_products_inst

  integer :: n, p, nv, nc, nz, np, nd
  integer, dimension(8) :: decomp_cpool_cncol_index = (/ 3, 4, 5, 2, 10, 11, 12, 13 /)
  integer, dimension(8) :: decomp_npool_cncol_index = (/ 18, 19, 20, 17,25, 26, 27, 28 /)
  !----------------

  n = 0
  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,num_zon    ! CN zone loop
      n = n + 1

      cncol(nc,nz, 1) = soilbiogeochem_carbonstate_inst%ctrunc_vr_col(n,1)

      do nd = 1,ndecomp_pools
         ! jkolassa: accounting for fact that pool order in CNCOL is different from CTSM
         cncol(nc,nz,decomp_cpool_cncol_index(nd)) = soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col   (n,1,nd)
         cncol(nc,nz,decomp_npool_cncol_index(nd)) = soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col (n,1,nd)
      end do

      cncol(nc,nz, 6) = bgc_vegetation_inst%cnveg_carbonstate_inst%totvegc_col              (n)
      ! jkolassa: variables below transitioned from being column-level to being gridcell-level in CLM;
      !           assuming here that quantities are spread over zones according to zone weight
      cncol(nc,nz, 7) = bgc_vegetation_inst%c_products_inst%prod100_grc(nc)*CN_zone_weight(nz) 
      cncol(nc,nz, 8) = bgc_vegetation_inst%c_products_inst%prod100_grc(nc)*CN_zone_weight(nz)
      cncol(nc,nz, 9) = bgc_vegetation_inst%cnveg_carbonstate_inst%seedc_grc(nc)*CN_zone_weight(nz)
      cncol(nc,nz,14) = bgc_vegetation_inst%cnveg_carbonstate_inst%totc_col                 (n)
      cncol(nc,nz,15) = soilbiogeochem_carbonstate_inst%totlitc_col     (n)
      cncol(nc,nz,16) = soilbiogeochem_nitrogenstate_inst%ntrunc_vr_col (n,1)


      cncol(nc,nz,21) = bgc_vegetation_inst%n_products_inst%prod100_grc(nc)*CN_zone_weight(nz) 
      cncol(nc,nz,22) = bgc_vegetation_inst%n_products_inst%prod100_grc(nc)*CN_zone_weight(nz)
      cncol(nc,nz,23) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%seedn_grc(nc)*CN_zone_weight(nz)
      cncol(nc,nz,24) = soilbiogeochem_nitrogenstate_inst%sminn_vr_col   (n,1)
      cncol(nc,nz,29) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%totn_col                (n)
      cncol(nc,nz,30) = soilbiogeochem_state_inst%fpg_col                (n)
      cncol(nc,nz,31) = bgc_vegetation_inst%cnveg_state_inst%annsum_counter_col              (n)
      cncol(nc,nz,32) = bgc_vegetation_inst%cnveg_state_inst%annavg_t2m_col                  (n)
      cncol(nc,nz,33) = bgc_vegetation_inst%cnveg_carbonflux_inst%annsum_npp_col             (n)
      cncol(nc,nz,34) = bgc_vegetation_inst%cnveg_state_inst%farea_burned_col                (n)
      cncol(nc,nz,35) = soilbiogeochem_state_inst%fpi_col                (n)
      cncol(nc,nz,36) = soilbiogeochem_nitrogenstate_inst%smin_no3_col   (n)
      cncol(nc,nz,37) = soilbiogeochem_nitrogenstate_inst%smin_nh4_col   (n)


      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,num_veg ! defined veg loop

          if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
            cnpft(nc,nz,nv,  1) = bgc_vegetation_inst%cnveg_carbonstate_inst%cpool_patch                 (np)
            cnpft(nc,nz,nv,  2) = bgc_vegetation_inst%cnveg_carbonstate_inst%deadcrootc_patch            (np)
            cnpft(nc,nz,nv,  3) = bgc_vegetation_inst%cnveg_carbonstate_inst%deadcrootc_storage_patch    (np)
            cnpft(nc,nz,nv,  4) = bgc_vegetation_inst%cnveg_carbonstate_inst%deadcrootc_xfer_patch       (np)
            cnpft(nc,nz,nv,  5) = bgc_vegetation_inst%cnveg_carbonstate_inst%deadstemc_patch             (np)
            cnpft(nc,nz,nv,  6) = bgc_vegetation_inst%cnveg_carbonstate_inst%deadstemc_storage_patch     (np)
            cnpft(nc,nz,nv,  7) = bgc_vegetation_inst%cnveg_carbonstate_inst%deadstemc_xfer_patch        (np)
            cnpft(nc,nz,nv,  8) = bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_patch                (np)
            cnpft(nc,nz,nv,  9) = bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_storage_patch        (np)
            cnpft(nc,nz,nv, 10) = bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_xfer_patch           (np)
            cnpft(nc,nz,nv, 11) = bgc_vegetation_inst%cnveg_carbonstate_inst%gresp_storage_patch         (np)
            cnpft(nc,nz,nv, 12) = bgc_vegetation_inst%cnveg_carbonstate_inst%gresp_xfer_patch            (np)
            cnpft(nc,nz,nv, 13) = bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch                 (np)
            cnpft(nc,nz,nv, 14) = bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_storage_patch         (np)
            cnpft(nc,nz,nv, 15) = bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_xfer_patch            (np)
            cnpft(nc,nz,nv, 16) = bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_patch            (np)
            cnpft(nc,nz,nv, 17) = bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_storage_patch    (np)
            cnpft(nc,nz,nv, 18) = bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_xfer_patch       (np)
            cnpft(nc,nz,nv, 19) = bgc_vegetation_inst%cnveg_carbonstate_inst%livestemc_patch             (np)
            cnpft(nc,nz,nv, 20) = bgc_vegetation_inst%cnveg_carbonstate_inst%livestemc_storage_patch     (np)
            cnpft(nc,nz,nv, 21) = bgc_vegetation_inst%cnveg_carbonstate_inst%livestemc_xfer_patch        (np)
            cnpft(nc,nz,nv, 22) = bgc_vegetation_inst%cnveg_carbonstate_inst%ctrunc_patch                (np)
            cnpft(nc,nz,nv, 23) = bgc_vegetation_inst%cnveg_carbonstate_inst%xsmrpool_patch              (np)
            cnpft(nc,nz,nv, 24) = bgc_vegetation_inst%cnveg_state_inst%annavg_t2m_patch                  (np)
            cnpft(nc,nz,nv, 25) = bgc_vegetation_inst%cnveg_state_inst%annmax_retransn_patch             (np)
            cnpft(nc,nz,nv, 26) = bgc_vegetation_inst%cnveg_carbonflux_inst%annsum_npp_patch             (np)
            cnpft(nc,nz,nv, 27) = bgc_vegetation_inst%cnveg_state_inst%annsum_potential_gpp_patch        (np)
            cnpft(nc,nz,nv, 28) = grc%dayl                                           (nc) ! jkolassa Dec 2023: dayl is a gridcell =-level variable in CLM, but is stored as patch-level variable in CatcCN restart
            cnpft(nc,nz,nv, 29) = bgc_vegetation_inst%cnveg_state_inst%days_active_patch                 (np)
            cnpft(nc,nz,nv, 30) = bgc_vegetation_inst%cnveg_state_inst%dormant_flag_patch                (np)
            cnpft(nc,nz,nv, 31) = bgc_vegetation_inst%cnveg_state_inst%offset_counter_patch              (np)
            cnpft(nc,nz,nv, 32) = bgc_vegetation_inst%cnveg_state_inst%offset_fdd_patch                  (np)
            cnpft(nc,nz,nv, 33) = bgc_vegetation_inst%cnveg_state_inst%offset_flag_patch                 (np)
            cnpft(nc,nz,nv, 34) = bgc_vegetation_inst%cnveg_state_inst%offset_swi_patch                  (np)
            cnpft(nc,nz,nv, 35) = bgc_vegetation_inst%cnveg_state_inst%onset_counter_patch               (np)
            cnpft(nc,nz,nv, 36) = bgc_vegetation_inst%cnveg_state_inst%onset_fdd_patch                   (np)
            cnpft(nc,nz,nv, 37) = bgc_vegetation_inst%cnveg_state_inst%onset_flag_patch                  (np)
            cnpft(nc,nz,nv, 38) = bgc_vegetation_inst%cnveg_state_inst%onset_gdd_patch                   (np)
            cnpft(nc,nz,nv, 39) = bgc_vegetation_inst%cnveg_state_inst%onset_gddflag_patch               (np)
            cnpft(nc,nz,nv, 40) = bgc_vegetation_inst%cnveg_state_inst%onset_swi_patch                   (np)
            cnpft(nc,nz,nv, 41) = bgc_vegetation_inst%cnveg_carbonflux_inst%prev_frootc_to_litter_patch  (np)
            cnpft(nc,nz,nv, 42) = bgc_vegetation_inst%cnveg_carbonflux_inst%prev_leafc_to_litter_patch   (np)
            cnpft(nc,nz,nv, 43) = bgc_vegetation_inst%cnveg_state_inst%tempavg_t2m_patch                 (np)
            cnpft(nc,nz,nv, 44) = bgc_vegetation_inst%cnveg_state_inst%tempmax_retransn_patch            (np)
            cnpft(nc,nz,nv, 45) = bgc_vegetation_inst%cnveg_carbonflux_inst%tempsum_npp_patch            (np)
            cnpft(nc,nz,nv, 46) = bgc_vegetation_inst%cnveg_state_inst%tempsum_potential_gpp_patch       (np)
            cnpft(nc,nz,nv, 47) = bgc_vegetation_inst%cnveg_carbonflux_inst%xsmrpool_recover_patch       (np)
            cnpft(nc,nz,nv, 48) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadcrootn_patch          (np)
            cnpft(nc,nz,nv, 49) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadcrootn_storage_patch  (np)
            cnpft(nc,nz,nv, 50) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadcrootn_xfer_patch     (np)
            cnpft(nc,nz,nv, 51) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadstemn_patch           (np)
            cnpft(nc,nz,nv, 52) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadstemn_storage_patch   (np)
            cnpft(nc,nz,nv, 53) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadstemn_xfer_patch      (np)
            cnpft(nc,nz,nv, 54) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%frootn_patch              (np)
            cnpft(nc,nz,nv, 55) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%frootn_storage_patch      (np)
            cnpft(nc,nz,nv, 56) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%frootn_xfer_patch         (np)
            cnpft(nc,nz,nv, 57) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch               (np)
            cnpft(nc,nz,nv, 58) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_storage_patch       (np)
            cnpft(nc,nz,nv, 59) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_xfer_patch          (np)
            cnpft(nc,nz,nv, 60) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%livecrootn_patch          (np)
            cnpft(nc,nz,nv, 61) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%livecrootn_storage_patch  (np)
            cnpft(nc,nz,nv, 62) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%livecrootn_xfer_patch     (np)
            cnpft(nc,nz,nv, 63) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%livestemn_patch           (np)
            cnpft(nc,nz,nv, 64) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%livestemn_storage_patch   (np)
            cnpft(nc,nz,nv, 65) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%livestemn_xfer_patch      (np)
            cnpft(nc,nz,nv, 66) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%npool_patch               (np)
            cnpft(nc,nz,nv, 67) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%ntrunc_patch              (np)
            cnpft(nc,nz,nv, 68) = bgc_vegetation_inst%cnveg_nitrogenstate_inst%retransn_patch            (np)
            cnpft(nc,nz,nv, 69) = canopystate_inst%elai_patch                        (np)
            cnpft(nc,nz,nv, 70) = canopystate_inst%esai_patch                        (np)
            cnpft(nc,nz,nv, 71) = canopystate_inst%hbot_patch                        (np)
            cnpft(nc,nz,nv, 72) = canopystate_inst%htop_patch                        (np)
            cnpft(nc,nz,nv, 73) = canopystate_inst%tlai_patch                        (np)
            cnpft(nc,nz,nv, 74) = canopystate_inst%tsai_patch                        (np)
            cnpft(nc,nz,nv, 75) = bgc_vegetation_inst%cnveg_nitrogenflux_inst%plant_ndemand_patch        (np)
            cnpft(nc,nz,nv, 76) = canopystate_inst%vegwp_patch                       (np,1)
            cnpft(nc,nz,nv, 77) = canopystate_inst%vegwp_patch                       (np,2)
            cnpft(nc,nz,nv, 78) = canopystate_inst%vegwp_patch                       (np,3)
            cnpft(nc,nz,nv, 79) = canopystate_inst%vegwp_patch                       (np,4)
            cnpft(nc,nz,nv, 80) = bgc_vegetation_inst%cnveg_carbonflux_inst%annsum_litfall_patch         (np)
            cnpft(nc,nz,nv, 81) = bgc_vegetation_inst%cnveg_carbonflux_inst%tempsum_litfall_patch        (np)
          endif

        end do ! defined veg loop
      end do   ! PFT index loop
    end do     ! CN zone loop
  end do       ! catchment tile loop

  return

 end subroutine CN_exit 

!--------------------------
  subroutine get_CN_LAI(nch,ityp,fveg,elai,esai,tlai,tsai)

 ! ARGUMENTS

  ! INPUT/OUTPUT
  integer, intent(in) :: nch ! number of tiles
  integer, dimension(nch,num_veg,num_zon), intent(in) :: ityp ! PFT index
  real, dimension(nch,num_veg,num_zon), intent(in) :: fveg    ! PFT fraction
  real, dimension(nch,num_veg,num_zon), intent(out)           :: elai   ! exposed leaf-area index
  real, dimension(nch,num_veg,num_zon), intent(out), optional :: esai   ! exposed stem-area index
  real, dimension(nch,num_veg,num_zon), intent(out), optional :: tlai   ! total leaf-area index
  real, dimension(nch,num_veg,num_zon), intent(out), optional :: tsai   ! total stem-area index

  ! LOCAL
  integer :: n, p, nv, nc, nz, np

!  real(r8), pointer :: elai_clm(:)
!  real(r8), pointer :: esai_clm(:)
!  real(r8), pointer :: tlai_clm(:)
!  real(r8), pointer :: tsai_clm(:)

  !------------------------------
 associate(&
  elai_clm => canopystate_inst%elai_patch , &
  esai_clm => canopystate_inst%esai_patch , &
  tlai_clm => canopystate_inst%tlai_patch , & 
  tsai_clm => canopystate_inst%tsai_patch   &
 )

                    elai = 0.
  if(present(esai)) esai = 0.
  if(present(tlai)) tlai = 0.
  if(present(tsai)) tsai = 0.

  n = 0
  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,num_zon  ! CN zone loop
      n = n + 1
      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,num_veg ! defined veg loop

! extract LAI & SAI from CN clmtype
! ---------------------------------
          if(ityp(nc,nv,nz)==p .and. ityp(nc,nv,nz)>0 .and. fveg(nc,nv,nz)>1.e-4) then
                              elai(nc,nv,nz) = elai_clm(np)
            if(present(esai)) esai(nc,nv,nz) = esai_clm(np)
            if(present(tlai)) tlai(nc,nv,nz) = tlai_clm(np)
            if(present(tsai)) tsai(nc,nv,nz) = tsai_clm(np)
          endif

        end do ! defined veg loop
      end do   ! PFT index loop
    end do     ! CN zone loop
  end do       ! catchment tile loop

  end associate

  end subroutine get_CN_LAI
!---------------------------

! subroutine FireMethodInit(bounds,paramfile)
!
!  use MAPL             , only : NetCDF4_FileFormatter
!
!
!  type(bounds_type), intent(in)    :: bounds
!  character(300), intent(in)     :: paramfile
!
!  type(Netcdf4_fileformatter) :: ncid
!  integer            :: rc, status
!  !--------------------------------
!
!  call create_cnfire_method(cnfire_method)
!  call cnfire_method%FireInit(bounds)
!
!  call ncid%open(trim(paramfile),pFIO_READ, __RC__)
!  call cnfire_method%CNFireReadParams( ncid )
!  call ncid%close(rc=status)
!
!  end subroutine FireMethodInit
end module CNCLM_DriverMod
