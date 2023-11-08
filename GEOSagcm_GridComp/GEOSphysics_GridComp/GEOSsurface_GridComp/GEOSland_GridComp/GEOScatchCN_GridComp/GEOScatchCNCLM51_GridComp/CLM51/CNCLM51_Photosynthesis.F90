 module CNCLM_Photosynthesis

 use MAPL_ConstantsMod
 use clm_varpar,         only : numpft, numrad, num_veg, num_zon, &
                                nlevcan
 use decompMod
 use PatchType
 use filterMod
 
 use CNVegNitrogenstateType
 use CNVegCarbonstateType
 use atm2lndType
 use TemperatureType
 use SoilStateType
 use pftconMod
 use WaterDiagnosticBulkType
 use SurfaceAlbedoType
 use SolarAbsorbedType
 use CanopyStateType
 use OzoneBaseMod
 use PhotosynthesisMod
 use WaterFluxBulkType
 use WaterStateType
 use WaterType
 use CNVegetationFacade

 implicit none

 private
 public catchcn_calc_rc

 contains

!---------------------------------------------------
 subroutine catchcn_calc_rc(nch,fveg,tc,qa,pbot,co2v,dayl_factor, &
            t10,tm,cond,psis,wet3,bee,capac,fwet,coszen,ityp,&
            pardir,pardif,albdir,albdif,dtc,dea,water_inst,bgc_vegetation_inst,rc,rc_dea,rc_dt,&
            laisun_out,laisha_out,psnsun_out,psnsha_out,lmrsun_out,&
            lmrsha_out,parabs,btran_out)

 use MAPL_SatVaporMod
 use QSatMod              , only : QSat
 use SurfaceAlbedoMod     , only : TwoStream
 use SurfaceRadiationMod  , only : CanopySunShadeFracs
 ! INPUTS
 integer, intent(in) :: nch               ! vector length

 real, dimension(nch,num_veg,num_zon), intent(in) :: fveg ! catchment vegetation fractions  
 real, intent(in) :: tc(nch,num_zon)              ! canopy temperature (K)
 real, intent(in) :: qa(nch,num_zon)              ! canopy air specific humidity (kg/kg)
 real, intent(in) :: pbot(nch)            ! surface pressure (Pa)
 real, intent(in) :: co2v(nch)            ! atmospheric carbon dioxide concentration
 real, intent(in) :: dayl_factor(nch)     ! daylength factor (0-1)
 real, intent(in) :: t10(nch)             ! 10-day "running mean" of the 2 m temperature (K)
 real, intent(in) :: tm(nch)              ! air temperature at agcm reference height (K)
 real, intent(in) :: cond(nch)            ! saturated hydraulic conductivity (m/s)
 real, intent(in) :: psis(nch)            ! saturated matric potential [m]
 real, intent(in) :: wet3(nch)            ! average soil profile wetness [-]
 real, intent(in) :: bee(nch)             ! Clapp-Hornberger 'b' [-]
 real, intent(in) :: capac(nch)           ! interception reservoir capacity [kg m^-2] 
 real, intent(in) :: fwet(nch)            ! fraction of canopy that is wet (0-1)
 real, intent(in) :: coszen(nch)          ! cosine solar zenith angle
 integer, intent(in) :: ityp(nch,num_veg,num_zon)    ! canopy vegetation index (PFT) 
 real, intent(in) :: pardir(nch)          ! direct PAR (W/m2)
 real, intent(in) :: pardif(nch)          ! diffuse PAR (W/m2)
 real, intent(in) :: albdir(nch,num_veg,num_zon,numrad)          ! direct albedo
 real, intent(in) :: albdif(nch,num_veg,num_zon,numrad)          ! diffuse albedo
 real, intent(in) :: dtc ! canopy temperature perturbation (K) [approx 1:10000]
 real, intent(in) :: dea ! vapor pressure perturbation (Pa) [approx 1:10000]
 type(water_type),intent(in) :: water_inst
 type(cn_vegetation_type), intent(in) :: bgc_vegetation_inst

 ! OUTPUTS
 real, dimension(nch,num_zon), intent(out) :: rc       ! unperturbed canopy stomatal resistance [s/m]
 real, dimension(nch,num_zon), intent(out) :: rc_dea   ! canopy stomatal resistance with vapor pressure pertubation [s/m]
 real, dimension(nch,num_zon), intent(out) :: rc_dt    ! canopy stomatal resistance with canopy temperature pertubation [s/m]
 real, dimension(nch,num_veg,num_zon), intent(out) :: laisun_out
 real, dimension(nch,num_veg,num_zon), intent(out) :: laisha_out
 real, dimension(nch,num_veg,num_zon), intent(out) :: psnsun_out
 real, dimension(nch,num_veg,num_zon), intent(out) :: psnsha_out
 real, dimension(nch,num_veg,num_zon), intent(out) :: lmrsun_out
 real, dimension(nch,num_veg,num_zon), intent(out) :: lmrsha_out
 real, dimension(nch,num_veg,num_zon), intent(out) :: parabs
 real, dimension(nch,num_veg,num_zon), intent(out) :: btran_out

! LOCAL

 ! CLM variables
! type(bounds_type)              :: bounds
! type(atm2lnd_type)             :: atm2lnd_inst
! type(temperature_type)         :: temperature_inst
! type(soilstate_type)           :: soilstate_inst
! type(waterdiagnosticbulk_type) :: waterdiagnosticbulk_inst
! type(surfalb_type)             :: surfalb_inst
! type(solarabs_type)            :: solarabs_inst
! type(canopystate_type)         :: canopystate_inst
! type(ozone_base_type)          :: ozone_inst
! type(photosyns_type)           :: photosyns_inst
!  type(waterfluxbulk_type)       :: waterfluxbulk_inst
! type(cnveg_nitrogenstate_type) :: cnveg_nitrogenstate_inst
! type(cnveg_carbonstate_type)   :: cnveg_carbonstate_inst
! type(waterstate_type)          :: waterstate_inst
! type(clumpfilter)              :: filter

 ! temporary and loop variables                                                                                        
 integer :: n, p, pft_num, nv, nc, nz, np, ib, nl, iv
 real    :: bare, tmp_albgrd_vis,tmp_albgrd_nir,&
            tmp_albgri_vis,tmp_albgri_nir, &
            tmp_parsun, tmp_parsha

 ! filter variables
 integer, allocatable, save :: filter_vegsol(:), filter_novegsol(:)
 integer                    :: num_vegsol, num_novegsol

 ! constants and parameters 
 real :: rair = MAPL_RDRY
 real :: extkn = 0.30_r8    ! nitrogen allocation coefficient
 integer, parameter :: npft = numpft+1

 ! local variables for stomatal resistance calculations
 real                                    :: rs, rs_dea, rs_dt, rcs, rcs_dea, rcs_dt
 real, dimension(nch*NUM_ZON*(numpft+1)) :: laisun, laisha, rssun, rssha
 real, dimension(nch*NUM_ZON*(numpft+1)) :: laisun_dea, laisha_dea, rssun_dea, rssha_dea
 real, dimension(nch*NUM_ZON*(numpft+1)) :: laisun_dt, laisha_dt, rssun_dt, rssha_dt

 ! local variables to compute Photosynthesis inputs
 real                           :: ws, wl
 real(r8), allocatable, dimension(:,:) :: rho, tau
 real, dimension (nch, NUM_ZON) :: esat_tv     ! vapor pressure inside leaf (sat vapor press at tc) (Pa)
 real, dimension (nch, NUM_ZON) :: eair        ! vapor pressure of canopy air
 real, dimension (nch)          :: oair        ! Atmospheric O2 partial pressure (Pa)
 real, dimension (nch)          :: deldT       ! d(es)/d(T)
 real, dimension (nch)          :: cair        ! compute CO2 partial pressure
 real(r8), dimension (nch)          :: rb          ! boundary layer resistance (s/m)
 real(r8), dimension (nch)          :: el          ! vapor pressure on leaf surface [pa]
 real(r8), dimension (nch, NUM_ZON) :: qsatl       ! leaf specific humidity [kg/kg]
 real(r8), dimension (nch, NUM_ZON) :: qsatldT     ! derivative of "qsatl" on "t_veg"
 real, dimension (nch, NUM_ZON) :: qaf         ! canopy air humidity [kg/kg]
 real(r8), dimension(nch,num_zon)   :: tc_in
 real(r8), dimension(nch)           :: pbot_in

 ! local inputs to Photosynthesis in CLM space
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: coszen_clm ! cosine solar zenith angle for next time step in CLM dimensions
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: esat_tv_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: eair_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: cair_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: oair_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: rb_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: dayl_factor_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: qsatl_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: qaf_clm
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: deldT_clm

 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: eair_pert
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: esat_tv_pert
 real(r8), dimension(nch*NUM_ZON*(numpft+1)) :: temp_unpert

 ! local pointers for Photosynthesis inputs
 real, pointer :: leafn(:)    ! leaf N (gN/m2)   
 real, pointer :: froot_carbon(:) ! fine root carbon (gC/m2) [pft]
 real, pointer :: croot_carbon(:) ! live coarse root carbon (gC/m2) [pft] 
 integer, pointer :: filter_nourbanp
 integer, pointer :: filter_num_nourbanp
 integer, pointer :: filter_exposedvegp
 integer, pointer :: filter_num_exposedvegp

 ! local outputs from Photosynthesis routine
  real(r8)  , allocatable, dimension(:)   :: bsun        ! sunlit canopy transpiration wetness factor (0 to 1)
  real(r8)  , allocatable, dimension(:)   :: bsha       ! shaded canopy transpiration wetness factor (0 to 1)
  real(r8)  , allocatable, dimension(:)   :: btran         ! transpiration wetness factor (0 to 1) [pft]

 ! associate variables

 associate(&
       rhol                    => pftcon%rhol                     , & ! Input:  leaf reflectance: 1=vis, 2=nir        
       rhos                    => pftcon%rhos                     , & ! Input:  stem reflectance: 1=vis, 2=nir        
       taul                    => pftcon%taul                     , & ! Input:  leaf transmittance: 1=vis, 2=nir      
       taus                    => pftcon%taus                     , & ! Input:  stem transmittance: 1=vis, 2=nir   
       xl                      => pftcon%xl                       , & 
       leafn                   => bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch , &
       froot_carbon            => bgc_vegetation_inst%cnveg_carbonstate_inst%frootc_patch  , &
       croot_carbon            => bgc_vegetation_inst%cnveg_carbonstate_inst%livecrootc_patch, &
       elai                    => canopystate_inst%elai_patch      , &
       esai                    => canopystate_inst%esai_patch      ,  &
       filter_nourbanp         => filter(1)%nourbanp               , &
       filter_num_nourbanp     => filter(1)%num_nourbanp           , &
       filter_exposedvegp      => filter(1)%exposedvegp            , &
       filter_num_exposedvegp  => filter(1)%num_exposedvegp          &
        )

! allocate filters
!-----------------------------

 allocate (filter_vegsol(bounds%endp-bounds%begp+1))
 allocate (filter_novegsol(bounds%endp-bounds%begp+1))
 num_vegsol   = 0
 num_novegsol = 0

! allocate variables for radiation calculations 
!---------------------------------

 allocate(rho(bounds%begp:bounds%endp,numrad))
 allocate(tau(bounds%begp:bounds%endp,numrad))

! allocate Photosynthesis outputs
!--------------------------------

 allocate(bsun(bounds%begp:bounds%endp))
 allocate(bsha(bounds%begp:bounds%endp))
 allocate(btran(bounds%begp:bounds%endp))

! compute saturation vapor pressure
! ---------------------------------
   do n = 1,nch
      do nz = 1,NUM_ZON
         esat_tv(n,nz) = MAPL_EQsat(tc(n,nz),DQ=deldT(n))
      end do
   end do

 ! compute canopy air vapor pressure
 !----------------------------------
   do n = 1,nch
      do nz = 1,NUM_ZON
         eair(n,nz) = pbot(n) * qa(n,nz) / (0.622 + qa(n,nz))  ! canopy air vapor pressure (Pa);  jk: this is different from the formulation in the CLM code, which is different from the formulation in the CLM documentation
      end do
   end do
 ! compute atmospheric O2 partial pressure
 !-----------------------------------------
  oair(:) = 0.20946*pbot

 ! compute CO2 partial pressure constant ratio [internal leaf CO2 partial pressure]
 !-------------------------------
  cair(:) = co2v(:)*pbot

 ! leaf boundary layer resistance
 !--------------------------------
  rb = 10.    ! jk: in the original Catchment-CN this was arbitrarily set to 10 by gkw, not sure why

 ! leaf specific humidity
 !------------------------

  tc_in      =  tc
  pbot_in    =  pbot

 do n = 1,nch
    do nz = 1,NUM_ZON
       call QSat(tc_in(n,nz), pbot_in(n), qsatl(n,nz), &
                 el(n), &
                 qsatldT(n,nz))
    end do
 end do

 ! canopy air humidity
 !--------------------

 do n = 1,nch
    do nz = 1,NUM_ZON
       qaf(n,nz) = qa(n,nz)
    end do
 end do

 ! atmospheric pressure and density downscaled to column level
 ! vegetation temperature, 2m 10-day running mean temperature, temperature at AGCM ref. height
 !------------------------------------------------
 p = 0
 n = 0

 do nc = 1,nch
    atm2lnd_inst%forc_solad_grc (nc,1)  = pardir(nc) 
    atm2lnd_inst%forc_solai_grc (nc,1)  = pardif(nc)
    do nz = 1,num_zon
       n = n + 1
       atm2lnd_inst%forc_pbot_downscaled_col (n)  = pbot(nc)
       atm2lnd_inst%forc_rho_downscaled_col  (n)  = (pbot(nc)-0.378*eair(nc,nz))/(rair*tc(nc,nz)) 

       soilstate_inst%hksat_col (n,1:nlevgrnd) = 1000.*COND(nc)                          ! saturated hydraulic conductivity mapped to CLM space
                                                                               ! and converted to [mm/s]
       soilstate_inst%hk_l_col   (n,1:nlevgrnd) = 1000.*COND(nc)*(wet3(nc)**(2*bee(nc)+3)) ! actual hydraulic conductivity mapped to CLM space
                                                                               ! and converted to [mm/s]
       soilstate_inst%smp_l_col  (n,1:nlevgrnd) = 1000.*PSIS(nc)*(max(1.e-06_r8,wet3(nc))**(-bee(nc)))    ! actual soil matric potential mapped to CLM space 
                                                                               ! and converted to [mm]
       soilstate_inst%bsw_col    (n,1:nlevgrnd) = bee(nc)                                 ! Clapp-Hornberger 'b'
       soilstate_inst%sucsat_col (n,1:nlevgrnd) = 1000.*psis(nc)*(-1)                     ! minimum soil suction [mm]

       ! compute column level direct and diffuse albedos (vis and nir) from pft level quantities
       tmp_albgrd_vis = 0.
       tmp_albgrd_nir = 0.
       tmp_albgri_vis = 0.
       tmp_albgri_nir = 0.

       do nv = 1,num_veg
        tmp_albgrd_vis = tmp_albgrd_vis + albdir(nc,nv,nz,1)*fveg(nc,nv,nz)
        tmp_albgrd_nir = tmp_albgrd_nir + albdir(nc,nv,nz,2)*fveg(nc,nv,nz)

        tmp_albgri_vis = tmp_albgri_vis + albdif(nc,nv,nz,1)*fveg(nc,nv,nz)
        tmp_albgri_nir = tmp_albgri_nir + albdif(nc,nv,nz,2)*fveg(nc,nv,nz)
       end do 

       surfalb_inst%albgrd_col  (n,1)    = tmp_albgrd_vis
       surfalb_inst%albgrd_col  (n,2)    = tmp_albgrd_nir
       surfalb_inst%albgri_col  (n,1)    = tmp_albgri_vis
       surfalb_inst%albgri_col  (n,2)    = tmp_albgri_nir

       do np = 0,numpft
          p = p + 1  
      
          ! initialize temperature_inst here and not in its own F90 file, because values of tc, t10, and tm are computed in GridComp
          temperature_inst%t_veg_patch(p) = tc(nc,nz)
          temperature_inst%t_a10_patch(p) = t10(nc)
          temperature_inst%thm_patch(p)   = tm(nc)

          ! map Photosynthesis inputs to CLM space
          esat_tv_clm    (p) = esat_tv(nc,nz)
          oair_clm       (p) = oair(nc)
          cair_clm       (p) = cair(nc)
          rb_clm         (p) = rb(nc)
          qsatl_clm      (p) = qsatl(nc,nz)
          qaf_clm        (p) = qaf(nc,nz)
          dayl_factor_clm(p) = dayl_factor(nc)
          coszen_clm     (p) = coszen(nc)
          deldT_clm      (p) = deldT(nc)

          ! compute canopy air vapor pressure (in CLM space)
          eair_clm      (p) = eair(nc,nz)

          if (coszen_clm(p)>0. .and. (elai(p) + esai(p))>0.) then 
              ! calculate solar vegetated filter
              num_vegsol = num_vegsol + 1
              filter_vegsol(num_vegsol) = p

              ! calculate rho (weighted reflectance) and tau (weighted transmittance) needed for call to TwoStream later
              wl = elai(p) / max( elai(p)+esai(p), 1.e-06_r8 )
              ws = esai(p) / max( elai(p)+esai(p), 1.e-06_r8 )
            
              do ib = 1, numrad
                 rho(p,ib) = max( rhol(np,ib)*wl + rhos(np,ib)*ws, 1.e-06_r8 )
                 tau(p,ib) = max( taul(np,ib)*wl + taus(np,ib)*ws, 1.e-06_r8 )
              end do
          else
              num_novegsol = num_novegsol + 1
              filter_novegsol(num_novegsol) = p
          end if

          if (nlevcan == 1) then   ! jk: currently only coded for one canopy layer
             surfalb_inst%tlai_z_patch(p,1) = elai(p)
             surfalb_inst%tsai_z_patch(p,1) = esai(p)
          end if
  
          do iv = 1, surfalb_inst%nrad_patch(p)
             surfalb_inst%fabd_sun_z_patch(p,iv) = 0._r8
             surfalb_inst%fabd_sha_z_patch(p,iv) = 0._r8
             surfalb_inst%fabi_sun_z_patch(p,iv) = 0._r8
             surfalb_inst%fabi_sha_z_patch(p,iv) = 0._r8
             surfalb_inst%fsun_z_patch(p,iv)     = 0._r8
          end do               
       
          if (nlevcan == 1) then
             surfalb_inst%vcmaxcintsun_patch(p) = 0._r8
             surfalb_inst%vcmaxcintsha_patch(p) = (1._r8 - exp(-extkn*elai(p))) / extkn
             if (elai(p) > 0._r8) then
                surfalb_inst%vcmaxcintsha_patch(p) = surfalb_inst%vcmaxcintsha_patch(p) / elai(p)
             else
                surfalb_inst%vcmaxcintsha_patch(p) = 0._r8
             end if            
          else if (nlevcan > 1) then
             surfalb_inst%vcmaxcintsun_patch(p) = 0._r8
             surfalb_inst%vcmaxcintsha_patch(p) = 0._r8
          end if

          water_inst%waterdiagnosticbulk_inst%fdry_patch(p)    = (1-fwet(nc))*elai(p)/max( elai(p)+esai(p), 1.e-06_r8 )
          water_inst%waterdiagnosticbulk_inst%fwet_patch(p)    = fwet(nc)
          water_inst%waterdiagnosticbulk_inst%fcansno_patch(p) = fwet(nc)   !jk: This is not a mistake, see notes on why we set fcansno = fwet
       end do 
    end do
 end do


 ! call TwoStream subroutine which computes surface albedo variables needed for the subsequent calls;
 ! jk Jan 2022:  In older versions of CatchCN, the calculations were copy and pasted prior to the Photsynthesis calls;
 ! In CLM this subroutine is called *after* the canopy flux calculations, but I decided to add it here to
 ! have all required inputs on first time step (similar to how it was done in older CatchCN versions)
 
 call TwoStream(bounds, &
        filter_vegsol, num_vegsol, &
        coszen_clm, rho, tau, &
        canopystate_inst, temperature_inst, water_inst%waterdiagnosticbulk_inst, surfalb_inst)

 ! compute canopy shaded and sunlit variables (jk: needed to fill solarabs_inst before PHS call)
 call CanopySunShadeFracs(filter_nourbanp, filter_num_nourbanp,  &
                          atm2lnd_inst, surfalb_inst,     &
                          canopystate_inst, solarabs_inst)

! jkolassa: Below are three calls to the photosynthesis subroutine, one unperturbed, 
!           one with perturbed vapor pressure and one with perturbed canopy temperature.
!           The unperturbed call is issued last, so that CLM objects have unperturbed values
!           going forward. 

!  compute resistance with small delta ea
 call photosyns_inst%TimeStepInit(bounds)
 eair_pert(:) = eair_clm(:) + dea

 call PhotosynthesisHydraulicStress ( bounds, filter(1)%num_exposedvegp, filter(1)%exposedvegp, &
       esat_tv_clm, eair_pert, oair_clm, cair_clm, rb_clm, bsun, bsha, btran, dayl_factor_clm, leafn, &
       qsatl_clm, qaf_clm, &
       atm2lnd_inst, temperature_inst, soilstate_inst, water_inst%waterdiagnosticbulk_inst, &
       surfalb_inst, solarabs_inst, canopystate_inst, ozone_inst, &
       photosyns_inst, water_inst%waterfluxbulk_inst, froot_carbon, croot_carbon)
 
  laisun_dea = canopystate_inst%laisun_patch
  laisha_dea = canopystate_inst%laisha_patch
  rssun_dea  = photosyns_inst%rssun_patch
  rssha_dea  = photosyns_inst%rssha_patch

!  compute resistance with small delta Tc

 call photosyns_inst%TimeStepInit(bounds)
   temp_unpert =  temperature_inst%t_veg_patch
   temperature_inst%t_veg_patch = temperature_inst%t_veg_patch + dtc
   esat_tv_pert(:) = esat_tv_clm(:) + deldT_clm(:)*dtc 

 call PhotosynthesisHydraulicStress ( bounds, filter(1)%num_exposedvegp, filter(1)%exposedvegp, &
       esat_tv_pert, eair_clm, oair_clm, cair_clm, rb_clm, bsun, bsha, btran, dayl_factor_clm, leafn, &
       qsatl_clm, qaf_clm, &
       atm2lnd_inst, temperature_inst, soilstate_inst, water_inst%waterdiagnosticbulk_inst, &
       surfalb_inst, solarabs_inst, canopystate_inst, ozone_inst, &
       photosyns_inst, water_inst%waterfluxbulk_inst, froot_carbon, croot_carbon)

  laisun_dt = canopystate_inst%laisun_patch
  laisha_dt = canopystate_inst%laisha_patch
  rssun_dt  = photosyns_inst%rssun_patch
  rssha_dt  = photosyns_inst%rssha_patch

!  compute unperturbed resistance

 call photosyns_inst%TimeStepInit(bounds)

   temperature_inst%t_veg_patch = temp_unpert ! reset canopy temperature to unperturbed value

 call PhotosynthesisHydraulicStress ( bounds, filter(1)%num_exposedvegp, filter(1)%exposedvegp, &
       esat_tv_clm, eair_clm, oair_clm, cair_clm, rb_clm, bsun, bsha, btran, dayl_factor_clm, leafn, &
       qsatl_clm, qaf_clm, &
       atm2lnd_inst, temperature_inst, soilstate_inst, water_inst%waterdiagnosticbulk_inst, &
       surfalb_inst, solarabs_inst, canopystate_inst, ozone_inst, &
       photosyns_inst, water_inst%waterfluxbulk_inst, froot_carbon, croot_carbon)
 
  laisun = canopystate_inst%laisun_patch
  laisha = canopystate_inst%laisha_patch
  rssun  = photosyns_inst%rssun_patch
  rssha  = photosyns_inst%rssha_patch   

 call PhotosynthesisTotal (filter(1)%num_exposedvegp, filter(1)%exposedvegp, &
       atm2lnd_inst, canopystate_inst, photosyns_inst)

 laisun_out = 0.
 laisha_out = 0.
 psnsun_out = 0.
 psnsha_out = 0.
 lmrsun_out = 0.
 lmrsha_out = 0.

   np = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop 

          rcs = 0.
          rcs_dea = 0.
          rcs_dt = 0.

          do p = 0,numpft  ! PFT index loop
             np = np + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then

                 ! stomatal resistances
                 rs = laisun(np)/max(rssun(np), 1.e-06_r8 ) + laisha(np)/max(rssha(np), 1.e-06_r8 )
                 rcs = rcs + fveg(nc,nv,nz)*rs 

                 rs_dea = laisun_dea(np)/max(rssun_dea(np), 1.e-06_r8 ) + laisha_dea(np)/max(rssha_dea(np), 1.e-06_r8 )
                 rcs_dea = rcs_dea + fveg(nc,nv,nz)*rs_dea

                 rs_dt = laisun_dt(np)/max(rssun_dt(np), 1.e-06_r8 ) + laisha_dt(np)/max(rssha_dt(np), 1.e-06_r8 )
                 rcs_dt = rcs_dt + fveg(nc,nv,nz)*rs_dt

                 ! LAI
                 laisun_out(nc,nv,nz) = laisun(np)
                 laisha_out(nc,nv,nz) = laisha(np)

                 if (isnan(laisun(np))) laisun_out(nc,nv,nz) = 0.
                 if (isnan(laisha(np))) laisha_out(nc,nv,nz) = 0.

                 ! Photosynthesis
                 psnsun_out(nc,nv,nz) = photosyns_inst%psnsun_patch(np)
                 psnsha_out(nc,nv,nz) = photosyns_inst%psnsha_patch(np)

                 if (isnan(psnsun_out(nc,nv,nz))) psnsun_out(nc,nv,nz) = 0.
                 if (isnan(psnsha_out(nc,nv,nz))) psnsha_out(nc,nv,nz) = 0.

                 ! Leaf maintenance respiration
                 lmrsun_out(nc,nv,nz) = photosyns_inst%lmrsun_patch(np)
                 lmrsha_out(nc,nv,nz) = photosyns_inst%lmrsha_patch(np)

                 if (isnan(lmrsun_out(nc,nv,nz))) lmrsun_out(nc,nv,nz) = 0.
                 if (isnan(lmrsha_out(nc,nv,nz))) lmrsha_out(nc,nv,nz) = 0.

                 ! total absorbed PAR
                 tmp_parsun = 0.
                 tmp_parsha = 0.
                 do nl = 1,nlevcan
                    tmp_parsun = tmp_parsun + solarabs_inst%parsun_z_patch(np,nl)
                    tmp_parsha = tmp_parsha + solarabs_inst%parsha_z_patch(np,nl)
                 end do

                 parabs(nc,nv,nz) =  tmp_parsun * laisun(np) + tmp_parsha * laisha(np)

                 ! transpiration wetness factor / water stress

                 btran_out(nc,nv,nz) = btran(np)

                end if ! ityp = p
          end do !nv
       end do ! p
       rc(nc,nz)     = 1./max(rcs,5.e-5) + rb(nc)       ! rc: unperturbed stomatal resistance (rs is stomatal conductance)
       rc_dea(nc,nz) = 1./max(rcs_dea,5.e-5) + rb(nc)   ! rc_dea: stomatal resistance with vapor pressure perturbation
       rc_dt(nc,nz)  = 1./max(rcs_dt,5.e-5) + rb(nc)    ! rc_dt: stomatal resistance with canopy temperature perturbation
     end do ! nz
  end do ! nc

 end associate 

 deallocate(filter_vegsol)
 deallocate(filter_novegsol)
 deallocate(rho)
 deallocate(tau)
 deallocate(bsun)
 deallocate(bsha)
 deallocate(btran)

 end subroutine catchcn_calc_rc

end module CNCLM_Photosynthesis

