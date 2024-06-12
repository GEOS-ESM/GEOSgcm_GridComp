module SoilStateType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use clm_varpar       , only : nlevsoi, nlevgrnd, nlevmaxurbgrnd, &
                                nlayer, nlevsno
  use clm_varcon       , only : spval
  use nanMod           , only : nan
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  !
  type, public :: soilstate_type

     ! sand/ clay/ organic matter
     real(r8), pointer :: sandfrac_patch       (:)   ! patch sand fraction
     real(r8), pointer :: clayfrac_patch       (:)   ! patch clay fraction
     real(r8), pointer :: mss_frc_cly_vld_col  (:)   ! col mass fraction clay limited to 0.20
     real(r8), pointer :: cellorg_col          (:,:) ! col organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellsand_col         (:,:) ! sand value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: cellclay_col         (:,:) ! clay value for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: bd_col               (:,:) ! col bulk density of dry soil material [kg/m^3] (CN)

     ! hydraulic properties
     real(r8), pointer :: hksat_col            (:,:) ! col hydraulic conductivity at saturation (mm H2O /s) 
     real(r8), pointer :: hksat_min_col        (:,:) ! col mineral hydraulic conductivity at saturation (hksat) (mm/s)
     real(r8), pointer :: hk_l_col             (:,:) ! col hydraulic conductivity (mm/s)
     real(r8), pointer :: smp_l_col            (:,:) ! col soil matric potential (mm)
     real(r8), pointer :: smpmin_col           (:)   ! col restriction for min of soil potential (mm) 
     real(r8), pointer :: bsw_col              (:,:) ! col Clapp and Hornberger "b" (nlevgrnd)  
     real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity) 
     real(r8), pointer :: watdry_col           (:,:) ! col btran parameter for btran = 0
     real(r8), pointer :: watopt_col           (:,:) ! col btran parameter for btran = 1
     real(r8), pointer :: watfc_col            (:,:) ! col volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd) 
     real(r8), pointer :: dsl_col              (:)   ! col dry surface layer thickness (mm)
     real(r8), pointer :: soilresis_col        (:)   ! col soil evaporative resistance S&L14 (s/m)
     real(r8), pointer :: soilbeta_col         (:)   ! col factor that reduces ground evaporation L&P1992(-)
     real(r8), pointer :: soilalpha_col        (:)   ! col factor that reduces ground saturated specific humidity (-)
     real(r8), pointer :: soilalpha_u_col      (:)   ! col urban factor that reduces ground saturated specific humidity (-) 
     real(r8), pointer :: soilpsi_col          (:,:) ! col soil water potential in each soil layer (MPa) (CN)
     real(r8), pointer :: psiwilt_col          (:,:) ! col soil water potential at wilting point (added by jkolassa to use for assessing water stress instead of globally constant value)
     real(r8), pointer :: wtfact_col           (:)   ! col maximum saturated fraction for a gridcell
     real(r8), pointer :: porosity_col         (:,:) ! col soil porisity (1-bulk_density/soil_density) (VIC)
     real(r8), pointer :: eff_porosity_col     (:,:) ! col effective porosity = porosity - vol_ice (nlevgrnd) 
     real(r8), pointer :: gwc_thr_col          (:)   ! col threshold soil moisture based on clay content
!scs: vangenuchten
     real(r8), pointer :: msw_col              (:,:) ! col vanGenuchtenClapp "m"
     real(r8), pointer :: nsw_col              (:,:) ! col vanGenuchtenClapp "n"
     real(r8), pointer :: alphasw_col          (:,:) ! col vanGenuchtenClapp "nalpha"
     real(r8), pointer :: watres_col           (:,:) ! residual soil water content
     ! thermal conductivity / heat capacity
     real(r8), pointer :: thk_col              (:,:) ! col thermal conductivity of each layer [W/m-K] 
     real(r8), pointer :: tkmg_col             (:,:) ! col thermal conductivity, soil minerals  [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: tkdry_col            (:,:) ! col thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd) 
     real(r8), pointer :: tksatu_col           (:,:) ! col thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd) 
     real(r8), pointer :: csol_col             (:,:) ! col heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd) 

     ! roots
     real(r8), pointer :: rootr_patch          (:,:) ! patch effective fraction of roots in each soil layer (SMS method only) (nlevgrnd)
     real(r8), pointer :: rootr_col            (:,:) ! col effective fraction of roots in each soil layer (SMS method only) (nlevgrnd)  
     real(r8), pointer :: rootfr_col           (:,:) ! col fraction of roots in each soil layer (nlevgrnd) 
     real(r8), pointer :: rootfr_patch         (:,:) ! patch fraction of roots for water in each soil layer (nlevgrnd)
     real(r8), pointer :: crootfr_patch        (:,:) ! patch fraction of roots for carbon in each soil layer (nlevgrnd)
     real(r8), pointer :: root_depth_patch     (:)   ! root depth
     real(r8), pointer :: rootr_road_perv_col  (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: rootfr_road_perv_col (:,:) ! col effective fraction of roots in each soil layer of urban pervious road
     real(r8), pointer :: k_soil_root_patch    (:,:) ! patch soil-root interface conductance [mm/s]
     real(r8), pointer :: root_conductance_patch(:,:) ! patch root conductance [mm/s]
     real(r8), pointer :: soil_conductance_patch(:,:) ! patch soil conductance [mm/s]

   contains 

     procedure, public :: Init

end type soilstate_type
type(soilstate_type), public, target, save :: soilstate_inst

contains 

!-----------------------------------------------------------
  subroutine Init(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM soil state type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT
    type(bounds_type), intent(in) :: bounds
    class(soilstate_type)         :: this

    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    !-----------------------

    begp = bounds%begp  ; endp = bounds%endp
    begc = bounds%begc  ; endc = bounds%endc

    allocate(this%mss_frc_cly_vld_col  (begc:endc))                     ; this%mss_frc_cly_vld_col  (:)   = nan
    allocate(this%sandfrac_patch       (begp:endp))                     ; this%sandfrac_patch       (:)   = nan
    allocate(this%clayfrac_patch       (begp:endp))                     ; this%clayfrac_patch       (:)   = nan
    allocate(this%cellorg_col          (begc:endc,nlevsoi))             ; this%cellorg_col          (:,:) = nan
    allocate(this%cellsand_col         (begc:endc,nlevsoi))             ; this%cellsand_col         (:,:) = nan
    allocate(this%cellclay_col         (begc:endc,nlevsoi))             ; this%cellclay_col         (:,:) = nan
    allocate(this%bd_col               (begc:endc,nlevgrnd))            ; this%bd_col               (:,:) = nan

    allocate(this%hksat_col            (begc:endc,nlevgrnd))            ; this%hksat_col            (:,:) = spval
    allocate(this%hksat_min_col        (begc:endc,nlevgrnd))            ; this%hksat_min_col        (:,:) = spval
    allocate(this%hk_l_col             (begc:endc,nlevgrnd))            ; this%hk_l_col             (:,:) = nan
    allocate(this%smp_l_col            (begc:endc,nlevgrnd))            ; this%smp_l_col            (:,:) = nan
    allocate(this%smpmin_col           (begc:endc))                     ; this%smpmin_col           (:)   = nan

    allocate(this%bsw_col              (begc:endc,nlevgrnd))            ; this%bsw_col              (:,:) = nan
    allocate(this%watsat_col           (begc:endc,nlevmaxurbgrnd))      ; this%watsat_col           (:,:) = nan
    allocate(this%watdry_col           (begc:endc,nlevgrnd))            ; this%watdry_col           (:,:) = spval
    allocate(this%watopt_col           (begc:endc,nlevgrnd))            ; this%watopt_col           (:,:) = spval
    allocate(this%watfc_col            (begc:endc,nlevgrnd))            ; this%watfc_col            (:,:) = nan
    allocate(this%sucsat_col           (begc:endc,nlevgrnd))            ; this%sucsat_col           (:,:) = spval
    allocate(this%dsl_col              (begc:endc))                     ; this%dsl_col         (:)   = spval!nan   
    allocate(this%soilresis_col        (begc:endc))                     ; this%soilresis_col         (:)   = spval!nan   
    allocate(this%soilbeta_col         (begc:endc))                     ; this%soilbeta_col         (:)   = nan
    allocate(this%soilalpha_col        (begc:endc))                     ; this%soilalpha_col        (:)   = nan
    allocate(this%soilalpha_u_col      (begc:endc))                     ; this%soilalpha_u_col      (:)   = nan
    allocate(this%soilpsi_col          (begc:endc,nlevgrnd))            ; this%soilpsi_col          (:,:) = nan
    allocate(this%psiwilt_col          (begc:endc,nlevgrnd))            ; this%psiwilt_col          (:,:) = nan
    allocate(this%wtfact_col           (begc:endc))                     ; this%wtfact_col           (:)   = nan
    allocate(this%porosity_col         (begc:endc,nlayer))              ; this%porosity_col         (:,:) = spval
    allocate(this%eff_porosity_col     (begc:endc,nlevgrnd))            ; this%eff_porosity_col     (:,:) = spval
    allocate(this%gwc_thr_col          (begc:endc))                     ; this%gwc_thr_col          (:)   = nan

    allocate(this%thk_col              (begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ; this%thk_col              (:,:) = nan
    allocate(this%tkmg_col             (begc:endc,nlevgrnd))            ; this%tkmg_col             (:,:) = nan
    allocate(this%tkdry_col            (begc:endc,nlevgrnd))            ; this%tkdry_col            (:,:) = nan
    allocate(this%tksatu_col           (begc:endc,nlevgrnd))            ; this%tksatu_col           (:,:) = nan
    allocate(this%csol_col             (begc:endc,nlevgrnd))            ; this%csol_col             (:,:) = nan

    allocate(this%rootr_patch          (begp:endp,1:nlevgrnd))          ; this%rootr_patch          (:,:) = nan
    allocate(this%root_depth_patch     (begp:endp))                     ; this%root_depth_patch     (:)   = nan
    allocate(this%rootr_col            (begc:endc,nlevgrnd))            ; this%rootr_col            (:,:) = nan
    allocate(this%rootr_road_perv_col  (begc:endc,1:nlevgrnd))          ; this%rootr_road_perv_col  (:,:) = nan
    allocate(this%rootfr_patch         (begp:endp,1:nlevgrnd))          ; this%rootfr_patch         (:,:) = nan
    allocate(this%crootfr_patch        (begp:endp,1:nlevgrnd))          ; this%crootfr_patch        (:,:) = nan
    allocate(this%rootfr_col           (begc:endc,1:nlevgrnd))          ; this%rootfr_col           (:,:) = nan
    allocate(this%rootfr_road_perv_col (begc:endc,1:nlevgrnd))          ; this%rootfr_road_perv_col (:,:) = nan
    allocate(this%k_soil_root_patch    (begp:endp,1:nlevsoi))           ; this%k_soil_root_patch (:,:) = nan
    allocate(this%root_conductance_patch(begp:endp,1:nlevsoi))          ; this%root_conductance_patch (:,:) = nan
    allocate(this%soil_conductance_patch(begp:endp,1:nlevsoi))          ; this%soil_conductance_patch (:,:) = nan
    allocate(this%msw_col              (begc:endc,1:nlevgrnd))          ; this%msw_col              (:,:) = nan
    allocate(this%nsw_col              (begc:endc,1:nlevgrnd))          ; this%nsw_col              (:,:) = nan
    allocate(this%alphasw_col          (begc:endc,1:nlevgrnd))          ; this%alphasw_col          (:,:) = nan
    allocate(this%watres_col           (begc:endc,1:nlevgrnd))          ; this%watres_col           (:,:) = nan

  end subroutine Init

end module SoilStateType
