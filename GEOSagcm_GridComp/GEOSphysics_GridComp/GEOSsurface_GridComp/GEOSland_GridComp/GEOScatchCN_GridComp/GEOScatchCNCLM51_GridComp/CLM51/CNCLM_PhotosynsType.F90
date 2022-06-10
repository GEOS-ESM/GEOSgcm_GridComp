module CNCLM_PhotosynsType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use clm_varctl           , only : use_luna
  use clm_varpar           , only : numpft, num_zon, num_veg, &
                                    var_col, var_pft
  use nanMod               , only : nan
  use CNCLM_pftconMod      , only : pftcon
  use CNCLM_decompMod      , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_photosyns_type

  type, public :: photosyns_type

     logical , pointer, private :: c3flag_patch      (:)   ! patch true if C3 and false if C4
     ! Plant hydraulic stress specific variables
     real(r8), pointer, private :: ac_phs_patch      (:,:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: aj_phs_patch      (:,:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ap_phs_patch      (:,:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ag_phs_patch      (:,:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_sun_patch      (:,:)   ! patch sunlit net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_sha_patch      (:,:)   ! patch shaded net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: vcmax_z_phs_patch (:,:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, private :: kp_z_phs_patch    (:,:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, private :: tpu_z_phs_patch   (:,:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, public  :: gs_mol_sun_patch  (:,:) ! patch sunlit leaf stomatal conductance (umol H2O/m**2/s)
     real(r8), pointer, public  :: gs_mol_sha_patch  (:,:) ! patch shaded leaf stomatal conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: gs_mol_sun_ln_patch (:,:) ! patch sunlit leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
     real(r8), pointer, private :: gs_mol_sha_ln_patch (:,:) ! patch shaded leaf stomatal conductance averaged over 1 hour before to 1 hour after local noon (umol H2O/m**2/s)
     real(r8), pointer, private :: ac_patch          (:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: aj_patch          (:,:) ! patch RuBP-limited gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ap_patch          (:,:) ! patch product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: ag_patch          (:,:) ! patch co-limited gross leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: an_patch          (:,:) ! patch net leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: vcmax_z_patch     (:,:) ! patch maximum rate of carboxylation (umol co2/m**2/s)
     real(r8), pointer, private :: cp_patch          (:)   ! patch CO2 compensation point (Pa)
     real(r8), pointer, private :: kc_patch          (:)   ! patch Michaelis-Menten constant for CO2 (Pa)
     real(r8), pointer, private :: ko_patch          (:)   ! patch Michaelis-Menten constant for O2 (Pa)
     real(r8), pointer, private :: qe_patch          (:)   ! patch quantum efficiency, used only for C4 (mol CO2 / mol photons)
     real(r8), pointer, private :: tpu_z_patch       (:,:) ! patch triose phosphate utilization rate (umol CO2/m**2/s)
     real(r8), pointer, private :: kp_z_patch        (:,:) ! patch initial slope of CO2 response curve (C4 plants)
     real(r8), pointer, private :: bbb_patch         (:)   ! patch Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: mbb_patch         (:)   ! patch Ball-Berry slope of conductance-photosynthesis relationship
     real(r8), pointer, private :: gs_mol_patch      (:,:) ! patch leaf stomatal conductance       (umol H2O/m**2/s)
     real(r8), pointer, private :: gb_mol_patch      (:)   ! patch leaf boundary layer conductance (umol H2O/m**2/s)
     real(r8), pointer, private :: rh_leaf_patch     (:)   ! patch fractional humidity at leaf surface (dimensionless)
     real(r8), pointer, private :: vpd_can_patch     (:)   ! patch canopy vapor pressure deficit (kPa)
     real(r8), pointer, private :: alphapsnsun_patch (:)   ! patch sunlit 13c fractionation ([])
     real(r8), pointer, private :: alphapsnsha_patch (:)   ! patch shaded 13c fractionation ([])

     real(r8), pointer, public  :: rc13_canair_patch (:)   ! patch C13O2/C12O2 in canopy air
     real(r8), pointer, public  :: rc13_psnsun_patch (:)   ! patch C13O2/C12O2 in sunlit canopy psn flux
     real(r8), pointer, public  :: rc13_psnsha_patch (:)   ! patch C13O2/C12O2 in shaded canopy psn flux

     real(r8), pointer, public  :: psnsun_patch      (:)   ! patch sunlit leaf photosynthesis     (umol CO2/m**2/s)
     real(r8), pointer, public  :: psnsha_patch      (:)   ! patch shaded leaf photosynthesis     (umol CO2/m**2/s)
     real(r8), pointer, public  :: c13_psnsun_patch  (:)   ! patch c13 sunlit leaf photosynthesis (umol 13CO2/m**2/s)
     real(r8), pointer, public  :: c13_psnsha_patch  (:)   ! patch c13 shaded leaf photosynthesis (umol 13CO2/m**2/s)
     real(r8), pointer, public  :: c14_psnsun_patch  (:)   ! patch c14 sunlit leaf photosynthesis (umol 14CO2/m**2/s)
     real(r8), pointer, public  :: c14_psnsha_patch  (:)   ! patch c14 shaded leaf photosynthesis (umol 14CO2/m**2/s)

     real(r8), pointer, private :: psnsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_z_patch    (:,:) ! patch canopy layer: shaded leaf photosynthesis   (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wc_patch   (:)   ! patch Rubsico-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wc_patch   (:)   ! patch Rubsico-limited shaded leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wj_patch   (:)   ! patch RuBP-limited sunlit leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wj_patch   (:)   ! patch RuBP-limited shaded leaf photosynthesis    (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsun_wp_patch   (:)   ! patch product-limited sunlit leaf photosynthesis (umol CO2/m**2/s)
     real(r8), pointer, private :: psnsha_wp_patch   (:)   ! patch product-limited shaded leaf photosynthesis (umol CO2/m**2/s)

     real(r8), pointer, public  :: fpsn_patch        (:)   ! patch photosynthesis                 (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wc_patch     (:)   ! patch Rubisco-limited photosynthesis (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wj_patch     (:)   ! patch RuBP-limited photosynthesis    (umol CO2/m**2 ground/s)
     real(r8), pointer, private :: fpsn_wp_patch     (:)   ! patch product-limited photosynthesis (umol CO2/m**2 ground/s)

     real(r8), pointer, public  :: lnca_patch        (:)   ! top leaf layer leaf N concentration (gN leaf/m^2)

     real(r8), pointer, public  :: lmrsun_patch      (:)   ! patch sunlit leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer, public  :: lmrsha_patch      (:)   ! patch shaded leaf maintenance respiration rate               (umol CO2/m**2/s)
     real(r8), pointer, private :: lmrsun_z_patch    (:,:) ! patch canopy layer: sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
     real(r8), pointer, private :: lmrsha_z_patch    (:,:) ! patch canopy layer: shaded leaf maintenance respiration rate (umol CO2/m**2/s)

     real(r8), pointer, public  :: cisun_z_patch     (:,:) ! patch intracellular sunlit leaf CO2 (Pa)
     real(r8), pointer, public  :: cisha_z_patch     (:,:) ! patch intracellular shaded leaf CO2 (Pa)

     real(r8), pointer, private :: rssun_z_patch     (:,:) ! patch canopy layer: sunlit leaf stomatal resistance (s/m)
     real(r8), pointer, private :: rssha_z_patch     (:,:) ! patch canopy layer: shaded leaf stomatal resistance (s/m)
     real(r8), pointer, public  :: rssun_patch       (:)   ! patch sunlit stomatal resistance (s/m)
     real(r8), pointer, public  :: rssha_patch       (:)   ! patch shaded stomatal resistance (s/m)
     real(r8), pointer, public  :: luvcmax25top_patch (:)   ! vcmax25 !     (umol/m2/s)
     real(r8), pointer, public  :: lujmax25top_patch  (:)   ! vcmax25 (umol/m2/s)
     real(r8), pointer, public  :: lutpu25top_patch   (:)   ! vcmax25 (umol/m2/s)
!!


     ! LUNA specific variables
     real(r8), pointer, public  :: vcmx25_z_patch    (:,:) ! patch  leaf Vc,max25 (umol CO2/m**2/s) for canopy layer 
     real(r8), pointer, public  :: jmx25_z_patch     (:,:) ! patch  leaf Jmax25 (umol electron/m**2/s) for canopy layer 
     real(r8), pointer, public  :: vcmx25_z_last_valid_patch  (:,:) ! patch  leaf Vc,max25 at the end of the growing season for the previous year
     real(r8), pointer, public  :: jmx25_z_last_valid_patch   (:,:) ! patch  leaf Jmax25 at the end of the growing season for the previous year
     real(r8), pointer, public  :: pnlc_z_patch      (:,:) ! patch proportion of leaf nitrogen allocated for light capture for canopy layer
     real(r8), pointer, public  :: enzs_z_patch      (:,:) ! enzyme decay status 1.0-fully active; 0-all decayed during stress
     real(r8), pointer, public  :: fpsn24_patch      (:)   ! 24 hour mean patch photosynthesis (umol CO2/m**2 ground/day)

     ! Logical switches for different options
     logical, public  :: rootstem_acc                      ! Respiratory acclimation for roots and stems
     logical, private :: light_inhibit                     ! If light should inhibit respiration
     integer, private :: leafresp_method                   ! leaf maintencence respiration at 25C for canopy top method to use
     integer, private :: stomatalcond_mtd                  ! Stomatal conduction method type
     logical, private :: modifyphoto_and_lmr_forcrop       ! Modify photosynthesis and LMR for crop

  end type photosyns_type
  type(photosyns_type), public, target, save :: photosyns_inst

contains 

!-------------------------------------------------------------
 subroutine init_photosyns_type(bounds, nch, ityp, fveg, cncol, cnpft, this, cn5_cold_start)

  ! !DESCRIPTION:
  ! Initialize CTSM photosynthesis type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    ! INPUT/OUTPUT
    type(bounds_type),                               intent(in) :: bounds
    integer,                                         intent(in) :: nch    ! number of Catchment tiles
    integer, dimension(nch,num_veg,num_zon),         intent(in) :: ityp   ! PFT index
    real,    dimension(nch,num_veg,num_zon),         intent(in) :: fveg   ! PFT fraction   
    real,    dimension(nch,num_zon,var_col),         intent(in) :: cncol  ! column-level restart variable array 
    real,    dimension(nch,num_zon,num_veg,var_pft), intent(in) :: cnpft  ! pft-level (patch-level) restart variable array
    logical, optional,                               intent(in) :: cn5_cold_start
    type(photosyns_type),                            intent(inout):: this

    ! LOCAL
    integer :: begp, endp  ! patch-level beginning and end index
    integer :: begc, endc  ! column-level beginning and end index 
    integer :: np, nc, nz, p, nv
    logical :: cold_start = .false.
    !------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
 
    ! check whether a cn5_cold_start option was set and change cold_start accordingly
    if (present(cn5_cold_start) .and. (cn5_cold_start==.true.)) then
       cold_start = .true.
    end if

    ! jkolassa: if cold_start is false, check that both CNCOL and CNPFT have the expected size for CNCLM50, else abort 
    if ((cold_start==.false.) .and. ((size(cncol,3).ne.var_col) .or. &
       (size(cnpft,3).ne.var_pft)))
       _ASSERT(.FALSE.,'option CNCLM50_cold_start = .FALSE. requires a CNCLM50 restart file')
    end if 


    allocate(this%c3flag_patch      (begp:endp))             ; this%c3flag_patch      (:)     =.false.
    allocate(this%ac_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ac_phs_patch      (:,:,:) = nan
    allocate(this%aj_phs_patch      (begp:endp,2,1:nlevcan)) ; this%aj_phs_patch      (:,:,:) = nan
    allocate(this%ap_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ap_phs_patch      (:,:,:) = nan
    allocate(this%ag_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ag_phs_patch      (:,:,:) = nan
    allocate(this%an_sun_patch      (begp:endp,1:nlevcan))   ; this%an_sun_patch      (:,:)   = nan
    allocate(this%an_sha_patch      (begp:endp,1:nlevcan))   ; this%an_sha_patch      (:,:)   = nan
    allocate(this%vcmax_z_phs_patch (begp:endp,2,1:nlevcan)) ; this%vcmax_z_phs_patch (:,:,:) = nan
    allocate(this%tpu_z_phs_patch   (begp:endp,2,1:nlevcan)) ; this%tpu_z_phs_patch   (:,:,:) = nan
    allocate(this%kp_z_phs_patch    (begp:endp,2,1:nlevcan)) ; this%kp_z_phs_patch    (:,:,:) = nan
    allocate(this%gs_mol_sun_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sun_patch  (:,:)   = nan
    allocate(this%gs_mol_sha_patch  (begp:endp,1:nlevcan))   ; this%gs_mol_sha_patch  (:,:)   = nan
    allocate(this%gs_mol_sun_ln_patch (begp:endp,1:nlevcan)) ; this%gs_mol_sun_ln_patch (:,:)   = nan
    allocate(this%gs_mol_sha_ln_patch (begp:endp,1:nlevcan)) ; this%gs_mol_sha_ln_patch (:,:)   = nan
    allocate(this%ac_patch          (begp:endp,1:nlevcan)) ; this%ac_patch          (:,:) = nan
    allocate(this%aj_patch          (begp:endp,1:nlevcan)) ; this%aj_patch          (:,:) = nan
    allocate(this%ap_patch          (begp:endp,1:nlevcan)) ; this%ap_patch          (:,:) = nan
    allocate(this%ag_patch          (begp:endp,1:nlevcan)) ; this%ag_patch          (:,:) = nan
    allocate(this%an_patch          (begp:endp,1:nlevcan)) ; this%an_patch          (:,:) = nan
    allocate(this%vcmax_z_patch     (begp:endp,1:nlevcan)) ; this%vcmax_z_patch     (:,:) = nan
    allocate(this%tpu_z_patch       (begp:endp,1:nlevcan)) ; this%tpu_z_patch       (:,:) = nan
    allocate(this%kp_z_patch        (begp:endp,1:nlevcan)) ; this%kp_z_patch        (:,:) = nan
    allocate(this%gs_mol_patch      (begp:endp,1:nlevcan)) ; this%gs_mol_patch      (:,:) = nan
    allocate(this%cp_patch          (begp:endp))           ; this%cp_patch          (:)   = nan
    allocate(this%kc_patch          (begp:endp))           ; this%kc_patch          (:)   = nan
    allocate(this%ko_patch          (begp:endp))           ; this%ko_patch          (:)   = nan
    allocate(this%qe_patch          (begp:endp))           ; this%qe_patch          (:)   = nan
    allocate(this%bbb_patch         (begp:endp))           ; this%bbb_patch         (:)   = nan
    allocate(this%mbb_patch         (begp:endp))           ; this%mbb_patch         (:)   = nan
    allocate(this%gb_mol_patch      (begp:endp))           ; this%gb_mol_patch      (:)   = nan
    allocate(this%rh_leaf_patch     (begp:endp))           ; this%rh_leaf_patch     (:)   = nan
    allocate(this%vpd_can_patch     (begp:endp))           ; this%vpd_can_patch     (:)   = nan
    allocate(this%psnsun_patch      (begp:endp))           ; this%psnsun_patch      (:)   = nan
    allocate(this%psnsha_patch      (begp:endp))           ; this%psnsha_patch      (:)   = nan
    allocate(this%c13_psnsun_patch  (begp:endp))           ; this%c13_psnsun_patch  (:)   = nan
    allocate(this%c13_psnsha_patch  (begp:endp))           ; this%c13_psnsha_patch  (:)   = nan
    allocate(this%c14_psnsun_patch  (begp:endp))           ; this%c14_psnsun_patch  (:)   = nan
    allocate(this%c14_psnsha_patch  (begp:endp))           ; this%c14_psnsha_patch  (:)   = nan

    allocate(this%psnsun_z_patch    (begp:endp,1:nlevcan)) ; this%psnsun_z_patch    (:,:) = nan
    allocate(this%psnsha_z_patch    (begp:endp,1:nlevcan)) ; this%psnsha_z_patch    (:,:) = nan
    allocate(this%psnsun_wc_patch   (begp:endp))           ; this%psnsun_wc_patch   (:)   = nan
    allocate(this%psnsha_wc_patch   (begp:endp))           ; this%psnsha_wc_patch   (:)   = nan
    allocate(this%psnsun_wj_patch   (begp:endp))           ; this%psnsun_wj_patch   (:)   = nan
    allocate(this%psnsha_wj_patch   (begp:endp))           ; this%psnsha_wj_patch   (:)   = nan
    allocate(this%psnsun_wp_patch   (begp:endp))           ; this%psnsun_wp_patch   (:)   = nan
    allocate(this%psnsha_wp_patch   (begp:endp))           ; this%psnsha_wp_patch   (:)   = nan
    allocate(this%fpsn_patch        (begp:endp))           ; this%fpsn_patch        (:)   = nan
    allocate(this%fpsn_wc_patch     (begp:endp))           ; this%fpsn_wc_patch     (:)   = nan
    allocate(this%fpsn_wj_patch     (begp:endp))           ; this%fpsn_wj_patch     (:)   = nan
    allocate(this%fpsn_wp_patch     (begp:endp))           ; this%fpsn_wp_patch     (:)   = nan

    allocate(this%lnca_patch        (begp:endp))           ; this%lnca_patch        (:)   = nan

    allocate(this%lmrsun_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsun_z_patch    (:,:) = nan
    allocate(this%lmrsha_z_patch    (begp:endp,1:nlevcan)) ; this%lmrsha_z_patch    (:,:) = nan
    allocate(this%lmrsun_patch      (begp:endp))           ; this%lmrsun_patch      (:)   = nan
    allocate(this%lmrsha_patch      (begp:endp))           ; this%lmrsha_patch      (:)   = nan

    allocate(this%alphapsnsun_patch (begp:endp))           ; this%alphapsnsun_patch (:)   = nan
    allocate(this%alphapsnsha_patch (begp:endp))           ; this%alphapsnsha_patch (:)   = nan
    allocate(this%rc13_canair_patch (begp:endp))           ; this%rc13_canair_patch (:)   = nan
    allocate(this%rc13_psnsun_patch (begp:endp))           ; this%rc13_psnsun_patch (:)   = nan
    allocate(this%rc13_psnsha_patch (begp:endp))           ; this%rc13_psnsha_patch (:)   = nan

    allocate(this%cisun_z_patch     (begp:endp,1:nlevcan)) ; this%cisun_z_patch     (:,:) = nan
    allocate(this%cisha_z_patch     (begp:endp,1:nlevcan)) ; this%cisha_z_patch     (:,:) = nan

    allocate(this%rssun_z_patch     (begp:endp,1:nlevcan)) ; this%rssun_z_patch     (:,:) = nan
    allocate(this%rssha_z_patch     (begp:endp,1:nlevcan)) ; this%rssha_z_patch     (:,:) = nan
    allocate(this%rssun_patch       (begp:endp))           ; this%rssun_patch       (:)   = nan
    allocate(this%rssha_patch       (begp:endp))           ; this%rssha_patch       (:)   = nan
    allocate(this%luvcmax25top_patch(begp:endp))           ; this%luvcmax25top_patch(:) = nan
    allocate(this%lujmax25top_patch (begp:endp))           ; this%lujmax25top_patch(:)  = nan
    allocate(this%lutpu25top_patch  (begp:endp))           ; this%lutpu25top_patch(:)   = nan
!!
!    allocate(this%psncanopy_patch   (begp:endp))           ; this%psncanopy_patch   (:)   = nan
!    allocate(this%lmrcanopy_patch   (begp:endp))           ; this%lmrcanopy_patch   (:)   = nan
    if(use_luna)then
      ! NOTE(bja, 2015-09) because these variables are only allocated
      ! when luna is turned on, they can not be placed into associate
      ! statements.
      allocate(this%vcmx25_z_patch  (begp:endp,1:nlevcan)) ; this%vcmx25_z_patch    (:,:) = 30._r8
      allocate(this%jmx25_z_patch   (begp:endp,1:nlevcan)) ; this%jmx25_z_patch     (:,:) = 60._r8
      allocate(this%vcmx25_z_last_valid_patch     (begp:endp,1:nlevcan)) ; this%vcmx25_z_last_valid_patch       (:,:) = 30._r8
      allocate(this%jmx25_z_last_valid_patch      (begp:endp,1:nlevcan)) ; this%jmx25_z_last_valid_patch        (:,:) = 60._r8
      allocate(this%pnlc_z_patch    (begp:endp,1:nlevcan)) ; this%pnlc_z_patch      (:,:) = 0.01_r8
      allocate(this%fpsn24_patch    (begp:endp))           ; this%fpsn24_patch      (:)   = nan
      allocate(this%enzs_z_patch    (begp:endp,1:nlevcan)) ; this%enzs_z_patch      (:,:) = 1._r8
    endif


     this%light_inhibit = .true.     ! jkolassa, Feb 2022: This is the default value for CTSM5.1; we could in the future control this through resource files

     this%leafresp_method = 2        ! jkolassa, Feb 2022: Default for CTSM5.1 if use_cn is true (2 corresponds to Atkin et al., 2015)

     this%stomatalcond_mtd = 2       ! jkolassa, Feb 2022: Default for CTSM5.1, corresponds to Medlyn et al., 2011

     this%modifyphoto_and_lmr_forcrop = .true. ! jkolassa, Feb 2022: Default for CLM50 and up


  ! initialize types from restart file or through cold start values

  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,num_zon    ! CN zone loop
       do p = 0,numpft  ! PFT index loop
          np = np + 1
          do nv = 1,num_veg ! defined veg loop
             if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
                if (cold_start) then
                   photosyns_inst%alphapsnsun_patch(np) = 0._r8
                   photosyns_inst%alphapsnsha_patch(np) = 0._r8
                else (cold_start=.false.) then
                    photosyns_inst%alphapsnsun_patch(np) = cnpft(nc,nz,nv, 75)
                    photosyns_inst%alphapsnsha_patch(np) = cnpft(nc,nz,nv, 76)
                end if 
              end if ! ityp =p  
          end do !nv
       end do ! p
     end do ! nz
  end do ! nc

  end subroutine init_photosyns_type

end module CNCLM_PhotosynsType
