module SolarAbsorbedType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use clm_varcon       , only : spval
  use clm_varpar       , only : nlevcan, numrad, nlevsno
  use clm_varctl       , only : use_luna
  use nanMod           , only : nan
  use decompMod        , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: solarabs_type

     ! Solar reflected
     real(r8), pointer :: fsr_patch              (:)   ! patch solar radiation reflected (W/m**2)
     real(r8), pointer :: fsrSF_patch            (:)   ! diagnostic snow-free patch solar radiation reflected (W/m**2)
     real(r8), pointer :: ssre_fsr_patch         (:)   ! snow radiative effect on patch solar radiation reflected (W/m**2)
     ! Solar Absorbed
     real(r8), pointer :: fsa_patch              (:)   ! patch solar radiation absorbed (total) (W/m**2)  
     real(r8), pointer :: fsa_u_patch            (:)   ! patch urban solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: fsa_r_patch            (:)   ! patch rural solar radiation absorbed (total) (W/m**2)
     real(r8), pointer :: parsun_z_patch         (:,:) ! patch absorbed PAR for sunlit leaves in canopy layer (W/m**2) 
     real(r8), pointer :: parsha_z_patch         (:,:) ! patch absorbed PAR for shaded leaves in canopy layer (W/m**2) 
     real(r8), pointer :: par240d_z_patch        (:,:) ! 10-day running mean of daytime patch absorbed PAR for leaves in canopy layer (W/m**2) 
     real(r8), pointer :: par240x_z_patch        (:,:) ! 10-day running mean of maximum patch absorbed PAR for leaves in canopy layer (W/m**2)
     real(r8), pointer :: par24d_z_patch         (:,:) ! daily accumulated  absorbed PAR for leaves in canopy layer from midnight to current step(J/m**2) 
     real(r8), pointer :: par24x_z_patch         (:,:) ! daily max of patch absorbed PAR for  leaves in canopy layer from midnight to current step(W/m**2)
     real(r8), pointer :: sabg_soil_patch        (:)   ! patch solar radiation absorbed by soil (W/m**2)       
     real(r8), pointer :: sabg_snow_patch        (:)   ! patch solar radiation absorbed by snow (W/m**2)       
     real(r8), pointer :: sabg_patch             (:)   ! patch solar radiation absorbed by ground (W/m**2)     
     real(r8), pointer :: sabg_chk_patch         (:)   ! patch fsno weighted sum (W/m**2)                                   
     real(r8), pointer :: sabg_lyr_patch         (:,:) ! patch absorbed radiation in each snow layer and top soil layer (pft,lyr) [W/m2]
     real(r8), pointer :: sabg_pen_patch         (:)   ! patch (rural) shortwave radiation penetrating top soisno layer [W/m2]

     real(r8), pointer :: sub_surf_abs_SW_patch  (:)   ! patch fraction of solar radiation absorbed below first snow layer
     real(r8), pointer :: sabv_patch             (:)   ! patch solar radiation absorbed by vegetation (W/m**2) 

     real(r8), pointer :: sabs_roof_dir_lun      (:,:) ! lun direct solar absorbed by roof per unit ground area per unit incident flux
     real(r8), pointer :: sabs_roof_dif_lun      (:,:) ! lun diffuse solar absorbed by roof per unit ground area per unit incident flux
     real(r8), pointer :: sabs_sunwall_dir_lun   (:,:) ! lun direct  solar absorbed by sunwall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_sunwall_dif_lun   (:,:) ! lun diffuse solar absorbed by sunwall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_shadewall_dir_lun (:,:) ! lun direct  solar absorbed by shadewall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_shadewall_dif_lun (:,:) ! lun diffuse solar absorbed by shadewall per unit wall area per unit incident flux
     real(r8), pointer :: sabs_improad_dir_lun   (:,:) ! lun direct  solar absorbed by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_improad_dif_lun   (:,:) ! lun diffuse solar absorbed by impervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_perroad_dir_lun   (:,:) ! lun direct  solar absorbed by pervious road per unit ground area per unit incident flux
     real(r8), pointer :: sabs_perroad_dif_lun   (:,:) ! lun diffuse solar absorbed by pervious road per unit ground area per unit incident flux

     ! Currently needed by lake code 
     ! TODO (MV 8/20/2014) should be moved in the future
     real(r8), pointer :: fsds_nir_d_patch       (:)   ! patch incident direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_i_patch       (:)   ! patch incident diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsds_nir_d_ln_patch    (:)   ! patch incident direct beam nir solar radiation at local noon (W/m**2)
     real(r8), pointer :: fsr_nir_d_patch        (:)   ! patch reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_i_patch        (:)   ! patch reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsr_nir_d_ln_patch     (:)   ! patch reflected direct beam nir solar radiation at local noon (W/m**2)
     ! optional diagnostic fluxes:
     real(r8), pointer :: fsrSF_nir_d_patch      (:)   ! snow-free patch reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: fsrSF_nir_i_patch      (:)   ! snow-free patch reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: fsrSF_nir_d_ln_patch   (:)   ! snow-free patch reflected direct beam nir solar radiation at local noon (W/m**2)
     real(r8), pointer :: ssre_fsr_nir_d_patch   (:)   ! snow-free patch reflected direct beam nir solar radiation (W/m**2)
     real(r8), pointer :: ssre_fsr_nir_i_patch   (:)   ! snow-free patch reflected diffuse nir solar radiation (W/m**2)
     real(r8), pointer :: ssre_fsr_nir_d_ln_patch(:)   ! snow-free patch reflected direct beam nir solar radiation at local noon (W/m**2)

   contains 

    procedure, public :: Init

 end type solarabs_type
 type(solarabs_type), public, target, save :: solarabs_inst

contains

!------------------------------------------------------
  subroutine Init(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM solar absorbed type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(solarabs_type)          :: this

    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    !---------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
    begl = bounds%begl ; endl = bounds%endl

    allocate(this%fsa_patch              (begp:endp))              ; this%fsa_patch              (:)   = nan
    allocate(this%fsa_u_patch            (begp:endp))              ; this%fsa_u_patch            (:)   = nan
    allocate(this%fsa_r_patch            (begp:endp))              ; this%fsa_r_patch            (:)   = nan
    allocate(this%parsun_z_patch         (begp:endp,1:nlevcan))    ; this%parsun_z_patch         (:,:) = nan
    allocate(this%parsha_z_patch         (begp:endp,1:nlevcan))    ; this%parsha_z_patch         (:,:) = nan
    if(use_luna)then
        allocate(this%par240d_z_patch    (begp:endp,1:nlevcan))    ; this%par240d_z_patch        (:,:) = spval
        allocate(this%par240x_z_patch    (begp:endp,1:nlevcan))    ; this%par240x_z_patch        (:,:) = spval
        allocate(this%par24d_z_patch     (begp:endp,1:nlevcan))    ; this%par24d_z_patch         (:,:) = spval
        allocate(this%par24x_z_patch     (begp:endp,1:nlevcan))    ; this%par24x_z_patch         (:,:) = spval
    endif
    allocate(this%sabv_patch             (begp:endp))              ; this%sabv_patch             (:)   = nan
    allocate(this%sabg_patch             (begp:endp))              ; this%sabg_patch             (:)   = nan
    allocate(this%sabg_lyr_patch         (begp:endp,-nlevsno+1:1)) ; this%sabg_lyr_patch         (:,:) = nan
    allocate(this%sabg_pen_patch         (begp:endp))              ; this%sabg_pen_patch         (:)   = nan
    allocate(this%sabg_soil_patch        (begp:endp))              ; this%sabg_soil_patch        (:)   = nan
    allocate(this%sabg_snow_patch        (begp:endp))              ; this%sabg_snow_patch        (:)   = nan
    allocate(this%sabg_chk_patch         (begp:endp))              ; this%sabg_chk_patch         (:)   = nan
    allocate(this%sabs_roof_dir_lun      (begl:endl,1:numrad))     ; this%sabs_roof_dir_lun      (:,:) = nan
    allocate(this%sabs_roof_dif_lun      (begl:endl,1:numrad))     ; this%sabs_roof_dif_lun      (:,:) = nan
    allocate(this%sabs_sunwall_dir_lun   (begl:endl,1:numrad))     ; this%sabs_sunwall_dir_lun   (:,:) = nan
    allocate(this%sabs_sunwall_dif_lun   (begl:endl,1:numrad))     ; this%sabs_sunwall_dif_lun   (:,:) = nan
    allocate(this%sabs_shadewall_dir_lun (begl:endl,1:numrad))     ; this%sabs_shadewall_dir_lun (:,:) = nan
    allocate(this%sabs_shadewall_dif_lun (begl:endl,1:numrad))     ; this%sabs_shadewall_dif_lun (:,:) = nan
    allocate(this%sabs_improad_dir_lun   (begl:endl,1:numrad))     ; this%sabs_improad_dir_lun   (:,:) = nan
    allocate(this%sabs_improad_dif_lun   (begl:endl,1:numrad))     ; this%sabs_improad_dif_lun   (:,:) = nan
    allocate(this%sabs_perroad_dir_lun   (begl:endl,1:numrad))     ; this%sabs_perroad_dir_lun   (:,:) = nan
    allocate(this%sabs_perroad_dif_lun   (begl:endl,1:numrad))     ; this%sabs_perroad_dif_lun   (:,:) = nan
    allocate(this%sub_surf_abs_SW_patch  (begp:endp))              ; this%sub_surf_abs_SW_patch  (:)   = nan
    allocate(this%fsr_patch              (begp:endp))              ; this%fsr_patch              (:)   = nan
    allocate(this%fsr_nir_d_patch        (begp:endp))              ; this%fsr_nir_d_patch        (:)   = nan
    allocate(this%fsr_nir_i_patch        (begp:endp))              ; this%fsr_nir_i_patch        (:)   = nan
    allocate(this%fsr_nir_d_ln_patch     (begp:endp))              ; this%fsr_nir_d_ln_patch     (:)   = nan
    allocate(this%fsrSF_patch            (begp:endp))              ; this%fsrSF_patch            (:)   = nan
    allocate(this%fsrSF_nir_d_patch      (begp:endp))              ; this%fsrSF_nir_d_patch      (:)   = nan
    allocate(this%fsrSF_nir_i_patch      (begp:endp))              ; this%fsrSF_nir_i_patch      (:)   = nan
    allocate(this%fsrSF_nir_d_ln_patch   (begp:endp))              ; this%fsrSF_nir_d_ln_patch   (:)   = nan
    allocate(this%ssre_fsr_patch         (begp:endp))              ; this%ssre_fsr_patch         (:)   = nan
    allocate(this%ssre_fsr_nir_d_patch   (begp:endp))              ; this%ssre_fsr_nir_d_patch   (:)   = nan
    allocate(this%ssre_fsr_nir_i_patch   (begp:endp))              ; this%ssre_fsr_nir_i_patch   (:)   = nan
    allocate(this%ssre_fsr_nir_d_ln_patch(begp:endp))              ; this%ssre_fsr_nir_d_ln_patch(:)   = nan
    allocate(this%fsds_nir_d_patch       (begp:endp))              ; this%fsds_nir_d_patch       (:)   = nan
    allocate(this%fsds_nir_i_patch       (begp:endp))              ; this%fsds_nir_i_patch       (:)   = nan
    allocate(this%fsds_nir_d_ln_patch    (begp:endp))              ; this%fsds_nir_d_ln_patch    (:)   = nan

  end subroutine Init

end module SolarAbsorbedType
