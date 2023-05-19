#include "MAPL_Generic.h"

module SoilBiogeochemStateType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan, &
                                nlevsno, nlevgrnd, nlevlak
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp, nlevsoi, &
                                VAR_COL, VAR_PFT, num_zon, num_veg, numpft
  use clm_varctl       , only : use_cn
  use clm_varcon       , only : spval
  use decompMod        , only : bounds_type
  use MAPL_ExceptionHandling

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

  public :: get_spinup_latitude_term

  ! !PUBLIC TYPES:
  type, public :: soilbiogeochem_state_type

     real(r8) , pointer :: leaf_prof_patch             (:,:)   ! (1/m) profile of leaves (vertical profiles for calculating fluxes)
     real(r8) , pointer :: froot_prof_patch            (:,:)   ! (1/m) profile of fine roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: croot_prof_patch            (:,:)   ! (1/m) profile of coarse roots (vertical profiles for calculating fluxes)
     real(r8) , pointer :: stem_prof_patch             (:,:)   ! (1/m) profile of stems (vertical profiles for calculating fluxes)
     real(r8) , pointer :: fpi_vr_col                  (:,:)   ! (no units) fraction of potential immobilization 
     real(r8) , pointer :: fpi_col                     (:)     ! (no units) fraction of potential immobilization 
     real(r8),  pointer :: fpg_col                     (:)     ! (no units) fraction of potential gpp 
     real(r8) , pointer :: rf_decomp_cascade_col       (:,:,:) ! (frac) respired fraction in decomposition step 
     real(r8) , pointer :: pathfrac_decomp_cascade_col (:,:,:) ! (frac) what fraction of C leaving a given pool passes through a given transition 
     real(r8) , pointer :: nfixation_prof_col          (:,:)   ! (1/m) profile for N fixation additions 
     real(r8) , pointer :: ndep_prof_col               (:,:)   ! (1/m) profile for N fixation additions 
     real(r8) , pointer :: som_adv_coef_col            (:,:)   ! (m2/s) SOM advective flux 
     real(r8) , pointer :: som_diffus_coef_col         (:,:)   ! (m2/s) SOM diffusivity due to bio/cryo-turbation 
     real(r8) , pointer :: plant_ndemand_col           (:)     ! column-level plant N demand

   contains

     procedure, public :: Init

  end type soilbiogeochem_state_type
  type(soilbiogeochem_state_type), public, target, save :: soilbiogeochem_state_inst

contains

!---------------------------------------
 subroutine Init(this, bounds, nch, cncol, cnpft, ityp, fveg, rc)

    !
    ! !ARGUMENTS:
    !INPUT/OUTPUT
    type(bounds_type),                     intent(in) :: bounds
    integer,                               intent(in) :: nch ! number of tiles
    integer, dimension(nch,NUM_VEG,NUM_ZON),intent(in) :: ityp ! PFT index
    real, dimension(nch,NUM_VEG,NUM_ZON),   intent(in) :: fveg    ! PFT fraction
    real, dimension(nch,NUM_ZON,VAR_COL),  intent(in) :: cncol ! gkw: column CN restart
    real, dimension(nch,NUM_ZON,NUM_VEG,VAR_PFT), intent(in) :: cnpft ! gkw: PFT CN restart
    class(soilbiogeochem_state_type)                  :: this
    integer, optional,                     intent(out) :: rc
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc,endc
    integer :: n, nc, nz, np, nv, p
    !-----------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc


    allocate(this%leaf_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%leaf_prof_patch     (:,:) = spval
    allocate(this%froot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%froot_prof_patch    (:,:) = spval
    allocate(this%croot_prof_patch    (begp:endp,1:nlevdecomp_full)) ; this%croot_prof_patch    (:,:) = spval
    allocate(this%stem_prof_patch     (begp:endp,1:nlevdecomp_full)) ; this%stem_prof_patch     (:,:) = spval
    allocate(this%fpi_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%fpi_vr_col          (:,:) = nan
    allocate(this%fpi_col             (begc:endc))                   ; this%fpi_col             (:)   = nan
    allocate(this%fpg_col             (begc:endc))                   ; this%fpg_col             (:)   = nan
    allocate(this%nfixation_prof_col  (begc:endc,1:nlevdecomp_full)) ; this%nfixation_prof_col  (:,:) = spval
    allocate(this%ndep_prof_col       (begc:endc,1:nlevdecomp_full)) ; this%ndep_prof_col       (:,:) = spval
    allocate(this%som_adv_coef_col    (begc:endc,1:nlevdecomp_full)) ; this%som_adv_coef_col    (:,:) = spval
    allocate(this%som_diffus_coef_col (begc:endc,1:nlevdecomp_full)) ; this%som_diffus_coef_col (:,:) = spval
    allocate(this%plant_ndemand_col   (begc:endc))                   ; this%plant_ndemand_col   (:)   = nan

    allocate(this%rf_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));
    this%rf_decomp_cascade_col(:,:,:) = nan

    allocate(this%pathfrac_decomp_cascade_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions));
    this%pathfrac_decomp_cascade_col(:,:,:) = nan

 ! initialize variables from restart file or set to cold start value
   n = 0
   do nc = 1,nch          ! catchment tile loop
      do nz = 1,num_zon   ! CN zone loop
          n = n + 1

          this%fpg_col(n) = cncol(nc,nz, 30)
          this%fpi_col(n) = cncol(nc,nz, 35)

          do np = 1,nlevdecomp_full
             this%fpi_vr_col(n,np) = cncol(nc,nz, 35)
          end do

          this%plant_ndemand_col(n) = 0._r8
          do p = 0,numpft  ! PFT index loop
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
                   this%plant_ndemand_col(n) = this%plant_ndemand_col(n) + cnpft(nc,nz,nv, 75)
                end if
             end do ! nv
          end do ! p
      end do !nz
   end do ! nc

 end  subroutine Init

!-----------------------------------------------
  function get_spinup_latitude_term(latitude) result(ans)

    !!DESCRIPTION:
    ! calculate a logistic function to scale spinup factors so that spinup is more accelerated in high latitude regions
    !
    ! !REVISION HISTORY
    ! charlie koven, nov. 2015
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: latitude
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans

    ans = 1._r8 + 50._r8 / ( 1._r8 + exp(-0.15_r8 * (abs(latitude) - 60._r8) ) )

    return
  end function get_spinup_latitude_term

end module SoilBiogeochemStateType
