module SoilBiogeochemStateType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use nanMod           , only : nan
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan, &
                                nlevsno, nlevgrnd, nlevlak, nlevsoifl
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp, nlevsoi
  use clm_varctl       , only : use_cn
  use clm_varcon       , only : spval
  use CNCLM_decompMod  , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_soilbiogeochem_state_type

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

  end type soilbiogeochem_state_type
  type(soilbiogeochem_state_type), public, target, save :: soilbiogeochem_state_inst

contains

!---------------------------------------
 subroutine init_soilbiogeochem_state_type(bounds, nch, cncol, cn5_cold_start,  this)

    !
    ! !ARGUMENTS:
    !INPUT/OUTPUT
    type(bounds_type),                     intent(in) :: bounds
    integer,                               intent(in) :: nch ! number of tiles
    real, dimension(nch,NUM_ZON,VAR_COL),  intent(in) :: cncol ! gkw: column CN restart
    logical, optional,                     intent(in) :: cn5_cold_start
    type(soilbiogeochem_state_type),       intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc,endc
    integer :: n, nc, nz, n, np
    logical :: cold_start = .false.
    !-----------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    ! check whether a cn5_cold_start option was set and change cold_start accordingly
    if (present(cn5_cold_start) .and. (cn5_cold_start==.true.)) then
       cold_start = .true.
    end if

    ! jkolassa: if cold_start is false, check that both CNCOL and CNPFT have the expected size for CNCLM50, else abort 
    if ((cold_start==.false.) .and. ((size(cncol,3).ne.var_col) .or. &
       (size(cnpft,3).ne.var_pft)))
       _ASSERT(.FALSE.,'option CNCLM50_cold_start = .FALSE. requires a CNCLM50 restart file')
    end if

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


          ! "new" variables: introduced in CNCLM50
          if (cold_start==.false.) then
             do nw = 1,nlevdecomp_full
                this%nfixation_prof_col(n,nw)    = cnpft(nc,nz,nv, XXX+(nw-1))
                this%ndep_prof_col(n,nw)         = cnpft(nc,nz,nv, XXX+(nw-1))
             end do
          elseif (cold_start) then
             this%nfixation_prof_col(n,1:nlevdecomp_full)    = 0._r8
             this%ndep_prof_col(n,1:nlevdecomp_full)    = 0._r8
          else
            _ASSERT(.FALSE.,'missing CNCLM50_cold_start setting')
          end if


          do np = 1,nlevdecomp_full
             this%fpi_vr_col(n,np) = cncol(nc,nz, 35)
          end do
      end do !nz
   end do ! nc

 end  subroutine init_soilbiogeochem_state_type
end module SoilBiogeochemStateType
