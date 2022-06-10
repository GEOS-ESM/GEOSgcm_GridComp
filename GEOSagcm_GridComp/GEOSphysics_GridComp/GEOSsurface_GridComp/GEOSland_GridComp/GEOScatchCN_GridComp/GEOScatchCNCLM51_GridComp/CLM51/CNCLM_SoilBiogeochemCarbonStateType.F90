module CNCLM_SoilBiogeochemCarbonStateType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use nanMod           , only : nan
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp, nlevsoi
  use clm_varctl       , only : use_soil_matrixcn
  use CNCLM_decompMod  , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_soilbiogeochem_carbonstate_type

  type, public :: soilbiogeochem_carbonstate_type

     ! all c pools involved in decomposition
     real(r8), pointer :: decomp_cpools_vr_col (:,:,:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: decomp0_cpools_vr_col(:,:,:) ! (gC/m3) vertically-resolved C baseline (initial value of this year) in decomposing (litter, cwd, soil) pools in dimension (col,nlev,npools)
     real(r8), pointer :: decomp_cpools_vr_SASUsave_col(:,:,:) ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: decomp_soilc_vr_col  (:,:)   ! (gC/m3) vertically-resolved decomposing total soil c pool
     real(r8), pointer :: ctrunc_vr_col        (:,:)   ! (gC/m3) vertically-resolved column-level sink for C truncation

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: ctrunc_col              (:)     ! (gC/m2) column-level sink for C truncation
     real(r8), pointer :: totlitc_col             (:)     ! (gC/m2) total litter carbon
     real(r8), pointer :: totlitc_1m_col          (:)     ! (gC/m2) total litter carbon to 1 meter
     real(r8), pointer :: totsomc_col             (:)     ! (gC/m2) total soil organic matter carbon
     real(r8), pointer :: totsomc_1m_col          (:)     ! (gC/m2) total soil organic matter carbon to 1 meter
     real(r8), pointer :: cwdc_col                (:)     ! (gC/m2) coarse woody debris C (diagnostic)
     real(r8), pointer :: decomp_cpools_1m_col    (:,:)   ! (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
     real(r8), pointer :: decomp_cpools_col       (:,:)   ! (gC/m2)  decomposing (litter, cwd, soil) c pools
     real(r8), pointer :: dyn_cbal_adjustments_col(:)     ! (gC/m2) adjustments to each column made in this timestep via dynamic column area adjustments (note: this variable only makes sense at the column-level: it is meaningless if averaged to the gridcell-level)
     integer           :: restart_file_spinup_state       ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
     real(r8)          :: totvegcthresh                   ! threshold for total vegetation carbon to zero out decomposition pools

     ! Matrix-cn
     real(r8), pointer :: matrix_cap_decomp_cpools_col    (:,:)   ! (gC/m2) C capacity in decomposing (litter, cwd, soil) N pools in dimension (col,npools)
     real(r8), pointer :: matrix_cap_decomp_cpools_vr_col (:,:,:) ! (gC/m3) vertically-resolved C capacity in decomposing (litter, cwd, soil) pools in dimension(col,nlev,npools)
     real(r8), pointer :: in_acc                          (:,:)   ! (gC/m3/yr) accumulated litter fall C input per year in dimension(col,nlev*npools)
     real(r8), pointer :: in_acc_2d                       (:,:,:) ! (gC/m3/yr) accumulated litter fall C input per year in dimension(col,nlev,npools)
     real(r8), pointer :: tran_acc                        (:,:,:) ! (gC/m3/yr) accumulated C transfers from j to i (col,i,j) per year in dimension(col,nlev*npools,nlev*npools)
     real(r8), pointer :: vert_up_tran_acc                (:,:,:) ! (gC/m3/yr) accumulated upward vertical C transport in dimension(col,nlev,npools)
     real(r8), pointer :: vert_down_tran_acc              (:,:,:) ! (gC/m3/yr) accumulated downward vertical C transport in dimension(col,nlev,npools)
     real(r8), pointer :: exit_acc                        (:,:,:) ! (gC/m3/yr) accumulated exit C in dimension(col,nlev,npools)
     real(r8), pointer :: hori_tran_acc                   (:,:,:) ! (gC/m3/yr) accumulated C transport between pools at the same level in dimension(col,nlev,ntransfers)
    ! type(sparse_matrix_type) :: AKXcacc                          ! (gC/m3/yr) accumulated N transfers from j to i (col,i,j) per year in dimension(col,nlev*npools,nlev*npools) in sparse matrix type
    ! type(vector_type) :: matrix_Cinter                           ! (gC/m3)    vertically-resolved decomposing (litter, cwd, soil) N pools in dimension(col,nlev*npools) in vector type

  end type soilbiogeochem_carbonstate_type
  type(soilbiogeochem_carbonstate_type), public, target, save :: soilbiogeochem_carbonstate_inst

contains

!-------------------------------------------
 subroutine init_soilbiogeochem_carbonstate_type(bounds, nch, cncol,  this)

    !
    ! !ARGUMENTS:
    !INPUT/OUTPUT
    type(bounds_type),                     intent(in) :: bounds
    integer,                               intent(in) :: nch ! number of tiles
    real, dimension(nch,NUM_ZON,VAR_COL),  intent(in) :: cncol ! gkw: column CN restart
    type(soilbiogeochem_carbonstate_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer               :: begc,endc
    integer :: n, nc, nz, n
    integer, dimension(8) :: decomp_cpool_cncol_index = (/ 3, 4, 5, 2, 10, 11, 12, 13 /)
    !-----------------------------------

    begc = bounds%begc ; endc = bounds%endc

    allocate( this%decomp_cpools_col    (begc :endc,1:ndecomp_pools))   ; this%decomp_cpools_col    (:,:) = nan
    allocate( this%decomp_cpools_1m_col (begc :endc,1:ndecomp_pools))   ; this%decomp_cpools_1m_col (:,:) = nan
    if(use_soil_matrixcn)then
       allocate( this%matrix_cap_decomp_cpools_col    (begc :endc,1:ndecomp_pools))   ; this%matrix_cap_decomp_cpools_col    (:,:) = nan
    end if

    allocate( this%ctrunc_vr_col(begc :endc,1:nlevdecomp_full)) ;
    this%ctrunc_vr_col        (:,:) = nan

    allocate(this%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_cpools_vr_col(:,:,:)= nan
    !matrix-spinup
    if(use_soil_matrixcn)then
       allocate(this%matrix_cap_decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%matrix_cap_decomp_cpools_vr_col(:,:,:)= nan
       allocate(this%decomp0_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%decomp0_cpools_vr_col(:,:,:)= nan
       allocate(this%decomp_cpools_vr_SASUsave_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%decomp_cpools_vr_SASUsave_col(:,:,:)= nan
       allocate(this%in_acc(begc:endc,1:nlevdecomp*ndecomp_pools))
       this%in_acc(:,:)= nan
       allocate(this%tran_acc(begc:endc,1:nlevdecomp*ndecomp_pools,1:nlevdecomp*ndecomp_pools))
       this%tran_acc(:,:,:)= nan

       allocate(this%in_acc_2d(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%in_acc_2d(:,:,:)= nan
       allocate(this%vert_up_tran_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%vert_up_tran_acc(:,:,:)= nan
       allocate(this%vert_down_tran_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%vert_down_tran_acc(:,:,:)= nan
       allocate(this%exit_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%exit_acc(:,:,:)= nan
       allocate(this%hori_tran_acc(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
       this%hori_tran_acc(:,:,:)= nan
       ! jkolassa May 2022: comment out the two functions below as currently use_soil_matrixcn = .false.
       !call this%AKXcacc%InitSM(ndecomp_pools*nlevdecomp,begc,endc,decomp_cascade_con%n_all_entries)
       !call this%matrix_Cinter%InitV        (ndecomp_pools*nlevdecomp,begc,endc)
    end if
    allocate(this%decomp_soilc_vr_col(begc:endc,1:nlevdecomp_full))
    this%decomp_soilc_vr_col(:,:)= nan

    allocate(this%ctrunc_col     (begc :endc)) ; this%ctrunc_col     (:) = nan
    if ( .not. use_fates ) then
       allocate(this%cwdc_col       (begc :endc)) ; this%cwdc_col       (:) = nan
    endif
    allocate(this%totlitc_col    (begc :endc)) ; this%totlitc_col    (:) = nan
    allocate(this%totsomc_col    (begc :endc)) ; this%totsomc_col    (:) = nan
    allocate(this%totlitc_1m_col (begc :endc)) ; this%totlitc_1m_col (:) = nan
    allocate(this%totsomc_1m_col (begc :endc)) ; this%totsomc_1m_col (:) = nan
    allocate(this%dyn_cbal_adjustments_col (begc:endc)) ; this%dyn_cbal_adjustments_col (:) = nan

    this%restart_file_spinup_state = huge(1)

 ! initialize variables from restart file or set to cold start value
   n = 0
   do nc = 1,nch           ! catchment tile loop
      do nz = 1,num_zon    ! CN zone loop
          n = n + 1

          this%ctrunc_vr_col (n) = cncol(nc,nz,1)
          this%totlitc_col   (n) = cncol(nc,nz,15)

          do np = 1,ndecomp_pools
             ! jkolassa May 2022: accounting for fact that pool order in CNCOL is different from CTSM
             this%decomp_cpools_col    (n,np) = cncol(nc,nz,decomp_cpool_cncol_index(np))
             this%decomp_cpools_col_1m (n,np) = cncol(nc,nz,decomp_cpool_cncol_index(np))
             ! jkolassa May 2022: loop has to be added below of we add more biogeochemical (or soil) layers
             this%decomp_cpools_vr_col (n,1,np) cncol(nc,nz,decomp_cpool_cncol_index(np))
          end do !np

          ! sum soil carbon pools
          this%totsomc_col  (n) = this%decomp_cpools_col(n,5) + this%decomp_cpools_col(n,6) &
                                + this%decomp_cpools_col(n,7) + this%decomp_cpools_col(n,8)
      end do !nz
   end do ! nc

 end init_soilbiogeochem_carbonstate_type

end module CNCLM_SoilBiogeochemCarbonStateType
