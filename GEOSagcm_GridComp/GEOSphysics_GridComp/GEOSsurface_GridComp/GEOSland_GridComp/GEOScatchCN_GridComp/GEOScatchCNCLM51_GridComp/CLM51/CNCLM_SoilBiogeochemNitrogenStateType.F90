module CNCLM_SoilBiogeochemNitrogenStateType

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
  public :: init_soilbiogeochem_nitrogenstate_type

  type, public :: soilbiogeochem_nitrogenstate_type

     real(r8), pointer :: decomp_npools_vr_col         (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: decomp0_npools_vr_col        (:,:,:) ! col (gN/m3) vertically-resolved N baseline (initial value of this year) in decomposing (litter, cwd, soil) pools in dimension (col,nlev,npools)
     real(r8), pointer :: decomp_npools_vr_SASUsave_col(:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools

     real(r8), pointer :: decomp_soiln_vr_col          (:,:)   ! col (gN/m3) vertically-resolved decomposing total soil N pool

     real(r8), pointer :: sminn_vr_col                 (:,:)   ! col (gN/m3) vertically-resolved soil mineral N
     real(r8), pointer :: ntrunc_vr_col                (:,:)   ! col (gN/m3) vertically-resolved column-level sink for N truncation

     ! nitrif_denitrif
     real(r8), pointer :: smin_no3_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_no3_col                 (:)     ! col (gN/m2) soil mineral NO3 pool
     real(r8), pointer :: smin_nh4_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4_col                 (:)     ! col (gN/m2) soil mineral NH4 pool

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: decomp_npools_col            (:,:)   ! col (gN/m2)  decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: decomp_npools_1m_col         (:,:)   ! col (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
     real(r8), pointer :: sminn_col                    (:)     ! col (gN/m2) soil mineral N
     real(r8), pointer :: ntrunc_col                   (:)     ! col (gN/m2) column-level sink for N truncation
     real(r8), pointer :: cwdn_col                     (:)     ! col (gN/m2) Diagnostic: coarse woody debris N
     real(r8), pointer :: totlitn_col                  (:)     ! col (gN/m2) total litter nitrogen
     real(r8), pointer :: totsomn_col                  (:)     ! col (gN/m2) total soil organic matter nitrogen
     real(r8), pointer :: totlitn_1m_col               (:)     ! col (gN/m2) total litter nitrogen to 1 meter
     real(r8), pointer :: totsomn_1m_col               (:)     ! col (gN/m2) total soil organic matter nitrogen to 1 meter
     real(r8), pointer :: dyn_nbal_adjustments_col (:) ! (gN/m2) adjustments to each column made in this timestep via dynamic column adjustments (note: this variable only makes sense at the column-level: it is meaningless if averaged to the gridcell-level)

     ! Track adjustments to no3 and nh4 pools separately, since those aren't included in
     ! the N balance check
     real(r8), pointer :: dyn_no3bal_adjustments_col (:) ! (gN/m2) NO3 adjustments to each column made in this timestep via dynamic column area adjustments (only makes sense at the column-level: meaningless if averaged to the gridcell-level)
     real(r8), pointer :: dyn_nh4bal_adjustments_col (:) ! (gN/m2) NH4 adjustments to each column made in this timestep via dynamic column adjustments (only makes sense at the column-level: meaningless if averaged to the gridcell-level)
     real(r8)          :: totvegcthresh                  ! threshold for total vegetation carbon to zero out decomposition pools

     ! Matrix-cn
     real(r8), pointer :: matrix_cap_decomp_npools_col    (:,:)   ! col (gN/m2) N capacity in decomposing (litter, cwd, soil) N pools in dimension (col,npools)
     real(r8), pointer :: matrix_cap_decomp_npools_vr_col (:,:,:) ! col (gN/m3) vertically-resolved N capacity in decomposing (litter, cwd, soil) pools in dimension(col,nlev,npools)
     real(r8), pointer :: in_nacc                         (:,:)   ! col (gN/m3/yr) accumulated litter fall N input per year in dimension(col,nlev*npools)
     real(r8), pointer :: in_nacc_2d                      (:,:,:) ! col (gN/m3/yr) accumulated litter fall N input per year in dimension(col,nlev,npools)
     real(r8), pointer :: tran_nacc                       (:,:,:) ! col (gN/m3/yr) accumulated N transfers from j to i (col,i,j) per year in dimension(col,nlev*npools,nlev*npools)
     real(r8), pointer :: vert_up_tran_nacc               (:,:,:) ! col (gN/m3/yr) accumulated upward vertical N transport in dimension(col,nlev,npools)
     real(r8), pointer :: vert_down_tran_nacc             (:,:,:) ! col (gN/m3/yr) accumulated downward vertical N transport in dimension(col,nlev,npools)
     real(r8), pointer :: exit_nacc                       (:,:,:) ! col (gN/m3/yr) accumulated exit N in dimension(col,nlev,npools)
     real(r8), pointer :: hori_tran_nacc                  (:,:,:) ! col (gN/m3/yr) accumulated N transport between pools at the same level in dimension(col,nlev,ntransfers)
    ! type(sparse_matrix_type) :: AKXnacc                          ! col (gN/m3/yr) accumulated N transfers from j to i (col,i,j) per year in dimension(col,nlev*npools,nlev*npools) in sparse matrix type
    ! type(vector_type) :: matrix_Ninter                           ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools in dimension(col,nlev*npools) in vector type

  end type soilbiogeochem_nitrogenstate_type
  type(soilbiogeochem_nitrogenstate_type), public, target, save :: soilbiogeochem_nitrogenstate_inst

contains

!-------------------------------------------
 subroutine init_soilbiogeochem_nitrogenstate_type(bounds, nch, cncol,  this)

    !
    ! !ARGUMENTS:
    !INPUT/OUTPUT
    type(bounds_type),                     intent(in) :: bounds
    integer,                               intent(in) :: nch ! number of tiles
    real, dimension(nch,NUM_ZON,VAR_COL),  intent(in) :: cncol ! gkw: column CN restart
    type(soilbiogeochem_nitrogenstate_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    integer               :: begc,endc
    integer, dimension(8) :: decomp_npool_cncol_index = (/ 18, 19, 20, 17,25, 26, 27, 28 /)
    !-----------------------------------

    begc = bounds%begc ; endc = bounds%endc

    allocate(this%sminn_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%sminn_vr_col         (:,:) = nan
    allocate(this%ntrunc_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%ntrunc_vr_col        (:,:) = nan
    allocate(this%smin_no3_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_vr_col      (:,:) = nan
    allocate(this%smin_nh4_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_vr_col      (:,:) = nan
    allocate(this%smin_no3_col         (begc:endc))                   ; this%smin_no3_col         (:)   = nan
    allocate(this%smin_nh4_col         (begc:endc))                   ; this%smin_nh4_col         (:)   = nan
    allocate(this%cwdn_col             (begc:endc))                   ; this%cwdn_col             (:)   = nan
    allocate(this%sminn_col            (begc:endc))                   ; this%sminn_col            (:)   = nan
    allocate(this%ntrunc_col           (begc:endc))                   ; this%ntrunc_col           (:)   = nan
    allocate(this%totlitn_col          (begc:endc))                   ; this%totlitn_col          (:)   = nan
    allocate(this%totsomn_col          (begc:endc))                   ; this%totsomn_col          (:)   = nan
    allocate(this%totlitn_1m_col       (begc:endc))                   ; this%totlitn_1m_col       (:)   = nan
    allocate(this%totsomn_1m_col       (begc:endc))                   ; this%totsomn_1m_col       (:)   = nan
    allocate(this%dyn_nbal_adjustments_col (begc:endc)) ; this%dyn_nbal_adjustments_col (:) = nan
    allocate(this%dyn_no3bal_adjustments_col (begc:endc)) ; this%dyn_no3bal_adjustments_col (:) = nan
    allocate(this%dyn_nh4bal_adjustments_col (begc:endc)) ; this%dyn_nh4bal_adjustments_col (:) = nan
    allocate(this%decomp_npools_col    (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_col    (:,:) = nan
    allocate(this%decomp_npools_1m_col (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_1m_col (:,:) = nan
    if(use_soil_matrixcn)then
       allocate(this%matrix_cap_decomp_npools_col    (begc:endc,1:ndecomp_pools))   ; this%matrix_cap_decomp_npools_col    (:,:) = nan
    end if

    allocate(this%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%decomp_npools_vr_col(:,:,:)= nan
    if(use_soil_matrixcn)then
       allocate(this%matrix_cap_decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
       this%matrix_cap_decomp_npools_vr_col(:,:,:)= nan
! for matrix-spinup
       allocate(this%decomp0_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
       this%decomp0_npools_vr_col(:,:,:)= nan
       allocate(this%decomp_npools_vr_SASUsave_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
       this%decomp_npools_vr_SASUsave_col(:,:,:)= nan
       allocate(this%in_nacc(begc:endc,1:nlevdecomp*ndecomp_pools))
       this%in_nacc(:,:)= nan
       allocate(this%tran_nacc(begc:endc,1:nlevdecomp*ndecomp_pools,1:nlevdecomp*ndecomp_pools))
       this%tran_nacc(:,:,:)= nan

       allocate(this%in_nacc_2d(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%in_nacc_2d(:,:,:)= nan
       allocate(this%vert_up_tran_nacc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%vert_up_tran_nacc(:,:,:)= nan
       allocate(this%vert_down_tran_nacc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%vert_down_tran_nacc(:,:,:)= nan
       allocate(this%exit_nacc(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
       this%exit_nacc(:,:,:)= nan
       allocate(this%hori_tran_nacc(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
       this%hori_tran_nacc(:,:,:)= nan
       call this%AKXnacc%InitSM(ndecomp_pools*nlevdecomp,begc,endc,decomp_cascade_con%n_all_entries)
       call this%matrix_Ninter%InitV (ndecomp_pools*nlevdecomp,begc,endc)
    end if
    allocate(this%decomp_soiln_vr_col(begc:endc,1:nlevdecomp_full))
    this%decomp_soiln_vr_col(:,:)= nan


 ! initialize variables from restart file or set to cold start value
   n = 0
   do nc = 1,nch        ! catchment tile loop
      do nz = 1,nzone    ! CN zone loop
          n = n + 1

          this%ntrunc_vr_col (n) = cncol(nc,nz,16)
          ! jkolassa May 2022: for now nlevdecomp_full = 1; will need to add loop if we introduce more soil layers
          this%sminn_vr_col  (n,1) = cncol(nc,nz,24)
          this%sminn_col  (n) = this%sminn_vr_col(n,1)

          do np = 1,ndecomp_pools
             ! jkolassa May 2022: accounting for fact that pool order in CNCOL is different from CTSM
             this%decomp_npools_col    (n,np) = cncol(nc,nz,decomp_npool_cncol_index(np))
             this%decomp_npools_col_1m (n,np) = cncol(nc,nz,decomp_npool_cncol_index(np))
             ! jkolassa May 2022: loop has to be added below of we add more biogeochemical (or soil) layers
             this%decomp_npools_vr_col (n,1,np) cncol(nc,nz,decomp_npool_cncol_index(np))
          end do !np
      end do !nz
   end do 

 end subroutine init_soilbiogeochem_nitrogenstate_type

end CNCLM_SoilBiogeochemNitrogenStateType
