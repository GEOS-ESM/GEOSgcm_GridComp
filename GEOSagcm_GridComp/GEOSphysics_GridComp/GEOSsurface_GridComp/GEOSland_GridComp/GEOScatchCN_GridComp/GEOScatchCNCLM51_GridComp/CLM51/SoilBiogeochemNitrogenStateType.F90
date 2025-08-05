 module SoilBiogeochemNitrogenStateType

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use abortutils       , only : endrun
  use nanMod           , only : nan
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp, nlevsoi, &
                                NUM_ZON, VAR_COL
  use clm_varcon       , only : spval, dzsoi_decomp, zisoi
  use clm_varctl       , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp, use_soil_matrixcn
  use decompMod        , only : bounds_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use LandunitType                       , only : lun
  use ColumnType                         , only : col
  use landunit_varcon                    , only : istcrop, istsoil


  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

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

  contains
   
     procedure , public :: Summary
     procedure , public  :: SetTotVgCThresh
     procedure , public :: Init

  end type soilbiogeochem_nitrogenstate_type
  type(soilbiogeochem_nitrogenstate_type), public, target, save :: soilbiogeochem_nitrogenstate_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

!-------------------------------------------
 subroutine Init(this, bounds, nch, cncol)

    !
    ! !ARGUMENTS:
    !INPUT/OUTPUT
    type(bounds_type),                     intent(in) :: bounds
    integer,                               intent(in) :: nch ! number of tiles
    real, dimension(nch,NUM_ZON,VAR_COL),  intent(in) :: cncol ! gkw: column CN restart
    class(soilbiogeochem_nitrogenstate_type)          :: this
    !
    ! !LOCAL VARIABLES:
    integer               :: begc,endc
    integer               :: n, nc, nz, np, l, c
    integer, dimension(8) :: decomp_npool_cncol_index = (/ 18, 19, 20, 17,25, 26, 27, 28 /)
    logical               :: no_cn51_rst = .false.
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
       !call this%AKXnacc%InitSM(ndecomp_pools*nlevdecomp,begc,endc,decomp_cascade_con%n_all_entries)
       !call this%matrix_Ninter%InitV (ndecomp_pools*nlevdecomp,begc,endc)
    end if
    allocate(this%decomp_soiln_vr_col(begc:endc,1:nlevdecomp_full))
    this%decomp_soiln_vr_col(:,:)= nan


 ! initialize variables from restart file or set to cold start value
   n = 0
   do nc = 1,nch        ! catchment tile loop
      do nz = 1,num_zon    ! CN zone loop
          n = n + 1

          this%ntrunc_vr_col (n,1:nlevdecomp_full) = cncol(nc,nz,16)
          ! jkolassa May 2022: for now nlevdecomp_full = 1; will need to add loop if we introduce more soil layers
          this%sminn_vr_col  (n,1:nlevdecomp_full) = cncol(nc,nz,24)
          this%sminn_col  (n) = this%sminn_vr_col(n,1)

          if (no_cn51_rst) then ! jkolassa Nov 2024: when no CN51 restart file is available compute NO3 and NH4 from N
             this%smin_no3_col(n) = (1.25/2.25)*this%sminn_col(n)
             this%smin_nh4_col(n) = this%sminn_col(n)/2.25
          else
             this%smin_no3_col(n) = cncol(nc,nz,36);
             this%smin_nh4_col(n) = cncol(nc,nz,37);
          end if

          this%smin_no3_vr_col(n,1:nlevdecomp_full) = this%smin_no3_col(n)
          this%smin_nh4_vr_col(n,1:nlevdecomp_full) = this%smin_nh4_col(n)

          do np = 1,ndecomp_pools
             ! jkolassa May 2022: accounting for fact that pool order in CNCOL is different from CTSM
             this%decomp_npools_col    (n,np) = cncol(nc,nz,decomp_npool_cncol_index(np))
             this%decomp_npools_1m_col (n,np) = cncol(nc,nz,decomp_npool_cncol_index(np))
             ! jkolassa May 2022: loop has to be added below if we add more biogeochemical (or soil) layers
             this%decomp_npools_vr_col (n,1,np) = cncol(nc,nz,decomp_npool_cncol_index(np))
          end do !np
      end do !nz
   end do 

    do c = begc, endc
       l = col%landunit(c)

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          this%totlitn_col(c)    = 0._r8
          this%totsomn_col(c)    = 0._r8
          this%totlitn_1m_col(c) = 0._r8
          this%totsomn_1m_col(c) = 0._r8
          this%cwdn_col(c)       = 0._r8

       end if
    end do




 end subroutine Init

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_allc, filter_allc)
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_allc       ! number of columns in allc filter
    integer           , intent(in) :: filter_allc(:) ! filter for all active columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l     ! indices
    integer  :: fc          ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

   ! vertically integrate NO3 NH4 N2O pools
   if (use_nitrif_denitrif) then
      do fc = 1,num_allc
         c = filter_allc(fc)
         this%smin_no3_col(c) = 0._r8
         this%smin_nh4_col(c) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%smin_no3_col(c) = &
                 this%smin_no3_col(c) + &
                 this%smin_no3_vr_col(c,j) * dzsoi_decomp(j)

            this%smin_nh4_col(c) = &
                 this%smin_nh4_col(c) + &
                 this%smin_nh4_vr_col(c,j) * dzsoi_decomp(j)
          end do
       end do

    end if

   ! vertically integrate each of the decomposing N pools
   do l = 1, ndecomp_pools
      do fc = 1,num_allc
         c = filter_allc(fc)
         this%decomp_npools_col(c,l) = 0._r8
         if(use_soil_matrixcn)then
            this%matrix_cap_decomp_npools_col(c,l) = 0._r8
         end if
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%decomp_npools_col(c,l) = &
                 this%decomp_npools_col(c,l) + &
                 this%decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
            if(use_soil_matrixcn)then
               this%matrix_cap_decomp_npools_col(c,l) = &
                    this%matrix_cap_decomp_npools_col(c,l) + &
                    this%matrix_cap_decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
            end if
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp > 1) then

      do l = 1, ndecomp_pools
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%decomp_npools_1m_col(c,l) = 0._r8
         end do
      end do


      ! vertically integrate each of the decomposing n pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) <= maxdepth ) then
               do fc = 1,num_allc
                  c = filter_allc(fc)
                  this%decomp_npools_1m_col(c,l) = &
                       this%decomp_npools_1m_col(c,l) + &
                       this%decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) < maxdepth ) then
               do fc = 1,num_allc
                  c = filter_allc(fc)
                  this%decomp_npools_1m_col(c,l) = &
                       this%decomp_npools_1m_col(c,l) + &
                       this%decomp_npools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do


      ! Add soil nitrogen pools together to produce vertically-resolved decomposing total soil N pool
      if ( nlevdecomp_full > 1 ) then
         do j = 1, nlevdecomp
            do fc = 1,num_allc
               c = filter_allc(fc)
               this%decomp_soiln_vr_col(c,j) = 0._r8
            end do
         end do
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do j = 1, nlevdecomp
                  do fc = 1,num_allc
                     c = filter_allc(fc)
                     this%decomp_soiln_vr_col(c,j) = this%decomp_soiln_vr_col(c,j) + &
                          this%decomp_npools_vr_col(c,j,l)
                  end do
               end do
            end if
         end do
      end if

      ! total litter nitrogen to 1 meter (TOTLITN_1m)
      do fc = 1,num_allc
         c = filter_allc(fc)
         this%totlitn_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_allc
               c = filter_allc(fc)
               this%totlitn_1m_col(c) = &
                    this%totlitn_1m_col(c) + &
                    this%decomp_npools_1m_col(c,l)
            end do
         end if
      end do


      ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      do fc = 1,num_allc
         c = filter_allc(fc)
         this%totsomn_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_allc
               c = filter_allc(fc)
               this%totsomn_1m_col(c) = this%totsomn_1m_col(c) + &
                    this%decomp_npools_1m_col(c,l)
            end do
         end if
      end do

   endif

   ! total litter nitrogen (TOTLITN)
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%totlitn_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%totlitn_col(c) = &
                 this%totlitn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do


   ! total soil organic matter nitrogen (TOTSOMN)
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%totsomn_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_soil(l) ) then
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%totsomn_col(c) = this%totsomn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do

   ! total cwdn
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%cwdn_col(c) = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_cwd(l) ) then
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%cwdn_col(c) = this%cwdn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do


   ! total sminn
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%sminn_col(c)      = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_allc
         c = filter_allc(fc)
         this%sminn_col(c) = this%sminn_col(c) + &
              this%sminn_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! total col_ntrunc
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%ntrunc_col(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_allc
         c = filter_allc(fc)
         this%ntrunc_col(c) = this%ntrunc_col(c) + &
              this%ntrunc_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

 end subroutine Summary

  !------------------------------------------------------------------------
  subroutine SetTotVgCThresh(this, totvegcthresh)

    class(soilbiogeochem_nitrogenstate_type)           :: this
    real(r8)                              , intent(in) :: totvegcthresh
    character(len=512) :: msg

    if ( totvegcthresh <= 0.0_r8 )then
        call endrun(msg=' Error totvegcthresh is zero or negative and should be > 0'//&
               errMsg(sourcefile, __LINE__))
    end if
    this%totvegcthresh = totvegcthresh

  end subroutine SetTotVgCThresh


end module SoilBiogeochemNitrogenStateType
