module SoilBiogeochemCarbonFluxType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan, ndecomp_cascade_outtransitions
  use clm_varpar       , only : nlevdecomp_full, nlevgrnd, nlevdecomp, nlevsoi, ndecomp_pools_vr
  use clm_varctl       , only : use_fates, use_soil_matrixcn, use_vertsoilc
  use clm_varcon       , only : spval, ispval, dzsoi_decomp
  use decompMod        , only : bounds_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: soilbiogeochem_carbonflux_type

     ! fire fluxes
     real(r8), pointer :: somc_fire_col                             (:)     ! (gC/m2/s) carbon emissions due to peat burning

     ! decomposition fluxes
     real(r8), pointer :: decomp_cpools_sourcesink_col              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)  
     real(r8), pointer :: decomp_cascade_hr_vr_col                  (:,:,:) ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     real(r8), pointer :: decomp_cascade_hr_col                     (:,:)   ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
     real(r8), pointer :: decomp_cascade_ctransfer_vr_col           (:,:,:) ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_ctransfer_col              (:,:)   ! vertically-integrated (diagnostic) C transferred along decomposition cascade (gC/m2/s)
     real(r8), pointer :: decomp_k_col                              (:,:,:) ! rate constant for decomposition (1./sec)
! for soil-matrix
     real(r8), pointer :: hr_vr_col                                 (:,:)   ! (gC/m3/s) total vertically-resolved het. resp. from decomposing C pools 
     real(r8), pointer :: o_scalar_col                              (:,:)   ! fraction by which decomposition is limited by anoxia
     real(r8), pointer :: w_scalar_col                              (:,:)   ! fraction by which decomposition is limited by moisture availability
     real(r8), pointer :: t_scalar_col                              (:,:)   ! fraction by which decomposition is limited by temperature
     real(r8), pointer :: som_c_leached_col                         (:)     ! (gC/m^2/s) total SOM C loss from vertical transport 
     real(r8), pointer :: decomp_cpools_leached_col                 (:,:)   ! (gC/m^2/s) C loss from vertical transport from each decomposing C pool 
     real(r8), pointer :: decomp_cpools_transport_tendency_col      (:,:,:) ! (gC/m^3/s) C tendency due to vertical transport in decomposing C pools 

     ! nitrif_denitrif
     real(r8), pointer :: phr_vr_col                                (:,:)   ! (gC/m3/s) potential hr (not N-limited) 
     real(r8), pointer :: fphr_col                                  (:,:)   ! fraction of potential heterotrophic respiration

     real(r8), pointer :: hr_col                                    (:)     ! (gC/m2/s) total heterotrophic respiration
     real(r8), pointer :: lithr_col                                 (:)     ! (gC/m2/s) litter heterotrophic respiration 
     real(r8), pointer :: somhr_col                                 (:)     ! (gC/m2/s) soil organic matter heterotrophic res   
     real(r8), pointer :: cwdhr_col                                 (:)     ! (gC/m2/s) coarse woody debris heterotrophic res
     real(r8), pointer :: soilc_change_col                          (:)     ! (gC/m2/s) FUN used soil C


     ! fluxes to receive carbon inputs from FATES
     real(r8), pointer :: FATES_c_to_litr_lab_c_col                 (:,:)   ! total labile    litter coming from ED. gC/m3/s
     real(r8), pointer :: FATES_c_to_litr_cel_c_col                 (:,:)   ! total cellulose    litter coming from ED. gC/m3/s
     real(r8), pointer :: FATES_c_to_litr_lig_c_col                 (:,:)   ! total lignin    litter coming from ED. gC/m3/s

     ! track tradiagonal matrix  
     real(r8), pointer :: matrix_decomp_fire_k_col                  (:,:)   ! decomposition rate due to fire (gC*m3)/(gC*m3*step))
     real(r8), pointer :: tri_ma_vr                                 (:,:)   ! vertical C transfer rate in sparse matrix format (gC*m3)/(gC*m3*step))


!     type(sparse_matrix_type)         :: AKsoilc                            ! A*K for C transfers between pools
!     type(sparse_matrix_type)         :: AVsoil                             ! V for C and N transfers between soil layers
!     type(sparse_matrix_type)         :: AKfiresoil                         ! Kfire for CN transfers from soil to atm due to fire
!     type(sparse_matrix_type)         :: AKallsoilc                         ! (A*K+V-Kfire) for soil C cycle
!     integer                          :: NE_AKallsoilc                      ! Number of entries in AKallsoilc, Automatically generated by functions SPMP_*
!     integer,pointer,dimension(:)     :: RI_AKallsoilc                      ! Row numbers of entries in AKallsoilc, Automatically generated by functions SPMP_*
!     integer,pointer,dimension(:)     :: CI_AKallsoilc                      ! Column numbers of entries in AKallsoilc, Automatically generated by functions SPMP_*
!     integer,pointer,dimension(:)     :: RI_a                               ! Row numbers of all entries from AKsoilc, Automatically generated by SetValueA
!     integer,pointer,dimension(:)     :: CI_a                               ! Column numbers of all entries from AKsoilc, Automatically generated by SetValueA
!
!     type(diag_matrix_type)           :: Ksoil                              ! CN turnover rate in different soil pools and layers
!     type(diag_matrix_type)           :: Xdiagsoil                          ! Temporary C and N state variable to calculate accumulation transfers
!
!     type(vector_type)                :: matrix_Cinput                      ! C input to different soil compartments (pools and layers) (gC/m3/step)

  contains

     procedure , public  :: SetValues
     procedure , public  :: Summary
     procedure , public  :: Init

  end type soilbiogeochem_carbonflux_type
  type(soilbiogeochem_carbonflux_type), public, target, save :: soilbiogeochem_carbonflux_inst

contains

!--------------------------------------------------------------
 subroutine Init(this, bounds)

     type(bounds_type),                    intent(in)    :: bounds
     class(soilbiogeochem_carbonflux_type)               :: this
     !
     ! !LOCAL VARIABLES:
     integer           :: begp,endp
     integer           :: begc,endc,Ntrans,Ntrans_diag

     !------------------------------------------------------------------------

     begp = bounds%begp; endp = bounds%endp
     begc = bounds%begc; endc = bounds%endc

     allocate(this%t_scalar_col      (begc:endc,1:nlevdecomp_full)); this%t_scalar_col      (:,:) =spval
     allocate(this%w_scalar_col      (begc:endc,1:nlevdecomp_full)); this%w_scalar_col      (:,:) =spval
     allocate(this%o_scalar_col      (begc:endc,1:nlevdecomp_full)); this%o_scalar_col      (:,:) =spval
     allocate(this%phr_vr_col        (begc:endc,1:nlevdecomp_full)); this%phr_vr_col        (:,:) =nan
     allocate(this%fphr_col          (begc:endc,1:nlevgrnd))       ; this%fphr_col          (:,:) =nan
     allocate(this%som_c_leached_col (begc:endc))                  ; this%som_c_leached_col (:)   =spval
     allocate(this%somc_fire_col     (begc:endc))                  ; this%somc_fire_col     (:)   =nan
     allocate(this%hr_vr_col         (begc:endc,1:nlevdecomp_full)); this%hr_vr_col         (:,:) =nan

     allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
     this%decomp_cpools_sourcesink_col(:,:,:)= nan

     allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
     this%decomp_cascade_hr_vr_col(:,:,:)= spval

     allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))
     this%decomp_cascade_hr_col(:,:)= nan

     allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
     this%decomp_cascade_ctransfer_vr_col(:,:,:)= nan

     allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))
     this%decomp_cascade_ctransfer_col(:,:)= nan

     allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
     this%decomp_k_col(:,:,:)= spval

     allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))
     this%decomp_cpools_leached_col(:,:)= nan

     allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
     this%decomp_cpools_transport_tendency_col(:,:,:)= nan

     allocate(this%hr_col                  (begc:endc)) ; this%hr_col                  (:) = nan
     allocate(this%lithr_col               (begc:endc)) ; this%lithr_col               (:) = nan
     allocate(this%somhr_col               (begc:endc)) ; this%somhr_col               (:) = nan
     allocate(this%cwdhr_col               (begc:endc)) ; this%cwdhr_col               (:) = nan
     allocate(this%soilc_change_col        (begc:endc)) ; this%soilc_change_col        (:) = nan

!     if(use_soil_matrixcn)then
!        allocate(this%matrix_decomp_fire_k_col(begc:endc,1:nlevdecomp*ndecomp_pools));  this%matrix_decomp_fire_k_col(:,:)= nan
!        Ntrans = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp
!        call this%AKsoilc%InitSM                (ndecomp_pools*nlevdecomp,begc,endc,Ntrans+ndecomp_pools*nlevdecomp)
!        call this%AVsoil%InitSM                 (ndecomp_pools*nlevdecomp,begc,endc,decomp_cascade_con%Ntri_setup)
!        call this%AKfiresoil%InitSM             (ndecomp_pools*nlevdecomp,begc,endc,ndecomp_pools*nlevdecomp)
!        call this%AKallsoilc%InitSM             (ndecomp_pools*nlevdecomp,begc,endc,Ntrans+decomp_cascade_con%Ntri_setup+nlevdecomp)
!        this%NE_AKallsoilc = Ntrans+ndecomp_pools*nlevdecomp+decomp_cascade_con%Ntri_setup+ndecomp_pools*nlevdecomp
!        allocate(this%RI_AKallsoilc(1:this%NE_AKallsoilc)); this%RI_AKallsoilc(1:this%NE_AKallsoilc)=-9999
!        allocate(this%CI_AKallsoilc(1:this%NE_AKallsoilc)); this%CI_AKallsoilc(1:this%NE_AKallsoilc)=-9999
!        Ntrans_diag = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp+ndecomp_pools_vr
!        allocate(this%RI_a(1:Ntrans_diag)); this%RI_a(1:Ntrans_diag) = -9999
!        allocate(this%CI_a(1:Ntrans_diag)); this%CI_a(1:Ntrans_diag) = -9999
!        call this%Ksoil%InitDM                  (ndecomp_pools*nlevdecomp,begc,endc)
!        call this%Xdiagsoil%InitDM              (ndecomp_pools*nlevdecomp,begc,endc)
!        call this%matrix_Cinput%InitV(ndecomp_pools*nlevdecomp,begc,endc)
!     end if
     if(use_soil_matrixcn .and. use_vertsoilc)then
        allocate(this%tri_ma_vr(begc:endc,1:decomp_cascade_con%Ntri_setup))
     else
        allocate(this%tri_ma_vr(1,1)); this%tri_ma_vr(:,:) = nan
     end if
     if ( use_fates ) then
        ! initialize these variables to be zero rather than a bad number since they are not zeroed every timestep (due to a need for them to persist)

        allocate(this%FATES_c_to_litr_lab_c_col(begc:endc,1:nlevdecomp_full))
        this%FATES_c_to_litr_lab_c_col(begc:endc,1:nlevdecomp_full) = 0._r8

        allocate(this%FATES_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full))
        this%FATES_c_to_litr_cel_c_col(begc:endc,1:nlevdecomp_full) = 0._r8

        allocate(this%FATES_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full))
        this%FATES_c_to_litr_lig_c_col(begc:endc,1:nlevdecomp_full) = 0._r8

     endif

 end subroutine Init

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon fluxes
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_carbonflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k,l    ! indices
    !------------------------------------------------------------------------

    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_hr_col(i,l)             = value_column
             this%decomp_cascade_hr_vr_col(i,j,l)        = value_column
             this%decomp_cascade_ctransfer_col(i,l)      = value_column
             this%decomp_cascade_ctransfer_vr_col(i,j,l) = value_column
             this%decomp_k_col(i,j,l)                    = value_column
          end do
       end do
    end do


    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cpools_leached_col(i,k) = value_column
       end do
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_transport_tendency_col(i,j,k) = value_column
             this%decomp_cpools_sourcesink_col(i,j,k)         = value_column
          end do
       end do
    end do

! for matrix 
    if(use_soil_matrixcn)then
       do k = 1, ndecomp_pools
          do j = 1, nlevdecomp
             do fi = 1,num_column
                i = filter_column(fi)
                this%matrix_decomp_fire_k_col(i,j+nlevdecomp*(k-1)) = value_column
             end do
          end do
       end do
      ! call this%matrix_Cinput%SetValueV_scaler(num_column,filter_column(1:num_column),value_column)
       ! IMPORTANT NOTE: Although it looks like the following if appears to be
       ! backwards (it should be 'if use_versoilc'), fixing it causes Carbon 
       ! balance checks to fail. EBK 10/21/2019
       ! Both use_vertsoilc and .not. use_vertsoilc should reset tri_ma_vr to 0. 
       ! Because single soil layer still add V matrix but as a zero matrix. CL 10/23/2019
        if(use_vertsoilc)then
          do k = 1,decomp_cascade_con%Ntri_setup
             do fi = 1,num_column
                i = filter_column(fi)
                this%tri_ma_vr(i,k) = value_column
             end do
          end do
        end if
    end if
    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%hr_vr_col(i,j) = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)
       this%hr_col(i)            = value_column
       this%somc_fire_col(i)     = value_column
       this%som_c_leached_col(i) = value_column
       this%somhr_col(i)         = value_column
       this%cwdhr_col(i)         = value_column
       this%lithr_col(i)         = value_column
       this%soilc_change_col(i)  = value_column
    end do

    ! NOTE: do not zero the fates to BGC C flux variables since they need to persist from the daily fates timestep s to the half-hourly BGC timesteps.  I.e. FATES_c_to_litr_lab_c_col, FATES_c_to_litr_cel_c_col, FATES_c_to_litr_lig_c_col

 end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, column-level carbon summary calculations
    !
    ! !USES:
    ! !ARGUMENTS:
    class(soilbiogeochem_carbonflux_type)           :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l
    integer  :: fc
    !-----------------------------------------------------------------------

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%som_c_leached_col(c) = 0._r8
    end do

    ! vertically integrate HR and decomposition cascade fluxes
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cascade_hr_col(c,k) = &
                  this%decomp_cascade_hr_col(c,k) + &
                  this%decomp_cascade_hr_vr_col(c,j,k) * dzsoi_decomp(j)

             this%decomp_cascade_ctransfer_col(c,k) = &
                  this%decomp_cascade_ctransfer_col(c,k) + &
                  this%decomp_cascade_ctransfer_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do

    ! total heterotrophic respiration, vertically resolved (HR)
    do j = 1,nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%hr_vr_col(c,j) = 0._r8
       end do
    end do
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%hr_vr_col(c,j) = &
                  this%hr_vr_col(c,j) + &
                  this%decomp_cascade_hr_vr_col(c,j,k)
          end do
       end do
    end do

    ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_cpools_leached_col(c,l) = 0._r8
       end do
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_leached_col(c,l) = this%decomp_cpools_leached_col(c,l) + &
                  this%decomp_cpools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%som_c_leached_col(c) = this%som_c_leached_col(c) + this%decomp_cpools_leached_col(c,l)
       end do
    end do


    ! soil organic matter heterotrophic respiration 
       associate(is_soil => decomp_cascade_con%is_soil) ! TRUE => pool is a soil pool  
         do k = 1, ndecomp_cascade_transitions
            if ( is_soil(decomp_cascade_con%cascade_donor_pool(k)) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%somhr_col(c) = this%somhr_col(c) + this%decomp_cascade_hr_col(c,k)
               end do
            end if
         end do
       end associate

    ! litter heterotrophic respiration (LITHR)
       associate(is_litter => decomp_cascade_con%is_litter) ! TRUE => pool is a litter pool
         do k = 1, ndecomp_cascade_transitions
            if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%lithr_col(c) = this%lithr_col(c) + this%decomp_cascade_hr_col(c,k)
               end do
            end if
         end do
       end associate

    ! coarse woody debris heterotrophic respiration (CWDHR)
    associate(is_cwd => decomp_cascade_con%is_cwd)  ! TRUE => pool is a cwd pool
      do k = 1, ndecomp_cascade_transitions
         if ( is_cwd(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%cwdhr_col(c) = this%cwdhr_col(c) + this%decomp_cascade_hr_col(c,k)
            end do
         end if
      end do
    end associate

    ! total heterotrophic respiration (HR)
       do fc = 1,num_soilc
          c = filter_soilc(fc)

          this%hr_col(c) = &
               this%lithr_col(c) + &
               this%cwdhr_col(c) + &
               this%somhr_col(c)

       end do

  end subroutine Summary

end module SoilBiogeochemCarbonFluxType

