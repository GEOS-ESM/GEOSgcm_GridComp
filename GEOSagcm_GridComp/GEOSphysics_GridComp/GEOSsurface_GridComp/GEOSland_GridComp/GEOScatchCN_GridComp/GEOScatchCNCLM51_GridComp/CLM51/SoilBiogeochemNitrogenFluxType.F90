module SoilBiogeochemNitrogenFluxType

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use nanMod           , only : nan
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, ndecomp_cascade_outtransitions
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp, ndecomp_pools_vr
  use clm_varctl       , only : use_nitrif_denitrif, use_vertsoilc, use_crop, use_soil_matrixcn
  use clm_varcon       , only : spval, dzsoi_decomp
  use decompMod        , only : bounds_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: SoilBiogeochem_nitrogenflux_type

     ! deposition fluxes
     real(r8), pointer :: ndep_to_sminn_col                         (:)     ! col atmospheric N deposition to soil mineral N (gN/m2/s)
     real(r8), pointer :: nfix_to_sminn_col                         (:)     ! col symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
     real(r8), pointer :: ffix_to_sminn_col                         (:)     ! col free living N fixation to soil mineral N (gN/m2/s)  
     real(r8), pointer :: fert_to_sminn_col                         (:)     ! col fertilizer N to soil mineral N (gN/m2/s)
     real(r8), pointer :: soyfixn_to_sminn_col                      (:)     ! col soybean fixation to soil mineral N (gN/m2/s)

     ! decomposition fluxes
     real(r8), pointer :: decomp_cascade_ntransfer_vr_col           (:,:,:) ! col vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_ntransfer_col              (:,:)   ! col vert-int (diagnostic) transfer of N from donor to receiver pool along decomp. cascade (gN/m2/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_vr_col          (:,:,:) ! col vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_col             (:,:)   ! col vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)

     ! Used to update concentrations concurrently with vertical transport
     ! vertically-resolved immobilization fluxes
     real(r8), pointer :: potential_immob_vr_col                    (:,:)   ! col vertically-resolved potential N immobilization (gN/m3/s) at each level
     real(r8), pointer :: potential_immob_col                       (:)     ! col vert-int (diagnostic) potential N immobilization (gN/m2/s)
     real(r8), pointer :: actual_immob_vr_col                       (:,:)   ! col vertically-resolved actual N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_col                          (:)     ! col vert-int (diagnostic) actual N immobilization (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral N (gN/m3/s)
     real(r8), pointer :: sminn_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
     real(r8), pointer :: supplement_to_sminn_vr_col                (:,:)   ! col vertically-resolved supplemental N supply (gN/m3/s)
     real(r8), pointer :: supplement_to_sminn_col                   (:)     ! col vert-int (diagnostic) supplemental N supply (gN/m2/s)
     real(r8), pointer :: gross_nmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of N mineralization (gN/m3/s)
     real(r8), pointer :: gross_nmin_col                            (:)     ! col vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
     real(r8), pointer :: net_nmin_vr_col                           (:,:)   ! col vertically-resolved net rate of N mineralization (gN/m3/s)
     real(r8), pointer :: net_nmin_col                              (:)     ! col vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_col                    (:)     ! col total soil N uptake of FUN        (gN/m2/s)
     ! ---------- NITRIF_DENITRIF  ---------------------

     ! nitrification / denitrification fluxes
     real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
     real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux
     real(r8), pointer :: f_nit_col                                 (:)     ! col (gN/m2/s) soil nitrification flux
     real(r8), pointer :: f_denit_col                               (:)     ! col (gN/m2/s) soil denitrification flux

     real(r8), pointer :: pot_f_nit_vr_col                          (:,:)   ! col (gN/m3/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_vr_col                        (:,:)   ! col (gN/m3/s) potential soil denitrification flux
     real(r8), pointer :: pot_f_nit_col                             (:)     ! col (gN/m2/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_col                           (:)     ! col (gN/m2/s) potential soil denitrification flux
     real(r8), pointer :: n2_n2o_ratio_denit_vr_col                 (:,:)   ! col ratio of N2 to N2O production by denitrification [gN/gN]
     real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_denit_col                           (:)     ! col flux of N2o from denitrification [gN/m^2/s]
     real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_nit_col                             (:)     ! col flux of N2o from nitrification [gN/m^2/s]

     ! immobilization / uptake fluxes
     real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
     real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
     real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
     real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)
     real(r8), pointer :: actual_immob_no3_col                      (:)     ! col actual immobilization of NO3 (gN/m2/s)
     real(r8), pointer :: actual_immob_nh4_col                      (:)     ! col actual immobilization of NH4 (gN/m2/s)
     real(r8), pointer :: smin_no3_to_plant_col                     (:)     ! col plant uptake of soil NO3 (gN/m2/s)
     real(r8), pointer :: smin_nh4_to_plant_col                     (:)     ! col plant uptake of soil Nh4 (gN/m2/s)

     ! leaching fluxes
     real(r8), pointer :: smin_no3_leached_vr_col                   (:,:)   ! col vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
     real(r8), pointer :: smin_no3_leached_col                      (:)     ! col soil mineral NO3 pool loss to leaching (gN/m2/s)
     real(r8), pointer :: smin_no3_runoff_vr_col                    (:,:)   ! col vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
     real(r8), pointer :: smin_no3_runoff_col                       (:)     ! col soil mineral NO3 pool loss to runoff (gN/m2/s)

     ! nitrification /denitrification diagnostic quantities
     real(r8), pointer :: smin_no3_massdens_vr_col                  (:,:)   ! col (ugN / g soil) soil nitrate concentration
     real(r8), pointer :: soil_bulkdensity_col                      (:,:)   ! col (kg soil / m3) bulk density of soil
     real(r8), pointer :: k_nitr_t_vr_col                           (:,:)
     real(r8), pointer :: k_nitr_ph_vr_col                          (:,:)
     real(r8), pointer :: k_nitr_h2o_vr_col                         (:,:)
     real(r8), pointer :: k_nitr_vr_col                             (:,:)
     real(r8), pointer :: wfps_vr_col                               (:,:)
     real(r8), pointer :: fmax_denit_carbonsubstrate_vr_col         (:,:)
     real(r8), pointer :: fmax_denit_nitrate_vr_col                 (:,:)
     real(r8), pointer :: f_denit_base_vr_col                       (:,:)   ! col nitrification and denitrification fluxes
     real(r8), pointer :: diffus_col                                (:,:)   ! col diffusivity (m2/s)
     real(r8), pointer :: ratio_k1_col                              (:,:)
     real(r8), pointer :: ratio_no3_co2_col                         (:,:)
     real(r8), pointer :: soil_co2_prod_col                         (:,:)
     real(r8), pointer :: fr_WFPS_col                               (:,:)

     real(r8), pointer :: r_psi_col                                 (:,:)
     real(r8), pointer :: anaerobic_frac_col                        (:,:)
     real(r8), pointer :: sminn_to_plant_fun_no3_vr_col             (:,:)   ! col total layer no3 uptake of FUN     (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_nh4_vr_col             (:,:)   ! col total layer nh4 uptake of FUN     (gN/m2/s)
     !----------- no NITRIF_DENITRIF--------------


     ! denitrification fluxes
     real(r8), pointer :: sminn_to_denit_decomp_cascade_vr_col      (:,:,:) ! col vertically-resolved denitrification along decomp cascade (gN/m3/s) 
     real(r8), pointer :: sminn_to_denit_decomp_cascade_col         (:,:)   ! col vertically-integrated (diagnostic) denitrification along decomp cascade (gN/m2/s) 
     real(r8), pointer :: sminn_to_denit_excess_vr_col              (:,:)   ! col vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
     real(r8), pointer :: sminn_to_denit_excess_col                 (:)     ! col vertically-integrated (diagnostic) denitrification from excess mineral N pool (gN/m2/s)

     ! leaching fluxes
     real(r8), pointer :: sminn_leached_vr_col                      (:,:)   ! col vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
     real(r8), pointer :: sminn_leached_col                         (:)     ! col soil mineral N pool loss to leaching (gN/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: denit_col                                 (:)     ! col total rate of denitrification (gN/m2/s)
     real(r8), pointer :: ninputs_col                               (:)     ! col column-level N inputs (gN/m2/s)
     real(r8), pointer :: noutputs_col                              (:)     ! col column-level N outputs (gN/m2/s)
     real(r8), pointer :: som_n_leached_col                         (:)     ! col total SOM N loss from vertical transport (gN/m^2/s)
     real(r8), pointer :: decomp_npools_leached_col                 (:,:)   ! col N loss from vertical transport from each decomposing N pool (gN/m^2/s)
     real(r8), pointer :: decomp_npools_transport_tendency_col      (:,:,:) ! col N tendency due to vertical transport in decomposing N pools (gN/m^3/s)

     ! all n pools involved in decomposition
     real(r8), pointer :: decomp_npools_sourcesink_col              (:,:,:) ! col (gN/m3) change in decomposing n pools 
                                                                            ! (sum of all additions and subtractions from stateupdate1).  
          real(r8), pointer :: sminn_to_plant_fun_vr_col                 (:,:)   ! col total layer soil N uptake of FUN  (gN/m2/s)

     ! track tradiagonal matrix  
!     type(sparse_matrix_type)     :: AKsoiln           ! A*K for N transfers between pools
!     type(sparse_matrix_type)     :: AKallsoiln        ! (A*K+V-Kfire) for soil N cycle
     integer                      :: NE_AKallsoiln     ! Number of non-zero entries in AKallsoiln. Automatically generated by functions SPMP_*
     integer,pointer,dimension(:) :: RI_AKallsoiln     ! Row numbers of entries in AKallsoiln. Automatically generated by functions in SPMP_*
     integer,pointer,dimension(:) :: CI_AKallsoiln     ! Column numbers of entries in AKallsoiln, Automatically generated by functions in SPMP_*
     integer,pointer,dimension(:) :: RI_na             ! Row numbers of all entries from AKsoiln. Automatically generated by SetValueA
     integer,pointer,dimension(:) :: CI_na             ! Column numbers of all entries from AKsoiln. Automatically generated by SetValueA
!     type(vector_type)            :: matrix_Ninput     ! N input to different soil compartments (pools and layers) (gN/m3/step)

  contains 

     procedure , public  :: SetValues
     procedure , public  :: Summary
     procedure , public  :: Init

  end type soilbiogeochem_nitrogenflux_type
  type(soilbiogeochem_nitrogenflux_type), public, target, save :: soilbiogeochem_nitrogenflux_inst

contains

!--------------------------------------------------------------
 subroutine Init(this, bounds)

     !ARGUMENTS
     implicit none
     !INPUT/OUTPUT
     type(bounds_type),                      intent(in)    :: bounds
     class(soilbiogeochem_nitrogenflux_type)               :: this
    !
    ! !LOCAL VARIABLES:
    integer           :: begc,endc,Ntrans,Ntrans_diag
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    allocate(this%ndep_to_sminn_col                 (begc:endc))                   ; this%ndep_to_sminn_col          (:)   = nan
    allocate(this%nfix_to_sminn_col                 (begc:endc))                   ; this%nfix_to_sminn_col          (:)   = nan
    allocate(this%ffix_to_sminn_col                 (begc:endc))                   ; this%ffix_to_sminn_col          (:)   = nan
    allocate(this%fert_to_sminn_col                 (begc:endc))                   ; this%fert_to_sminn_col          (:)   = nan
    allocate(this%soyfixn_to_sminn_col              (begc:endc))                   ; this%soyfixn_to_sminn_col       (:)   = nan
    allocate(this%sminn_to_plant_col                (begc:endc))                   ; this%sminn_to_plant_col         (:)   = nan
    allocate(this%potential_immob_col               (begc:endc))                   ; this%potential_immob_col        (:)   = nan
    allocate(this%actual_immob_col                  (begc:endc))                   ; this%actual_immob_col           (:)   = nan
    allocate(this%gross_nmin_col                    (begc:endc))                   ; this%gross_nmin_col             (:)   = nan
    allocate(this%net_nmin_col                      (begc:endc))                   ; this%net_nmin_col               (:)   = nan
    allocate(this%denit_col                         (begc:endc))                   ; this%denit_col                  (:)   = nan
    allocate(this%supplement_to_sminn_col           (begc:endc))                   ; this%supplement_to_sminn_col    (:)   = nan
    allocate(this%ninputs_col                       (begc:endc))                   ; this%ninputs_col                (:)   = nan
    allocate(this%noutputs_col                      (begc:endc))                   ; this%noutputs_col               (:)   = nan
    allocate(this%som_n_leached_col                 (begc:endc))                   ; this%som_n_leached_col          (:)   = nan


    allocate(this%r_psi_col                         (begc:endc,1:nlevdecomp_full)) ; this%r_psi_col                  (:,:) = spval
    allocate(this%anaerobic_frac_col                (begc:endc,1:nlevdecomp_full)) ; this%anaerobic_frac_col         (:,:) = spval
    allocate(this%potential_immob_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_vr_col     (:,:) = nan
    allocate(this%actual_immob_vr_col               (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_vr_col        (:,:) = nan
    allocate(this%sminn_to_plant_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_vr_col      (:,:) = nan
    allocate(this%supplement_to_sminn_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminn_vr_col (:,:) = nan
    allocate(this%gross_nmin_vr_col                 (begc:endc,1:nlevdecomp_full)) ; this%gross_nmin_vr_col          (:,:) = nan
    allocate(this%net_nmin_vr_col                   (begc:endc,1:nlevdecomp_full)) ; this%net_nmin_vr_col            (:,:) = nan
    allocate(this%sminn_to_plant_fun_col            (begc:endc))                   ; this%sminn_to_plant_fun_col     (:)     = nan
    allocate(this%sminn_to_plant_fun_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_fun_vr_col  (:,:)   = nan
    allocate(this%sminn_to_plant_fun_no3_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_fun_no3_vr_col(:,:) = nan
    allocate(this%sminn_to_plant_fun_nh4_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_fun_nh4_vr_col(:,:) = nan
    allocate(this%f_nit_vr_col                      (begc:endc,1:nlevdecomp_full)) ; this%f_nit_vr_col               (:,:) = nan
    allocate(this%f_denit_vr_col                    (begc:endc,1:nlevdecomp_full)) ; this%f_denit_vr_col             (:,:) = nan
    allocate(this%smin_no3_leached_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_leached_vr_col    (:,:) = nan
    allocate(this%smin_no3_leached_col              (begc:endc))                   ; this%smin_no3_leached_col       (:)   = nan
    allocate(this%smin_no3_runoff_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_runoff_vr_col     (:,:) = nan
    allocate(this%smin_no3_runoff_col               (begc:endc))                   ; this%smin_no3_runoff_col        (:)   = nan
    allocate(this%pot_f_nit_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%pot_f_nit_vr_col           (:,:) = nan
    allocate(this%pot_f_nit_col                     (begc:endc))                   ; this%pot_f_nit_col              (:)   = nan
    allocate(this%pot_f_denit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%pot_f_denit_vr_col         (:,:) = nan
    allocate(this%pot_f_denit_col                   (begc:endc))                   ; this%pot_f_denit_col            (:)   = nan
    allocate(this%actual_immob_no3_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_no3_vr_col    (:,:) = nan
    allocate(this%actual_immob_nh4_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_nh4_vr_col    (:,:) = nan
    allocate(this%smin_no3_to_plant_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_to_plant_vr_col   (:,:) = nan
    allocate(this%smin_nh4_to_plant_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_to_plant_vr_col   (:,:) = nan
    allocate(this%f_nit_col                         (begc:endc))                   ; this%f_nit_col                  (:)   = nan
    allocate(this%f_denit_col                       (begc:endc))                   ; this%f_denit_col                (:)   = nan
    allocate(this%n2_n2o_ratio_denit_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%n2_n2o_ratio_denit_vr_col  (:,:) = nan
    allocate(this%f_n2o_denit_col                   (begc:endc))                   ; this%f_n2o_denit_col            (:)   = nan
    allocate(this%f_n2o_denit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_denit_vr_col         (:,:) = nan
    allocate(this%f_n2o_nit_col                     (begc:endc))                   ; this%f_n2o_nit_col              (:)   = nan
    allocate(this%f_n2o_nit_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_nit_vr_col           (:,:) = nan


    allocate(this%smin_no3_massdens_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_massdens_vr_col   (:,:) = nan
    allocate(this%soil_bulkdensity_col              (begc:endc,1:nlevdecomp_full)) ; this%soil_bulkdensity_col       (:,:) = nan
    allocate(this%k_nitr_t_vr_col                   (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_t_vr_col            (:,:) = nan
    allocate(this%k_nitr_ph_vr_col                  (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_ph_vr_col           (:,:) = nan
    allocate(this%k_nitr_h2o_vr_col                 (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_h2o_vr_col          (:,:) = nan
    allocate(this%k_nitr_vr_col                     (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_vr_col              (:,:) = nan
    allocate(this%wfps_vr_col                       (begc:endc,1:nlevdecomp_full)) ; this%wfps_vr_col                (:,:) = nan
    allocate(this%f_denit_base_vr_col               (begc:endc,1:nlevdecomp_full)) ; this%f_denit_base_vr_col        (:,:) = nan
    allocate(this%diffus_col                        (begc:endc,1:nlevdecomp_full)) ; this%diffus_col                 (:,:) = spval
    allocate(this%ratio_k1_col                      (begc:endc,1:nlevdecomp_full)) ; this%ratio_k1_col               (:,:) = nan
    allocate(this%ratio_no3_co2_col                 (begc:endc,1:nlevdecomp_full)) ; this%ratio_no3_co2_col          (:,:) = spval
    allocate(this%soil_co2_prod_col                 (begc:endc,1:nlevdecomp_full)) ; this%soil_co2_prod_col          (:,:) = nan
    allocate(this%fr_WFPS_col                       (begc:endc,1:nlevdecomp_full)) ; this%fr_WFPS_col                (:,:) = spval

    allocate(this%fmax_denit_carbonsubstrate_vr_col (begc:endc,1:nlevdecomp_full)) ;
    this%fmax_denit_carbonsubstrate_vr_col (:,:) = nan
    allocate(this%fmax_denit_nitrate_vr_col         (begc:endc,1:nlevdecomp_full)) ;
    this%fmax_denit_nitrate_vr_col         (:,:) = nan

    allocate(this%decomp_cascade_ntransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_sminn_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_ntransfer_col      (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%decomp_cascade_sminn_flux_col     (begc:endc,1:ndecomp_cascade_transitions                   ))

    this%decomp_cascade_ntransfer_vr_col  (:,:,:) = nan
    this%decomp_cascade_sminn_flux_vr_col (:,:,:) = nan
    this%decomp_cascade_ntransfer_col     (:,:)   = nan
    this%decomp_cascade_sminn_flux_col    (:,:)   = nan

    allocate(this%sminn_to_denit_decomp_cascade_vr_col (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%sminn_to_denit_decomp_cascade_col    (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%sminn_to_denit_excess_vr_col         (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminn_to_denit_excess_col            (begc:endc                                                 ))
    allocate(this%sminn_leached_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminn_leached_col                    (begc:endc                                                 ))
    allocate(this%decomp_npools_leached_col            (begc:endc,1:ndecomp_pools                                 ))
    allocate(this%decomp_npools_transport_tendency_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))

    this%sminn_to_denit_decomp_cascade_vr_col (:,:,:) = nan
    this%sminn_to_denit_decomp_cascade_col    (:,:)   = nan
    this%sminn_to_denit_excess_vr_col         (:,:)   = nan
    this%sminn_to_denit_excess_col            (:)     = nan
    this%sminn_leached_vr_col                 (:,:)   = nan
    this%sminn_leached_col                    (:)     = nan
    this%decomp_npools_leached_col            (:,:)   = nan
    this%decomp_npools_transport_tendency_col (:,:,:) = nan

    allocate(this%decomp_npools_sourcesink_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_npools_sourcesink_col (:,:,:) = nan
    if(use_soil_matrixcn)then

       Ntrans = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp
 !      call this%AKsoiln%InitSM                (ndecomp_pools*nlevdecomp,begc,endc,Ntrans+ndecomp_pools*nlevdecomp)
 !      call this%AKallsoiln%InitSM             (ndecomp_pools*nlevdecomp,begc,endc,Ntrans+decomp_cascade_con%Ntri_setup+nlevdecomp)
        this%NE_AKallsoiln = (Ntrans+nlevdecomp*ndecomp_pools) + (Ntrans+decomp_cascade_con%Ntri_setup + nlevdecomp) + (ndecomp_pools*nlevdecomp)
        allocate(this%RI_AKallsoiln(1:this%NE_AKallsoiln)); this%RI_AKallsoiln(1:this%NE_AKallsoiln)=-9999
        allocate(this%CI_AKallsoiln(1:this%NE_AKallsoiln)); this%CI_AKallsoiln(1:this%NE_AKallsoiln)=-9999
        Ntrans_diag = (ndecomp_cascade_transitions-ndecomp_cascade_outtransitions)*nlevdecomp+ndecomp_pools_vr
        allocate(this%RI_na(1:Ntrans_diag)); this%RI_na(1:Ntrans_diag) = -9999
        allocate(this%CI_na(1:Ntrans_diag)); this%CI_na(1:Ntrans_diag) = -9999
 !      call this%matrix_Ninput%InitV (ndecomp_pools*nlevdecomp,begc,endc)
    end if

 end subroutine Init

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class(soilbiogeochem_nitrogenflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          if (.not. use_nitrif_denitrif) then
             this%sminn_to_denit_excess_vr_col(i,j)      = value_column
             this%sminn_leached_vr_col(i,j)              = value_column
             this%sminn_to_plant_fun_vr_col(i,j)         = value_column
          else
             this%f_nit_vr_col(i,j)                      = value_column
             this%f_denit_vr_col(i,j)                    = value_column
             this%smin_no3_leached_vr_col(i,j)           = value_column
             this%smin_no3_runoff_vr_col(i,j)            = value_column
             this%n2_n2o_ratio_denit_vr_col(i,j)         = value_column
             this%pot_f_nit_vr_col(i,j)                  = value_column
             this%pot_f_denit_vr_col(i,j)                = value_column
             this%actual_immob_no3_vr_col(i,j)           = value_column
             this%actual_immob_nh4_vr_col(i,j)           = value_column
             this%smin_no3_to_plant_vr_col(i,j)          = value_column
             this%smin_nh4_to_plant_vr_col(i,j)          = value_column
             this%f_n2o_denit_vr_col(i,j)                = value_column
             this%f_n2o_nit_vr_col(i,j)                  = value_column

             this%smin_no3_massdens_vr_col(i,j)          = value_column
             this%k_nitr_t_vr_col(i,j)                   = value_column
             this%k_nitr_ph_vr_col(i,j)                  = value_column
             this%k_nitr_h2o_vr_col(i,j)                 = value_column
             this%k_nitr_vr_col(i,j)                     = value_column
             this%wfps_vr_col(i,j)                       = value_column
             this%fmax_denit_carbonsubstrate_vr_col(i,j) = value_column
             this%fmax_denit_nitrate_vr_col(i,j)         = value_column
             this%f_denit_base_vr_col(i,j)               = value_column

             this%diffus_col(i,j)                        = value_column
             this%ratio_k1_col(i,j)                      = value_column
             this%ratio_no3_co2_col(i,j)                 = value_column
             this%soil_co2_prod_col(i,j)                 = value_column
             this%fr_WFPS_col(i,j)                       = value_column
             this%soil_bulkdensity_col(i,j)              = value_column

             this%r_psi_col(i,j)                         = value_column
             this%anaerobic_frac_col(i,j)                = value_column
          end if
          this%potential_immob_vr_col(i,j)               = value_column
          this%actual_immob_vr_col(i,j)                  = value_column
          this%sminn_to_plant_vr_col(i,j)                = value_column
          this%supplement_to_sminn_vr_col(i,j)           = value_column
          this%gross_nmin_vr_col(i,j)                    = value_column
          this%net_nmin_vr_col(i,j)                      = value_column
          this%sminn_to_plant_fun_no3_vr_col(i,j)        = value_column
          this%sminn_to_plant_fun_nh4_vr_col(i,j)        = value_column
       end do
    end do


    do fi = 1,num_column
       i = filter_column(fi)

       this%ndep_to_sminn_col(i)            = value_column
       this%nfix_to_sminn_col(i)             = value_column
       this%ffix_to_sminn_col(i)             = value_column
       this%fert_to_sminn_col(i)             = value_column
       this%soyfixn_to_sminn_col(i)          = value_column
       this%potential_immob_col(i)           = value_column
       this%actual_immob_col(i)              = value_column
       this%sminn_to_plant_col(i)            = value_column
       this%supplement_to_sminn_col(i)       = value_column
       this%gross_nmin_col(i)                = value_column
       this%net_nmin_col(i)                  = value_column
       this%denit_col(i)                     = value_column
       this%sminn_to_plant_fun_col(i)        = value_column
       if (use_nitrif_denitrif) then
          this%f_nit_col(i)                  = value_column
          this%pot_f_nit_col(i)              = value_column
          this%f_denit_col(i)                = value_column
          this%pot_f_denit_col(i)            = value_column
          this%f_n2o_denit_col(i)            = value_column
          this%f_n2o_nit_col(i)              = value_column
          this%smin_no3_leached_col(i)       = value_column
          this%smin_no3_runoff_col(i)        = value_column
       else
          this%sminn_to_denit_excess_col(i)  = value_column
          this%sminn_leached_col(i)          = value_column
       end if
       this%ninputs_col(i)                   = value_column
       this%noutputs_col(i)                  = value_column
       this%som_n_leached_col(i)             = value_column
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_npools_leached_col(i,k) = value_column
       end do
    end do

    if(use_soil_matrixcn)then
!       call this%matrix_Ninput%SetValueV_scaler(num_column,filter_column(1:num_column),value_column)
    end if


    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_transport_tendency_col(i,j,k) = value_column
          end do
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cascade_ntransfer_col(i,l) = value_column
          this%decomp_cascade_sminn_flux_col(i,l) = value_column
          if (.not. use_nitrif_denitrif) then
             this%sminn_to_denit_decomp_cascade_col(i,l) = value_column
          end if
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_ntransfer_vr_col(i,j,l) = value_column
             this%decomp_cascade_sminn_flux_vr_col(i,j,l) = value_column
             if (.not. use_nitrif_denitrif) then
                this%sminn_to_denit_decomp_cascade_vr_col(i,j,l) = value_column
             end if
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_sourcesink_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc)
    !
    ! !USES:
    use clm_varpar , only: nlevdecomp, ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl , only: use_nitrif_denitrif
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l   ! indices
    integer  :: fc        ! filter indices
    !-----------------------------------------------------------------------

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%denit_col(c) = 0._r8
       this%supplement_to_sminn_col(c) = 0._r8
       this%som_n_leached_col(c)       = 0._r8
    end do

    ! vertically integrate decomposing N cascade fluxes and soil mineral N fluxes associated with decomposition cascade
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             this%decomp_cascade_ntransfer_col(c,k) = &
                  this%decomp_cascade_ntransfer_col(c,k) + &
                  this%decomp_cascade_ntransfer_vr_col(c,j,k) * dzsoi_decomp(j)

             this%decomp_cascade_sminn_flux_col(c,k) = &
                  this%decomp_cascade_sminn_flux_col(c,k) + &
                  this%decomp_cascade_sminn_flux_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do

    if (.not. use_nitrif_denitrif) then

       ! vertically integrate each denitrification flux
       do l = 1, ndecomp_cascade_transitions
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%sminn_to_denit_decomp_cascade_col(c,l) = &
                     this%sminn_to_denit_decomp_cascade_col(c,l) + &
                     this%sminn_to_denit_decomp_cascade_vr_col(c,j,l) * dzsoi_decomp(j)
             end do
          end do
       end do

       ! vertically integrate bulk denitrification and  leaching flux
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%sminn_to_denit_excess_col(c) = &
                  this%sminn_to_denit_excess_col(c) + &
                  this%sminn_to_denit_excess_vr_col(c,j) * dzsoi_decomp(j)

             this%sminn_leached_col(c) = &
                  this%sminn_leached_col(c) + &
                  this%sminn_leached_vr_col(c,j) * dzsoi_decomp(j)
          end do
       end do


       ! total N denitrification (DENIT)
       do l = 1, ndecomp_cascade_transitions
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%denit_col(c) = &
                  this%denit_col(c) + &
                  this%sminn_to_denit_decomp_cascade_col(c,l)
          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%denit_col(c) =  &
               this%denit_col(c) + &
               this%sminn_to_denit_excess_col(c)
       end do

    else


       ! vertically integrate NO3 NH4 N2O fluxes and pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             ! nitrification and denitrification fluxes
             this%f_nit_col(c) = &
                  this%f_nit_col(c) + &
                  this%f_nit_vr_col(c,j) * dzsoi_decomp(j)

             this%f_denit_col(c) = &
                  this%f_denit_col(c) + &
                  this%f_denit_vr_col(c,j) * dzsoi_decomp(j)

             this%pot_f_nit_col(c) = &
                  this%pot_f_nit_col(c) + &
                  this%pot_f_nit_vr_col(c,j) * dzsoi_decomp(j)

             this%pot_f_denit_col(c) = &
                  this%pot_f_denit_col(c) + &
                  this%pot_f_denit_vr_col(c,j) * dzsoi_decomp(j)

             this%f_n2o_nit_col(c) = &
                  this%f_n2o_nit_col(c) + &
                  this%f_n2o_nit_vr_col(c,j) * dzsoi_decomp(j)

             this%f_n2o_denit_col(c) = &
                  this%f_n2o_denit_col(c) + &
                  this%f_n2o_denit_vr_col(c,j) * dzsoi_decomp(j)

             ! leaching/runoff flux
             this%smin_no3_leached_col(c) = &
                  this%smin_no3_leached_col(c) + &
                  this%smin_no3_leached_vr_col(c,j) * dzsoi_decomp(j)

             this%smin_no3_runoff_col(c) = &
                  this%smin_no3_runoff_col(c) + &
                  this%smin_no3_runoff_vr_col(c,j) * dzsoi_decomp(j)

          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%denit_col(c) = this%f_denit_col(c)
       end do

    end if

    ! supplementary N supplement_to_sminn
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%supplement_to_sminn_col(c) = &
               this%supplement_to_sminn_col(c) + &
               this%supplement_to_sminn_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_npools_leached_col(c,l) = 0._r8
       end do

       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_npools_leached_col(c,l) = &
                  this%decomp_npools_leached_col(c,l) + &
                  this%decomp_npools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%som_n_leached_col(c) = &
               this%som_n_leached_col(c) + &
               this%decomp_npools_leached_col(c,l)
       end do
    end do

  end subroutine Summary

end module SoilBiogeochemNitrogenFluxType
