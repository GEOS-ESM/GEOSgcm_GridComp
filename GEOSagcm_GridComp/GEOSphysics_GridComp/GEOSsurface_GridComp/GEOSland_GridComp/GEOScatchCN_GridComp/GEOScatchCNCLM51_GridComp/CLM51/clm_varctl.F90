module clm_varctl

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varctl
!
! !DESCRIPTION:
! Module containing run control variables
!
! !USES:
  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none

  public :: init_clm_varctl          ! set parameters
  public :: cnallocate_carbon_only

  logical, public :: use_nguardrail         = .true.  ! true => use precision control

  logical, public :: use_luna = .false.            ! true => use  LUNA
  logical, public :: use_fates = .false.           ! true => use fates
  logical, public :: use_hydrstress = .true.       ! true => use plant hydraulic stress calculation


  ! If prognostic crops are turned on
  logical, public :: use_crop = .false.

  logical, public :: use_lch4            = .false.
  logical, public :: use_nitrif_denitrif = .false.
  logical, public :: use_vertsoilc       = .false.
  logical, public :: use_century_decomp  = .false.
  logical, public :: use_cn              = .true.
  logical, public :: use_cndv            = .false.
  logical, public :: use_grainproduct    = .false.
  logical, public :: use_dynroot         = .false.
  logical, public :: use_bedrock = .false. ! true => use spatially variable soil depth
  logical, public :: use_extralakelayers = .false.
  logical, public :: use_biomass_heat_storage = .false.

  logical, public :: use_c13 = .false.                  ! true => use C-13 model
  logical, public :: use_c14 = .false.                  ! true => use C-14 model
  
  ! use subgrid fluxes
  logical,  public :: use_subgrid_fluxes = .true.

  !----------------------------------------------------------
  ! SSRE diagnostic
  !----------------------------------------------------------
  logical, public :: use_SSRE = .false.   ! flag for SSRE diagnostic

  !----------------------------------------------------------
  ! CN matrix
  !----------------------------------------------------------  
  logical, public :: use_matrixcn = .false. !.false.              ! true => use cn matrix
  logical, public :: use_soil_matrixcn = .false.! true => use cn matrix

  real(r8), public :: nfix_timeconst = -1.2345_r8

  !----------------------------------------------------------
  ! Unit Numbers
  !----------------------------------------------------------
  !
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6; jkolassa: This is following CTSM, iulog is not set to output_unit

  !----------------------------------------------------------
  !  flexibleCN
  !----------------------------------------------------------
  logical, public :: use_flexibleCN = .false.
  logical, public :: CNratio_floating = .false.
  integer, public :: CN_evergreen_phenology_opt = 0
  logical, public :: lnc_opt = .false.
  logical, public :: reduce_dayl_factor = .false.
  integer, public :: vcmax_opt = 0

  !----------------------------------------------------------
  ! BGC logic and datasets
  !----------------------------------------------------------

  ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
  ! default value reset in controlMod
  logical, public :: anoxia  = .true.

  ! State of the model for the accelerated decomposition (AD) spinup. 
  ! 0 (default) = normal model; 1 = AD SPINUP
  integer, public :: spinup_state = 0

  logical, public :: use_snicar_frc      = .false.

  integer, public :: carbon_resp_opt = 0

  ! Set in CNAllocationInit (TODO - had to move it here to avoid circular dependency)
  logical, private:: carbon_only
contains

!---------------------------------------
 subroutine init_clm_varctl()

 !---
  if (nfix_timeconst == -1.2345_r8) then
     if (use_nitrif_denitrif) then
        nfix_timeconst = 10._r8
     else
        nfix_timeconst = 0._r8
     end if
  end if

 end subroutine init_clm_varctl

  ! Get module carbon_only flag
  logical function CNAllocate_Carbon_only()
    cnallocate_carbon_only = carbon_only
  end function CNAllocate_Carbon_only

end module clm_varctl
