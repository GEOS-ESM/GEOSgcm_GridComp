module clm_varcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varcon
!
! !DESCRIPTION:
! Module containing various model constants
!
! !USES:
  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use shr_const_mod, only: SHR_CONST_G,     &
                           SHR_CONST_RHOFW, &
                           SHR_CONST_TKFRZ, &
                           SHR_CONST_CDAY,  &
                           SHR_CONST_RGAS,  &
                           SHR_CONST_PI,    &
                           SHR_CONST_PDB
  use clm_varpar   , only: nlevgrnd, nlevdecomp_full, numrad, ngases

! !PUBLIC TYPES:
  implicit none
  save
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 27 February 2008: Keith Oleson; Add forcing height and aerodynamic parameters
!
!EOP
!-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(r8) :: rpi    = SHR_CONST_PI

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real(r8) :: grav   = SHR_CONST_G      !gravity constant [m/s2]
  real(r8) :: denh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real(r8) :: rgas   = SHR_CONST_RGAS   !universal gas constant [J/K/kmole]
  real(r8) :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8), public, parameter ::  secspday= SHR_CONST_CDAY  ! Seconds per day
  real(r8), public, parameter ::  secsphr = 3600._r8        ! Seconds in an hour
  real(r8), public, parameter ::  spval = 1.e36_r8          ! special value for real data
  integer , public, parameter ::  ispval = -9999            ! special value for int data
  integer, public, parameter  ::  fun_period  = 1            ! A FUN parameter, and probably needs to be changed for testing
  real(r8),public, parameter  ::  smallValue  = 1.e-12_r8    ! A small values used by FUN



  !------------------------------------------------------------------
  ! Soil depths
  !------------------------------------------------------------------

  real(r8), pointer :: zsoi(:)         !soil z  (layers)
  real(r8), pointer :: dzsoi(:)        !soil dz (thickness)
  real(r8), pointer :: zisoi(:)        !soil zi (interfaces)
  real(r8), pointer :: dzsoi_decomp(:) !soil dz (thickness)
  real(r8), public, parameter :: zmin_bedrock = 0.4_r8 ! minimum soil depth [m]

  !------------------------------------------------------------------
  ! Set subgrid names
  !------------------------------------------------------------------

  character(len=16), public, parameter :: grlnd  = 'lndgrid'      ! name of lndgrid
  character(len=16), public, parameter :: namea  = 'gridcellatm'  ! name of atmgrid
  character(len=16), public, parameter :: nameg  = 'gridcell'     ! name of gridcells
  character(len=16), public, parameter :: namel  = 'landunit'     ! name of landunits
  character(len=16), public, parameter :: namec  = 'column'       ! name of columns
  character(len=16), public, parameter :: namep  = 'pft'          ! name of patches
  character(len=16), public, parameter :: nameCohort = 'cohort'   ! name of cohorts (ED specific)

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------
  
  real(r8), public :: betads  = 0.5_r8            ! two-stream parameter betad for snow
  real(r8), public :: betais  = 0.5_r8            ! two-stream parameter betai for snow
  real(r8), public :: omegas(numrad) = (/0.8_r8, 0.4_r8/)  ! two-stream parameter omega for snow by band

  integer, parameter, public :: max_lunit  = 9  !maximum value that lun%itype can have

   ! typical del13C for C3 photosynthesis (permil, relative to PDB)
  real(r8), public, parameter :: c3_del13c = -28._r8

  ! typical del13C for C4 photosynthesis (permil, relative to PDB)
  real(r8), public, parameter :: c4_del13c = -13._r8

  ! isotope ratio (13c/12c) for C3 photosynthesis
  real(r8), public, parameter :: c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)

  ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
  real(r8), public, parameter :: c3_r2 = c3_r1/(1._r8 + c3_r1)

  real(r8), public :: c13ratio = 1. !jkolassa Jan 2023: dummy value since this is only needed to compile the code, but not used
  real(r8), public :: c14ratio = 1. !jkolassa Jan 2023: dummy value since this is only needed to compile the code, but not used


  real(r8), public, parameter :: nitrif_n2o_loss_frac = 6.e-4_r8 ! fraction of N lost as N2O in nitrification (Li et al., 2000)

  real(r8), public, parameter :: degpsec = 15._r8/3600.0_r8 ! Degree's earth rotates per second
  integer,  public, parameter :: isecspday= secspday        ! Integer seconds per day

  real(r8), public, parameter :: c_to_b = 2.0_r8         ! conversion between mass carbon and total biomass (g biomass /g C)

  !------------------------------------------------------------------
  ! (Non-tunable) Constants for the CH4 submodel (Tuneable constants in ch4varcon)
  !------------------------------------------------------------------
  ! Note some of these constants are also used in CNNitrifDenitrifMod

  integer, private :: i  ! loop index

  real(r8), public :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)
  data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
  data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

  real(r8), public :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  data (d_con_g(1,i),i=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759_r8, 0.00117_r8/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325_r8, 0.0009_r8/ ! CO2


! !PUBLIC MEMBER FUNCTIONS:
  public clm_varcon_init          ! Initialze constants that need to be initialized

! !REVISION HISTORY:
! Created by Mariana Vertenstein

!EOP
!-----------------------------------------------------------------------
contains
!-------------------------------
  subroutine clm_varcon_init()
!
! !DESCRIPTION:
! This subroutine initializes constants in clm_varcon. MUST be called 
! after the clm_varpar_init.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
!
!EOP
!------------------------------------------------------------------------------
  allocate( zsoi(1:nlevgrnd) )
  allocate( dzsoi(1:nlevgrnd) )
  allocate( zisoi(0:nlevgrnd) )
  allocate( dzsoi_decomp(1:nlevdecomp_full) )

  ! jkolassa Aug 2022: This follows previous implementations of Catchment-CN and works as long as we use a single soil layer (for CN); we will have to update this if we increase the number of soil layers.
  zsoi(1)  = 0.5
  dzsoi(1) = 1.
  zisoi(0) = 0.
  zisoi(1) = 1.
  dzsoi_decomp(1) = dzsoi(1)

  end subroutine clm_varcon_init
end module clm_varcon
