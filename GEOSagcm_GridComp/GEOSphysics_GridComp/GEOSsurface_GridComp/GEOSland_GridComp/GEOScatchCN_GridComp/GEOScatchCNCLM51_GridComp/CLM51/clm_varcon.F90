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
  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use shr_const_mod, only: SHR_CONST_G,     &
                           SHR_CONST_RHOFW, &
                           SHR_CONST_TKFRZ, &
                           SHR_CONST_CDAY,  &
                           SHR_CONST_RGAS,  &
                           SHR_CONST_PI,    &
                           SHR_CONST_PDB
  use clm_varpar   , only: nlevgrnd, nlevdecomp_full

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
  real(r8), public, parameter ::  spval = 1.e36_r8  ! special value for real data
  integer , public, parameter :: ispval = -9999     ! special value for int data

  !------------------------------------------------------------------
  ! Soil depths
  !------------------------------------------------------------------

  real(r8), pointer :: zsoi(:)         !soil z  (layers)
  real(r8), pointer :: dzsoi(:)        !soil dz (thickness)
  real(r8), pointer :: zisoi(:)        !soil zi (interfaces)
  real(r8), pointer :: dzsoi_decomp(:) !soil dz (thickness)


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

  integer, parameter, public :: max_lunit  = 9  !maximum value that lun%itype can have

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
