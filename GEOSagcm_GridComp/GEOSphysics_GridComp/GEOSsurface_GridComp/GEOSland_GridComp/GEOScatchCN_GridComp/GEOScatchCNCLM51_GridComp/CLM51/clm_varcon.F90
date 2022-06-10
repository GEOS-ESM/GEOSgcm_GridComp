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


end module clm_varcon
