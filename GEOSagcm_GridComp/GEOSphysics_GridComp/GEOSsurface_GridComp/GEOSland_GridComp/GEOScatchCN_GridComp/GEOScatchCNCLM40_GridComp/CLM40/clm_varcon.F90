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
  use shr_const_mod, only: SHR_CONST_G,     &
                           SHR_CONST_RHOFW, &
                           SHR_CONST_TKFRZ, &
                           SHR_CONST_CDAY,  &
                           SHR_CONST_RGAS,  &
                           SHR_CONST_PI
!
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

  real :: rpi    = SHR_CONST_PI

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real :: grav   = SHR_CONST_G      !gravity constant [m/s2]
  real :: denh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real :: rgas   = SHR_CONST_RGAS   !universal gas constant [J/K/kmole]
  real :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]

  real, public, parameter ::  secspday= SHR_CONST_CDAY  ! Seconds per day
  real, public, parameter ::  spval = 1.e36  ! special value for real data
  integer , public, parameter :: ispval = -9999     ! special value for int data

  ! These are tunable constants from clm2_3

  !------------------------------------------------------------------
  ! Initialize water type constants
  !------------------------------------------------------------------

  ! "land unit " types
  !   1     soil (includes vegetated landunits)
  !   2     land ice (glacier)
  !   3     deep lake
  !   4     shallow lake
  !   5     wetland (swamp, marsh, etc.)
  !   6     urban

  integer :: istsoil = 1  !soil         landunit type
  integer :: istice  = 2  !land ice     landunit type
  integer :: istdlak = 3  !deep lake    landunit type
  integer :: istslak = 4  !shallow lake landunit type
  integer :: istwet  = 5  !wetland      landunit type
  integer :: isturb  = 6  !urban        landunit type

end module clm_varcon
