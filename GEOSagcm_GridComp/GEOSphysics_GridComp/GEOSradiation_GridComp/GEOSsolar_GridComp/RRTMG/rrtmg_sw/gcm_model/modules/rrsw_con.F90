module rrsw_con

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw constants

! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! grav   :  real   : acceleration of gravity
! planck :  real   : planck constant
! boltz  :  real   : boltzmann constant
! clight :  real   : speed of light
! avogad :  real   : avogadro constant 
! alosmt :  real   : loschmidt constant
! gascon :  real   : molar gas constant
! radcn1 :  real   : first radiation constant
! radcn2 :  real   : second radiation constant
! sbcnst :  real   : stefan-boltzmann constant
!------------------------------------------------------------------

      real, parameter :: pi = 3.14159265358979323846d0
      real, parameter :: oneminus = 1. - 1.e-06

      real :: grav
      real :: planck, boltz, clight
      real :: avogad, alosmt, gascon
      real :: radcn1, radcn2
      real :: sbcnst

end module rrsw_con

