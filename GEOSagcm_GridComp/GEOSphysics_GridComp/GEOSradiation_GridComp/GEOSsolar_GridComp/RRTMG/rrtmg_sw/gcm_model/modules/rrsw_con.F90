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
!oneminus:  real   : 1.-1.e-6
! pi     :  real   : pi
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

      real  :: oneminus, pi, grav
      real  :: planck, boltz, clight
      real  :: avogad, alosmt, gascon
      real  :: radcn1, radcn2
      real  :: sbcnst

end module rrsw_con

