module rrlw_con

   implicit none
   save

!------------------------------------------------------------------
! rrtmg_lw constants

! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! fluxfac:  real   : radiance to flux conversion factor 
! heatfac:  real   : flux to heating rate conversion factor
!oneminus:  real   : 1.-1.e-6
! pi     :  real   : pi
! d2r    :  real   : multiply by to convert degrees to radians
! r2d    :  real   : multiply by to convert radians to degrees
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
!  secdy :  real   : seconds per day  
!------------------------------------------------------------------

   real, parameter :: pi = 3.14159265358979323846d0
   real, parameter :: d2r = pi / 180.
   real, parameter :: r2d = 180. / pi
   real, parameter :: oneminus = 1. - 1.e-6
   real, parameter :: fluxfac = pi * 2.e4 

   real  :: heatfac
   real  :: grav
   real  :: planck, boltz, clight
   real  :: avogad, alosmt, gascon
   real  :: radcn1, radcn2
   real  :: sbcnst, secdy

end module rrlw_con

