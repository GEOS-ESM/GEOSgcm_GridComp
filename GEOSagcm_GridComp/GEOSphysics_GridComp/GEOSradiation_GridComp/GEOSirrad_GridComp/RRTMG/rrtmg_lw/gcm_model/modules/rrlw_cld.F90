      module rrlw_cld

       !use parkind, only : rb => kind 

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw cloud property coefficients

! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! abscld1:  real   : 
! absice0:  real   : 
! absice1:  real   : 
! absice2:  real   : 
! absice3:  real   : 
! absliq0:  real   : 
! absliq1:  real   : 
!------------------------------------------------------------------

      real  :: abscld1
      real  , dimension(2) :: absice0
      real  , dimension(2,5) :: absice1
      real  , dimension(43,16) :: absice2
      real  , dimension(46,16) :: absice3
      real  , dimension(200,16) :: absice4
      real :: absliq0
      real  , dimension(58,16) :: absliq1

      end module rrlw_cld

