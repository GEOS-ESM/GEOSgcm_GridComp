      module rrsw_kg25

      !use parkind ,only : im => kind , rb => kind 
      use parrrsw, only : ng25

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 25
! band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
! Revised: MJIacono, AER, nov2015, solar variability
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real     
!sfluxrefo: real     
!irradnceo: real     
!facbrghto: real     
!snsptdrko: real     
! abso3ao : real     
! abso3bo : real     
! raylo   : real     
!-----------------------------------------------------------------

      integer, parameter :: no25 = 16

      real :: kao(5,13,no25)
      real :: sfluxrefo(no25)
      real :: irradnceo(no25)
      real :: facbrghto(no25),snsptdrko(no25)
      real :: abso3ao(no25), abso3bo(no25)
      real :: raylo(no25)

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 25
! band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
! Revised: MJIacono, AER, nov2015, solar variability
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real     
! absa    : real
! sfluxref: real     
! irradnce: real     
! facbrght: real     
! snsptdrk: real     
! abso3a  : real     
! abso3b  : real     
! rayl    : real     
!-----------------------------------------------------------------

      real :: ka(5,13,ng25), absa(65,ng25)
      real :: sfluxref(ng25)
      real :: irradnce(ng25)
      real :: facbrght(ng25),snsptdrk(ng25)
      real :: abso3a(ng25), abso3b(ng25)
      real :: rayl(ng25)

      equivalence (ka(1,1,1),absa(1,1))

      end module rrsw_kg25

