      module rrsw_kg19

      !use parkind ,only : im => kind , rb => kind_rb
      use parrrsw, only : ng19

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 19
! band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
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
! kbo     : real     
! selfrefo: real     
! forrefo : real
!sfluxrefo: real     
!irradnceo: real     
!facbrghto: real     
!snsptdrko: real     
!-----------------------------------------------------------------

      integer, parameter :: no19 = 16

      real :: kao(9,5,13,no19)
      real :: kbo(5,13:59,no19)
      real :: selfrefo(10,no19), forrefo(3,no19)
      real :: sfluxrefo(no19,9)
      real :: irradnceo(no19,9)
      real :: facbrghto(no19,9),snsptdrko(no19,9)

      real :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 19
! band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
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
! kb      : real     
! absa    : real
! absb    : real
! selfref : real     
! forref  : real
! sfluxref: real     
! irradnce: real     
! facbrght: real     
! snsptdrk: real     
!-----------------------------------------------------------------

      real :: ka(9,5,13,ng19), absa(585,ng19)
      real :: kb(5,13:59,ng19), absb(235,ng19)
      real :: selfref(10,ng19), forref(3,ng19)
      real :: sfluxref(ng19,9)
      real :: irradnce(ng19,9)
      real :: facbrght(ng19,9),snsptdrk(ng19,9)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg19

