      module rrsw_kg17

      !use parkind ,only : im => kind , rb => kind 
      use parrrsw, only : ng17

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 17
! band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
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

      integer, parameter :: no17 = 16

      real :: kao(9,5,13,no17)
      real :: kbo(5,5,13:59,no17)
      real :: selfrefo(10,no17), forrefo(4,no17)
      real :: sfluxrefo(no17,5)
      real :: irradnceo(no17,5)
      real :: facbrghto(no17,5),snsptdrko(no17,5)

      real :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 17
! band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)
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

      real :: ka(9,5,13,ng17) , absa(585,ng17)
      real :: kb(5,5,13:59,ng17), absb(1175,ng17)
      real :: selfref(10,ng17), forref(4,ng17)
      real :: sfluxref(ng17,5)
      real :: irradnce(ng17,5)
      real :: facbrght(ng17,5),snsptdrk(ng17,5)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

      end module rrsw_kg17

