      module rrsw_kg18

      !use parkind ,only : im => kind , rb => kind 
      use parrrsw, only : ng18

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 18
! band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
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

      integer, parameter :: no18 = 16

      real :: kao(9,5,13,no18)
      real :: kbo(5,13:59,no18)
      real :: selfrefo(10,no18), forrefo(3,no18)
      real :: sfluxrefo(no18,9)
      real :: irradnceo(no18,9)
      real :: facbrghto(no18,9),snsptdrko(no18,9)

      real :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 18
! band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)
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

      real :: ka(9,5,13,ng18), absa(585,ng18)
      real :: kb(5,13:59,ng18), absb(235,ng18)
      real :: selfref(10,ng18), forref(3,ng18)
      real :: sfluxref(ng18,9)
      real :: irradnce(ng18,9)
      real :: facbrght(ng18,9),snsptdrk(ng18,9)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg18

