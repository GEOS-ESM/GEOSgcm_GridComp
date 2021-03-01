      module rrsw_kg16

      !use parkind ,only : im => kind , rb => kind 
      use parrrsw, only : ng16

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
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

      integer, parameter :: no16 = 16

      real :: kao(9,5,13,no16)
      real :: kbo(5,13:59,no16)
      real :: selfrefo(10,no16), forrefo(3,no16)
      real :: sfluxrefo(no16)
      real :: irradnceo(no16)
      real :: facbrghto(no16),snsptdrko(no16)

      real :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 16
! band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)
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

      real :: ka(9,5,13,ng16) , absa(585,ng16)
      real :: kb(5,13:59,ng16), absb(235,ng16)
      real :: selfref(10,ng16), forref(3,ng16)
      real :: sfluxref(ng16)
      real :: irradnce(ng16)
      real :: facbrght(ng16),snsptdrk(ng16)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg16

