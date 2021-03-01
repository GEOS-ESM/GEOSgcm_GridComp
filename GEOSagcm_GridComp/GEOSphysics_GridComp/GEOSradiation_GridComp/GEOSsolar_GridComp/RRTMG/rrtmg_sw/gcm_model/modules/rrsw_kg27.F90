      module rrsw_kg27

      !use parkind ,only : im => kind , rb => kind 
      use parrrsw, only : ng27

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 27
! band 27: 29000-38000 cm-1 (low - o3; high - o3)
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
!sfluxrefo: real     
!irradnceo: real     
!facbrghto: real     
!snsptdrko: real     
! raylo   : real     
!-----------------------------------------------------------------

      integer, parameter :: no27 = 16

      real :: kao(5,13,no27)
      real :: kbo(5,13:59,no27)
      real :: sfluxrefo(no27)
      real :: irradnceo(no27)
      real :: facbrghto(no27),snsptdrko(no27)
      real :: raylo(no27)

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 27
! band 27: 29000-38000 cm-1 (low - o3; high - o3)
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
! sfluxref: real     
! irradnce: real     
! facbrght: real     
! snsptdrk: real     
! rayl    : real     
!-----------------------------------------------------------------

      real :: ka(5,13,ng27), absa(65,ng27)
      real :: kb(5,13:59,ng27), absb(235,ng27)
      real :: sfluxref(ng27)
      real :: irradnce(ng27)
      real :: facbrght(ng27),snsptdrk(ng27)
      real :: rayl(ng27)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg27

