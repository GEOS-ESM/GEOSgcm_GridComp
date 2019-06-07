      module rrsw_kg26

      !use parkind ,only : im => kind , rb => kind 
      use parrrsw, only : ng26

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 26
! band 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
! Revised: MJIacono, AER, nov2015, solar variability
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!sfluxrefo: real     
!irradnceo: real     
!facbrghto: real     
!snsptdrko: real     
! raylo   : real     
!-----------------------------------------------------------------

      integer, parameter :: no26 = 16

      real :: sfluxrefo(no26)
      real :: irradnceo(no26)
      real :: facbrghto(no26),snsptdrko(no26)
      real :: raylo(no26)

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 26
! band 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
! Revised: MJIacono, AER, nov2015, solar variability
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! sfluxref: real     
! irradnce: real     
! facbrght: real     
! snsptdrk: real     
! rayl    : real     
!-----------------------------------------------------------------

      real :: sfluxref(ng26)
      real :: irradnce(ng26)
      real :: facbrght(ng26),snsptdrk(ng26)
      real :: rayl(ng26)

      end module rrsw_kg26

