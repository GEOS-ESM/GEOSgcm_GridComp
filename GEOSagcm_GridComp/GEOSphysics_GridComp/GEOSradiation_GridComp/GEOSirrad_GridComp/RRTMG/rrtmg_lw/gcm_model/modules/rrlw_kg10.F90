module rrlw_kg10

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no10 = 16

   real :: fracrefao(no10), fracrefbo(no10)
   real :: kao(5,13,no10), kbo(5,13:59,no10)
   real :: selfrefo(10,no10), forrefo(4,no10)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng10 = 6

   real :: ka(5,13,   ng10), absa( 65,ng10)
   real :: kb(5,13:59,ng10), absb(235,ng10)
   equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

   real :: fracrefa(ng10), fracrefb(ng10)
   real :: selfref(10,ng10), forref(4,ng10)

end module rrlw_kg10
