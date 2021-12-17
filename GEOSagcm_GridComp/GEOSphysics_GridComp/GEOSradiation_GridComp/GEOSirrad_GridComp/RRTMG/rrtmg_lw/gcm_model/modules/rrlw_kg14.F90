module rrlw_kg14

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 14
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no14 = 16

   real :: fracrefao(no14), fracrefbo(no14)
   real :: kao(5,13,no14), kbo(5,13:59,no14)
   real :: selfrefo(10,no14), forrefo(4,no14)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 14
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng14 = 2

   real :: ka(5,13,   ng14), absa( 65,ng14)
   real :: kb(5,13:59,ng14), absb(235,ng14)
   equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng14), fracrefb(ng14)
   real :: selfref(10,ng14), forref(4,ng14)

end module rrlw_kg14
