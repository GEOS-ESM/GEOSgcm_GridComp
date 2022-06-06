module rrlw_kg11

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no11 = 16

   real :: fracrefao(no11), fracrefbo(no11)
   real :: kao(5,13,no11), kbo(5,13:59,no11)
   real :: kao_mo2(19,no11), kbo_mo2(19,no11)
   real :: selfrefo(10,no11), forrefo(4,no11)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 11
! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng11 = 8

   real :: ka(5,13,   ng11), absa( 65,ng11)
   real :: kb(5,13:59,ng11), absb(235,ng11)
   equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng11), fracrefb(ng11)
   real :: ka_mo2(19,ng11), kb_mo2(19,ng11)
   real :: selfref(10,ng11), forref(4,ng11)

end module rrlw_kg11
