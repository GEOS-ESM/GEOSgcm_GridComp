module rrlw_kg07

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no7 = 16

   real :: fracrefao(no7,9), fracrefbo(no7)
   real :: kao(9,5,13,no7), kbo(5,13:59,no7)
   real :: kao_mco2(9,19,no7), kbo_mco2(19,no7)
   real :: selfrefo(10,no7), forrefo(4,no7)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng7 = 12

   real :: ka(9,5,13, ng7), absa(585,ng7)
   real :: kb(5,13:59,ng7), absb(235,ng7)
   equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng7,9), fracrefb(ng7)
   real :: ka_mco2(9,19,ng7), kb_mco2(19,ng7)
   real :: selfref(10,ng7), forref(4,ng7)

end module rrlw_kg07
