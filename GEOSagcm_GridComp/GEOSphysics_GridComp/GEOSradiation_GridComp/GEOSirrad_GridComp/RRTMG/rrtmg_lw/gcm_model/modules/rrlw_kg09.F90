module rrlw_kg09

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no9 = 16

   real :: fracrefao(no9,9), fracrefbo(no9)
   real :: kao(9,5,13,no9), kbo(5,13:59,no9)
   real :: kao_mn2o(9,19,no9), kbo_mn2o(19,no9)
   real :: selfrefo(10,no9), forrefo(4,no9)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng9 = 12

   real :: ka(9,5,13, ng9), absa(585,ng9)
   real :: kb(5,13:59,ng9), absb(235,ng9)
   equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng9,9), fracrefb(ng9)
   real :: ka_mn2o(9,19,ng9), kb_mn2o(19,ng9)
   real :: selfref(10,ng9), forref(4,ng9)

end module rrlw_kg09
