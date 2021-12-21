module rrlw_kg03

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no3 = 16

   real :: fracrefao(no3,9), fracrefbo(no3,5)
   real :: kao(9,5,13,no3), kbo(5,5,13:59,no3)
   real :: kao_mn2o(9,19,no3), kbo_mn2o(5,19,no3)
   real :: selfrefo(10,no3), forrefo(4,no3)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 3
! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng3 = 16

   real :: ka(9,5,13,   ng3), absa( 585,ng3)
   real :: kb(5,5,13:59,ng3), absb(1175,ng3)
   equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

   real :: fracrefa(ng3,9), fracrefb(ng3,5)
   real :: ka_mn2o(9,19,ng3), kb_mn2o(5,19,ng3)
   real :: selfref(10,ng3), forref(4,ng3)

end module rrlw_kg03
