module rrlw_kg04

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no4 = 16

   real :: fracrefao(no4,9), fracrefbo(no4,5)
   real :: kao(9,5,13,no4), kbo(5,5,13:59,no4)
   real :: selfrefo(10,no4), forrefo(4,no4)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng4 = 14

   real :: ka(9,5,13,   ng4), absa( 585,ng4)
   real :: kb(5,5,13:59,ng4), absb(1175,ng4)
   equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

   real :: fracrefa(ng4,9), fracrefb(ng4,5)
   real :: selfref(10,ng4), forref(4,ng4)

end module rrlw_kg04
