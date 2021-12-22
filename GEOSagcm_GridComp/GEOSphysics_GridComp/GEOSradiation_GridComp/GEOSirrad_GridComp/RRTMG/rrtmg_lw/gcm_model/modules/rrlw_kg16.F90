module rrlw_kg16

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no16 = 16

   real :: fracrefao(no16,9), fracrefbo(no16)
   real :: kao(9,5,13,no16), kbo(5,13:59,no16)
   real :: selfrefo(10,no16), forrefo(4,no16)
   
!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 16
! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng16 = 2

   real :: ka(9,5,13, ng16), absa(585,ng16)
   real :: kb(5,13:59,ng16), absb(235,ng16)
   equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng16,9), fracrefb(ng16)
   real :: selfref(10,ng16), forref(4,ng16)

end module rrlw_kg16
