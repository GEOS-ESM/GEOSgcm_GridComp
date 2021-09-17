module rrlw_kg05

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no5 = 16

   real :: fracrefao(no5,9), fracrefbo(no5,5) 
   real :: kao(9,5,13,no5), kbo(5,5,13:59,no5)
   real :: kao_mo3(9,19,no5)
   real :: selfrefo(10,no5), forrefo(4,no5)
   real :: ccl4o(no5)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
 
   integer, parameter :: ng5 = 16

   real :: ka(9,5,13,   ng5), absa( 585,ng5)
   real :: kb(5,5,13:59,ng5), absb(1175,ng5)
   equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))

   real :: fracrefa(ng5,9), fracrefb(ng5,5)
   real :: ka_mo3(9,19,ng5)
   real :: selfref(10,ng5), forref(4,ng5)
   real :: ccl4(ng5)

end module rrlw_kg05
