module rrlw_kg02

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no2 = 16

   real :: fracrefao(no2), fracrefbo(no2)
   real :: kao(5,13,no2), kbo(5,13:59,no2)
   real :: selfrefo(10,no2), forrefo(4,no2)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng2 = 12

   real :: ka(5,13,   ng2), absa( 65,ng2)
   real :: kb(5,13:59,ng2), absb(235,ng2)
   equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng2), fracrefb(ng2)
   real :: selfref(10,ng2), forref(4,ng2)

end module rrlw_kg02
