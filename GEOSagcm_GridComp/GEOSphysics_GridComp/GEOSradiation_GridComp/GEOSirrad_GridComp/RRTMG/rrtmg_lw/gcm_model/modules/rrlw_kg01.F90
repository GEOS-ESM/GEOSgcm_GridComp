module rrlw_kg01

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no1 = 16

   real :: fracrefao(no1), fracrefbo(no1)
   real :: kao(5,13,no1), kbo(5,13:59,no1)
   real :: kao_mn2(19,no1), kbo_mn2(19,no1)
   real :: selfrefo(10,no1), forrefo(4,no1)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng1 = 10
      
   real :: ka(5,13,   ng1), absa( 65,ng1)
   real :: kb(5,13:59,ng1), absb(235,ng1)
   equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

   real :: fracrefa(ng1), fracrefb(ng1)
   real :: ka_mn2(19,ng1), kb_mn2(19,ng1)
   real :: selfref(10,ng1), forref(4,ng1)
     
end module rrlw_kg01
