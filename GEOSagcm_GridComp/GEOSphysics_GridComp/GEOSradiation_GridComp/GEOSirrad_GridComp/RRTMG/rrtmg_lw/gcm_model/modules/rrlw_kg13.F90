module rrlw_kg13

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no13 = 16

   real :: fracrefao(no13,9), fracrefbo(no13)
   real :: kao(9,5,13,no13)
   real :: kao_mco2(9,19,no13)
   real :: kbo_mo3(19,no13)
   real :: kao_mco(9,19,no13)
   real :: selfrefo(10,no13), forrefo(4,no13)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng13 = 4

   real :: ka(9,5,13,ng13), absa(585,ng13)
   equivalence (ka(1,1,1,1),absa(1,1))

   real :: fracrefa(ng13,9), fracrefb(ng13)
   real :: ka_mco2(9,19,ng13)
   real :: kb_mo3(19,ng13)
   real :: ka_mco(9,19,ng13)
   real :: selfref(10,ng13), forref(4,ng13)

end module rrlw_kg13
