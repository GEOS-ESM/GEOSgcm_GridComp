module rrlw_kg08

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no8 = 16

   real :: fracrefao(no8), fracrefbo(no8)
   real :: kao(5,13,no8), kbo(5,13:59,no8)
   real :: kao_mco2(19,no8), kbo_mco2(19,no8)
   real :: kao_mo3(19,no8)
   real :: kao_mn2o(19,no8), kbo_mn2o(19,no8)
   real :: selfrefo(10,no8), forrefo(4,no8)
   real :: cfc12o(no8), cfc22adjo(no8)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng8 = 8

   real :: ka(5,13,   ng8), absa( 65,ng8)
   real :: kb(5,13:59,ng8), absb(235,ng8)
   equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

   real :: fracrefa(ng8), fracrefb(ng8)
   real :: ka_mco2(19,ng8), kb_mco2(19,ng8)
   real :: ka_mo3(19,ng8)
   real :: ka_mn2o(19,ng8), kb_mn2o(19,ng8)
   real :: selfref(10,ng8), forref(4,ng8)
   real :: cfc12(ng8), cfc22adj(ng8)

end module rrlw_kg08
