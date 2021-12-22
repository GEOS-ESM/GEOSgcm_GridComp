module rrlw_kg06

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 6
! band 6:  820-980 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no6 = 16

   real :: fracrefao(no6)
   real :: kao(5,13,no6)
   real :: kao_mco2(19,no6)
   real :: selfrefo(10,no6), forrefo(4,no6)
   real :: cfc11adjo(no6), cfc12o(no6)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 6
! band 6:  820-980 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng6 = 8

   real :: ka(5,13,ng6), absa(65,ng6)
   equivalence (ka(1,1,1),absa(1,1))

   real :: fracrefa(ng6)
   real :: ka_mco2(19,ng6)
   real :: selfref(10,ng6), forref(4,ng6)
   real :: cfc11adj(ng6), cfc12(ng6)

end module rrlw_kg06
