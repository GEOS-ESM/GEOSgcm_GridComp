module rrlw_kg12

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 12
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no12 = 16

   real :: fracrefao(no12,9)
   real :: kao(9,5,13,no12)
   real :: selfrefo(10,no12), forrefo(4,no12)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 12
! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng12 = 8

   real :: ka(9,5,13,ng12), absa(585,ng12)
   equivalence (ka(1,1,1,1),absa(1,1))

   real :: fracrefa(ng12,9)
   real :: selfref(10,ng12), forref(4,ng12)

end module rrlw_kg12
