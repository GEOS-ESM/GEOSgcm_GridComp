module rrlw_kg15

   implicit none
   save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: no15 = 16

   real :: fracrefao(no15,9)
   real :: kao(9,5,13,no15)
   real :: kao_mn2(9,19,no15)
   real :: selfrefo(10,no15), forrefo(4,no15)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------

   integer, parameter :: ng15 = 2

   real :: ka(9,5,13,ng15), absa(585,ng15)
   equivalence (ka(1,1,1,1),absa(1,1))

   real :: fracrefa(ng15,9)
   real :: ka_mn2(9,19,ng15)
   real :: selfref(10,ng15), forref(4,ng15)

end module rrlw_kg15
