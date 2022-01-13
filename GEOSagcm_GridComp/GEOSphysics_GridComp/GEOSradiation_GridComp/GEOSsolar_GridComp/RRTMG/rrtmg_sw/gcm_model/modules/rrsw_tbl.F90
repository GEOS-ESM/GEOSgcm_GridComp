      module rrsw_tbl

      !use parkind, only : im => kind , rb => kind 

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw lookup table arrays

! Initial version: MJIacono, AER, may2007
! Revised: MJIacono, AER, aug2007
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ntbl   :  integer: Lookup table dimension
! tblint :  real   : Lookup table conversion factor
! tau_tbl:  real   : Clear-sky optical depth 
! exp_tbl:  real   : Exponential lookup table for transmittance
! od_lo  :  real   : Value of tau below which expansion is used
!                  : in place of lookup table
! pade   :  real   : Pade approximation constant
! bpade  :  real   : Inverse of Pade constant
!------------------------------------------------------------------

      integer, parameter :: ntbl = 10000

      real, parameter :: tblint = 10000.0 

      real, parameter :: od_lo = 0.06 

      real :: tau_tbl
      real, dimension(0:ntbl) :: exp_tbl

      real, parameter :: pade = 0.278 
      real :: bpade

      end module rrsw_tbl

