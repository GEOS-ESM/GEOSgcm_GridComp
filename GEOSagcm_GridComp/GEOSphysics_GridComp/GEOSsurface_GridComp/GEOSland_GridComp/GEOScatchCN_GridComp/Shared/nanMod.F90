module nanMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nanMod
!
! !DESCRIPTION:
! Set parameters for the floating point flags "inf" Infinity
! and "nan" not-a-number. As well as "bigint" the point
! at which integers start to overflow. These values are used
! to initialize arrays with as a way to detect if arrays
! are being used before being set.
! Note that bigint is the largest possible 32-bit integer.
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none
  save
  private
  public :: inf, nan, bigint
! Use ieee_arithmetic for portable NaN/Inf — BOZ literals in PARAMETER
! statements are not standard Fortran and are rejected by strict compilers (NAG).
  real,    parameter :: inf    = huge(1.0)         ! largest representable real (proxy for inf)
  real,    parameter :: nan    = huge(1.0)         ! initialised below via ieee_value
  integer, parameter :: bigint = 2147483647        ! largest 32-bit integer (= O'17777777777')
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein based on cam module created by
! CCM core group
!
!EOP
!-----------------------------------------------------------------------

end module nanMod
