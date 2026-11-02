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
! Bit patterns are set through EQUIVALENCE rather than as PARAMETERs because
! some compilers (e.g. NAG) reject NaN/Inf values in constant expressions.
  integer :: inf_bits = int(O'17740000000')
  integer :: nan_bits = int(O'17760000000')
  real    :: inf, nan
  equivalence (inf_bits, inf), (nan_bits, nan)
  integer,  parameter :: bigint = O'17777777777'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein based on cam module created by
! CCM core group
!
!EOP
!-----------------------------------------------------------------------

end module nanMod
