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
#ifdef __NAG_COMPILER_RELEASE
  use, intrinsic :: iso_fortran_env, only: INT32, REAL32
#endif
  implicit none
  save
  private
  public :: inf, nan, bigint
#ifdef __NAG_COMPILER_RELEASE
! NAG rejects NaN/Inf named constants, so overlay the same bit patterns
! on integer variables instead.
  integer(INT32) :: inf4_bits = int(O'17740000000', INT32)
  integer(INT32) :: nan4_bits = int(O'17760000000', INT32)
  real(REAL32)   :: inf, nan
  equivalence (inf, inf4_bits), (nan, nan4_bits)
#else
! signaling nan
  real*8, parameter :: inf8 = O'0777600000000000000000'
  real*8, parameter :: nan8 = O'0777610000000000000000'
  real*4, parameter :: inf4 = O'17740000000'
  real*4, parameter :: nan4 = O'17760000000'
  real,   parameter :: inf = inf4
  real,   parameter :: nan = nan4
#endif
  integer,  parameter :: bigint = O'17777777777'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein based on cam module created by
! CCM core group
!
!EOP
!-----------------------------------------------------------------------

end module nanMod
