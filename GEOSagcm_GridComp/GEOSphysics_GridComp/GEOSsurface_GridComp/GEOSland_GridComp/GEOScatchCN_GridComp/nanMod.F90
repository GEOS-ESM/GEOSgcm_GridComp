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

  use iso_fortran_env
  use ieee_arithmetic

! !PUBLIC TYPES:
  implicit none
  save
  private
  public :: inf, nan, inf64, nan64, bigint
! signaling nan
  integer,  parameter :: bigint = O'17777777777'

  contains

   function nan() result(nan_32)
      real(REAL32) :: nan_32

      nan_32 = ieee_value(nan_32,  ieee_quiet_nan)

   end function nan

   function nan64() result(nan_64)
      real(REAL64) :: nan_64

      nan_64 = ieee_value(nan_64,  ieee_quiet_nan)

   end function nan64

   function inf() result(inf_32)
      real(REAL32) :: inf_32

      inf_32 = ieee_value(inf_32,  ieee_positive_inf)

   end function inf

   function inf64() result(inf_64)
      real(REAL64) :: inf_64

      inf_64 = ieee_value(inf_64,  ieee_positive_inf)

   end function inf64


!
! !REVISION HISTORY:
! Created by Mariana Vertenstein based on cam module created by
! CCM core group
!
!EOP
!-----------------------------------------------------------------------

end module nanMod
