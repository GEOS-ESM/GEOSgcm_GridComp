! Module: mo_rng_mklvsl

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

!
! Implementation of random number interface using Intel Vector Statistics Library
!

include 'mkl_vsl.f90' ! Needed to generate MKL_VSL modules below

module mo_rng_mklvsl

  use mo_rte_kind, only : dp
  use mo_rng
  USE MKL_VSL_TYPE
  USE MKL_VSL
  implicit none

  ! Mersenne Twister
  ! Alternatives are VSL_BRNG_SFMT19937, maybe VSL_BRNG_MT2203?
  integer, parameter :: rng_type = VSL_BRNG_MT19937

  ! Uniform distribution
  integer, parameter :: rng_method = VSL_RNG_METHOD_UNIFORM_STD

  ! extended random number generator class
  type, extends(ty_rng), public :: ty_rng_mklvsl
    TYPE (VSL_STREAM_STATE) :: stream
  contains
    procedure :: get_random_vec      => get_rng_mkl_vec
    procedure :: get_random_vec_mask => get_rng_mkl_mask
    procedure :: init_rng            => init_rng_mkl
    procedure :: end_rng             => end_rng_mkl
  end type

contains

  subroutine stop_on_err(msg)

    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if (msg /= "") then
      write(error_unit, *) msg
      stop
    end if

  end subroutine

  ! -------------------------------------------------------------------------------------
  ! Provide num random numbers following a uniform distribution between 0 and 1
  !
  function get_rng_mkl_vec(this, num) result(u)

    class(ty_rng_mklvsl) :: this
    integer, intent(in) :: num
    real(dp), dimension(num) :: u

    integer :: status
    status = vdrnguniform(rng_method, this%stream, num, u, 0._dp, 1._dp)
    if (status /= VSL_STATUS_OK) call stop_on_err("Error getting random numbers")

  end function get_rng_mkl_vec

  ! -------------------------------------------------------------------------------------
  ! Provide random numbers for the TRUE elements of MASK
  !
  function get_rng_mkl_mask(this, mask) result(u)

    class(ty_rng_mklvsl) :: this
    logical, dimension(:), intent(in) :: mask
    real(dp), dimension(size(mask)) :: u

    u(:) = unpack(get_rng_mkl_vec(this, count(mask)), &
      MASK = mask, FIELD = 0._dp)

  end function get_rng_mkl_mask

  ! -------------------------------------------------------------------------------------
  ! Initialize the random number state.
  !
  subroutine init_rng_mkl(this, seeds)

    class(ty_rng_mklvsl) :: this
    integer, dimension(:), intent(in) :: seeds

    integer :: status
    status = vslnewstream(this%stream, rng_type, seeds(1))
    if (status /= VSL_STATUS_OK) &
      call stop_on_err("Error initializing random number stream")

  end subroutine init_rng_mkl

  ! -------------------------------------------------------------------------------------
  ! Release any resources associated with the RNG.
  !
  subroutine end_rng_mkl(this)

    class(ty_rng_mklvsl) :: this

    integer :: status
    status = vsldeletestream(this%stream)
    if (status /= VSL_STATUS_OK) &
      call stop_on_err("Error finalizing random number stream")

  end subroutine end_rng_mkl
  ! -------------------------------------------------------------------------------------

end module mo_rng_mklvsl
