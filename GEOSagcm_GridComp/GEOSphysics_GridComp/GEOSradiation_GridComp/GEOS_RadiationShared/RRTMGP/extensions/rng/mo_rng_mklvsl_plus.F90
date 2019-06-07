! Module: mo_rng_mklvsl_plus

! An extension by Peter Norris, of GMAO, NASA/GSFC, of mo_rng_mklvsl 
! from mo_rng_mklvsl.F90 of the RTE+RRTMGP distribution, the original 
! copyright message included here:
! -------------------------------------------------------------------
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
! -------------------------------------------------------------------

!
! Revised Implementation of random number interface using Intel Vector
! Statistics Library. Extends mo_rng_mklvsl to provide flexibility thru
! the initializers of which MKLVSL BRNG to use. Also uses ACCURATE mode
! of the uniform [0.,1.) generator.
!
! Note on accurate (_ACCURATE) mode uniform distribution:
! If only need fast (_STD) mode, can remove the get_rng_mkl_*_accurate()
! and just inherit the get_rng_mkl_*() from mo_rng_mklvsl
!

module mo_rng_mklvsl_plus

  use MKL_VSL_TYPE
  use MKL_VSL
  use mo_rte_kind, only: dp
  !use mo_rrtmgp_errors, only: write_message ! pmn: replaced by stop_on_err
  use mo_rng_mklvsl, only: ty_rng_mklvsl
  implicit none
  private
  public :: ty_rng_mklvsl_plus

  ! extended random number generator class
  type, extends(ty_rng_mklvsl) :: ty_rng_mklvsl_plus
  contains
    procedure :: get_random_vec      => get_rng_mkl_vec_accurate
    procedure :: get_random_vec_mask => get_rng_mkl_mask_accurate
    procedure :: init_rng_mkl_plus_scalar
    procedure :: init_rng_mkl_plus_vector
    generic, public :: init          => &
      init_rng_mkl_plus_scalar, init_rng_mkl_plus_vector
  end type

contains

  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg

    if(msg /= "") then
      write(error_unit, *) msg
      stop
    end if
  end subroutine

  ! -------------------------------------------------------------------------------------
  ! Provide num random numbers following a uniform distribution between 0 and 1
  !
  function get_rng_mkl_vec_accurate(this, num) result(u)

    class(ty_rng_mklvsl_plus) :: this
    integer,  intent(in) :: num
    real(DP), dimension(num) :: u

    integer :: status
    status = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, &
                          this%stream, num, u, 0._dp, 1._dp)
    if (status /= VSL_STATUS_OK) &
      call stop_on_err("Error getting random numbers")

  end function get_rng_mkl_vec_accurate

  ! -------------------------------------------------------------------------------------
  ! Provide random numbers for the TRUE elements of MASK
  !
  function get_rng_mkl_mask_accurate(this, mask) result(u)

    class(ty_rng_mklvsl_plus) :: this
    logical, dimension(:), intent(in) :: mask
    real(DP), dimension(size(mask)) :: u

    u(:) = unpack(get_rng_mkl_vec_accurate(this, COUNT(mask)), &
                  MASK = mask, FIELD = 0._dp)

  end function get_rng_mkl_mask_accurate

  ! -------------------------------------------------------------------------------------
  ! Initialize the random number state.
  !
  subroutine init_rng_mkl_plus_scalar(this, brng, seed)

    class(ty_rng_mklvsl_plus) :: this
    integer, intent(in) :: brng
    integer, intent(in) :: seed

    integer :: status
    status = vslNewStream(this%stream, brng, seed)
    if (status /= VSL_STATUS_OK) &
      call stop_on_err("Error initializing random number stream: scalar")

  end subroutine init_rng_mkl_plus_scalar
  !
  subroutine init_rng_mkl_plus_vector(this, brng, seeds)

    class(ty_rng_mklvsl_plus) :: this
    integer, intent(in) :: brng
    integer, dimension(:), intent(in) :: seeds

    integer :: status
    status = vslNewStreamEx(this%stream, brng, size(seeds), seeds)
    if (status /= VSL_STATUS_OK) &
      call stop_on_err("Error initializing random number stream: vector")

  end subroutine init_rng_mkl_plus_vector
  ! -------------------------------------------------------------------------------------

end module mo_rng_mklvsl_plus
