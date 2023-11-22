module output_mod

  implicit none

  private

  public OutputArrays_T, write_arrays, write_difference

  character(len=*), parameter :: fmt_diff = '(1x, a10, 1x, a1, 1x, e15.9, 1x, a1, 1x, e15.9)'
  character(len=*), parameter :: fmt_out = '(1x, a10, 3x, e18.10, 3x, e18.10, 3x, e18.10)'

  type OutputArrays_T
     real, allocatable, dimension (:, :) :: rain, snow, ice, graupel
     real, allocatable, dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
     real, allocatable, dimension (:, :, :) :: revap ! Rain evaporation
     real, allocatable, dimension (:, :, :) :: isubl ! Ice sublimation
   contains
     procedure, public :: write_arrays
  end type OutputArrays_T

  interface OutputArrays_T
     procedure :: initialize_
  end interface OutputArrays_T

contains

  function initialize_(iis, iie, jjs, jje, kks, kke) result (arr)

    ! Arguments
    integer, intent(in) :: iis, iie, jjs, jje, kks, kke
    type(OutputArrays_T) :: arr ! output

    ! Start
    ! print *, 'Initializing output arrays:'
    ! print *, iis, iie, jjs, jje, kks, kke
    allocate(arr%rain(iis:iie, jjs:jje), source = 0.)
    ! print *, 'shape(rain): ', shape(arr%rain)
    allocate(arr%snow, arr%ice, arr%graupel, mold = arr%rain)
    arr%snow = 0.
    arr%ice = 0.
    arr%graupel = 0.
    
    allocate(arr%m2_rain(iis:iie, jjs:jje, kks:kke), source = 0.)
    ! print *, 'shape(m2_rain): ', shape(arr%m2_rain)
    allocate(arr%m2_sol, arr%revap, arr%isubl, mold = arr%m2_rain)
    arr%m2_sol = 0.
    arr%revap = 0.
    arr%isubl = 0.

  end function initialize_

  subroutine write_arrays(self)

    ! Arguments
    class(OutputArrays_T), intent(in) :: self

    ! Start
    write(*, fmt_out) 'revap', minval(self%revap), maxval(self%revap), sum(self%revap)
    write(*, fmt_out) 'isubl', minval(self%isubl), maxval(self%isubl), sum(self%isubl)
    write(*, fmt_out) 'rain', minval(self%rain), maxval(self%rain), sum(self%rain)
    write(*, fmt_out) 'snow', minval(self%snow), maxval(self%snow), sum(self%snow)
    write(*, fmt_out) 'ice', minval(self%ice), maxval(self%ice), sum(self%ice)
    write(*, fmt_out) 'graupel', minval(self%graupel), maxval(self%graupel), sum(self%graupel)
    write(*, fmt_out) 'm2_rain', minval(self%m2_rain), maxval(self%m2_rain), sum(self%m2_rain)
    write(*, fmt_out) 'm2_sol', minval(self%m2_sol), maxval(self%m2_sol), sum(self%m2_sol)

  end subroutine write_arrays

  subroutine write_difference(arr1, arr2)

    ! Arguments
    type(OutputArrays_T), intent(in) :: arr1, arr2

    ! Start
    print *, ''
    print *, '-----------|-----------------|-----------------'
    print *, '   out var |     abs error   |     rel error'
    print *, '-----------|-----------------|-----------------'
    
    write(*, fmt_diff) 'revap', '|', norm2(arr2%revap-arr1%revap), '|', norm2(arr2%revap-arr1%revap)/norm2(arr1%revap)
    write(*, fmt_diff) 'isubl', '|', norm2(arr2%isubl-arr1%isubl), '|', norm2(arr2%isubl-arr1%isubl)/norm2(arr1%isubl)
    write(*, fmt_diff) 'rain', '|', norm2(arr2%rain-arr1%rain), '|', norm2(arr2%rain-arr1%rain)/norm2(arr1%rain)
    write(*, fmt_diff) 'snow', '|', norm2(arr2%snow-arr1%snow), '|', norm2(arr2%snow-arr1%snow)/norm2(arr1%snow)
    write(*, fmt_diff) 'ice', '|', norm2(arr2%ice-arr1%ice), '|', norm2(arr2%ice-arr1%ice)/norm2(arr1%ice)
    write(*, fmt_diff) 'graupel', '|', norm2(arr2%graupel-arr1%graupel), '|', norm2(arr2%graupel-arr1%graupel)/norm2(arr1%graupel)
    write(*, fmt_diff) 'm2_rain', '|', norm2(arr2%m2_rain-arr1%m2_rain), '|', norm2(arr2%m2_rain-arr1%m2_rain)/norm2(arr1%m2_rain)
    write(*, fmt_diff) 'm2_sol', '|', norm2(arr2%m2_sol-arr1%m2_sol), '|', norm2(arr2%m2_sol-arr1%m2_sol)/norm2(arr1%m2_sol)
    print *, '-----------|-----------------|-----------------'

  end subroutine write_difference

end module output_mod
