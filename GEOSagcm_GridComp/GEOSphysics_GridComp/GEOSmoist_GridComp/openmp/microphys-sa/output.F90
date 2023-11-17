module output_mod

  use constants_mod, only: fmt_e

  implicit none

  private

  public OutputArrays_T, write_difference

  type OutputArrays_T
     real, allocatable, dimension (:, :) :: rain, snow, ice, graupel
     real, allocatable, dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
     real, allocatable, dimension (:, :, :) :: revap ! Rain evaporation
     real, allocatable, dimension (:, :, :) :: isubl ! Ice sublimation
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
    allocate(arr%rain(iis:iie, jjs:jje), source = 0.)
    allocate(arr%snow, arr%ice, arr%graupel, mold = arr%rain)
    arr%snow = 0.
    arr%ice = 0.
    arr%graupel = 0.
    arr%rain = 0.
    
    allocate(arr%m2_rain(iis:iie, jjs:jje, kks:kke), source = 0.)
    allocate(arr%m2_sol, arr%revap, arr%isubl, mold = arr%m2_rain)
    arr%m2_sol = 0.
    arr%revap = 0.
    arr%isubl = 0.

  end function initialize_

  subroutine write_difference(arr1, arr2)

    ! Arguments
    type(OutputArrays_T), intent(in) :: arr1, arr2

    ! Start
    print *, ''
    print *, '-----------|-----------------'
    print *, '   out var |  abs error'
    print *, '-----------|-----------------'
    
    write(*, fmt_e) 'revap', '|', norm2(arr2%revap-arr1%revap) ! /norm2(arr1%revap)
    write(*, fmt_e) 'isubl', '|', norm2(arr2%isubl-arr1%isubl) ! /norm2(arr1%isubl)
    write(*, fmt_e) 'rain', '|', norm2(arr2%rain-arr1%rain) ! /norm2(arr1%rain)
    write(*, fmt_e) 'snow', '|', norm2(arr2%snow-arr1%snow) ! /norm2(arr1%snow)
    write(*, fmt_e) 'graupel', '|', norm2(arr2%graupel-arr1%graupel) ! /norm2(arr1%graupel)
    write(*, fmt_e) 'm2_rain', '|', norm2(arr2%m2_rain-arr1%m2_rain) ! /norm2(arr1%m2_rain)
    write(*, fmt_e) 'm2_sol', '|', norm2(arr2%m2_sol-arr1%m2_sol) ! /norm2(arr1%m2_sol)
    print *, '-----------|-----------------'

  end subroutine write_difference

end module output_mod
