module igrf_input_module

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  public :: get_igrf_field

  real, parameter :: GAUSS_TO_TESLA = 1.0e-4

  interface
    subroutine FELDCOF(YEAR)
      real, intent(in) :: YEAR
    end subroutine FELDCOF

    subroutine FELDG(GLAT, GLON, ALT, BNORTH, BEAST, BDOWN, BABS)
      real, intent(in) :: GLAT, GLON, ALT
      real, intent(out) :: BNORTH, BEAST, BDOWN, BABS
    end subroutine FELDG
  end interface

contains

  subroutine get_igrf_field(latitudes_deg, longitudes_deg, iyear, doy, alt_km, &
                            b_north_t, b_east_t, b_down_t, b_mag_t)
    ! Evaluate the IGRF field bundled with IRI at each GEOS model point.
    ! FELDG returns local north/east/down components in Gauss. This wrapper
    ! converts all magnetic-field outputs to Tesla.

    real, intent(in) :: latitudes_deg(:,:), longitudes_deg(:,:)
    integer, intent(in) :: iyear, doy
    real, intent(in) :: alt_km(:,:,:)

    real, intent(out) :: b_north_t(:,:,:), b_east_t(:,:,:)
    real, intent(out) :: b_down_t(:,:,:), b_mag_t(:,:,:)

    integer :: i, j, k
    integer :: days_in_year
    real :: decimal_year
    real :: b_north_g, b_east_g, b_down_g, b_mag_g

    call validate_shapes(latitudes_deg, longitudes_deg, alt_km, &
                         b_north_t, b_east_t, b_down_t, b_mag_t)

    days_in_year = 365
    if (mod(iyear, 400) == 0 .or. &
        (mod(iyear, 4) == 0 .and. mod(iyear, 100) /= 0)) then
      days_in_year = 366
    end if

    decimal_year = real(iyear) + real(doy - 1) / real(days_in_year)

    ! Initialize/interpolate the IGRF coefficients for the requested date.
    ! Use FELDCOF/FELDG, not the legacy POGO FIELDG routine in igrf.for.
    call FELDCOF(decimal_year)

    b_north_t = 0.0
    b_east_t = 0.0
    b_down_t = 0.0
    b_mag_t = 0.0

    do k = 1, size(alt_km, 3)
      do j = 1, size(alt_km, 2)
        do i = 1, size(alt_km, 1)
          if (.not. ieee_is_finite(latitudes_deg(i,j)) .or. &
              .not. ieee_is_finite(longitudes_deg(i,j)) .or. &
              .not. ieee_is_finite(alt_km(i,j,k))) cycle

          call FELDG(latitudes_deg(i,j), longitudes_deg(i,j), alt_km(i,j,k), &
                     b_north_g, b_east_g, b_down_g, b_mag_g)

          if (.not. ieee_is_finite(b_north_g) .or. &
              .not. ieee_is_finite(b_east_g) .or. &
              .not. ieee_is_finite(b_down_g) .or. &
              .not. ieee_is_finite(b_mag_g)) cycle

          b_north_t(i,j,k) = b_north_g * GAUSS_TO_TESLA
          b_east_t(i,j,k)  = b_east_g  * GAUSS_TO_TESLA
          b_down_t(i,j,k)  = b_down_g  * GAUSS_TO_TESLA
          b_mag_t(i,j,k)   = b_mag_g   * GAUSS_TO_TESLA
        end do
      end do
    end do

  end subroutine get_igrf_field


  subroutine validate_shapes(latitudes_deg, longitudes_deg, alt_km, &
                             b_north_t, b_east_t, b_down_t, b_mag_t)
    real, intent(in) :: latitudes_deg(:,:), longitudes_deg(:,:)
    real, intent(in) :: alt_km(:,:,:)
    real, intent(out) :: b_north_t(:,:,:), b_east_t(:,:,:)
    real, intent(out) :: b_down_t(:,:,:), b_mag_t(:,:,:)

    integer :: nx, ny, nz

    nx = size(alt_km, 1)
    ny = size(alt_km, 2)
    nz = size(alt_km, 3)

    if (size(latitudes_deg,1) /= nx .or. size(latitudes_deg,2) /= ny) then
      error stop 'IGRF: latitude shape does not match alt_km'
    end if
    if (size(longitudes_deg,1) /= nx .or. size(longitudes_deg,2) /= ny) then
      error stop 'IGRF: longitude shape does not match alt_km'
    end if
    if (.not. same_shape_3d(b_north_t, nx, ny, nz)) &
      error stop 'IGRF: b_north_t shape does not match alt_km'
    if (.not. same_shape_3d(b_east_t, nx, ny, nz)) &
      error stop 'IGRF: b_east_t shape does not match alt_km'
    if (.not. same_shape_3d(b_down_t, nx, ny, nz)) &
      error stop 'IGRF: b_down_t shape does not match alt_km'
    if (.not. same_shape_3d(b_mag_t, nx, ny, nz)) &
      error stop 'IGRF: b_mag_t shape does not match alt_km'

  end subroutine validate_shapes


  logical function same_shape_3d(field, nx, ny, nz)
    real, intent(in) :: field(:,:,:)
    integer, intent(in) :: nx, ny, nz

    same_shape_3d = size(field,1) == nx .and. &
                    size(field,2) == ny .and. &
                    size(field,3) == nz
  end function same_shape_3d

end module igrf_input_module
