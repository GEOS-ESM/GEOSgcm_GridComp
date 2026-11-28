! LIS geometry helpers. No NetCDF or MAPL dependencies in this module.

module LISRasterGeometryMod

  implicit none
  private

  integer, parameter, public :: LIS_OCEAN = 0
  integer, parameter, public :: LIS_LAKE  = 190000000
  integer, parameter, public :: LIS_ICE   = 200000000

  public :: reference_mask_mode
  public :: mask_surface_type
  public :: lis_id_at_pixel

contains

  pure integer function mask_surface_type(value, max_land_id) result(surface)

    integer, intent(in) :: value
    integer, intent(in) :: max_land_id

    surface = -1

    if (value == LIS_OCEAN) then
       surface = 0
    else if (value >= 1 .and. value <= max_land_id) then
       surface = 100
    else if (value == LIS_LAKE) then
       surface = 19
    else if (value == LIS_ICE) then
       surface = 20
    end if

  end function mask_surface_type


  pure integer function reference_mask_mode(block, max_land_id) result(winner)

    ! Reproduce the choice made for GEOS5_10arcsec_mask in mkLandRaster:
    ! exclude ocean, count identifiers rather than surface classes, and
    ! choose the lowest encoded identifier on ties.
    !
    ! Lake and land-ice values are temporarily placed above the land-ID
    ! range so that the original identifier ordering is preserved.

    integer, intent(in) :: block(:,:)
    integer, intent(in) :: max_land_id

    integer :: values(size(block))
    integer :: n
    integer :: i
    integer :: j
    integer :: value
    integer :: run_count
    integer :: best_count

    winner = block(1,1)

    if (all(block == winner)) return

    winner = LIS_OCEAN
    n = 0

    do j = 1, size(block, 2)
       do i = 1, size(block, 1)

          value = block(i,j)

          if (value == LIS_OCEAN) cycle

          if (value == LIS_LAKE) value = max_land_id + 1
          if (value == LIS_ICE)  value = max_land_id + 2

          n = n + 1
          values(n) = value

       end do
    end do

    if (n == 0) return

    ! Production reduction blocks are small (3x3). Use an allocation-free
    ! insertion sort before counting repeated identifiers.

    do i = 2, n

       value = values(i)
       j = i - 1

       do while (j >= 1)
          if (values(j) <= value) exit

          values(j + 1) = values(j)
          j = j - 1
       end do

       values(j + 1) = value

    end do

    winner = values(1)
    best_count = 0
    i = 1

    do while (i <= n)

       value = values(i)
       run_count = 0

       do while (i <= n)
          if (values(i) /= value) exit

          run_count = run_count + 1
          i = i + 1
       end do

       if (run_count > best_count) then
          winner = value
          best_count = run_count
       end if

    end do

    if (winner == max_land_id + 1) winner = LIS_LAKE
    if (winner == max_land_id + 2) winner = LIS_ICE

  end function reference_mask_mode


  pure integer function lis_id_at_pixel(i, j, nx, ny, lis_im, lis_jm) result(id)

    ! One-based indices. The caller verifies that the working raster nests
    ! exactly within the global regular LIS grid.

    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: lis_im
    integer, intent(in) :: lis_jm

    integer :: lis_i
    integer :: lis_j

    lis_i = (i - 1) / (nx / lis_im) + 1
    lis_j = (j - 1) / (ny / lis_jm) + 1

    id = (lis_j - 1) * lis_im + lis_i

  end function lis_id_at_pixel

end module LISRasterGeometryMod
