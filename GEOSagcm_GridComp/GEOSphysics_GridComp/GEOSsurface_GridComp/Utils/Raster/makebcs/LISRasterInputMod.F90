! Read a LIS-indexed GEOS surface mask without Catchment land parameters.
! NetCDF variables are accessed by name rather than by variable order.

module LISRasterInputMod

  use iso_fortran_env, only: real32, real64, int64, error_unit
  use netcdf
  use LISRasterGeometryMod, only: reference_mask_mode, &
                                  mask_surface_type,   &
                                  lis_id_at_pixel

  implicit none
  private

  integer, parameter :: LIS_CACHE_SIZE   = 512 * 1024 * 1024
  integer, parameter :: LIS_CACHE_NELEMS = 10007

  public :: read_lis_raster

contains

  subroutine require(condition, message)

    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not. condition) then
       write(error_unit, '(A)') 'LIS raster: '//trim(message)
       error stop 1
    end if

  end subroutine require


  subroutine nc_check(status, context)

    integer, intent(in) :: status
    character(*), intent(in) :: context

    if (status /= nf90_noerr) then
       write(error_unit, '(A)') &
            'LIS raster: '//trim(context)//': '//trim(nf90_strerror(status))
       error stop 1
    end if

  end subroutine nc_check


  integer function dim_length(ncid, name) result(length)

    integer, intent(in) :: ncid
    character(*), intent(in) :: name

    integer :: dimid

    call nc_check(nf90_inq_dimid(ncid, name, dimid), &
                  'dimension '//name)

    call nc_check(nf90_inquire_dimension(ncid, dimid, len=length), &
                  'length of '//name)

  end function dim_length


  integer function variable(ncid, name, dimensions, xtype) result(varid)

    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    character(*), intent(in) :: dimensions(:)
    integer, optional, intent(in) :: xtype

    integer :: ndims
    integer :: vartype
    integer :: dimids(nf90_max_var_dims)
    integer :: k
    integer :: status

    character(nf90_max_name) :: dimname

    real(real64) :: scale
    real(real64) :: offset

    call nc_check(nf90_inq_varid(ncid, name, varid), &
                  'variable '//name)

    call nc_check(nf90_inquire_variable(ncid, varid, &
                                        xtype=vartype, &
                                        ndims=ndims,   &
                                        dimids=dimids), &
                  'schema '//name)

    call require(ndims == size(dimensions), &
                 'wrong rank for '//name)

    if (present(xtype)) then
       call require(vartype == xtype, &
                    'wrong data type for '//name)
    end if

    do k = 1, ndims

       call nc_check(nf90_inquire_dimension(ncid, dimids(k), &
                                            name=dimname), &
                     'dimension name for '//name)

       call require(trim(dimname) == trim(dimensions(k)), &
                    'wrong dimension order for '//name)

    end do

    ! nf90_get_var does not apply packed-data scale/offset. Do not
    ! interpret scaled IDs or coordinates as physical values.
    !
    ! nf90_get_att has INTENT(OUT), so a missing attribute does not
    ! preserve a value assigned before the call.

    status = nf90_get_att(ncid, varid, 'scale_factor', scale)

    if (status == nf90_enotatt) then
       scale = 1.0_real64
    else
       call nc_check(status, 'scale_factor '//name)
    end if

    status = nf90_get_att(ncid, varid, 'add_offset', offset)

    if (status == nf90_enotatt) then
       offset = 0.0_real64
    else
       call nc_check(status, 'add_offset '//name)
    end if

    call require(scale == 1.0_real64 .and. offset == 0.0_real64, &
                 'packed/scaled variables unsupported: '//name)

  end function variable


  subroutine read_axis(ncid, name, dimname, axis)

    integer, intent(in) :: ncid
    character(*), intent(in) :: name
    character(*), intent(in) :: dimname
    real(real64), intent(out) :: axis(:)

    integer :: varid

    varid = variable(ncid, name, [dimname])

    call nc_check(nf90_get_var(ncid, varid, axis), &
                  'reading '//name)

  end subroutine read_axis


  subroutine check_global_centers(axis, west, span, name)

    real(real64), intent(in) :: axis(:)
    real(real64), intent(in) :: west
    real(real64), intent(in) :: span
    character(*), intent(in) :: name

    real(real64) :: tolerance
    real(real64) :: expected

    integer :: k

    tolerance = max(1.e-8_real64, &
                    2.0_real64 * &
                    real(spacing(real(maxval(abs(axis)), real32)), real64))

    do k = 1, size(axis)

       expected = west + &
                  (real(k, real64) - 0.5_real64) * span / size(axis)

       ! This also rejects NaNs. The high-resolution reference
       ! coordinates contain float32 quantization despite being
       ! stored as float64.

       call require(abs(axis(k) - expected) <= tolerance, &
                    'non-global or misaligned '//name)

    end do

  end subroutine check_global_centers


  subroutine read_lis_raster(filename, reference, nx, ny, raster, &
                             n_lis, expected_ref_count, verbose)

    character(*), intent(in) :: filename
    character(*), intent(in) :: reference

    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: expected_ref_count

    integer, pointer, intent(out) :: raster(:,:)
    integer, intent(out) :: n_lis

    logical, intent(in) :: verbose

    integer :: ncid
    integer :: refid
    integer :: lis_varid
    integer :: reference_varid
    integer :: varid

    integer :: reference_count
    integer :: hr_im
    integer :: hr_jm
    integer :: lis_im
    integer :: lis_jm

    integer :: x_factor
    integer :: y_factor

    integer :: i
    integer :: j
    integer :: k
    integer :: block_i
    integer :: block_j

    integer :: lis_id
    integer :: surface
    integer :: reference_surface
    integer :: selected_id
    integer :: x_start

    integer, allocatable :: lis_values(:,:)
    integer, allocatable :: reference_values(:,:)
    integer, allocatable :: lookup(:)

    real(real64), allocatable :: lat(:)
    real(real64), allocatable :: lon(:)
    real(real64), allocatable :: reference_lat(:)
    real(real64), allocatable :: reference_lon(:)
    real(real64), allocatable :: native_lat(:,:)
    real(real64), allocatable :: native_lon(:,:)

    real(real64) :: expected
    real(real64) :: tolerance

    integer(int64) :: surface_counts(0:3)

    call require(nx > 0 .and. ny > 0, &
                 'working raster dimensions must be positive')

    call require(int(nx, int64) * int(ny, int64) <= &
                 int(huge(0), int64), &
                 'working raster exceeds current 32-bit tile pixel '// &
                 'counters; use 30 arcsec or coarser')

    call require(len_trim(reference) > 0, &
                 'LIS mode requires the original GEOS reference mask (-r)')

    ! Keep one full zonal band of the supplied LISIndex chunks in
    ! memory. The production file uses 90 chunks of
    ! 1440 x 720 x 4 bytes, approximately 356 MiB.
    !
    ! Three-row reads revisit that band 240 times. A smaller NetCDF
    ! cache can therefore cause repeated decompression. This affects
    ! only I/O performance, not the input or reduction algorithm.

    call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid, &
                            cache_size=LIS_CACHE_SIZE, &
                            cache_nelems=LIS_CACHE_NELEMS), &
                  'opening '//filename)

    call nc_check(nf90_open(trim(reference), nf90_nowrite, refid), &
                  'opening '//reference)

    ! ------------------------------------------------------------------
    ! Input dimensions
    ! ------------------------------------------------------------------

    hr_im = dim_length(ncid, 'N_lon')
    hr_jm = dim_length(ncid, 'N_lat')

    call require(hr_im == dim_length(refid, 'N_lon') .and. &
                 hr_jm == dim_length(refid, 'N_lat'), &
                 'LIS and reference high-resolution dimensions differ')

    lis_im = dim_length(ncid, 'east_west')
    lis_jm = dim_length(ncid, 'north_south')

    call require(lis_im > 0 .and. lis_jm > 0, &
                 'native dimensions must be positive')

    n_lis = dim_length(ncid, 'N_LIS_IDs')

    call require(int(lis_im, int64) * int(lis_jm, int64) == &
                 int(n_lis, int64), &
                 'N_LIS_IDs differs from native grid size')

    call require(n_lis > 0 .and. n_lis < 10000000, &
                 'LIS source IDs exceed the supported land/ocean '// &
                 'ordering key range')

    reference_count = dim_length(refid, 'N_catchments')

    call require(reference_count == expected_ref_count, &
                 'reference catchment ID bound differs from this GEOS build')

    call require(mod(hr_im, nx) == 0 .and. mod(hr_jm, ny) == 0, &
                 'working raster must divide the source raster')

    call require(mod(nx, lis_im) == 0 .and. mod(ny, lis_jm) == 0, &
                 'working pixels must nest inside LIS cells; choose an '// &
                 'aligned, sufficiently fine raster')

    x_factor = hr_im / nx
    y_factor = hr_jm / ny

    call require(x_factor >= 1 .and. y_factor >= 1, &
                 'working raster is finer than source')

    ! ------------------------------------------------------------------
    ! Input variables and high-resolution coordinates
    ! ------------------------------------------------------------------

    lis_varid = variable(ncid, 'LISIndex', &
                         [character(5) :: 'N_lon', 'N_lat'], nf90_int)

    reference_varid = variable(refid, 'CatchIndex', &
                               [character(5) :: 'N_lon', 'N_lat'], &
                               nf90_int)

    allocate(lat(hr_jm))
    allocate(lon(hr_im))
    allocate(reference_lat(hr_jm))
    allocate(reference_lon(hr_im))

    call read_axis(ncid, 'latitude', 'N_lat', lat)
    call read_axis(ncid, 'longitude', 'N_lon', lon)

    call read_axis(refid, 'latitude', 'N_lat', reference_lat)
    call read_axis(refid, 'longitude', 'N_lon', reference_lon)

    call require(all(lat == reference_lat) .and. &
                 all(lon == reference_lon), &
                 'LIS/reference coordinate values differ')

    call check_global_centers(lat, -90.0_real64, 180.0_real64, &
                              'high-resolution latitude')

    call check_global_centers(lon, -180.0_real64, 360.0_real64, &
                              'high-resolution longitude')

    deallocate(lat)
    deallocate(lon)
    deallocate(reference_lat)
    deallocate(reference_lon)

    ! ------------------------------------------------------------------
    ! Native LIS grid and ID ordering
    ! ------------------------------------------------------------------

    allocate(lookup(n_lis))
    allocate(native_lat(lis_im, lis_jm))
    allocate(native_lon(lis_im, lis_jm))

    varid = variable(ncid, 'LISgridID', &
                     [character(9) :: 'N_LIS_IDs'], nf90_int)

    call nc_check(nf90_get_var(ncid, varid, lookup), &
                  'reading LISgridID')

    do k = 1, n_lis
       call require(lookup(k) == k, &
                    'LISgridID must be sequential 1..N '// &
                    '(row-major native grid)')
    end do

    varid = variable(ncid, 'LIS_lat', &
                     [character(11) :: 'east_west', 'north_south'])

    call nc_check(nf90_get_var(ncid, varid, native_lat), &
                  'reading LIS_lat')

    varid = variable(ncid, 'LIS_lon', &
                     [character(11) :: 'east_west', 'north_south'])

    call nc_check(nf90_get_var(ncid, varid, native_lon), &
                  'reading LIS_lon')

    tolerance = 1.e-5_real64

    do j = 1, lis_jm

       expected = -90.0_real64 + &
                  (real(j, real64) - 0.5_real64) * &
                  180.0_real64 / lis_jm

       call require(all(abs(native_lat(:,j) - expected) <= tolerance), &
                    'LIS_lat must use global DE/PE cell centers')

    end do

    do i = 1, lis_im

       expected = -180.0_real64 + &
                  (real(i, real64) - 0.5_real64) * &
                  360.0_real64 / lis_im

       call require(all(abs(native_lon(i,:) - expected) <= tolerance), &
                    'LIS_lon must use global DE/PE cell centers')

    end do

    deallocate(lookup)
    deallocate(native_lat)
    deallocate(native_lon)

    ! ------------------------------------------------------------------
    ! Reduce the high-resolution masks to the working raster
    ! ------------------------------------------------------------------

    allocate(raster(nx, ny))
    allocate(lis_values(hr_im, y_factor))
    allocate(reference_values(hr_im, y_factor))

    raster = 0
    surface_counts = 0_int64

    do j = 1, ny

       call nc_check( &
            nf90_get_var(ncid, lis_varid, lis_values, &
                         start=[1, (j - 1) * y_factor + 1], &
                         count=[hr_im, y_factor]), &
            'LISIndex row block')

       call nc_check( &
            nf90_get_var(refid, reference_varid, reference_values, &
                         start=[1, (j - 1) * y_factor + 1], &
                         count=[hr_im, y_factor]), &
            'CatchIndex row block')

       do i = 1, nx

          x_start = (i - 1) * x_factor

          lis_id = lis_id_at_pixel(i, j, nx, ny, lis_im, lis_jm)

          do block_j = 1, y_factor
             do block_i = x_start + 1, x_start + x_factor

                surface = mask_surface_type( &
                     lis_values(block_i, block_j), n_lis)

                reference_surface = mask_surface_type( &
                     reference_values(block_i, block_j), reference_count)

                if (surface < 0 .or. reference_surface < 0 .or. &
                    surface /= reference_surface) then

                   write(error_unit, *) &
                        'native pixel (i,j), reference, LIS: ', &
                        block_i, (j - 1) * y_factor + block_j, &
                        reference_values(block_i, block_j), &
                        lis_values(block_i, block_j)

                   call require(.false., &
                                'invalid value or changed high-resolution '// &
                                'surface class')

                end if

                if (surface == 100 .and. &
                    lis_values(block_i, block_j) /= lis_id) then

                   write(error_unit, *) &
                        'native pixel (i,j), expected LIS ID, actual: ', &
                        block_i, (j - 1) * y_factor + block_j, &
                        lis_id, lis_values(block_i, block_j)

                   call require(.false., &
                                'LIS ID does not contain its '// &
                                'high-resolution land pixel')

                end if

             end do
          end do

          ! Reproduce the standard GEOS mask reduction using the original
          ! CatchIndex values. If the selected identifier is land, replace
          ! it with the native LIS cell ID.

          selected_id = reference_mask_mode( &
               reference_values(x_start + 1:x_start + x_factor, :), &
               reference_count)

          if (selected_id >= 1 .and. selected_id <= reference_count) then
             selected_id = lis_id
          end if

          raster(i,j) = selected_id

          select case (mask_surface_type(selected_id, n_lis))

          case (0)
             surface_counts(0) = surface_counts(0) + 1

          case (100)
             surface_counts(1) = surface_counts(1) + 1

          case (19)
             surface_counts(2) = surface_counts(2) + 1

          case (20)
             surface_counts(3) = surface_counts(3) + 1

          end select

       end do

       if (verbose .and. (mod(j, 200) == 0 .or. j == ny)) then
          write(*, *) 'LIS mask rows: ', j, ' / ', ny
       end if

    end do

    call nc_check(nf90_close(ncid), &
                  'closing LIS mask')

    call nc_check(nf90_close(refid), &
                  'closing reference mask')

    write(*, *) 'LIS native dimensions (i,j), ID bound: ', &
                lis_im, lis_jm, n_lis

    write(*, *) 'Working-raster pixels: ocean, land, lake, ice: ', &
                surface_counts

    write(*, '(A)') &
         'Surface selection: original GEOS identifier mode; '// &
         'LIS IDs only on resulting land.'

  end subroutine read_lis_raster

end module LISRasterInputMod
