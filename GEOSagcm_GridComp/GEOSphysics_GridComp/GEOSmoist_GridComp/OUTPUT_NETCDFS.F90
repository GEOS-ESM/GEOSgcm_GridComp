module OUTPUT_NETCDFS
  use netcdf
  implicit none
  private
  public :: create_netcdf_file, add_variable

  ! Name of the unlimited time dimension in every file
  character(len=*), parameter :: TIME_DIM = "time"

  interface add_variable
    ! real(4)
    module procedure add_r4_0d, add_r4_1d, add_r4_2d, add_r4_3d
    ! real(8)
    module procedure add_r8_0d, add_r8_1d, add_r8_2d, add_r8_3d
    ! integer(4)
    module procedure add_i4_0d, add_i4_1d, add_i4_2d, add_i4_3d
    ! integer(8)
    module procedure add_i8_0d, add_i8_1d, add_i8_2d, add_i8_3d
    ! logical
    module procedure add_l_0d,  add_l_1d,  add_l_2d,  add_l_3d
  end interface add_variable

contains

  ! ---------------------------------------------------------------------------
  ! File creation — creates an empty file with an unlimited time dimension
  ! ---------------------------------------------------------------------------

  subroutine create_netcdf_file(file_name)
    character(len=*), intent(in) :: file_name
    integer :: ncid, tdimid
    integer :: slash_pos
    character(len=512) :: cwd

    ! Print the directory where the file will be saved
    slash_pos = scan(file_name, "/", back=.true.)
    if (slash_pos > 0) then
      print *, "NetCDF saving to directory: ", file_name(1:slash_pos)
    else
      call get_environment_variable("PWD", cwd)
      print *, "NetCDF saving to directory: ", trim(cwd)
    end if

    call check(nf90_create(file_name, NF90_NETCDF4, ncid))
    call check(nf90_def_dim(ncid, TIME_DIM, NF90_UNLIMITED, tdimid))
    call check(nf90_close(ncid))
  end subroutine create_netcdf_file


  ! ---------------------------------------------------------------------------
  ! Internal helpers: open/define/close for first call, then write per type
  ! ---------------------------------------------------------------------------

  ! Open the file, define the variable if needed, return ncid+varid+start+cnt
  !
  ! Time index is tracked per-variable using a hidden scalar counter variable
  ! named "__t__VARNAME" stored in the same file. This means multiple variables
  ! can be written independently without interfering with each other's time index.
  subroutine prepare_var(file_name, var_name, nc_type, dims, ncid, varid, start, cnt)
    character(len=*), intent(in)  :: file_name, var_name
    integer,          intent(in)  :: nc_type, dims(:)
    integer,          intent(out) :: ncid, varid
    integer, allocatable, intent(out) :: start(:), cnt(:)

    integer :: tdimid, ndims, t_index, status
    integer :: xdimid, ydimid, zdimid
    integer :: ctr_varid, ctr_dimid
    character(len=256) :: ctr_name

    ndims    = size(dims)
    ctr_name = '__t__'//trim(var_name)

    call check(nf90_open(file_name, NF90_WRITE, ncid))
    call check(nf90_inq_dimid(ncid, TIME_DIM, tdimid))

    status = nf90_inq_varid(ncid, var_name, varid)
    if (status /= NF90_NOERR) then
      ! First call: define the data variable and its counter
      call check(nf90_redef(ncid))

      if (ndims == 0) then
        call check(nf90_def_var(ncid, var_name, nc_type, [tdimid], varid))
      else if (ndims == 1) then
        call check(nf90_def_dim(ncid, 'dim_x_'//trim(var_name), dims(1), xdimid))
        call check(nf90_def_var(ncid, var_name, nc_type, [xdimid, tdimid], varid))
      else if (ndims == 2) then
        call check(nf90_def_dim(ncid, 'dim_x_'//trim(var_name), dims(1), xdimid))
        call check(nf90_def_dim(ncid, 'dim_y_'//trim(var_name), dims(2), ydimid))
        call check(nf90_def_var(ncid, var_name, nc_type, [xdimid, ydimid, tdimid], varid))
      else if (ndims == 3) then
        call check(nf90_def_dim(ncid, 'dim_x_'//trim(var_name), dims(1), xdimid))
        call check(nf90_def_dim(ncid, 'dim_y_'//trim(var_name), dims(2), ydimid))
        call check(nf90_def_dim(ncid, 'dim_z_'//trim(var_name), dims(3), zdimid))
        call check(nf90_def_var(ncid, var_name, nc_type, [xdimid, ydimid, zdimid, tdimid], varid))
      end if

      ! Define a scalar counter to track how many timesteps have been written
      call check(nf90_def_var(ncid, trim(ctr_name), NF90_INT, varid=ctr_varid))
      call check(nf90_enddef(ncid))

      t_index = 1
      call check(nf90_put_var(ncid, ctr_varid, t_index))
    else
      ! Subsequent call: read and increment the per-variable counter
      call check(nf90_inq_varid(ncid, trim(ctr_name), ctr_varid))
      call check(nf90_get_var(ncid, ctr_varid, t_index))
      t_index = t_index + 1
      call check(nf90_put_var(ncid, ctr_varid, t_index))
    end if

    if (ndims == 0) then
      allocate(start(1), cnt(1))
      start = [t_index];  cnt = [1]
    else
      allocate(start(ndims+1), cnt(ndims+1))
      start(1:ndims) = 1;     start(ndims+1) = t_index
      cnt(1:ndims)   = dims;  cnt(ndims+1)   = 1
    end if
  end subroutine prepare_var


  ! ---------------------------------------------------------------------------
  ! real(4)
  ! ---------------------------------------------------------------------------

  subroutine add_r4_0d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(4),          intent(in) :: array
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer :: dims(0)
    call prepare_var(file_name, var_name, NF90_FLOAT, dims, ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_r4_1d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(4),          intent(in) :: array(:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_FLOAT, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_r4_2d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(4),          intent(in) :: array(:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_FLOAT, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_r4_3d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(4),          intent(in) :: array(:,:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_FLOAT, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine


  ! ---------------------------------------------------------------------------
  ! real(8)
  ! ---------------------------------------------------------------------------

  subroutine add_r8_0d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(8),          intent(in) :: array
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer :: dims(0)
    call prepare_var(file_name, var_name, NF90_DOUBLE, dims, ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_r8_1d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(8),          intent(in) :: array(:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_DOUBLE, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_r8_2d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(8),          intent(in) :: array(:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_DOUBLE, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_r8_3d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    real(8),          intent(in) :: array(:,:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_DOUBLE, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine


  ! ---------------------------------------------------------------------------
  ! integer(4)
  ! ---------------------------------------------------------------------------

  subroutine add_i4_0d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(4),          intent(in) :: array
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer :: dims(0)
    call prepare_var(file_name, var_name, NF90_INT, dims, ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_i4_1d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(4),       intent(in) :: array(:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_INT, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_i4_2d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(4),       intent(in) :: array(:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_INT, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_i4_3d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(4),       intent(in) :: array(:,:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_INT, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine


  ! ---------------------------------------------------------------------------
  ! integer(8)
  ! ---------------------------------------------------------------------------

  subroutine add_i8_0d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(8),          intent(in) :: array
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer :: dims(0)
    call prepare_var(file_name, var_name, NF90_INT64, dims, ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_i8_1d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(8),       intent(in) :: array(:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_INT64, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_i8_2d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(8),       intent(in) :: array(:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_INT64, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_i8_3d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    integer(8),       intent(in) :: array(:,:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    call prepare_var(file_name, var_name, NF90_INT64, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, array, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine


  ! ---------------------------------------------------------------------------
  ! logical (stored as NF90_BYTE, 0/1)
  ! ---------------------------------------------------------------------------

  subroutine add_l_0d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    logical,          intent(in) :: array
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer :: dims(0)
    integer(1) :: b
    b = merge(1_1, 0_1, array)
    call prepare_var(file_name, var_name, NF90_BYTE, dims, ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, b))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_l_1d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    logical,          intent(in) :: array(:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer(1), allocatable :: b(:)
    allocate(b(size(array)))
    where (array); b = 1_1; elsewhere; b = 0_1; end where
    call prepare_var(file_name, var_name, NF90_BYTE, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, b, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_l_2d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    logical,          intent(in) :: array(:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer(1), allocatable :: b(:,:)
    allocate(b(size(array,1), size(array,2)))
    where (array); b = 1_1; elsewhere; b = 0_1; end where
    call prepare_var(file_name, var_name, NF90_BYTE, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, b, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine

  subroutine add_l_3d(file_name, array, var_name)
    character(len=*), intent(in) :: file_name, var_name
    logical,          intent(in) :: array(:,:,:)
    integer :: ncid, varid
    integer, allocatable :: start(:), cnt(:)
    integer(1), allocatable :: b(:,:,:)
    allocate(b(size(array,1), size(array,2), size(array,3)))
    where (array); b = 1_1; elsewhere; b = 0_1; end where
    call prepare_var(file_name, var_name, NF90_BYTE, shape(array), ncid, varid, start, cnt)
    call check(nf90_put_var(ncid, varid, b, start=start, count=cnt))
    call check(nf90_close(ncid))
  end subroutine


  ! ---------------------------------------------------------------------------
  ! Utilities

  ! ---------------------------------------------------------------------------

  subroutine check(status)
    integer, intent(in) :: status
    if (status /= NF90_NOERR) then
      print *, "NetCDF error: ", trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine check

  pure function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=20)   :: s
    write(s, '(i0)') i
    s = adjustl(s)
  end function itoa

end module OUTPUT_NETCDFS