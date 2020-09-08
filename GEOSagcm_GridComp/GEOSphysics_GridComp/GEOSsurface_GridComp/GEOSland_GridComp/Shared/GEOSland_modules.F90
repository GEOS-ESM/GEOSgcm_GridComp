#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

MODULE GEOSland_modules

  use ESMF

  implicit none

  PRIVATE

  PUBLIC :: modis_date, read_modis_data

    ! ---------------------------------------------------------------------------
    
  integer function modis_date (DOY, interval) result (MOD_DOY)
    
    implicit none
    integer, intent(in) :: DOY, interval
    integer, parameter  :: N_MODIS_DAYS8 = 46, N_MODIS_DAYS5 = 74
    integer, dimension (N_MODIS_DAYS8), target ::     &
         MODIS_DOYS8 = (/                             &
         1  ,  9, 17, 25, 33, 41, 49, 57, 65,         &
         73 , 81, 89, 97,105,113,121,129,137,         &
         145,153,161,169,177,185,193,201,209,         &
         217,225,233,241,249,257,265,273,281,         &
         289,297,305,313,321,329,337,345,353,361/)
    
    integer, dimension (N_MODIS_DAYS5), target ::     &
         MODIS_DOYS5 = (/                             &
         1  ,  6, 11, 16, 21, 26, 31, 36, 41, 46,     &
         51 , 56, 61, 66, 71, 76, 81, 86, 91, 96,     &
         101,106,111,116,121,126,131,136,141,146,     &
         151,156,161,166,171,176,181,186,191,196,     &
         201,206,211,216,221,226,231,236,241,246,     &
         251,256,261,266,271,276,281,286,291,296,     &
         301,306,311,316,321,326,331,336,341,346,     &
         351,356,361,366/)
    integer, dimension(:), pointer :: MODIS_DOYS
    integer :: i,N_MODIS_DATES

    select case (interval)
    case (8)
       MODIS_DOYS => MODIS_DOYS8
       N_MODIS_DATES = N_MODIS_DAYS8
    case (5)
       MODIS_DOYS => MODIS_DOYS5
       N_MODIS_DATES = N_MODIS_DAYS5
    end select
    
    if (DOY < MODIS_DOYS(N_MODIS_DATES)) then
       do i = 1, N_MODIS_DATES
          if (MODIS_DOYS(i) > DOY) exit
       end do
       MOD_DOY = MODIS_DOYS(i-1)
    else
       MOD_DOY = MODIS_DOYS(N_MODIS_DATES)
    endif
    
  end function modis_date
    
  ! ---------------------------------------------------------------------------
  
  subroutine read_modis_data (MAPL,CUR_YY, MOD_DOY, b4_modis_date, &
       GRIDNAME, MODIS_PATH, label, MODIS_DATA, MODIS_NIR) 

    implicit none
    type(MAPL_MetaComp), intent(in),pointer  :: MAPL
    character (*), intent (in)               :: GRIDNAME, MODIS_PATH, label
    integer, intent (in)                     :: CUR_YY, MOD_DOY
    logical, intent (in)                     :: b4_modis_date
    real, dimension(:), intent(out)          :: MODIS_DATA
    real, dimension(:), intent(out),optional :: MODIS_NIR
    type(ESMF_Grid)                          :: TILEGRID
    type(MAPL_LocStream)                     :: LOCSTREAM
    integer, pointer                         :: mask(:)
    integer                                  :: status, unit
    character*300                            :: filename
    CHARACTER(len=7)                         :: YYYYDoY
    logical                                  :: file_exists
    
    if(b4_modis_date) then
       WRITE (YYYYDoY,'(a4,i3.3)') 'YYYY',MOD_DOY
    else
       WRITE (YYYYDoY,'(i4.4,i3.3)') CUR_YY,MOD_DOY
    endif
    
    filename = trim(MODIS_PATH)//'/'//trim(GRIDNAME)//'/'//trim(label)//'_data.'//YYYYDoY
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
       ASSERT_(.FALSE.)
    endif
    call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)            ; VERIFY_(STATUS)
    call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS); VERIFY_(STATUS)
    call MAPL_TileMaskGet(tilegrid,  mask, rc=status)              ; VERIFY_(STATUS)
    unit = GETFILE(trim(filename), form="unformatted", RC=STATUS)  ; VERIFY_(STATUS)
    call MAPL_VarRead(unit,tilegrid,MODIS_DATA,mask=mask,RC=STATUS); VERIFY_(STATUS)
    if(trim(label) == 'alb') &
         call MAPL_VarRead(unit,tilegrid,MODIS_NIR,mask=mask,RC=STATUS); VERIFY_(STATUS)
    call FREE_FILE(unit, RC=STATUS)                                ; VERIFY_(STATUS)
    
  end subroutine read_modis_data

END MODULE GEOSland_modules

! ********************************************************************

module GEOSland_io_hdf5

  use hdf5

  implicit none

  private

  integer, parameter :: UNINIT_INT = -99999
  character(len=*), parameter :: UNINIT_STR = ""

  type, public :: hdf5read
     private
     character(len=256) :: file_name = UNINIT_STR
     integer(hid_t) :: file_id = UNINIT_INT
     character(len=256) :: dset_name = UNINIT_STR
     integer(hid_t) :: dset_id = UNINIT_INT, dspace_id = UNINIT_INT, dtype_id = UNINIT_INT
     integer :: dset_rank = UNINIT_INT
     ! 7 is the max dimension of a fortran array
     integer(hsize_t) :: dset_size(7) = UNINIT_INT, dset_max_size(7) = UNINIT_INT
   contains
     ! public
     procedure, public  :: openFile
     procedure, public  :: closeFile
     procedure, public  :: queryDataset
     generic,   public  :: readDataset => readDataset1DReal, readDataset1DReal8, readDataset1DInt, readDataset1DChar24, readDataset2DReal
     ! private
     procedure, private :: readDataset1DReal
     procedure, private :: readDataset1DReal8
     procedure, private :: readDataset1DInt
     procedure, private :: readDataset1DChar24
     procedure, private :: readDataset2DReal
     procedure, private :: uninitDataset
  end type hdf5read

contains

  ! open file
  subroutine openFile(this, filename)

    ! input/output variables
    ! NEED class(hdf5read) instead of type(hdf5read)
    class (hdf5read), intent(inout) :: this
    character(len=*), intent(in) :: filename

    ! local variable
    integer :: hdf5err

    ! set obj param val
    this%file_name = filename

    ! initialize fortran interface
    call h5open_f(hdf5err)
    call checkErrCode_('h5open_f', hdf5err)

    ! open existing file
    call h5fopen_f(this%file_name, H5F_ACC_RDONLY_F, this%file_id, hdf5err)
    call checkErrCode_('h5fopen_f', hdf5err)

  end subroutine openFile

  ! close already opened file
  subroutine closeFile(this)

    ! input/output variables
    class (hdf5read), intent(inout) :: this

    ! local variable
    integer :: hdf5err

    ! ensure that dataset has been closed
    if (this%dset_name/=UNINIT_STR) stop "ERROR: Open dataset needs to be closed first. Stopping!"

    ! close file
    call h5fclose_f(this%file_id, hdf5err)
    call checkErrCode_('h5fclose_f', hdf5err)
    this%file_name = UNINIT_STR
    this%file_id = UNINIT_INT

    ! close fortran interface
    call h5close_f(hdf5err)
    call checkErrCode_('h5close_f', hdf5err)
    
  end subroutine closeFile

  ! query dataset for number of dims and its shape
  subroutine queryDataset(this, dsetName, dsetRank, dsetSize)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    character(len=*), intent(in) :: dsetName
    integer, intent(out) :: dsetRank
    integer, intent(out) :: dsetSize(7)

    ! local variable
    integer :: hdf5err

    ! ensure that file_name is set i.e. openFile
    ! must have been called prior to this routine
    if (this%file_name==UNINIT_STR) stop "ERROR: No open file available. Stopping!"

    ! set obj param val
    this%dset_name = dsetname

    ! open datset from already opened file
    call h5dopen_f(this%file_id, this%dset_name, this%dset_id, hdf5err)
    call checkErrCode_('h5dopen_f', hdf5err)

    ! get dataspace id
    call h5dget_space_f(this%dset_id, this%dspace_id, hdf5err)
    call checkErrCode_('h5dget_space_f', hdf5err)

    ! get num of dimensions
    call h5sget_simple_extent_ndims_f(this%dspace_id, this%dset_rank, hdf5err)
    call checkErrCode_('h5sget_simple_extent_ndims_f', hdf5err)
    dsetRank = this%dset_rank

    ! get size of array
    call h5sget_simple_extent_dims_f(this%dspace_id, this%dset_size, this%dset_max_size, hdf5err)
    call checkErrCode_('h5sget_simple_extent_dims_f', hdf5err)
    dsetSize = this%dset_size

  end subroutine queryDataset


  ! uninitalize dataset
  subroutine uninitDataset(this)

    ! input/output variables
    class (hdf5read), intent(inout) :: this

    ! un-initialize everything related to
    ! the dataset queried/read
    this%dset_name = UNINIT_STR
    this%dset_id = UNINIT_INT
    this%dspace_id = UNINIT_INT
    this%dset_rank = UNINIT_INT
    this%dset_size = UNINIT_INT
    this%dset_max_size = UNINIT_INT
    this%dtype_id = UNINIT_INT

  end subroutine uninitDataset


  ! read the dataset that was queried earlier
  subroutine readDataset1DChar24(this, dataChar)
    
    ! input/output variables
    class (hdf5read), intent(inout) :: this
    character(len=24), intent(out) :: dataChar(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, dataChar, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DChar24


  ! read the dataset that was queried earlier
  subroutine readDataset1DReal(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real, intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, data1D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DReal


  ! read the dataset that was queried earlier
  subroutine readDataset1DReal8(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real(8), intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, data1D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DReal8


  subroutine readDataset1DInt(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    integer, intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       !call h5dread_f(this%dset_id, this%dtype_id, data1D, this%dset_size, hdf5err)
       call h5dread_f(this%dset_id, H5T_NATIVE_INTEGER, data1D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DInt


  subroutine readDataset2DReal(this, data2D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real, intent(out) :: data2D(:,:)

    ! local variable
    integer :: hdf5err

    ! ensure that dset_name is set i.e. openDataset
    ! must have been called prior to this routine
    if (this%dset_name==UNINIT_STR) stop "ERROR: No open dataset available. Stopping!"

    if (this%dset_size(1)==0) then
       print *, 'Datset ', trim(this%dset_name), ' in file ', trim(this%file_name), ' is empty'
    else
       ! get data type
       call h5dget_type_f(this%dset_id, this%dtype_id, hdf5err)
       
       ! read data
       call h5dread_f(this%dset_id, this%dtype_id, data2D, this%dset_size, hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! close dataset
    call h5dclose_f(this%dset_id, hdf5err)
    call checkErrCode_('h5dclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset2DReal

  ! check return code
  ! (not part of class hdf5read)
  subroutine checkErrCode_(routineName, hdf5errCode)

    ! input/output variables
    character(len=*), intent(in) :: routineName
    integer, intent(in) :: hdf5errCode

    if (hdf5errCode<0) then
       write(*,*) 'ERROR: ', routineName, ' returned NEGATIVE err code. Stopping!'
       stop
    end if

  end subroutine checkErrCode_

end module GEOSland_io_hdf5

 
