#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

! This is good to force the river routing model offline using MERRA-2/SMAP-L4 time-averaged 
! Total runoff ( Surface RUNOFF + BASEFLOW) diagnostics fields.

module route_module

  implicit none
  private
  
  public :: tile_coord_type, optimize_latlon, PFAF_FILE, is_indom
  character*300, parameter :: PFAF_FILE = '/discover/nobackup/rreichle/l_data/LandBCs_files_for_mkCatchParam/V001/SRTM-TopoData/Pfafcatch-routing.dat'
  
  type :: tile_coord_type
          
     integer :: tile_id    ! unique tile ID
     integer :: pfaf       ! Pfafstetter number (for land tiles, NOT unique)
     real    :: com_lon    ! center-of-mass longitude
     real    :: com_lat    ! center-of-mass latitude
     integer :: i_indg     ! i index (w.r.t. *global* grid that cuts tiles) 
     integer :: j_indg     ! j index (w.r.t. *global* grid that cuts tiles)
     
  end type tile_coord_type
   
  
contains
  
  ! ------------------------------------------------------------------
  
  subroutine optimize_latlon(TILFILE, N_proc, my_pid, domain_lonlat, tile_coord)
    implicit none
    character(*), intent(in) :: TILFILE
    integer, intent (in)     :: N_proc, my_pid
    real, dimension (4), intent (in) :: domain_lonlat
    type(tile_coord_type), dimension(:), pointer, intent (inout) :: tile_coord
    integer :: N_tile,N_lon,N_lat,N_grid
    integer,allocatable :: landPosition(:)
    integer,allocatable :: IMS(:),JMS(:)
    integer,allocatable :: local_land(:)
    integer :: total_land
    integer :: n,typ,tmpint
    real ::  tmpreal
    integer :: avg_land,n0,local
    integer :: i,s,e,j,k,n1,n2
    logical :: file_exist
    character(len=100):: tmpLine
    character(len=100):: gridname
    real :: rate,rates(60),maxf(60)
    integer :: IMGLOB, JMGLOB
    integer :: face(6),face_land(6)
    integer :: tile_count, tmp_pfaf, n_local
    real    :: tmp_lon, tmp_lat
    integer, allocatable, dimension (:) :: tileid, i_ind, j_ind, pfaf_id
    integer, allocatable, dimension(:,:):: lon_tileid, proc_tileid
    real,    allocatable, dimension (:) :: lon, lat
    
    inquire(file=trim(TILFILE),exist=file_exist)
    if( .not. file_exist) then
       print *, trim (TILFILE)
       stop ( "tile file not exist")
    endif
    !     read (arg,*) N_proc
    
    open (10, file=trim(TILFILE), form='formatted', action='read')
    read (10,*) N_tile
    read (10,*) N_grid         ! some number (?)
    read (10,*) gridname       ! some string describing tile definition grid (?)
    read (10,*) N_lon
    read (10,*) n_lat
    
    allocate(IMS(N_Proc))
    allocate(local_land(N_Proc))
    allocate(tileid    (N_TILE))
    allocate(i_ind     (N_TILE))
    allocate(j_ind     (N_TILE))
    allocate(pfaf_id   (N_TILE))
    allocate(lon       (N_TILE))
    allocate(lat       (N_TILE))
    
    IMS=0
    local_land = 0
    tile_count = 0
    
    if(N_grid == 2) then
       read (10,*)          ! some string describing ocean grid                   (?)
       read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
       read (10,*)   
       read(10,'(A)') tmpLine
    else
       read(10,'(A)') tmpLine
       if (index(tmpLine,"OCEAN") /=0) then
          read (10,*)          ! # ocean grid cells in longitude direction (N_i_ocn) (?)
          read (10,*)   
          read(10,'(A)') tmpLine
       endif
    endif
    
    if (index(gridname,'EASE') /=0) then
       s=0
       e=N_lon-1
    else
       s=1
       e=N_lon
    endif
    
    allocate(landPosition(s:e))
    allocate(lon_tileid  (s:e, N_TILE))
    allocate(proc_tileid (N_PROC, 2*N_TILE/N_PROC))
    
    landPosition= 0
    total_land  = 0
    lon_tileid  = 0
    proc_tileid = 0
    
    ! 1) read through tile file, put the land tile into the N_lon of bucket
    
    read (tmpLine,*)  &
         typ,         &   !  1
         tmp_pfaf,    &   !  2  *
         tmp_lon,     &   !  3
         tmp_lat,     &   !  4
         i,j  !,             &   !  5
    
    if(typ==100) then
       tile_count = tile_count + 1
       tileid (tile_count) = tile_count
       pfaf_id(tile_count) = tmp_pfaf
       i_ind  (tile_count) = i
       j_ind  (tile_count) = j
       lon    (tile_count) = tmp_lon
       lat    (tile_count) = tmp_lat
       if (is_indom (domain_lonlat,tmp_lon,tmp_lat)) then
          total_land=total_land+1
          landPosition(i) = landPosition(i)+1
          lon_tileid (i,landPosition(i)) = tile_count
       endif
    endif
    
    do n = 2,N_tile
       read (10,*)       &
            typ,         &   !  1
            tmp_pfaf,    &   !  2  *
            tmp_lon,     &   !  3
            tmp_lat,     &   !  4
            i,j !
       !tmpint,        &   !  6
       !tmpreal,       &   !  7
       !tmpint,        &   !  8
       !tmpreal,       &   !  9  *
       !tmpint,        &   ! 10
       !tmpreal,       &   ! 11
       !tmpint       ! 12  * (previously "tile_id")
       if(typ==100) then
          tile_count = tile_count + 1
          tileid (tile_count) = tile_count
          pfaf_id(tile_count) = tmp_pfaf
          i_ind  (tile_count) = i
          j_ind  (tile_count) = j
          lon    (tile_count) = tmp_lon
          lat    (tile_count) = tmp_lat
          if (is_indom (domain_lonlat,tmp_lon,tmp_lat)) then
             total_land=total_land+1
             landPosition(i) = landPosition(i)+1
             lon_tileid (i,landPosition(i)) = tile_count   
          endif
       endif
       ! assume all lands are at the beginning 
       if(typ /= 100 .and. typ /=1100 ) exit
    enddo
    
    close(10)
    
    if(sum(landPosition) /= total_land) print*, "wrong counting of land"
    
    do n=1,60
       rates(n) = (n-1)*0.15
    enddo
    
    maxf=rms(rates)
    n=minloc(maxf,DIM=1)
    rate = rates(n)
    
    ! 2) each process should have average land tiles
    
    avg_land = ceiling(1.0*total_land/N_proc)
    ! print*,"avg_land",avg_land
    
    ! rate is used to readjust the avg_land
    ! in case that the last processors don't have any land tiles,
    ! we can increase ther rates
    
    avg_land = avg_land - nint(rate*avg_land/N_proc)
    !      print*,"re adjust the avg_land",avg_land
    tmpint = 0
    local  = 1
    n0     = s-1
    n_local= 1
    
    do n=s,e
       if(landPosition(n) > 0) then
          proc_tileid(local,n_local : n_local + landPosition(n) -1) = lon_tileid (n,1:landPosition(n))
          n_local = n_local + landPosition(n)
       endif
       tmpint=tmpint+landPosition(n)
       if(local == N_proc .and. n < e) cycle ! all lefteover goes to the last process
       if((tmpint .ge. avg_land) .or. (n==e)) then
          local_land(local)=tmpint
          IMS(local)=n-n0
          tmpint=0
          n0=n
          local = local + 1
          n_local = 1
       endif
    enddo
    
    n_local = count (mask = (proc_tileid (my_pid,:)> 0))
    allocate (tile_coord (1: n_local))
    tile_coord%tile_id = proc_tileid (my_pid,1:n_local)
    tile_coord%pfaf    = pfaf_id (tile_coord%tile_id)
    tile_coord%com_lon = lon  (tile_coord%tile_id)
    tile_coord%com_lat = lat  (tile_coord%tile_id)
    tile_coord%i_indg  = i_ind(tile_coord%tile_id)
    tile_coord%j_indg  = j_ind(tile_coord%tile_id)
    
    deallocate (tileid,i_ind,j_ind,pfaf_id,lon,lat)
                
  contains 
        
    elemental function rms(rates) result (f)
      real :: f
      real,intent(in) :: rates
      integer :: tmpint,local
      integer :: n0,proc,n
      integer :: avg_land
      integer,allocatable :: local_land(:)
      
      allocate (local_land(N_proc))
      local_land = 0
      avg_land = ceiling(1.0*total_land/N_proc)
      avg_land = avg_land -nint(rates*avg_land/N_proc)
      
      tmpint = 0
      local = 1
      n0 = s-1
      do n=s,e
         tmpint=tmpint+landPosition(n)
         if(local == N_proc .and. n < e) cycle ! all lefteover goes to the last process
         if((tmpint .ge. avg_land) .or. (n==e)) then
            local_land(local)=tmpint
            tmpint=0
            n0=n
            local = local + 1
         endif
      enddo
      f = 0.0
      do proc = 1, N_proc
         f =max(f,1.0*abs(local_land(proc)-avg_land))
      enddo
      deallocate(local_land)
    end function rms
    
    elemental function rms_cs(rates) result (f)
      real :: f
      real,intent(in) :: rates
      integer :: tmpint,local
      integer :: proc,n
      integer :: avg_land
      integer,allocatable :: local_land(:)
      integer :: n1,n2
      
      allocate (local_land(face(k)))
      local_land = 0
      avg_land = ceiling(1.0*face_land(k)/face(k))
      avg_land = avg_land -nint(rates*avg_land/face(k))
      if (avg_land <=0) then
         f = face_land(k)
         return
      endif
      
      tmpint = 0
      local = 1
      
      n1 = (k-1)*IMGLOB+1
      n2 = k*IMGLOB
      tmpint = 0
      do n = n1,n2
         tmpint=tmpint+landPosition(n)
         if(local == face(k) .and. n < n2) cycle ! all lefteover goes to the last process
         if((tmpint .ge. avg_land) .or. (n==n2)) then
            local_land(local)= tmpint
            tmpint=0
            local = local + 1
         endif
      enddo
      
      f = 0.0
      do proc = 1, face(k)
         ! punish for no land tiles
         f =max(f,1.0*abs(local_land(proc)-avg_land))
      enddo
      deallocate(local_land)
    end function rms_cs
  end subroutine optimize_latlon


  logical function is_indom (domain_box,x0,y0) result (f)
    real, intent (in) :: domain_box (4),x0,y0
    real :: x1,x2,y1,y2
    
    x1 = domain_box (1)
    x2 = domain_box (3)
    y1 = domain_box (2)
    y2 = domain_box (4)
    
    f = .false.
    if (((x0 >= x1).and.(x0 <= x2)) .and. ((y0 >= y1).and.(y0 <= y2))) &
         f = .true.
    
  end function is_indom

  
end module route_module

! --------------------------------------------------------------------------------------

module io_hdf5

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

end module io_hdf5

! ======================================================================
!                                MAIN DRIVER
! ======================================================================

PROGRAM GEOSroute_driver

  use io_hdf5,                          ONLY:     &
       hdf5read 
  use routing_model,                    ONLY:     &
       river_routing, SEARCH_DNST, RRM_TIMESTEP
  USE catch_constants, ONLY:          &
       N_Pfaf_Catchs, N_Pfaf_LandCatchs
  use ESMF
  use MAPL
  use route_module, ONLY: tile_coord_type, optimize_latlon, PFAF_FILE, is_indom     
  
  implicit none

  INCLUDE 'netcdf.inc'
  INCLUDE 'mpif.h'

  TYPE :: route_progn_type

     REAL    :: WSTREAM        ! AMOUNT OF WATER IN "LOCAL STREAM"
     REAL    :: WRIVER         ! AMOUNT OF WATER IN RIVER
     
  END TYPE route_progn_type

  TYPE :: route_param_type
     
     REAL    :: AREACAT     ! AREA OF CATCHMENT (km^2)
     REAL    :: LENGSC      ! LENGTHSCALE OF CATCHMENT FOR UPSTREAM-TO-DOWNSTREAM CALC (km)  
     INTEGER :: DNSTR       ! CATCHMENTE INDEX OF THE DOWNTREAM CATCHMENT 

  END TYPE route_param_type

  TYPE :: route_diagn_type

     REAL    :: QSFLOW         ! TRANSFER OF MOISTURE FROM STREAM VARIABLE TO RIVER VARIABLE
     REAL    :: QINFLOW        ! TRANSFER OF RIVER WATER FROM UPSTREAM CATCHMENTS
     REAL    :: QOUTFLOW       ! TRANSFER OF RIVER WATER TO DOWNSTREAM CATCHMENT
     REAL    :: RUNOFF         ! Input runoff from the land surface model (m3/s)
     REAL    :: BFLOW

  END TYPE route_diagn_type
  
  integer                 :: comm_rank, comm_size, status, mpistatus(MPI_STATUS_SIZE)  
  logical                 :: root_proc=.true.
  integer, parameter      :: logunit = 999
  type(ESMF_VM)           :: VM
  type(ESMF_Time)         :: CURRENT_TIME,FORCE_RING, DAY_RING,MON_RING, FIRST_TIME, END_TIME
  type(ESMF_TimeInterval) :: MODEL_DT, FORCE_DT, DAY_DT, MON_DT
  type(ESMF_Calendar)     :: gregorianCalendar
  integer                 :: FORCE_TIMESTEP, FIRST_DATE, LAST_DATE, &
       YY, MO, DA, HO, NC, NR
  real, dimension(4)      :: domain_lonlat
  character*20            :: DATA_SOURCE,data_label
  character*300           :: TILFILE, BCS_PATH, filename, WORK_PATH, MET_PATH, RUNOFF_FILE
  REAL,    DIMENSION(:,:), POINTER         :: CAT_AREA, SRUNOFF, BASEFLOW
  INTEGER, DIMENSION(:,:), POINTER         :: CAT_INDEX, AllActive,dstCatchID 
  INTEGER, DIMENSION(:)  , POINTER         :: NCats_in_GRID
  INTEGER, DIMENSION(:), ALLOCATABLE       :: LocalActive, srcProcsID, LocDstCatchID 
  integer, dimension (:),allocatable       :: MyActive, DomActive
  REAL,    DIMENSION(:), ALLOCATABLE       :: RUNOFF_DOMAIN, runoff_tile
  REAL,    DIMENSION(:), ALLOCATABLE       :: AREACAT, LENGSC, WSTREAM, WRIVER, QSFLOW, QOUTFLOW, RUNOFF
  integer   :: INFO, req, mpierr  
  type(tile_coord_type), dimension(:), pointer :: tile_coord => null()
  INTEGER                                  :: Local_Min, Local_Max, Pfaf_Min, Pfaf_Max,   &
       N_PfafD,N_PfafL, NCFOutID, N_Active, ThisCycle, N_catl, i, N, FirstCatch,Check_Min,&
       Check_Max, FileID , k,dset_rank,  dset_size(7),fcst_day
  logical                                  :: file_exists, ease_grid
  integer, dimension (:), allocatable      :: tmp_index
  CHARACTER*10                             :: YYYYMMDDHH
  CHARACTER*8                              :: YYYYMMDD
  CHARACTER*6                              :: YYYYMM, smap_label
  CHARACTER*4                              :: YYYY, HHMM
  CHARACTER*3                              :: SER
  CHARACTER*2                              :: MM, DD
  character*40                             :: MYNAME
  real                                     :: tmp_real (5), rbuff
  type(hdf5read)                           :: h5r
  
  ! These are local to the processor

  TYPE (route_progn_type), DIMENSION (:), allocatable :: rt_progn
  TYPE (route_param_type), DIMENSION (:), allocatable :: rt_param
  TYPE (route_diagn_type), DIMENSION (:), allocatable :: rt_diagn_daily
  
  ! Allocated only inside root_proc for the full domain
  
  TYPE (route_diagn_type), DIMENSION (:), allocatable :: rt_diagn
   
  call MPI_Init(status)                                ; VERIFY_(STATUS)
  call MPI_COMM_Size(MPI_COMM_WORLD,comm_size,status)  ; VERIFY_(STATUS)
  call MPI_COMM_Rank(MPI_COMM_WORLD,comm_rank,status)  ; VERIFY_(STATUS)
  if (comm_rank /= 0) root_proc = .false.

  call ESMF_Initialize (vm=vm, logKindFlag=ESMF_LOGKIND_NONE, rc=status) ; VERIFY_(STATUS)
  call ESMF_CalendarSetDefault ( ESMF_CALKIND_GREGORIAN) 

  ! Process input data
  
  call getenv ("MYNAME"        ,MYNAME        )
  call Read_Resource (MET_PATH, WORK_PATH, domain_lonlat,FIRST_DATE, LAST_DATE)

  write (YYYYMM,'(i6.6)')FIRST_DATE/100
  write (YYYY  ,'(i4.4)')FIRST_DATE/10000
  MM = YYYYMM(5:6)
  
  if (root_proc)  open (logunit, file = trim(WORK_PATH)//'river_routing.log', form ='formatted',status='unknown', action = 'write')
  
  if (index(met_path, 'merra2')  /= 0) then
     BCS_PATH       = '/discover/nobackup/smahanam/bcs/Icarus-NL/Icarus-NL_Reynolds/DC0576xPC0361_MERRA2-GRID/'
     TILFILE        = trim(BCS_PATH)//'til/DC_00576x00360_PC_0576x0360.til'
     FORCE_TIMESTEP = 1
     NC             = 576 
     NR             = 361
     data_label = 'MERRA2'
     data_source= 'MERRA2'     
     if (root_proc) write (logunit,*)'Total Runoff source  : MERRA2'
  endif

  if (index(met_path, 'L4_SM') /= 0) then
     BCS_PATH       = '/discover/nobackup/ltakacs/bcs/Icarus-NLv3/Icarus-NLv3_EASE/SMAP_EASEv2_M09/'
     data_label = 'SMAP_L4_SM_gph'
     data_source= 'SMAP_L4'
     if(index(met_path, 'Vv4030') /=0) smap_label = 'Vv4030'
     if(index(met_path, 'OL4001') /=0) smap_label = 'OL4001'
     TILFILE        = trim(BCS_PATH)//'SMAP_EASEv2_M09_3856x1624.til'     
     FORCE_TIMESTEP = 3
     NC             = 3856
     NR             = 1624
     if (root_proc) write (logunit,*)'Total Runoff source  : SMAP_L4/'//smap_label
  endif

  if (index(met_path, 'GEOS5_s2s')     /= 0) then
     BCS_PATH       = '/discover/nobackup/smahanam/bcs/Icarus-NL/Icarus-NL_Reynolds/DC0720xPC0361_FCAST-GRID/'
     i = index(MET_PATH,'/',back=.true.)
     data_label = trim(MET_PATH(i+1:i+5))
     data_source= 'S2S'
     TILFILE        = trim(BCS_PATH)//'til/DC_00720x00360_PC_0720x0360.til'     
     FORCE_TIMESTEP = 24
     NC             = 720
     NR             = 361
     if (root_proc) write (logunit,*)'Total Runoff source  : GEOS5 '//trim(data_source)
  endif

  if (root_proc) then
     write (logunit,*) 'ROUTING MODEL SIMULATION BY : ',trim (MYNAME)
     write (logunit,*) 'RUNOFF DATA SOURCE          : ',trim (DATA_SOURCE)
     write (logunit,*) ' '
  endif

    ! Domain decomposition

  call optimize_latlon(TILFILE,comm_size,comm_rank+1, domain_lonlat, tile_coord)
  N_catl = size (tile_coord%tile_id)

  if (index(met_path, 'L4_SM') /= 0) then
     do n = 1, N_catl
        tile_coord(n)%i_indg = tile_coord(n)%i_indg + 1
        tile_coord(n)%j_indg = tile_coord(n)%j_indg + 1
     end do
  endif
  
  ! read Grid2Catch

  call MPI_Info_create(info, STATUS)
  call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS)
    
  filename = trim(BCS_PATH)//'/clsm/Grid2Catch_TransferData.nc'
  if (root_proc) write (logunit,*) 'reading '//trim(filename)// ' and setting up the domain'
       
  call Read_Grid2Catch (N_catl, INFO, tile_coord, filename,  &
       Local_Min, Local_Max, CAT_AREA, CAT_INDEX, NCats_in_GRID)
       
  ! compute min and max Pfafstetter catchment indices across processors within the domain
    
  call MPI_Reduce(Local_Min,Pfaf_Min,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_WORLD,mpierr)
  call MPI_Reduce(Local_Max,Pfaf_Max,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
    
  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    
  call MPI_BCAST (Pfaf_Min , 1, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
  call MPI_BCAST (Pfaf_Max , 1, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
    
  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    
  N_PfafL = Local_Max - Local_Min + 1
  N_PfafD = Pfaf_Max  - Pfaf_Min  + 1

  allocate (rt_progn       (1:N_PfafL))
  allocate (rt_param       (1:N_PfafL))
  allocate (rt_diagn       (1:N_PfafD))
  allocate (rt_diagn_daily (1:N_PfafD))
  
  RT_DIAGN_DAILY%QSFLOW   = 0. 
  RT_DIAGN_DAILY%QINFLOW  = 0. 
  RT_DIAGN_DAILY%QOUTFLOW = 0. 
  RT_DIAGN_DAILY%RUNOFF   = 0.

  ! Reading restarts
    
  inquire(file=TRIM(work_path)//'/Pfaf_domain.dat', exist=file_exists)
    
  if (file_exists) then
     FirstCatch = Local_Min - Pfaf_Min + 1
     filename = TRIM(work_path)//'/rs/'//YYYY//'/route_internal_rst.'//YYYYMM//'01'
     call Read_Restarts (FirstCatch, N_PfafL, info, RT_PARAM, RT_PROGN, filename, .false.)
  else
     FirstCatch =  Local_Min 
     filename = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/route/' &
          //'route_internal_rst.YYYY'//YYYYMM(5:6)//'01'
     call Read_Restarts (FirstCatch, N_PfafL, info, RT_PARAM, RT_PROGN, filename, .true.)
  endif

  ! Pfafstetter catchment Domain Decomposition :         
  ! --------------------------------------------
  
  ! AllActive      : Processor(s) where the catchment is active  (identical in any processor). 
  ! srcProcsID     : For all active catchments anywhere tells which processor is the principal owner of the catchment (identical in any processor).
  ! DstCatchID     : 2-D array contains downstream catchID and downstream processor (identical in any processor)
  ! LocDstCatchID  : Downstream catchID when for catchments that are local to the processor.
  ! LocalActive    : This is local to the processor, tells whether the catchment is active.
  
  allocate (AllActive    (1:N_PfafD, 1: comm_size))
  allocate (dstCatchID   (1:N_PfafD, 1: comm_size)) 
  allocate (LocalActive  (1:N_PfafD))
  allocate (srcProcsID   (1:N_PfafD))
  allocate (LocDstCatchID(1:N_PfafD))
  
  AllActive    = -9999
  LocalActive  = -9999
  srcProcsID   = -9999
  dstCatchID   = -9999
  LocDstCatchID= -9999
  
  LocDstCatchID (Local_Min - Pfaf_Min + 1 :Local_Max - Pfaf_Min + 1)  = RT_PARAM%DNSTR
  
  CALL DOMAIN_DECOMP (N_PfafD, Local_Min, Local_Max,                 &
       Pfaf_Min, Pfaf_Max, NCats_in_GRID, CAT_INDEX, LocalActive,    &
       LocDstCatchID, srcProcsID,dstCatchID, AllActive)
  
  N_Active = count (LocalActive > 0)
  
  allocate (MyActive (1 : N_Active))
  allocate (DomActive(1 : N_Active))
  allocate (tmp_index(1 : N_PfafD ))
  
  forall (N=1:N_PfafD) tmp_index(N) = N
  
  DomActive = pack (tmp_index, mask = (LocalActive > 0))
  MyActive  = DomActive - (Local_Min - Pfaf_Min)
  deallocate (tmp_index)
  
  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  
  HAVERESTARTS : IF (file_exists) THEN
     
     if (root_proc) then
        ! restarts files are available from LDASsa spin up
          
        OPEN (10, file = TRIM(work_path)//'/Pfaf_domain.dat', form = 'formatted', &
             status= 'old', action = 'read')
        
        read (10,*) Check_Min, Check_Max
        
        CLOSE (10, STATUS = 'KEEP')
        
        IF((Check_Min /= Pfaf_Min).OR.(Check_Max /= Pfaf_Max)) then
           
           print *, 'Pfaf_Min Pfaf_Max mismatch !',Check_Min, Check_Max,Pfaf_Min, Pfaf_Max
           stop
           
        ENDIF
               
        write (logunit,*) 'reading exisiting restarts for PfafIndex range : ', Pfaf_Min, Pfaf_Max
        write (logunit,*) 'restart file   : ', trim(filename)
                    
     endif
       
  ELSE 
       
     if (root_proc) then          
        write (logunit,*) '# of Pfaf catchs and the catchment index range :',N_PfafD, Pfaf_Min, Pfaf_Max
        write (logunit,*)'   '
        write (logunit,*) 'Reading global cold restarts ', trim(filename)
          
        OPEN (10, file = TRIM(work_path)//'/Pfaf_domain.dat', form = 'formatted', &
             status= 'unknown', action = 'write')
        
        write (10,*) Pfaf_Min, Pfaf_Max
        
        ! mask out catchments that lie outside of the simulation window.
        
        DO I = Pfaf_Min,Pfaf_Max
           WRITE (10,*) I 
        END DO
        
        CLOSE (10, STATUS = 'KEEP')
     endif
     
     call Write_Restarts (N_PfafL, N_PfafD, Local_Min, Pfaf_Min,  &
          srcProcsID, rt_progn, rt_param, work_path, YYYYMM)  
     
  endif HAVERESTARTS
 
  ! runoff input file

  allocate (SRUNOFF  (1:NC, 1: NR))
  allocate (BASEFLOW (1:NC, 1: NR))
  ALLOCATE (RUNOFF_DOMAIN(1:N_PfafD))
  ALLOCATE (RUNOFF_TILE(1:N_CATL))
        
  ! -----------  
  !  Time loop
  ! -----------

  YY = FIRST_DATE / 10000
  MO = (FIRST_DATE - YY*10000) / 100
  DA = FIRST_DATE - (YY*10000 + MO*100)
  call ESMF_TimeSet (FIRST_TIME, yy=YY, mm=MO, dd=DA, rc=status)  ; VERIFY_(STATUS)
    
  YY = LAST_DATE / 10000
  MO = (LAST_DATE - YY*10000) / 100
  DA = LAST_DATE - (YY*10000 + MO*100)
  call ESMF_TimeSet (END_TIME, yy=YY, mm=MO, dd=DA, rc=status)      ; VERIFY_(STATUS)
   
  call ESMF_TimeIntervalSet(MODEL_DT, s=RRM_TIMESTEP,   rc=status ) ; VERIFY_(STATUS)
  call ESMF_TimeIntervalSet(FORCE_DT, h=FORCE_TIMESTEP, rc=status ) ; VERIFY_(STATUS)
  call ESMF_TimeIntervalSet(DAY_DT,   h=24,             rc=status ) ; VERIFY_(STATUS)

  CURRENT_TIME = FIRST_TIME
  DAY_RING     = CURRENT_TIME
  MON_RING     = CURRENT_TIME
  FORCE_RING   = CURRENT_TIME

  if(trim(data_source) == 'S2S') then
     STATUS      = NF_OPEN_PAR (trim(MET_PATH),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,FileID) ; VERIFY_(STATUS)
     fcst_day = 1
     if (root_proc) write (logunit,*)trim(MET_PATH)
  endif
  
  do while (CURRENT_TIME  <= END_TIME)

     NEW_DAY : if(CURRENT_TIME == DAY_RING) then

        ! 1) write out daily averages in monthly output file
        if (CURRENT_TIME /= FIRST_TIME) then
           DAY_RING = CURRENT_TIME - DAY_DT
           call ESMF_TimeGet (DAY_RING, yy=YY, mm=MO, dd=DA, rc=status)      ; VERIFY_(STATUS)
           filename = TRIM(work_path)//'/output/'//YYYY//'/route_daily_avg_'//YYYYMM//'.nc4'
           if (root_proc) write (logunit,*) 'writing output file   : ', trim(filename)
           I = 10000 * YY + 100 * MO + DA
           call write_output (NCFOutID,N_PfafD,I, DA,srcProcsID, RT_DIAGN_DAILY)           
        endif
        
        ! 2) open MERRA2 daily file
        if(trim(data_source) == 'MERRA2') then 
           call ESMF_TimeGet (CURRENT_TIME, yy=YY, mm=MO, dd=DA, rc=status)  ; VERIFY_(STATUS)
           write (YYYYMMDD,    '(i4.4,i2.2,i2.2)') YY, MO, DA
           write (YYYY    ,    '(i4.4)')           YY
           write (MM      ,    '(i2.2)')           MO
           write (DD      ,    '(i2.2)')           DA
           ! hourly data in Daily files 
           SER = '100'
           if (YY > 1991) SER = '200'
           if (YY > 2000) SER = '300'
           if (YY > 2010) SER = '400'
           
           STATUS      = NF_CLOSE (FileID)
           if (root_proc) write (logunit,*) '/MERRA2_'//SER//'/Y'//YYYY//'/M'//MM//'/MERRA2_'//SER//'.tavg1_2d_lnd_Nx.'//YYYYMMDD//'.nc4'
           RUNOFF_FILE = trim(MET_PATH)//'/data/products//MERRA2_'//SER//'/Y'//YYYY//'/M'//MM//'/MERRA2_'//SER//'.tavg1_2d_lnd_Nx.'//YYYYMMDD//'.nc4'                
           STATUS      = NF_OPEN_PAR (trim(RUNOFF_FILE),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,FileID) ; VERIFY_(STATUS)  
           
        endif
        
        ! 3) reset daily alarm
        DAY_RING = CURRENT_TIME + DAY_DT

     endif NEW_DAY
    
     NEW_MONTH : if(CURRENT_TIME == MON_RING) then        

        call ESMF_TimeGet (CURRENT_TIME, yy=YY, mm=MO, dd=DA, rc=status)  ; VERIFY_(STATUS)
        write (YYYYMM,    '(i4.4,i2.2)') YY, MO
        write (YYYY,      '(i4.4)') YY
        
        ! 1) write restarts
        if (CURRENT_TIME /= FIRST_TIME)                                   &
             call Write_Restarts (N_PfafL, N_PfafD, Local_Min, Pfaf_Min,  &
             srcProcsID, rt_progn, rt_param, work_path, YYYYMM)
        
        ! 2) open monthly output file
        if (root_proc) then
           STATUS      = NF_CLOSE (NCFOutID) 
           filename = TRIM(work_path)//'/output/'//YYYY//'/route_daily_avg_'//YYYYMM//'.nc4'
           call open_monthly_outfile (NCFOutID, N_PfafD, filename, MyName)
        endif
        
        ! 3) reset monthly alarm
        if(MO == 12) then
           call ESMF_TimeSet (MON_RING,     yy=YY+1, mm=1,dd=1, rc=status)  ; VERIFY_(STATUS)
        else
           call ESMF_TimeSet (MON_RING,     yy=YY, mm=MO+1,dd=1, rc=status)  ; VERIFY_(STATUS)
        endif        

     endif NEW_MONTH
     
     READ_FORCING : if(CURRENT_TIME == FORCE_RING) then
        
        call ESMF_TimeGet (CURRENT_TIME, yy=YY, mm=MO, dd=DA, H= HO, rc=status)  ; VERIFY_(STATUS)
        select case ( trim(data_source)) 
        case ('MERRA2')
           ! MERRA 2 read from the opened daily file
           STATUS = NF_GET_VARA_REAL(FileID, VarID(FileID,'BASEFLOW' ), (/1,1,HO +1/), (/NC,NR,1/),BASEFLOW) ; VERIFY_(STATUS)
           STATUS = NF_GET_VARA_REAL(FileID, VarID(FileID,'RUNOFF'   ), (/1,1,HO +1/), (/NC,NR,1/),SRUNOFF ) ; VERIFY_(STATUS)

        case ('SMAP_L4')
           write (YYYYMMDD,    '(i4.4,i2.2,i2.2)') YY, MO, DA
           write (YYYY    ,    '(i4.4)')           YY
           write (MM      ,    '(i2.2)')           MO
           write (DD      ,    '(i2.2)')           DA
           write (HHMM , '(i4.4)') (HO+1)*100 + 30                    
           RUNOFF_FILE = trim(MET_PATH)//'/Y'//YYYY//'/M'//MM//'/D'//DD//'/'//trim(data_label)//   &
                '_'//YYYYMMDD//'T'//HHMM//'00_'//smap_label//'_001.h5'

           call h5r%openFile(RUNOFF_FILE)
           call h5r%queryDataset('/Geophysical_Data/baseflow_flux', dset_rank, dset_size)
           call h5r%readDataset(baseflow)
           call h5r%queryDataset('/Geophysical_Data/overland_runoff_flux', dset_rank, dset_size)
           call h5r%readDataset(srunoff)
           call h5r%closeFile
           if (root_proc) write (logunit,*) '/Y'//YYYY//'/M'//MM//'/D'//DD//'/'//trim(data_label)//   &
                '_'//YYYYMMDD//'T'//HHMM//'00_'//smap_label//'_001.h5'
           
        case ('S2S')
           STATUS = NF_GET_VARA_REAL(FileID, VarID(FileID,'RUNOFF'   ), (/1,1,fcst_day/), (/NC,NR,1/),SRUNOFF ) ; VERIFY_(STATUS)
           BASEFLOW = 0.
           if (root_proc) write (logunit,*) 'READING FORECAST TIME STEP : ',fcst_day
           fcst_day = fcst_day + 1
           
        case default
           print *, 'UNKNOWN DATA SOURCE : '//trim(data_source)
           stop
        end select
        
        ! compute total runoff at each grid cell within the local domain
        runoff_tile = 0. 
        DO N = 1, N_catl
           if((baseflow (tile_coord(n)%i_indg, tile_coord(n)%j_indg) >= 0.).AND.        &
                (baseflow (tile_coord(n)%i_indg, tile_coord(n)%j_indg) <= 1000.)) then 
              
              runoff_tile (n) =  baseflow (tile_coord(n)%i_indg, tile_coord(n)%j_indg) + &
                   srunoff  (tile_coord(n)%i_indg, tile_coord(n)%j_indg)              
           endif
        ENDDO
        
        ! 2) reset forcing alarm
        FORCE_RING = FORCE_RING + FORCE_DT
        
     endif READ_FORCING

     !-------------------------------------------------------------------
     !    Route tile-space runoff_tile through Pfafstetter watersheds
     !-------------------------------------------------------------------
     
     RUNOFF_DOMAIN = 0.
                   
     DO N = 1, N_catl
        if(NCats_in_GRID(N) > 0) then
           DO I = 1, NCats_in_GRID(N)                         
              K = CAT_INDEX(I,N) - Pfaf_Min + 1                         
              RUNOFF_DOMAIN (K) = RUNOFF_DOMAIN (K) + 1000. * runoff_tile(N)*Cat_Area(I,N) ! (m3/s)              
           END DO
        endif
     END DO
       
     ! For catchment-tiles that contribute to the main catchment in some other processor, send runoff to the corresponding srcProcsID(N)           
     do N = 1, N_PfafD                      
        if ((LocalActive (N) > 0).and.(srcProcsID(N) /= comm_rank)) then           
           ! Send this catchment's contribution to the owner                          
           tmp_real (1) = RUNOFF_DOMAIN (N)                         
           call MPI_ISend(tmp_real(1),1,MPI_real,srcProcsID(N),999,MPI_COMM_WORLD,req,status)
           call MPI_WAIT (req,MPI_STATUS_IGNORE,status)                         
           RUNOFF_DOMAIN (N) = 0.                         
        else                         
           if(srcProcsID(N) == comm_rank) then                            
              do i = 1,comm_size                               
                 if((i-1 /= comm_rank).and.(AllActive (N,i) > 0))  then                                  
                    call MPI_RECV(tmp_real(1),1,MPI_real,i-1,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)
                    RUNOFF_DOMAIN (N) = RUNOFF_DOMAIN (N) + tmp_real (1)                    
                 endif
              end do
           endif
        endif                      
     end do
                   
     call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
                   
     !  call river routing model
                   
     if(allocated(rt_diagn) .eqv. .false. ) allocate (rt_diagn (1:N_PfafD))
                   
     RT_DIAGN%QSFLOW   = 0. 
     RT_DIAGN%QINFLOW  = 0. 
     RT_DIAGN%QOUTFLOW = 0. 
     RT_DIAGN%RUNOFF   = 0. 
                   
     QSFLOW   = 0. 
     QOUTFLOW = 0. 
                   
     if(allocated (LENGSC  ) .eqv. .false.) allocate (LENGSC  (1:N_Active))
     if(allocated (AREACAT ) .eqv. .false.) allocate (AREACAT (1:N_Active))
     if(allocated (WSTREAM ) .eqv. .false.) allocate (WSTREAM (1:N_Active))
     if(allocated (WRIVER  ) .eqv. .false.) allocate (WRIVER  (1:N_Active))
     if(allocated (QSFLOW  ) .eqv. .false.) allocate (QSFLOW  (1:N_Active))
     if(allocated (QOUTFLOW) .eqv. .false.) allocate (QOUTFLOW(1:N_Active))  
     if(allocated (RUNOFF  ) .eqv. .false.) allocate (RUNOFF  (1:N_Active))  
     
     DO N = 1, size (DomActive)
        I = MyActive (N)
        WSTREAM (N) = RT_PROGN(I)%WSTREAM
        WRIVER  (N) = RT_PROGN(I)%WRIVER
        LENGSC(N)   = RT_PARAM(I)%LENGSC
        AREACAT(N)  = RT_PARAM(I)%AREACAT
        I = DomActive(N)
        RUNOFF (N)  = RUNOFF_DOMAIN(I)
     END DO
                   
     CALL RIVER_ROUTING  (N_Active,REAL(RRM_TIMESTEP),RUNOFF,AREACAT,LENGSC, &
          WSTREAM,WRIVER, QSFLOW,QOUTFLOW) 
     
     DO N = 1, size (DomActive)
        
        I = DomActive(N)
        RT_DIAGN(I)%QSFLOW   = QSFLOW  (N)
        RT_DIAGN(I)%QOUTFLOW = QOUTFLOW(N)
        RT_DIAGN(I)%RUNOFF   = RUNOFF  (N)
        
        I = MyActive (N)
        RT_PROGN(I)%WSTREAM  = WSTREAM (N) 
        RT_PROGN(I)%WRIVER   = WRIVER  (N) 
        
     END DO
                   
     ! update downstream catchments
                   
     do N = 1,N_PfafD 
        
        if ((srcProcsID (N) == comm_rank).and.(srcProcsID (LocDstCatchID (N)) == comm_rank)) then ! destination is local
           
           RT_DIAGN(LocDstCatchID (N))%QINFLOW = RT_DIAGN(LocDstCatchID (N))%QINFLOW + RT_DIAGN(N)%QOUTFLOW
           I = LocDstCatchID (N) - Local_Min + Pfaf_Min
           if(LocDstCatchID (N) /= N) RT_PROGN(I)%WRIVER  = RT_PROGN(I)%WRIVER + & 
                RT_DIAGN(N)%QOUTFLOW * real(RRM_TIMESTEP) ! LocDstCatchID (N) /= N : dont want to keep adding to outlet catchments
           
        elseif ((srcProcsID (N) == comm_rank).and.(srcProcsID (LocDstCatchID (N)) /= comm_rank)) then 
                         
           if(srcProcsID (LocDstCatchID (N)) >= 0) then
              
              ! Send to downstream processor
              
              call MPI_ISend(RT_DIAGN(N)%QOUTFLOW,1,MPI_real,srcProcsID (LocDstCatchID(N)),999,MPI_COMM_WORLD,req,status)
              call MPI_WAIT(req,MPI_STATUS_IGNORE,status) 
              
           endif
           
        elseif ((srcProcsID (N) /= comm_rank).and.(srcProcsID (N) >= 0)) then
           
           K = srcProcsID (dstCatchID(N,srcProcsID (N)+1))
           
           if (k == comm_rank) then
              
              do i = 1,comm_size
                 
                 if(comm_rank /= i-1) then 
                    if((srcProcsID  (n) == i-1).and.(srcProcsID (dstCatchID(N, i)) == comm_rank))then 
                       
                       call MPI_RECV(rbuff,1,MPI_real, srcProcsID (N),999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)
                       K = dstCatchID(N,i) - Local_Min + Pfaf_Min
                       RT_DIAGN(dstCatchID(N,i))%QINFLOW = RT_DIAGN(dstCatchID(N,i))%QINFLOW + rbuff
                       RT_PROGN(K)%WRIVER  = RT_PROGN(K)%WRIVER  + rbuff * real(RRM_TIMESTEP)
                       
                    endif
                 endif
              end do
           endif           
        endif        
     end do
                   
     call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) 
     
     ! compute daily averaged post processing
     
     DO N = 1, size (DomActive)
        
        I = DomActive(N)
        RT_DIAGN_DAILY(I)%QSFLOW   = RT_DIAGN_DAILY(I)%QSFLOW   + RT_DIAGN(I)%QSFLOW  * real (RRM_TIMESTEP) / 86400.
        RT_DIAGN_DAILY(I)%QINFLOW  = RT_DIAGN_DAILY(I)%QINFLOW  + RT_DIAGN(I)%QINFLOW * real (RRM_TIMESTEP) / 86400.
        RT_DIAGN_DAILY(I)%QOUTFLOW = RT_DIAGN_DAILY(I)%QOUTFLOW + RT_DIAGN(I)%QOUTFLOW* real (RRM_TIMESTEP) / 86400.
        RT_DIAGN_DAILY(I)%RUNOFF   = RT_DIAGN_DAILY(I)%RUNOFF   + RT_DIAGN(I)%RUNOFF  * real (RRM_TIMESTEP) / 86400.
        
     END DO

     ! next time step
          
     CURRENT_TIME = CURRENT_TIME + MODEL_DT

  end do
  
  contains

    ! ------------------------------------------------------------------

    SUBROUTINE Read_Grid2Catch (N_catl, INFO, tilecoord, filename,      &
         Local_Min, Local_Max, CAT_AREA, CAT_INDEX, NCats_in_GRID)
      
      implicit none
      type(tile_coord_type), dimension(:), INTENT (IN)  :: tilecoord
      character (*),                    INTENT (IN)  :: filename
      INTEGER,                          INTENT (IN)  :: N_catl, INFO
      REAL,    DIMENSION(:,:), POINTER, INTENT (OUT) :: CAT_AREA
      INTEGER, DIMENSION(:,:), POINTER, INTENT (OUT) :: CAT_INDEX
      INTEGER, DIMENSION(:)  , POINTER, INTENT (OUT) :: NCats_in_GRID
      REAL,    DIMENSION(:), ALLOCATABLE :: TMP_AREA, PF_LON_MN, PF_LON_MX, PF_LAT_MN, PF_LAT_MX
      INTEGER,  DIMENSION(:), ALLOCATABLE :: TMP_INDEX 
      INTEGER,                          INTENT (OUT) :: Local_Min, Local_Max    
      INTEGER                                        :: NCIG, STATUS, NCFID, VID, MAX_CAT_PER_CELL, N, J, K
      INTEGER*8  :: PFAF_CODE
      
      open (10, file = PFAF_FILE,  &
           form = 'FORMATTED', status= 'old', action = 'read')
      READ (10, *) K

      allocate (PF_LON_MN (1:K))
      allocate (PF_LON_MX (1:K))
      allocate (PF_LAT_MN (1:K))
      allocate (PF_LAT_MX (1:K))
      
      DO N = 1, K
         READ (10, *) J,PFAF_CODE, PF_LON_MN(N), PF_LON_MX(N), PF_LAT_MN(N), PF_LAT_MX(N)
      END DO
      close (10, status = 'keep')
      
      Local_Max   = 0
      Local_Min   = N_Pfaf_Catchs + 1   
      
      STATUS = NF_OPEN_PAR   (trim(filename),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,NCFID)                      
      STATUS = NF_INQ_DIMID (NCFID, 'MAX_CAT_PER_CELL', vid)
      STATUS = NF_INQ_DIMLEN(NCFID, vid, MAX_CAT_PER_CELL) 
      
      allocate (CAT_INDEX     (1:MAX_CAT_PER_CELL, N_catl))
      allocate (CAT_AREA      (1:MAX_CAT_PER_CELL, N_catl))
      allocate (NCats_in_GRID (1:N_catl))
      allocate (TMP_INDEX     (1:MAX_CAT_PER_CELL))
      allocate (TMP_AREA      (1:MAX_CAT_PER_CELL))
      
      CAT_INDEX = -9999
      CAT_AREA  = 0.
      TMP_INDEX = -9999
      TMP_AREA  = 0.
      NCats_in_GRID = 0
      
      DO N = 1, N_catl
         
         STATUS = NF_GET_VARA_INT (NCFID, VarID(NCFID,'NCats_in_GRID'), (/tilecoord(n)%tile_id/), (/1/), NCIG)
         STATUS = NF_GET_VARA_INT (NCFID, VarID(NCFID,'Pfaf_Index'),    (/1,tilecoord(n)%tile_id/),(/NCIG,1/), TMP_INDEX(1:NCIG))
         STATUS = NF_GET_VARA_REAL(NCFID, VarID(NCFID,'Pfaf_Area' ),    (/1,tilecoord(n)%tile_id/),(/NCIG,1/), TMP_Area (1:NCIG))

         k = 0
         DO J = 1, NCIG
            IF(TMP_INDEX(J) > 0) THEN
               if(is_indom ((/                 &
                    PF_LON_MN(TMP_INDEX(J)),PF_LAT_MN(TMP_INDEX(J)),   &
                    PF_LON_MX(TMP_INDEX(J)),PF_LAT_MX(TMP_INDEX(J))/), &
                    tilecoord(n)%com_lon,tilecoord(n)%com_lat)) then
                  k = k + 1
                  NCats_in_GRID(N) = k
                  CAT_INDEX(K, N) = TMP_INDEX(J)
                  CAT_AREA (K, N) = TMP_AREA(J)
               end if
            ENDIF
         END DO

         TMP_INDEX = -9999
         TMP_AREA  = 0.
         
         
         IF ((MINVAL (CAT_INDEX(1:NCats_in_GRID(N), N)) > 0 ).AND.(MINVAL (CAT_INDEX(1:NCats_in_GRID(N), N)) < Local_Min )) Local_Min = MINVAL (CAT_INDEX (1:NCats_in_GRID(N), N))
         IF ( MAXVAL (CAT_INDEX(1:NCats_in_GRID(N), N)) > Local_Max ) Local_Max = MAXVAL (CAT_INDEX (1:NCats_in_GRID(N), N))
         
      END DO
      
      STATUS   = NF_CLOSE (NCFID)
      
    END SUBROUTINE Read_Grid2Catch

    ! ------------------------------------------------------------------

    SUBROUTINE Read_Restarts ( &
         FirstCatch, N_PfafL, info, RT_PARAM, RT_PROGN, filename, bypass_submerged )
      
      implicit none
      
      logical,       intent (in)                                 :: bypass_submerged
      character (*), intent (in)                                 :: filename
      integer,       intent (in)                                 :: FirstCatch, N_PfafL, info
      TYPE (route_progn_type), DIMENSION (N_PfafL), intent (out) :: rt_progn
      TYPE (route_param_type), DIMENSION (N_PfafL), intent (out) :: rt_param
      integer                                                    :: STATUS, NCFID, N, Pfaf_land,K
      real                                                       :: VDUM
      real, allocatable, dimension (:)                           :: temp_var
      integer, allocatable, dimension (:)                        :: DNST,  DNST2, Pfaf_all
      
      
      allocate (temp_var (1:N_PfafL))
      
      STATUS = NF_OPEN_PAR   (trim(filename),IOR(NF_NOWRITE,NF_MPIIO),MPI_COMM_WORLD, info,NCFID)
      
      STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'AREACAT'   ), (/FirstCatch/), (/N_PfafL/),RT_PARAM(:)%AREACAT )
      STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'LENGSC'    ), (/FirstCatch/), (/N_PfafL/),RT_PARAM(:)%LENGSC  )
      STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'DNSTR'     ), (/FirstCatch/), (/N_PfafL/),temp_var            )
      STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WSTREAM'   ), (/FirstCatch/), (/N_PfafL/),RT_PROGN(:)%WSTREAM )
      STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WRIVER'    ), (/FirstCatch/), (/N_PfafL/),RT_PROGN(:)%WRIVER  )
      
      RT_PARAM(:)%DNSTR = NINT (temp_var(:))
      
      if (bypass_submerged) then
         
         ! Out of 291284 (N_Pfaf_Catchs in catch_constants.f90) only 290188 (N_Pfaf_LandCatchs in catch_constants.f90)
         ! are overland Pfafstetetter catchments. The rest is submerged under lakes or large river segments.
         ! Here, we check if the assigned downstream catchment is submerged and if so we bypass the submerged catchments(s)
         ! and link to the first available over land catchment downstream of the catchment in question.
         
         deallocate (temp_var)
         allocate   (temp_var (1:N_Pfaf_Catchs))
         allocate   (DNST (1:N_Pfaf_Catchs), DNST2 (1:N_Pfaf_Catchs), Pfaf_all (1:N_Pfaf_Catchs))
         STATUS  = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'DNSTR'     ), (/1/), (/N_Pfaf_Catchs/),temp_var      )
         DNST    = NINT (temp_var(:))
         
         Pfaf_all(:) = -9999
         DNST2       = -1
         
         ! Read land only catchment information
         
         OPEN (10, FILE = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/route/' &
              //'Pfafstetter.til', FORM='FORMATTED', STATUS = 'OLD', ACTION = 'READ')
         
         DO N = 1,5            
            READ (10,'(A)')            
         END DO
         
         DO N = 1, N_Pfaf_LandCatchs
            READ (10,*)K,VDUM,VDUM, VDUM, Pfaf_land, K
            Pfaf_all (Pfaf_land) = Pfaf_land
         END DO
         
         close (10, status = 'keep')
         
         DO N = 1, N_Pfaf_Catchs
            
            if(Pfaf_all(N) >= 1) call SEARCH_DNST (N, N_Pfaf_Catchs, DNST, Pfaf_all, DNST2(N))
            
         END DO
         
         RT_PARAM(:)%DNSTR = DNST2 (FirstCatch : FirstCatch + N_PfafL -1)
         
         deallocate (DNST, DNST2, Pfaf_all)
         
      endif
      
      STATUS   = NF_CLOSE (NCFID)    
      
      deallocate (temp_var)

    END SUBROUTINE Read_Restarts

    ! ------------------------------------------------------------------

    SUBROUTINE DOMAIN_DECOMP (N_PfafD, Local_Min, Local_Max,          &
         Pfaf_Min, Pfaf_Max, NCats_in_GRID, CAT_INDEX, LocalActive,   &
         LocDstCatchID, srcProcsID,dstCatchID,Allactive)
      
      implicit none
      
      INTEGER, INTENT (IN)                         :: N_PfafD, Local_Min, Local_Max
      INTEGER, INTENT (IN)                         :: Pfaf_Min, Pfaf_Max
      INTEGER, DIMENSION (:),  INTENT (IN)         :: NCats_in_GRID
      INTEGER, DIMENSION (:,:),INTENT (IN)         :: CAT_INDEX
      INTEGER, DIMENSION (N_PfafD), INTENT (INOUT) :: LocalActive, srcProcsID, LocDstCatchID
      INTEGER, DIMENSION (N_PfafD,comm_size), INTENT (INOUT) :: Allactive,dstCatchID
      
      INTEGER, DIMENSION(:)  ,ALLOCATABLE  :: global_buff, scounts, rdispls, rcounts
      INTEGER                              :: N_active, I,J,K,N,i1,i2,NProcs, req
      
      ! STEP 1: Identify active catchments within the local processor. If the catchment is active in 
      !         more than 1 processor, choose an owner.
      ! --------------------------------------------------------------------------------------------
      
      do N = 1, size (CAT_INDEX,2)
         if(NCats_in_GRID(N) > 0) then
            DO I = 1, NCats_in_GRID(N)
               LocalActive(CAT_INDEX(I,N) - Pfaf_Min + 1) =  CAT_INDEX(I,N) - Pfaf_Min + 1
            END DO
         endif
      end do
      
      allocate (global_buff (N_PfafD * comm_size))
      allocate (scounts(comm_size),rdispls(comm_size),rcounts(comm_size))  
      
      scounts = N_PfafD
      rcounts = N_PfafD
      
      rdispls(1) = 0
      global_buff= 0
      
      do i=2,comm_size
         rdispls(i)=rdispls(i-1)+rcounts(i-1)
      enddo

      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      
      ! call MPI_allgatherv  (                          &
      !     LocalActive, scounts         ,MPI_INTEGER, &
      !     global_buff, rcounts, rdispls,MPI_INTEGER, &
      !     MPI_COMM_WORLD, mpierr)

      do i = 1, comm_size
         
         if ((i == 1).and.(comm_rank == 0)) then
            global_buff((i-1)*N_PfafD+1:i*N_PfafD) = LocalActive(:)
         elseif (i > 1) then
            if(I-1 == comm_rank) then
              ! send to root
               call MPI_ISend(LocalActive,N_PfafD,MPI_INTEGER,0,993,MPI_COMM_WORLD,req,mpierr)
               call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
            else if (comm_rank == 0) then
               ! root receives
               call MPI_RECV(global_buff((i-1)*N_PfafD+1:i*N_PfafD),N_PfafD , MPI_INTEGER, i-1,993,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
            endif
         endif
      end do
      
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call MPI_BCAST (global_buff, N_PfafD*comm_size, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
      
      do i=1,comm_size
         Allactive (:,i) = global_buff((i-1)*N_PfafD+1:i*N_PfafD)
      enddo
      
      if (root_proc) then
         
         DO N = 1, N_PfafD
            NPROCS = count(Allactive(N,:) >= 1)
            if(NPROCS > 0)then
               if (NPROCS == 1) then
                  srcProcsID (N) = maxloc(Allactive(N,:),dim=1) - 1
               else
                  i1 = MAX(N - 5,1)
                  i2 = MIN(N + 5, N_PfafD)
                  N_active = 0
                  do I = 1,comm_size
                     if(Allactive (N,I) >= 1) then
                        if(count (Allactive(I1:I2,I) > 0) > N_active) then
                           N_active = count (Allactive(I1:I2,I) > 0)
                           J        = I
                        endif
                     endif
                  end do
                  srcProcsID (N) = J - 1
               endif
            endif
         END DO
         
      endif
      
      call MPI_BCAST (srcProcsID, N_PfafD, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
      
      ! STEP 2: convert downstream catchment indeces (from -1 OR 1:291284) of catchments that are
      !            in the local processor to full domain indeces.
      ! ------------------------------------------------------------------------------------------
      
      do N = Local_Min, Local_Max
         
         I = N - Pfaf_Min + 1
         
         if(LocalActive (I) >=1) then 
            
            if ((LocDstCatchID (I) == -1).OR.(LocDstCatchID (I) < Pfaf_Min) &
                 .OR. (LocDstCatchID (I) > Pfaf_Max)) then
               ! (a) DNST Catch lies outside of the full domain or a sink hole, set to drain to self 
               LocDstCatchID (I)    = I 
               
            else
               ! (b) DNST catch is an active catchment somewhere 
               LocDstCatchID  (I) =  LocDstCatchID (I) - Pfaf_Min + 1  ! changed to domain index           
               
            endif
            
         else
            
            LocDstCatchID  (I) = -9999 ! is inactive
            
         endif
      end do
      
      global_buff= 0
      
!      call MPI_allgatherv  (                          &
!           LocDstCatchID , scounts  ,MPI_INTEGER,     &
!           global_buff, rcounts, rdispls,MPI_INTEGER, &
!           MPI_COMM_WORLD, mpierr)

      do i = 1, comm_size         
         if ((i == 1).and.(comm_rank == 0)) then
            global_buff((i-1)*N_PfafD+1:i*N_PfafD) = LocDstCatchID(:)
         elseif (i > 1) then
            if(I-1 == comm_rank) then
              ! send to root
               call MPI_ISend(LocDstCatchID,N_PfafD,MPI_INTEGER,0,994,MPI_COMM_WORLD,req,mpierr)
               call MPI_WAIT (req,MPI_STATUS_IGNORE,mpierr)    
            else if (comm_rank == 0) then
               ! root receives
               call MPI_RECV(global_buff((i-1)*N_PfafD+1:i*N_PfafD),N_PfafD , MPI_INTEGER, i-1,994,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpierr)
            endif
         endif
      end do
      call MPI_Barrier(MPI_COMM_WORLD, mpierr)
      call MPI_BCAST (global_buff, N_PfafD*comm_size, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)
      
      do i=1,comm_size
         dstCatchID (:,i) = global_buff((i-1)*N_PfafD+1:i*N_PfafD)
      enddo
      
      deallocate (global_buff, scounts, rdispls, rcounts)
      
    END SUBROUTINE DOMAIN_DECOMP

    ! ------------------------------------------------------------------

    SUBROUTINE write_output (NCFOutID, N_PfafD, I, day, srcProcsID, RT_DIAGN_DAILY)

      implicit none
      integer, intent (in)                                         :: N_PfafD, NCFOutID,  I, day
      TYPE (route_diagn_type), DIMENSION (:), allocatable          :: rt_diagn_daily_f
      INTEGER,                 DIMENSION (N_PfafD), intent (in)    :: srcProcsID
      TYPE (route_diagn_type), DIMENSION (N_PfafD), intent (inout) :: RT_DIAGN_DAILY
      integer                                                      :: n,status,req
      real                                                         :: tmp_real (5), rbuff
      
      if(root_proc) then
         allocate (rt_diagn_daily_f (1:N_PfafD))
         rt_diagn_daily_f%QSFLOW   = MAPL_UNDEF
         rt_diagn_daily_f%QINFLOW  = MAPL_UNDEF
         rt_diagn_daily_f%QOUTFLOW = MAPL_UNDEF
         rt_diagn_daily_f%RUNOFF   = MAPL_UNDEF
      endif
      
      DO N = 1, N_PfafD
         
         if((srcProcsID(N) == 0).and.(comm_rank == 0)) then
            
            RT_DIAGN_DAILY_F(N)%QSFLOW   = RT_DIAGN_DAILY(N)%QSFLOW  
            RT_DIAGN_DAILY_F(N)%QINFLOW  = RT_DIAGN_DAILY(N)%QINFLOW 
            RT_DIAGN_DAILY_F(N)%QOUTFLOW = RT_DIAGN_DAILY(N)%QOUTFLOW
            RT_DIAGN_DAILY_F(N)%RUNOFF   = RT_DIAGN_DAILY(N)%RUNOFF
            
         else if (srcProcsID(N) > 0) then
            
            if (srcProcsID(N) == comm_rank) then 
               
               tmp_real (1) = RT_DIAGN_DAILY(N)%QSFLOW  
               tmp_real (2) = RT_DIAGN_DAILY(N)%QINFLOW 
               tmp_real (3) = RT_DIAGN_DAILY(N)%QOUTFLOW
               tmp_real (4) = RT_DIAGN_DAILY(N)%RUNOFF
               
               call MPI_ISend(tmp_real(1:4) ,4,MPI_real,0,999,MPI_COMM_WORLD,req,status)
               call MPI_WAIT (req,MPI_STATUS_IGNORE,status)
               
            else if (comm_rank == 0) then
               
               call MPI_RECV(tmp_real(1:4),4,MPI_real,srcProcsID(N),999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)
               
               RT_DIAGN_DAILY_F(N)%QSFLOW   = tmp_real (1)
               RT_DIAGN_DAILY_F(N)%QINFLOW  = tmp_real (2)
               RT_DIAGN_DAILY_F(N)%QOUTFLOW = tmp_real (3)  
               RT_DIAGN_DAILY_F(N)%RUNOFF   = tmp_real (4)  
               
            endif
         endif
      END DO
      
      call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
      
      if(comm_rank == 0) then
         
         status = NF_PUT_VARA_INT (NCFOutID, 1,(/day/),(/1/),I)
         status = NF_PUT_VARA_REAL(NCFOutID, 2,(/1,day/),(/N_PfafD,1/),RT_DIAGN_DAILY_F(:)%QSFLOW  )
         status = NF_PUT_VARA_REAL(NCFOutID, 3,(/1,day/),(/N_PfafD,1/),RT_DIAGN_DAILY_F(:)%QINFLOW )
         status = NF_PUT_VARA_REAL(NCFOutID, 4,(/1,day/),(/N_PfafD,1/),RT_DIAGN_DAILY_F(:)%QOUTFLOW)
         status = NF_PUT_VARA_REAL(NCFOutID, 5,(/1,day/),(/N_PfafD,1/),RT_DIAGN_DAILY_F(:)%RUNOFF  )
         deallocate (rt_diagn_daily_f)
      ENDIF

      RT_DIAGN_DAILY%QSFLOW   = 0. 
      RT_DIAGN_DAILY%QINFLOW  = 0. 
      RT_DIAGN_DAILY%QOUTFLOW = 0. 
      RT_DIAGN_DAILY%RUNOFF   = 0.
      
    END SUBROUTINE write_output
      
    ! ------------------------------------------------------------------
  
    SUBROUTINE Write_Restarts (N_PfafL, N_PfafD, Local_Min, Pfaf_Min,  &
         srcProcsID, rt_progn, rt_param, work_path, YYYYMM)
      
      implicit none
      
      integer, intent (in) :: N_PfafL, N_PfafD, Local_Min, Pfaf_Min
      character(*),                                 intent (in) :: work_path
      character*6,                                  intent (in) :: YYYYMM
      INTEGER,                 DIMENSION (N_PfafD), intent (in) :: srcProcsID
      TYPE (route_progn_type), DIMENSION (N_PfafL), intent (in) :: rt_progn
      TYPE (route_param_type), DIMENSION (N_PfafL), intent (in) :: rt_param
      
      type(Netcdf4_FileFormatter)                         :: InNCIO, OutNCIO
      type(FileMetadata)                                  :: meta_data
      integer                                             :: RC, nVars, N, I, req, status
      character*300                                       :: filename
      real                                                :: tmp_real (5)
      TYPE (route_progn_type), allocatable, DIMENSION (:) :: rts_f
      TYPE (route_param_type), allocatable, DIMENSION (:) :: rtp_f
      
      if (root_proc) then
         
         allocate (rts_f            (1:N_PfafD))
         allocate (rtp_f            (1:N_PfafD))
         
         rts_f%WSTREAM = MAPL_UNDEF
         rts_f%WRIVER  = MAPL_UNDEF
         rtp_f%AREACAT = MAPL_UNDEF
         rtp_f%LENGSC  = MAPL_UNDEF
         rtp_f%DNSTR   = -9999
         
      endif
      
      call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
      
      DO N = 1, N_PfafD
         
         I = N - Local_Min + Pfaf_Min
         
         if((srcProcsID(N) == 0).and.(comm_rank == 0)) then
            
            rts_f(N)%WSTREAM = RT_PROGN(I)%WSTREAM  
            rts_f(N)%WRIVER  = RT_PROGN(I)%WRIVER 
            Rtp_f(N)%AREACAT = RT_PARAM(I)%AREACAT
            Rtp_f(N)%LENGSC  = RT_PARAM(I)%LENGSC
            Rtp_f(N)%DNSTR   = RT_PARAM(I)%DNSTR
          
         else if (srcProcsID(N) > 0) then
            
            if (srcProcsID(N) == comm_rank) then 
               
               tmp_real (1) = RT_PROGN(I)%wstream
               tmp_real (2) = RT_PROGN(I)%WRIVER
               tmp_real (3) = RT_PARAM(I)%AREACAT
               tmp_real (4) = RT_PARAM(I)%LENGSC
               tmp_real (5) = real(RT_PARAM(I)%DNSTR)
               
               call MPI_ISend(tmp_real ,5,MPI_real,0,999,MPI_COMM_WORLD,req,status)
               call MPI_WAIT (req,MPI_STATUS_IGNORE,status)
               
            else if (comm_rank == 0) then
               
               call MPI_RECV(tmp_real,5,MPI_real,srcProcsID(N),999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)
               
               RTS_F(N)%wstream = tmp_real (1)
               RTS_F(N)%WRIVER  = tmp_real (2)
               RTP_F(N)%AREACAT = tmp_real (3)
               RTP_F(N)%LENGSC  = tmp_real (4)
               RTP_F(N)%DNSTR   = nint(tmp_real (5))
               
            endif
         endif
      END DO
      
      call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
      
      if (root_proc) then
         
         filename = TRIM(work_path)//'/rs/'//YYYYMM(1:4)//'/route_internal_rst.'//YYYYMM//'01'
         
         write (logunit,*) 'writing restart file   : ', trim(filename)
         
         call InNCIO%open('/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/route/' &
              //'route_internal_rst.YYYY0101', pFIO_READ, rc=rc)
         meta_data = InNCIO%read(rc=rc)
         call InNCIO%close(rc=rc)

         call meta_data%modify_dimension('tile',N_PfafD,rc=rc)
         call OutNCIO%create(trim(filename),rc=rc)
         call OutNCIO%write(meta_data, rc=rc)
                  
         call MAPL_VarWrite (OutNCIO,'AREACAT',RTP_F(:)%AREACAT )    
         call MAPL_VarWrite (OutNCIO,'LENGSC', RTP_F(:)%LENGSC  )    
         call MAPL_VarWrite (OutNCIO,'DNSTR',  REAL(RTP_F(:)%DNSTR))
         call MAPL_VarWrite (OutNCIO,'WSTREAM',RTS_F(:)%WSTREAM )    
         call MAPL_VarWrite (OutNCIO,'WRIVER', RTS_F(:)%WRIVER  )
         call OutNCIO%close(rc=rc)

         deallocate (RTP_F, RTS_F)
         
      endif
      
    END SUBROUTINE Write_Restarts

    ! ------------------------------------------------------------------

    SUBROUTINE open_monthly_outfile (NCFOutID,N_PfafD, filename, MyName)

      implicit none
      integer, intent (in)       :: N_PfafD
      integer, intent (inout)    :: NCFOutID
      character (*), intent (in) :: filename, MyName
      integer, dimension(8)      :: date_time_values
      character (22)             :: time_stamp
      integer                    :: status, d2(2),CellID,TimID, VID
      
      
      status = NF_CLOSE (NCFOutID)  
      status = NF_CREATE (filename, NF_NETCDF4, NCFOutID)
      status = NF_DEF_DIM(NCFOutID, 'N_PfafD' , N_PfafD,CellID)
      status = NF_DEF_DIM(NCFOutID, 'TSTEP'   , NF_UNLIMITED, TimID)
      
      d2(1)  = CellID
      d2(2)  = TimID
      
      status = NF_DEF_VAR(NCFOutID, 'DATE'    , NF_INT  , 1 ,TimID, vid)
      status = NF_DEF_VAR(NCFOutID, 'QSFLOW'  , NF_FLOAT, 2 ,d2   , vid)
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',&
           LEN_TRIM('TRANSFER_OF_MOISTURE_FROM_STREAM_VARIABLE_TO_RIVER_VARIABLE'), &
           trim('TRANSFER_OF_MOISTURE_FROM_STREAM_VARIABLE_TO_RIVER_VARIABLE')) 
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'UNITS',    &
           LEN_TRIM('m+3 s-1'), trim('m+3 s-1')) 
      status = NF_DEF_VAR(NCFOutID, 'QINFLOW' , NF_FLOAT, 2 ,d2   , vid)
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',&
           LEN_TRIM('TRANSFER_OF_RIVER_WATER_FROM_UPSTREAM_CATCHMENTS'), &
           trim('TRANSFER_OF_RIVER_WATER_FROM_UPSTREAM_CATCHMENTS')) 
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'UNITS',    &
           LEN_TRIM('m+3 s-1'), trim('m+3 s-1')) 
      status = NF_DEF_VAR(NCFOutID, 'QOUTFLOW' , NF_FLOAT, 2 ,d2   , vid)
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',&
           LEN_TRIM('TRANSFER_OF_RIVER_WATER_TO_DOWNSTREAM_CATCHMENT'), &
           trim('TRANSFER_OF_RIVER_WATER_TO_DOWNSTREAM_CATCHMENT')) 
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'UNITS',    &
           LEN_TRIM('m+3 s-1'), trim('m+3 s-1')) 
      status = NF_DEF_VAR(NCFOutID, 'RUNOFF' , NF_FLOAT, 2 ,d2   , vid)
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',&
           LEN_TRIM('RUNOFF_INPUT_FROM_LAND_SURFACE_MODEL'), &
           trim('RUNOFF_INPUT_FROM_LAND_SURFACE_MODEL')) 
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'UNITS',    &
           LEN_TRIM('m+3 s-1'), trim('m+3 s-1')) 
      status = NF_DEF_VAR(NCFOutID, 'BFLOW' , NF_FLOAT, 2 ,d2   , vid)
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',&
           LEN_TRIM('BFLOW_INPUT_FROM_LAND_SURFACE_MODEL'), &
           trim('BFLOW_INPUT_FROM_LAND_SURFACE_MODEL')) 
      status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'UNITS',    &
           LEN_TRIM('m+3 s-1'), trim('m+3 s-1')) 
      !
      !  Global attributes
      !
      call date_and_time(VALUES=date_time_values)
      
      write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
           date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
           date_time_values(5),':',date_time_values(6),':',date_time_values(7)
      
      status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'CreatedBy', LEN_TRIM(MYNAME),  &
           trim(MYNAME))
      status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))
      status = NF_ENDDEF(NCFOutID)    

    END SUBROUTINE open_monthly_outfile
    
    ! ----------------------------------------------------------------------

    integer function VarID (NCFID, VNAME) 
      
      integer, intent (in)      :: NCFID
      character(*), intent (in) :: VNAME
      integer                   :: status
      
      STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID) ; VERIFY_(STATUS)
      
    end function VarID
 
   ! -----------------------------------------------------------------------------
   
    SUBROUTINE Read_Resource (MET_PATH, WORK_PATH, domain_lonlat,BEG_DATE, END_DATE)
    
      implicit none
      
      integer :: n,status
      character(*), intent (inout)      :: WORK_PATH, MET_PATH
      real, dimension (4), intent(inout):: domain_lonlat
      integer, intent(inout)            :: BEG_DATE, END_DATE
      character*120, dimension (15)     :: route_rc
      character*120                     :: tmpstr,tmpstr2
      
      open (10,file = 'run_rrm.x',form = 'formatted', status='old',action = 'read')
      
      do n = 1,size (route_rc)
         read(10,'(a)',IOSTAT=status)tmpstr
         if(status == 0) route_rc(n) = tmpstr
      end do
      
      close (10,status='keep')

      domain_lonlat = (/-180.,-90.,180.,90./)

      do n = 1,size (route_rc)
         
         tmpstr = route_rc(n)

         if(.not.index(tmpstr,'#')) then 
            
            if(index(tmpstr,'WORK_PATH'  )) WORK_PATH   = tmpstr(index(tmpstr,':')+2:index(tmpstr,' ',back=.true.)-1)
            if(index(tmpstr,'MET_PATH'   )) MET_PATH    = tmpstr(index(tmpstr,':')+2:index(tmpstr,' ',back=.true.)-1)
            if(index(tmpstr,'BEG_DATE'   )) read (tmpstr(index(tmpstr,':')+1:index(tmpstr,' ',back=.true.)-1), *) BEG_DATE
            if(index(tmpstr,'END_DATE'   )) read (tmpstr(index(tmpstr,':')+1:index(tmpstr,' ',back=.true.)-1), *) END_DATE
            if(index(tmpstr,'MIN_LON'    )) read (tmpstr(index(tmpstr,':')+1:index(tmpstr,' ',back=.true.)-1), *) domain_lonlat (1)
            if(index(tmpstr,'MAX_LON'    )) read (tmpstr(index(tmpstr,':')+1:index(tmpstr,' ',back=.true.)-1), *) domain_lonlat (3)
            if(index(tmpstr,'MIN_LAT'    )) read (tmpstr(index(tmpstr,':')+1:index(tmpstr,' ',back=.true.)-1), *) domain_lonlat (2)
            if(index(tmpstr,'MAX_LAT'    )) read (tmpstr(index(tmpstr,':')+1:index(tmpstr,' ',back=.true.)-1), *) domain_lonlat (4)
            
         end if
      end do
      
    END SUBROUTINE Read_Resource
    
    FUNCTION UpperCase ( Input_String ) RESULT ( Output_String ) 
      CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
      CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      CHARACTER( * ), INTENT( IN ) :: Input_String 
      CHARACTER( LEN( Input_String ) ) :: Output_String 
      INTEGER :: i, n 
      
      Output_String = Input_String 
      
      DO i = 1, LEN( Output_String ) 
         n = INDEX( LOWER_CASE, Output_String( i:i ) ) 
         IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n ) 
      END DO
    END FUNCTION UpperCase

  END PROGRAM GEOSroute_driver

    
    


