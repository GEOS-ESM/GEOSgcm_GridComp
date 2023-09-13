#include "MAPL_Generic.h"

module GEOS_GigatrajGridCompMod
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated
  use, intrinsic :: iso_c_binding, only : c_loc
  use ESMF
  use MAPL
  use MAPL_VerticalDataMod
  use mpi
  use GEOS_Giga_interOpMod
  implicit none

  public :: SetServices

  private
  type horde
    integer :: num_parcels
    integer, allocatable :: IDS(:)
    real, allocatable    :: lats(:), lons(:), zs(:)
  end type

  type GigaTrajInternal
    integer :: npes
    integer :: npz ! number of pressure levels
    type (ESMF_Grid) :: LatLonGrid
    type (ESMF_Grid) :: CubedGrid
    class (AbstractRegridder), pointer :: cube2latlon => null()
    integer, allocatable :: CellToRank(:,:)
    type(horde) :: parcels
    type(c_ptr) :: metSrc
    type(ESMF_Time) :: startTime
    character(len=ESMF_MAXSTR), allocatable :: ExtraFieldNames(:)
    type(VerticalData) :: vdata
    logical :: regrid_to_latlon
  end type

  type GigatrajInternalWrap
    type (GigaTrajInternal), pointer :: PTR
  end type


contains

  subroutine SetServices ( GC, RC )
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    type (GigaTrajInternal), pointer   :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)        :: wrap

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // 'SetServices'

    ! Register services for this component
    ! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, _RC )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  GetInitVars , _RC )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run         , _RC )


    ! Internal state

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter, _RC)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter, _RC)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'OMEGA',                                     &
         LONG_NAME  = 'vertical_pressure_velocity',                &
         UNITS      = 'Pa s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter, _RC)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PL',                                        &
         LONG_NAME  = 'mid_level_pressure',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,  _RC)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'edge_pressure',                             &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,    _RC)

    allocate(GigaTrajInternalPtr)
    wrap%ptr => GigaTrajInternalPtr
    call ESMF_UserCompSetInternalState(GC, 'GigaTrajInternal', wrap, _RC) 

    call MAPL_GenericSetServices(GC, _RC )

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,_RC)
    call MAPL_TimerAdd(GC, name="RUN"           ,_RC)

    RETURN_(ESMF_SUCCESS)  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code
  
    character(len=ESMF_MAXSTR)           :: IAm 
    integer                              :: STATUS
    character(len=ESMF_MAXSTR)           :: COMP_NAME
  
    ! Local derived type aliases
    type (MAPL_MetaComp),  pointer  :: MPL 
    type (ESMF_VM)                  :: vm
    integer :: I1, I2, J1, J2, comm, npes, my_rank, rank, ierror, NX, NY, NPZ
    type(ESMF_Grid) :: CubedGrid
    integer, allocatable :: I1s(:), J1s(:), I2s(:),J2s(:)
    integer :: DIMS(3), counts(3), i, j, k, l
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
    type (GigatrajInternalWrap)   :: wrap
    real :: dlat, dlon, lon_start
    real, allocatable, target :: lats_center(:), lons_center(:), levs_center(:)
    real, allocatable, target :: cube_lats_center(:, :), cube_lons_center(:,:)
    real, pointer :: levs_ptr(:), ptr(:,:)
    type (ESMF_TIME) :: CurrentTime
    character(len=20), target :: ctime 
    character(len=:), allocatable, target :: name_, unit_
    type(ESMF_Alarm)  :: GigaTrajOutAlarm, GigaTrajRebalanceAlarm
    type(ESMF_TimeInterval)      :: parcels_DT, Rebalance_DT
    type(ESMF_TimeInterval)      :: ModelTimeStep
    type(ESMF_State)             :: INTERNAL
    integer :: minutes_, imc, jmc
    character(len=ESMF_MAXSTR) :: parcels_file
    character(len=ESMF_MAXSTR) :: grid_name
    character(len=ESMF_MAXSTR) :: regrid_to_latlon
    real(ESMF_KIND_R8), pointer     :: centerX(:,:)
    real(ESMF_KIND_R8), pointer     :: centerY(:,:)
    type(ESMF_Field)   :: field
    type(ESMF_RouteHandle) :: rh

    call ESMF_GridCompGet ( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // "Initialize"

    call MAPL_GetObjectFromGC ( GC, MPL, _RC)

    call MAPL_TimerOn(MPL,"TOTAL")
    call MAPL_TimerOn(MPL,"INITIALIZE")

    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, _RC)

    call MAPL_GetResource(MPL, minutes_, "GIGATRAJ_REBALANCE_MINUTES:", default=30, _RC) 
    call ESMF_TimeIntervalSet(Rebalance_DT, m= minutes_, _RC)
    call MAPL_GetResource(MPL, minutes_, "GIGATRAJ_OUTPUT_MINUTES:", default=30, _RC)
    call ESMF_TimeIntervalSet(parcels_DT,   m= minutes_, _RC)

    GigaTrajOutAlarm = ESMF_AlarmCreate(                   &
         clock,                                            &
         name='GigatrajOut',                               &
         ringTime= CurrentTime + parcels_DT-ModelTimeStep, &
         ringInterval=parcels_DT,                          &
         ringTimeStepCount=1,                              &
         sticky=.false., _RC)

    GigaTrajRebalanceAlarm = ESMF_AlarmCreate(                   &
         clock,                                            &
         name='GigatrajRebalance',                               &
         ringTime= CurrentTime + parcels_DT-ModelTimeStep, &
         ringInterval=parcels_DT,                          &
         ringTimeStepCount=1,                              &
         sticky=.false., _RC)


    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  _RC)

    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
    call MPI_Comm_size(comm, npes, ierror); _VERIFY(ierror)
    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)
    
    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, status); _VERIFY(STATUS)
    GigaTrajInternalPtr => wrap%ptr
    GigaTrajInternalPtr%npes = npes

    call ESMF_GridCompGet(GC, grid=CubedGrid, _RC)
    call MAPL_GridGet(CubedGrid, globalCellCountPerDim=DIMS, _RC)

    levs_center = [1000.,  975.,  950.,  925.,  900., 875., 850., 825., 800., 775., 750., 725.,700., 650., 600., 550., 500., &
                   450., 400., 350., 300., 250., 200., 150., 100.,  70.,  50.,  40.,  30., 20., 10., 7., 5., 4., 3., 2., &
                   1., 0.7, 0.5, 0.4,  0.3, 0.2, 0.1, 0.07, 0.05, 0.04, 0.03, 0.02]
    npz = size(levs_center, 1)
    GigaTrajInternalPtr%npz = npz

    call MAPL_GetResource(MPL, NX, "NX:", _RC)
    call MAPL_GetResource(MPL, grid_name, "AGCM_GRIDNAME:", _RC)

    GigaTrajInternalPtr%CubedGrid = grid_manager%make_grid(&
                 CubedSphereGridFactory(grid_name=trim(grid_name),im_world = DIMS(1), lm=npz, nx=NX, ny=NX, rc=status)); _VERIFY(status)

    call MAPL_GetResource(MPL, regrid_to_latlon, "GIGATRAJ_REGRID_TO_LATLON:", default='YES', _RC)
    GigaTrajInternalPtr%regrid_to_latlon = .true.

    if (trim(regrid_to_latlon) == "NO") GigaTrajInternalPtr%regrid_to_latlon = .false.
    
    if ( GigaTrajInternalPtr%regrid_to_latlon ) then

       call MAPL_MakeDecomposition(NX,NY,_RC)

       GigaTrajInternalPtr%LatLonGrid = grid_manager%make_grid( &
                 LatLonGridFactory(im_world=DIMS(1)*4, jm_world=DIMS(1)*2+1, lm=npz,  &
                 nx=NX, ny=NY, pole='PC', dateline= 'DC', rc=status) ); _VERIFY(status) 

       GigaTrajInternalPtr%cube2latlon => new_regridder_manager%make_regridder(GigaTrajInternalPtr%CubedGrid,  GigaTrajInternalPtr%LatLonGrid, REGRID_METHOD_CONSERVE, _RC)

       call MAPL_Grid_interior(GigaTrajInternalPtr%LatLonGrid ,i1,i2,j1,j2)

       call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

       ! lat and lon centers need to hold the halo with width 1

       dlon = 360.0/dims(1)
       ! DE
       !lons_center = [(dlon*(i-1)+dlon/2., i= i1-1, i2+1)]
       ! DC    
       lons_center = [(dlon*(i-1), i= i1-1, i2+1)]

       !where(lons_center < 0. ) lons_center  = 0.
       !where(lons_center >360.) lons_center = 360.

       !PE
       !dlat = 180.0/dims(2) 
       !lats_center = [(-dlat/2. + dlat*j-90.0, j= j1-1, j2+1)] 
       !PC
       dlat = 180.0/(dims(2)-1)  ! PC
       lats_center = [(-90.0 + (j-1)*dlat, j= j1-1, j2+1)] 
       !where(lats_center <-90.) lats_center = -90.
       !where(lats_center >90. ) lats_center =  90.

    else

      call MAPL_Grid_interior(CubedGrid ,i1,i2,j1,j2)
      imc = i2-i1 + 1
      jmc = j2-j1 + 1
      allocate(cube_lats_center(imc+2, jmc+2))
      allocate(cube_lons_center(imc+2, jmc+2))
      allocate(ptr(0:imc+1, 0:jmc+1))

      call ESMF_GridGetCoord(CubedGrid, coordDim=1, localDE=0, &
             staggerloc=ESMF_STAGGERLOC_CENTER, &
             farrayPtr=centerX, rc=status)
      VERIFY_(STATUS)

      field = ESMF_FieldCreate(CubedGrid,ptr,staggerLoc=ESMF_STAGGERLOC_CENTER,totalLWidth=[1,1],totalUWidth=[1,1],rc=status)
      _VERIFY(status)
      call ESMF_FieldHaloStore(field,rh,rc=status)
      _VERIFY(status)

      ptr(1:imc,1:jmc)=centerX
      call ESMF_FieldHalo(field,rh,rc=status)
       _VERIFY(status)
      cube_lons_center = ptr

      call ESMF_GridGetCoord(CubedGrid , coordDim=2, localDE=0, &
              staggerloc=ESMF_STAGGERLOC_CENTER, &
              farrayPtr=centerY, rc=status)
      VERIFY_(STATUS)
      ptr(1:imc,1:jmc)=centerY
      call ESMF_FieldHalo(field,rh,rc=status)
       _VERIFY(status)
      cube_lats_center = ptr

      deallocate(ptr)
      call ESMF_FieldDestroy(field,rc=status)
      _VERIFY(status)
      call ESMF_FieldHaloRelease(rh,rc=status)
      _VERIFY(status)

      cube_lons_center = cube_lons_center*180.0/MAPL_PI
      cube_lats_center = cube_lats_center*180.0/MAPL_PI

    endif

    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime)
    ctime(20:20) = c_null_char

    allocate(I1s(npes),J1s(npes))
    allocate(I2s(npes),J2s(npes))

    call MPI_Allgather(i1, 1, MPI_INTEGER, I1s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)
    call MPI_Allgather(i2, 1, MPI_INTEGER, I2s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)
    call MPI_Allgather(j1, 1, MPI_INTEGER, J1s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)
    call MPI_Allgather(j2, 1, MPI_INTEGER, J2s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)

    allocate(GigaTrajInternalPtr%CellToRank(DIMS(1),DIMS(2)))

    do rank = 0, npes -1
       I1 = I1s(rank+1)
       I2 = I2s(rank+1)      
       J1 = J1s(rank+1)
       J2 = J2s(rank+1)
       GigaTrajInternalPtr%CellToRank(I1:I2,J1:J2) = rank
    enddo

    ! WJiang notes: the vcoord should be consistent with the HISTORY.rc
    levs_ptr=>levs_center
    GigaTrajInternalPtr%vdata = VerticalData(levs_ptr, vcoord = 'log(PLE)', vscale = 100.0, vunit = 'hPa',_RC)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
      call MAPL_Grid_interior(GigaTrajInternalPtr%LatLonGrid,i1,i2,j1,j2)
      GigaTrajInternalPtr%metSrc = initMetGEOSDistributedLatLonData(comm, c_loc(GigaTrajInternalPtr%CellToRank), DIMS(1), DIMS(2),   &
                                   npz, counts(1)+2, counts(2)+2, npz, &
                                   c_loc(lons_center), c_loc(lats_center), c_loc(levs_center), c_loc(ctime))
      deallocate(lons_center, lats_center,levs_center)
    else
      call MAPL_Grid_interior(CubedGrid ,i1,i2,j1,j2)
      GigaTrajInternalPtr%metSrc = initMetGEOSDistributedCubedData(comm, c_loc(GigaTrajInternalPtr%CellToRank), DIMS(1),   &
                                   npz, i1, i2, j1, j2, npz, &
                                   c_loc(cube_lons_center), c_loc(cube_lats_center), c_loc(levs_center), c_loc(ctime))
      deallocate(cube_lons_center, cube_lats_center,levs_center)

    endif

    call read_parcels(GC, GigaTrajInternalPtr, _RC) 

    call MAPL_TimerOff(MPL,"INITIALIZE")
    call MAPL_TimerOff(MPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine GetInitVars ( GC, IMPORT, EXPORT, CLOCK, RC )
   type(ESMF_GridComp), intent(inout) :: GC      
   type(ESMF_State),    intent(inout) :: IMPORT 
   type(ESMF_State),    intent(inout) :: EXPORT
   type(ESMF_Clock),    intent(inout) :: CLOCK  
   integer, optional,   intent(  out) :: RC     

   integer :: i, status, k, FC, l, lm
   character(len=ESMF_MAXSTR) :: IAm 
   type (ESMF_State)          :: INTERNAL, leaf_export
   type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
   type (GigatrajInternalWrap)   :: wrap
   character(len=ESMF_MAXSTR)    :: GigaRstFile
   character(len=ESMF_MAXSTR)    :: other_fields
   character(len=ESMF_MAXSTR)    :: NAME
   type(ESMF_Field) :: tmp_field
   type (ESMF_FieldBundle)       :: TRI
   type (ESMF_TIME) :: CurrentTime
   character(len=20), target :: ctime 
   type (MAPL_MetaComp),  pointer  :: MPL
   logical, save :: init = .false.
   real, dimension(:,:,:) , pointer :: ptr3d
   type(Netcdf4_fileformatter) :: formatter
   type(FileMetadata) :: meta
   character(len=ESMF_MAXSTR) :: parcels_file
   character(len=:), allocatable :: fieldname
   character(len=20) :: diffusions(4)
   character (len=ESMF_MAXSTR), allocatable  :: itemNameList(:)
   character (len=ESMF_MAXSTR), allocatable  :: fieldnames(:)
   integer :: nitems
   logical :: file_exists

   Iam = "getInitVars"

   if (init) then
     RETURN_(ESMF_SUCCESS)
   endif

   call MAPL_GetObjectFromGC ( GC, MPL, _RC)

   call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
   call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime) 
   ctime(20:20) = c_null_char

   call MAPL_Get(MPL, INTERNAL_ESMF_STATE=INTERNAL, _RC)

   call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
   GigaTrajInternalPtr => wrap%ptr

   call MAPL_GetResource(MPL, GigaRstFile, 'GIGATRAJ_INTERNAL_RESTART_FILE:', default="NONE", RC=STATUS )

   if (trim(GigaRstFile) == 'NONE') then
      ! without restart file, get value from import
      call init_metsrc_field0(GC,  IMPORT,  ctime, 'PLE', _RC)
   else
      INQUIRE(FILE= GigaRstFile, EXIST=file_exists) 
      _ASSERT( file_exists, " GIGATRAJ_INTERNAL_RESTART_FILE does not exist")
      call init_metsrc_field0(GC, INTERNAL, ctime, 'PL',  _RC)
   endif

   call MAPL_GetResource(MPL, other_fields, "GIGATRAJ_EXTRA_FIELDS:", default='NONE', _RC)

   if (other_fields /= 'NONE') then
      call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC)
      if (MAPL_AM_I_ROOT()) then
        call formatter%open(trim(parcels_file), pFIO_WRITE, _RC)
        meta = formatter%read(_RC)
      endif

      call getExtraFieldNames(other_fields,GigaTrajInternalPtr%ExtraFieldNames)

      do i = 1, size(GigaTrajInternalPtr%ExtraFieldNames)
         if (index(GigaTrajInternalPtr%ExtraFieldNames(i), 'CA.bc') /=0) then
           call MAPL_ExportStateGet([import], 'CA.bc', leaf_export, _RC)
           call ESMF_StateGet(leaf_export, trim(GigaTrajInternalPtr%ExtraFieldNames(i)), tmp_field, _RC)
           call ESMF_FieldGet(tmp_field,farrayPtr=ptr3d,rc=status)
           if( .not. associated(ptr3d)) then
              call MAPL_AllocateCoupling(tmp_field, _RC)
           endif
           call ESMF_StateAddReplace(import, [tmp_field], _RC)
           do k = 1, size(ptr3d, 3)
              fieldname = trim(GigaTrajInternalPtr%ExtraFieldNames(i))//'00'//i_to_string(k)
              call create_extra_var(fieldname)
           enddo
           cycle
         endif

         if (index(GigaTrajInternalPtr%ExtraFieldNames(i), 'CA.oc') /=0) then
           call MAPL_ExportStateGet([import], 'CA.oc', leaf_export, _RC)
           call ESMF_StateGet(leaf_export, trim(GigaTrajInternalPtr%ExtraFieldNames(i)), tmp_field, _RC)
           call ESMF_FieldGet(tmp_field,farrayPtr=ptr3d,rc=status)
           if( .not. associated(ptr3d)) then
              call MAPL_AllocateCoupling(tmp_field, _RC)
           endif
           call ESMF_StateAddReplace(import, [tmp_field], _RC)
           do k = 1, size(ptr3d, 3)
              fieldname = trim(GigaTrajInternalPtr%ExtraFieldNames(i))//'00'//i_to_string(k)
              call create_extra_var(fieldname)
           enddo
           cycle
         endif

         if (index(GigaTrajInternalPtr%ExtraFieldNames(i), 'TRI') /=0) then
           call ESMF_StateGet(import, 'PHYSICS_Exports/TURBULENCE_Exports/TRI', TRI, _RC)
           !call ESMF_FieldBundleGet(TRI, fieldCount=nitems, _RC)
           !allocate(itemNameList(nitems))
           !call ESMF_FieldBundleGet(TRI,fieldnamelist=itemNameList,rc=status)
           allocate(fieldnames(4))
           fieldnames(1) = "CA.bc::CA.bcphilicIT"
           fieldnames(2) = "CA.bc::CA.bcphobicIT"
           fieldnames(3) = "CA.oc::CA.ocphilicIT"
           fieldnames(4) = "CA.oc::CA.ocphobicIT"
           do k = 1, 4
             !call ESMF_FieldBundleGet(TRI, trim(itemNameList(k)), field=tmp_field, _RC)
             !call ESMF_FieldBundleGet(TRI, trim(fieldnames(k)), field=tmp_field, _RC)
             fieldname = trim(fieldnames(k))
             call ESMF_FieldBundleGet(TRI, fieldname, field=tmp_field, _RC)

             !call MAPL_AllocateCoupling(tmp_field, _RC)

             call ESMF_FieldGet(tmp_field,farrayPtr=ptr3d, rc=status)
             !if (MAPL_AM_I_Root()) then
             if ( status /=0 )then
               print*, 'status, associated(ptr3d)', status, associated(ptr3d)
               print*, "not created: ", fieldname 
               print*, "Add to historty.rc to triger the allocation" 
               _FAIL(" Not allocated diffusion tendency")
             !else
             !   print*, "createded: ", itemNameList(k)
             !   print*, "assocated: ", associated(ptr3d), itemNameList(k)    
             endif
             !endif
             call create_extra_var(fieldname(8:))
           enddo
           deallocate(fieldnames)
           cycle
         endif
         call create_extra_var(trim(GigaTrajInternalPtr%ExtraFieldNames(i)))
      enddo
   endif

   if (MAPL_AM_I_Root()) then
      call formatter%close()
   endif

   init = .true.

   RETURN_(ESMF_SUCCESS)
   contains
     subroutine getExtraFieldNames(other_fields, Fieldnames)
       character(*), intent(in) :: other_fields
       character(len=ESMF_MAXSTR), allocatable, intent(out) :: Fieldnames(:)
       integer :: num_field, i, k, endl, num_
       num_field = 1
       k = 1
       do
         i = index(other_fields(k:),';')
         if (i == 0) exit
         k = k+i
         num_field = num_field+1
       enddo

       allocate(Fieldnames(num_field))
       k    = 1
       num_ = 1
       do
         i = index(other_fields(k:),';')
         if (i == 0) then
            endl = len(other_fields)
         else
            endl = (k-1)+i-1
         endif
         FieldNames(num_) = trim(adjustl(other_fields(k:endl)))
         num_ = num_ + 1
         k = endl + 2
         if (num_ > num_field) exit
       enddo
     end subroutine getExtraFieldNames

     ! if the field name is not in original file
     subroutine create_extra_var(fieldname)
       character(*), intent(in) :: fieldname
       type(Variable) :: var
       character(len=:), allocatable :: long_name, units, var_name
       real, allocatable :: tmp(:)
       select case (fieldname)
       case('TH')
         var_name = 'theta'
         long_name = "air_potential_temperature"
         units = "K"
       case('T')         
         var_name = 't' 
         long_name = "air_temperature"
         units = "K"
       case('PALT')            
         var_name = 'palt'
         long_name = "pressure_altitude"
         units = "km"
       case default
         var_name = trim(fieldname)
         long_name = "unkown"
         units = "1"
         !print*, "Not yet define attribute of "//var_name
       end select

       if (index(var_name, 'CA.bcSD') /=0) then
          long_name = "Carbonaceous_Aerosol_Sedimentation_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'CA.bcSV') /=0) then
          long_name = "Carbonaceous_Aerosol_Convective_Scavenging_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'CA.bcDP') /=0) then
          long_name = "Carbonaceous_Aerosol_Dry_Deposition_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'CA.bcWT') /=0) then
          long_name = "Carbonaceous_Aerosol_Wet_Deposition_bin" 
          units =  "kg m-2 s-1"
       endif

       if (index(var_name, 'CA.ocSD') /=0) then
          long_name = "Carbonaceous_Aerosol_Sedimentation_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'CA.ocSV') /=0) then
          long_name = "Carbonaceous_Aerosol_Convective_Scavenging_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'CA.ocDP') /=0) then
          long_name = "Carbonaceous_Aerosol_Dry_Deposition_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'CA.ocWT') /=0) then
          long_name = "Carbonaceous_Aerosol_Wet_Deposition_bin" 
          units =  "kg m-2 s-1"
       endif
       if (index(var_name, 'phobic') /=0 .or. index(var_name, 'philic')/=0) then
          long_name = "Carbonaceous_Aerosol_Mixing_Ratio" 
          units =  "kg kg-1"
       endif

       if( meta%has_variable(var_name)) return
       if (MAPL_AM_I_Root()) then
         var = variable(type=pFIO_REAL32, dimensions='id,time')
         call var%add_attribute('long_name', long_name)
         call var%add_attribute('units', units)
         call var%add_attribute('positive', "up")
         call var%add_attribute('_FillValue', -999.99)
         call var%add_attribute('missing_value', -999.99)
         call meta%add_variable(var_name, var)
         call formatter%add_variable(meta, var_name)
       endif
     end subroutine create_extra_var
 end subroutine GetInitVars

 subroutine Init_metsrc_field0 (GC, state, ctime, PL, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state
    character(*), target,   intent(in) :: ctime
    character(*),           intent(in) :: PL
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm 
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap
    real, dimension(:,:,:), pointer     :: U, V, W, P, PL0, PLE
    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable :: U_inter, V_inter, W_inter, P_inter
    real, dimension(:,:,:), allocatable, target  :: haloU, haloV, haloW, haloP
    integer :: counts(3), dims(3), d1,d2,km,lm,l

    Iam = "init_metsrc_field0"

    call MAPL_GetPointer(state, U, "U", _RC)
    call MAPL_GetPointer(state, V, "V", _RC)
    call MAPL_GetPointer(state, W, "OMEGA", _RC)


    PL0=>null()
    if (PL == 'PL') then
       call MAPL_GetPointer(state, P, "PL", _RC)
    else if (PL == 'PLE') then
       call MAPL_GetPointer(state, PLE, "PLE", _RC)
       d1 =size(PLE,1)
       d2 =size(PLE,2)
       km =size(PLE,3)-1
       allocate(PL0(d1,d2,km))
       ! WJ notes, PLE's lower bound is (1,1,0)
       PL0 = (PLE(:,:,1:km)+PLE(:,:,0:km-1))*0.5
       P => PL0
       W = 0.0
    endif

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr
    call ESMF_StateGet(state, 'PLE', field=GigaTrajInternalPtr%vdata%interp_var, rc=status)
    call GigaTrajInternalPtr%vdata%setup_eta_to_pressure(_RC)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)
    else
       call MAPL_GridGet(GigaTrajInternalPtr%CubedGrid,  localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)
    endif
   !#lm = counts(3)
    lm = size(GigaTrajInternalPtr%vdata%levs)

    d1 =size(U,1)
    d2 =size(U,2)

    allocate(U_latlon(counts(1),counts(2),lm))
    allocate(V_latlon(counts(1),counts(2),lm))
    allocate(W_latlon(counts(1),counts(2),lm))
    allocate(P_latlon(counts(1),counts(2),lm))

    allocate(U_inter(d1,d2,lm))
    allocate(V_inter(d1,d2,lm))
    allocate(W_inter(d1,d2,lm))
    allocate(P_inter(d1,d2,lm))

    allocate(haloU(counts(1)+2, counts(2)+2,lm), source = 0.0)
    allocate(haloV(counts(1)+2, counts(2)+2,lm), source = 0.0)
    allocate(haloW(counts(1)+2, counts(2)+2,lm), source = 0.0)
    allocate(haloP(counts(1)+2, counts(2)+2,lm), source = 0.0)

    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(U, U_inter,rc=status)
    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(V, V_inter,rc=status)
    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(W, W_inter,rc=status)
    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(P, P_inter,rc=status)
    if ( GigaTrajInternalPtr%regrid_to_latlon) then
       call GigaTrajInternalPtr%cube2latlon%regrid(U_inter, V_inter, U_latlon, V_latlon, _RC)
       call GigaTrajInternalPtr%cube2latlon%regrid(W_inter, W_latlon, _RC)
       call GigaTrajInternalPtr%cube2latlon%regrid(P_inter, P_latlon, _RC)

       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, U_Latlon, haloU, _RC)
       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, V_Latlon, haloV, _RC)
       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, W_Latlon, haloW, _RC)
       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, P_Latlon, haloP, _RC)
    else
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, U_inter, haloU, _RC)
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, V_inter, haloV, _RC)
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, W_inter, haloW, _RC)
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, P_inter, haloP, _RC)
    endif

    call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(haloU), c_loc(haloV), c_loc(haloW), c_loc(haloP))

    if(associated(PL0)) deallocate(PL0)
    RETURN_(ESMF_SUCCESS)
    
 end subroutine init_metsrc_field0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)       :: IAm 
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME
    integer        :: CSTAT, ESTAT, YY, MM, HH, DD, H, M,S
    character(512) :: CMSG
    character(256) :: command_line
    character(19)  :: begdate, enddate
    character(64)  :: format_string
    type(ESMF_TimeInterval) :: ModelTimeStep
    type(ESMF_Time)         :: CurrentTime
    type(ESMF_Grid)         :: grid_
  
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap
  
    integer :: num_parcels, my_rank,lm, l, d1,d2
    real, allocatable, target :: lats(:), lons(:), zs(:)
    real, allocatable, target :: U(:), V(:), W(:)
    integer ::counts(3), DIMS(3), rank, comm, ierror
    type (ESMF_VM)   :: vm
  
    real, dimension(:,:,:), pointer     :: U_cube, V_cube, W_cube, P_cube, PLE_Cube, with_halo
    real, dimension(:,:,:), pointer     :: U_internal, V_internal, W_internal, P_internal, PLE_internal

    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable :: U_inter, V_inter, W_inter, P_inter
    real, dimension(:,:,:), allocatable, target  :: haloU, haloV, haloW, haloP

    real(ESMF_KIND_R8) :: DT

    type(ESMF_Field)   :: halo_field
    type(ESMF_RouteHandle) :: rh
    character(len=20), target :: ctime, ctime0
    type(ESMF_State)            :: INTERNAL
    type(MAPL_MetaComp),pointer :: MPL  
    character(len=ESMF_MAXSTR) :: parcels_file

    !RETURN_(ESMF_SUCCESS)
    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)

    call MAPL_GetObjectFromGC ( GC, MPL, _RC)
    call MAPL_Get (MPL, INTERNAL_ESMF_STATE=INTERNAL, _RC)

    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, _RC)
 
    ! W.J note: this run is after agcm's run. The clock is not yet ticked
    !           So the values we are using are at (CurrentTime + ModelTimeStep)
    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime0) 
    ctime0(20:20) = c_null_char
    CurrentTime = CurrentTime + ModelTimeStep
    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime) 
    ctime(20:20) = c_null_char

    call ESMF_TimeIntervalGet(ModelTimeStep,d_r8=DT, _RC)

!---------------
! Step 1) Regrid the metData field from cubed to lat-lon
!---------------

    call MAPL_GetPointer(Import, U_cube, "U", _RC)
    call MAPL_GetPointer(Import, V_cube, "V", _RC)
    call MAPL_GetPointer(Import, W_cube, "OMEGA", _RC)
    call MAPL_GetPointer(Import, P_cube, "PL", _RC)
    call MAPL_GetPointer(Import, PLE_cube, "PLE", _RC)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    lm = size(GigaTrajInternalPtr%vdata%levs)
    d1 = size(u_cube,1)
    d2 = size(u_cube,2)
    allocate(U_inter(d1, d2,lm), source = 0.0)
    allocate(V_inter(d1, d2,lm), source = 0.0)
    allocate(W_inter(d1, d2,lm), source = 0.0)
    allocate(P_inter(d1, d2,lm), source = 0.0)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       grid_ = GigaTrajInternalPtr%LatLonGrid
       call MAPL_GridGet(grid_, localCellCountPerDim=counts, &
                      globalCellCountPerDim=DIMS, _RC)

       allocate(U_latlon(counts(1), counts(2),lm), source = 0.0)
       allocate(V_latlon(counts(1), counts(2),lm), source = 0.0)
       allocate(W_latlon(counts(1), counts(2),lm), source = 0.0)
       allocate(P_latlon(counts(1), counts(2),lm), source = 0.0)
    else
       grid_ = GigaTrajInternalPtr%CubedGrid
       call MAPL_GridGet(grid_, localCellCountPerDim=counts, &
                      globalCellCountPerDim=DIMS, _RC)
    endif

    allocate(haloU(counts(1)+2, counts(2)+2,lm), source = 0.0)
    allocate(haloV(counts(1)+2, counts(2)+2,lm), source = 0.0)
    allocate(haloW(counts(1)+2, counts(2)+2,lm), source = 0.0)
    allocate(haloP(counts(1)+2, counts(2)+2,lm), source = 0.0)


    call ESMF_StateGet(import, 'PLE', field=GigaTrajInternalPtr%vdata%interp_var, rc=status)
    call GigaTrajInternalPtr%vdata%setup_eta_to_pressure(_RC)

    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(U_cube, U_inter,rc=status)
    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(V_cube, V_inter,rc=status)
    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(W_cube, W_inter,rc=status)
    call  GigaTrajInternalPtr%vdata%regrid_eta_to_pressure(P_cube, P_inter,rc=status)

    if (GigaTrajInternalPtr%regrid_to_latlon) then

       call GigaTrajInternalPtr%cube2latlon%regrid(U_inter,V_inter, U_latlon, V_latlon, _RC)
       call GigaTrajInternalPtr%cube2latlon%regrid(W_inter, W_latlon, _RC)
       call GigaTrajInternalPtr%cube2latlon%regrid(P_inter, P_latlon, _RC)
    endif

!---------------
! Step 2) Get halo
!---------------

    if (GigaTrajInternalPtr%regrid_to_latlon) then
      call esmf_halo(grid_, U_Latlon, haloU, _RC)
      call esmf_halo(grid_, V_Latlon, haloV, _RC)
      call esmf_halo(grid_, W_Latlon, haloW, _RC)
      call esmf_halo(grid_, P_Latlon, haloP, _RC)
    else
      call esmf_halo(grid_, U_inter, haloU, _RC)
      call esmf_halo(grid_, V_inter, haloV, _RC)
      call esmf_halo(grid_, W_inter, haloW, _RC)
      call esmf_halo(grid_, P_inter, haloP, _RC)

    endif

!---------------
! Step 3) Update
!---------------
   
    call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(haloU), c_loc(haloV), c_loc(haloW), c_loc(haloP))

!---------------
! Step 3) Time advance 
!---------------
    call rk4a_advance( GigaTrajInternalPtr%metSrc, c_loc(ctime0), DT, GigaTrajInternalPtr%parcels%num_parcels,  &
                       c_loc(GigaTrajInternalPtr%parcels%lons), &
                       c_loc(GigaTrajInternalPtr%parcels%lats), &
                       c_loc(GigaTrajInternalPtr%parcels%zs))
 
!---------------
! Step 4) Update internal 
!---------------
    call MAPL_GetPointer(INTERNAL, U_internal, "U", _RC)
    call MAPL_GetPointer(INTERNAL, V_internal, "V", _RC)
    call MAPL_GetPointer(INTERNAL, W_internal, "OMEGA", _RC)
    call MAPL_GetPointer(INTERNAL, P_internal, "PL", _RC)
    call MAPL_GetPointer(INTERNAL, PLE_internal, "PLE", _RC)
 
    U_internal = U_cube
    V_internal = V_cube
    W_internal = W_cube
    P_internal = P_cube
    PLE_internal = PLE_cube
   
    if (allocated(U_Latlon)) deallocate( U_Latlon, V_latlon, W_latlon, P_latlon) 
    deallocate(haloU, haloV, haloW, haloP)
  
!---------------
! Step 5) rebalance parcels among processors ( configurable with alarm) 
!---------------
    call rebalance_parcels(clock, GigaTrajInternalPtr%parcels, GigaTrajInternalPtr%CellToRank, comm, grid_, _RC)
!---------------
! Step 6) write out parcel positions and related fields ( configurable with alarm) 
!---------------

    call write_parcels(GC, import, clock, currentTime, _RC)

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

  subroutine mapl_halo(grid, U, V, W, P,  haloU, haloV, haloW, haloP, rc)
    type(ESMF_Grid), intent(in) :: grid
    real, dimension(:,:,:), intent(in) :: U, V, W, P
    real, dimension(:,:,:), intent(inout) :: haloU, haloV, haloW, haloP
    integer, optional,   intent(  out) :: RC

    character(len=ESMF_MAXSTR)              :: IAm
    integer :: counts(3), dims(3), k
    integer :: status
    class(AbstractGridFactory), pointer  :: factory
    
    Iam = "Gigatraj Halo"
    call MAPL_GridGet(grid, localCellCountPerDim=counts, &
                      globalCellCountPerDim=DIMS, _RC)

    haloU(2:counts(1)+1, 2:counts(2)+1, :) = U
    haloV(2:counts(1)+1, 2:counts(2)+1, :) = V
    haloW(2:counts(1)+1, 2:counts(2)+1, :) = W
    haloP(2:counts(1)+1, 2:counts(2)+1, :) = P

    factory =>grid_manager%get_factory(grid)
    do k =1, counts(3)
      call factory%halo(haloU(:,:,k), halo_width=1, _RC)
      call factory%halo(haloV(:,:,k), halo_width=1, _RC)
      call factory%halo(haloW(:,:,k), halo_width=1, _RC)
      call factory%halo(haloP(:,:,k), halo_width=1, _RC)
    enddo

    RETURN_(ESMF_SUCCESS)
  end subroutine

  subroutine esmf_halo(grid, Field,haloField, rc)
    type(ESMF_Grid), intent(in) :: grid
    real, dimension(:,:,:), intent(in) :: Field
    real, dimension(:,:,:), intent(inout) :: haloField
    integer, optional,   intent(  out) :: RC

    character(len=ESMF_MAXSTR)              :: IAm
    integer :: counts(3), k, count3
    integer :: status
    type(ESMF_Field)   :: halo_field
    type(ESMF_RouteHandle) :: rh
    real, dimension(:,:,:), pointer  :: with_halo

    Iam = "Gigatraj ESMF Halo"
    call MAPL_GridGet(grid, localCellCountPerDim=counts, _RC)

    count3 = size(field,3) ! may be nbins

    halo_field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R4, name='halo_field',  &
                                ungriddedLBound=[1],ungriddedUBound=[count3], &
                                totalLWidth=[1,1],totalUWidth=[1,1])
    call ESMF_FieldHaloStore(halo_field, rh, _RC)

    call ESMF_FieldGet(halo_field, farrayPtr=with_halo, _RC)
    !
    ! W.Y note, the pointer with_halo's lbound is 0
    !
    with_halo(1:counts(1), 1:counts(2), :) = Field
    call ESMF_FieldHalo(halo_field, rh, _RC)
    haloField = with_halo

    call ESMF_FieldDestroy( halo_field)

    RETURN_(ESMF_SUCCESS)
  end subroutine esmf_halo

  ! move the parcels to the PE where they belong to
  subroutine rebalance_parcels(clock, parcels, CellToRank, comm, grid, rc)
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    type(horde), intent(inout) :: parcels
    integer, dimension(:,:), intent(in) :: CellToRank
    integer :: comm
    type(ESMF_Grid), intent(inout) :: grid
    integer, optional, intent(out) :: rc

    integer :: status
    integer :: DIMS(3)
    character(len=:), allocatable :: Iam
    integer :: num_parcels0, num_parcels
    real, dimension(:), allocatable    :: lons0, lats0, zs0
    integer, dimension(:), allocatable :: IDs0

    integer :: i, npes, ierror, rank, my_rank, pos
    real :: dlon, dlat

    real, allocatable :: lons_send(:), lats_send(:), zs_send(:)
    integer, allocatable :: ids_send(:)
    type(ESMF_Alarm)  :: GigaTrajRebalanceAlarm
    integer, allocatable :: counts_send(:),counts_recv(:), II(:), JJ(:), ranks(:)
    integer, allocatable :: disp_send(:), disp_recv(:), tmp_position(:)

    Iam = "rebalance_parcels"

    call ESMF_ClockGetAlarm(clock, 'GigatrajRebalance', GigaTrajRebalanceAlarm, _RC)

    if ( .not. ESMF_AlarmIsRinging(GigaTrajRebalanceAlarm)) then
       RETURN_(ESMF_SUCCESS)
    endif

    call move_alloc( parcels%lons, lons0)
    call move_alloc( parcels%lats, lats0)
    call move_alloc( parcels%zs,   zs0)
    call move_alloc( parcels%IDs,  IDs0)
    num_parcels0 = parcels%num_parcels

    where (lons0 < 0)  lons0 =lons0 + 360.0
    where (lons0 >360) lons0 =lons0 - 360.0

    allocate(II(num_parcels0), JJ(num_parcels0))
    call MAPL_GridGet(Grid, globalCellCountPerDim=DIMS)
    if (DIMS(2) == 6*DIMS(1)) then
       call MAPL_GetGlobalHorzIJIndex(num_parcels0, II, JJ, lons0/180.0*MAPL_PI, lats0/180.0*MAPL_PI, Grid=Grid, rc=status)
    else

       dlon = 360.0 / DIMS(1)
       dlat = 180.0 / (DIMS(2)-1)
       ! DC
       II = min( max(ceiling ((lons0+dlon/2.)/dlon),1), DIMS(1))
       JJ = min( max(ceiling ((lats0+dlat/2.+ 90.0)/dlat),1), DIMS(2))

    endif

    call MPI_Comm_size(comm, npes, ierror)
    call MPI_Comm_rank(comm, my_rank, ierror)
    
    allocate(ranks    (num_parcels0))
    allocate(lons_send(num_parcels0))
    allocate(lats_send(num_parcels0))
    allocate(zs_send  (num_parcels0))
    allocate(IDs_send (num_parcels0))

    allocate(counts_send(npes))
    allocate(counts_recv(npes))
    allocate(disp_send(npes))
    allocate(disp_recv(npes))

    do i = 1, num_parcels0
       ranks(i) = CellToRank(II(i), JJ(i))
    enddo

    do rank = 0, npes-1
       counts_send(rank+1) = count(ranks == rank)
    enddo

    call MPI_AllToALL(counts_send, 1, MPI_INTEGER, counts_recv, 1, MPI_INTEGER, comm, ierror)

    disp_send = 0
    do rank = 1, npes-1
       disp_send(rank+1) = disp_send(rank)+ counts_send(rank)
    enddo
    disp_recv = 0
    do rank = 1, npes-1
       disp_recv(rank+1) = disp_recv(rank)+ counts_recv(rank)
    enddo

    ! re-arranged lats lons, and ids
    tmp_position = disp_send
    parcels%num_parcels  = sum(counts_recv)
    num_parcels =  parcels%num_parcels 
    allocate(parcels%lons(num_parcels ))
    allocate(parcels%lats(num_parcels ))
    allocate(parcels%zs  (num_parcels ))
    allocate(parcels%IDs (num_parcels ))

    do i = 1, num_parcels0
       rank   = ranks(i)
       pos    = tmp_position(rank+1) +1
       lons_send(pos) = lons0(i)
       lats_send(pos) = lats0(i)
       zs_send(pos)   =   zs0(i)
       IDs_send(pos)  =  IDs0(i)
       tmp_position(rank+1) = tmp_position(rank+1) + 1
    enddo

    call MPI_AllToALLv(lons_send, counts_send, disp_send, MPI_REAL, parcels%lons, counts_recv, disp_recv, MPI_REAL, comm, ierror)
    call MPI_AllToALLv(lats_send, counts_send, disp_send, MPI_REAL, parcels%lats, counts_recv, disp_recv, MPI_REAL, comm, ierror)
    call MPI_AllToALLv(zs_send,   counts_send, disp_send, MPI_REAL, parcels%zs,   counts_recv, disp_recv, MPI_REAL, comm, ierror)
    call MPI_AllToALLv(ids_send,  counts_send, disp_send, MPI_INTEGER, parcels%IDs,  counts_recv, disp_recv, MPI_INTEGER, comm, ierror)

     RETURN_(ESMF_SUCCESS)
  end subroutine rebalance_parcels

  ! Scatter parcels from root after reading parcels file
  subroutine scatter_parcels(num_parcels0, lons0, lats0, zs0, IDs0, CellToRank, Grid, comm, lons, lats, zs, IDs, num_parcels)
    integer :: num_parcels0
    real, dimension(:), intent(inout) :: lons0
    real, dimension(:), intent(in)    :: lats0, zs0
    integer, dimension(:), intent(in) :: IDs0
    integer, dimension(:,:), intent(in) :: CellToRank
    type(ESMF_GRID), intent(inout) :: Grid 
    integer, intent(in) :: comm
    real, dimension(:), allocatable, intent(out) :: lons, lats, zs
    integer, dimension(:), allocatable, intent(out) :: IDs
    integer, intent(out) :: num_parcels

    integer :: DIMS(3)
    integer :: i, npes, ierror, rank, my_rank, counts_recv, pos, status
    real :: dlon, dlat

    real, allocatable :: lons_send(:), lats_send(:), zs_send(:)
    integer, allocatable :: ids_send(:)

    integer, allocatable :: counts_send(:), II(:), JJ(:), ranks(:)
    integer, allocatable :: disp_send(:), tmp_position(:)

    call MPI_Comm_size(comm, npes, ierror)
    call MPI_Comm_rank(comm, my_rank, ierror)

    call MAPL_GridGet(Grid, globalCellCountPerDim=DIMS)  
  
    allocate(II(num_parcels0), JJ(num_parcels0))

    allocate(counts_send(npes), source = 0)
    allocate(disp_send(npes), source = 0)
    if (my_rank == 0) then
       if (DIMS(2) == 6*DIMS(1)) then
         call MAPL_GetGlobalHorzIJIndex(num_parcels0, II, JJ, lons0/180.0*MAPL_PI, lats0/180.0*MAPL_PI, Grid=Grid, rc=status)
       else
          dlon = 360.0 / DIMS(1)
          dlat = 180.0 / (DIMS(2)-1) !PC

          where (lons0 < 0)   lons0 =lons0 + 360.0
          where (lons0 > 360) lons0 =lons0 - 360.0
          II = min( max(ceiling ((lons0+dlon/2.0)         /dlon),1), DIMS(1))
          JJ = min( max(ceiling ((lats0+90.0+dlat/2.0)/dlat),1), DIMS(2))
       endif

       allocate(ranks(num_parcels0))
       do i = 1, num_parcels0
          ranks(i) = CellToRank(II(i), JJ(i))
       enddo

       do rank = 0, npes-1
          counts_send(rank+1) = count(ranks == rank)
       enddo

       do rank = 1, npes-1
          disp_send(rank+1) = disp_send(rank)+ counts_send(rank)
       enddo
    endif

    call MPI_Scatter(counts_send, 1, MPI_INTEGER, counts_recv, 1, MPI_INTEGER, 0, comm, ierror)

    ! re-arranged lats lons, and ids
    tmp_position = disp_send
    num_parcels  = counts_recv

    allocate(lons_send(num_parcels0))
    allocate(lons     (num_parcels ))
    allocate(lats_send(num_parcels0))
    allocate(lats     (num_parcels ))
    allocate(zs_send  (num_parcels0))
    allocate(zs       (num_parcels ))
    allocate(IDs_send (num_parcels0))
    allocate(IDs      (num_parcels ))

    do i = 1, num_parcels0
       rank   = ranks(i)
       pos = tmp_position(rank+1) +1
       lons_send(pos) = lons0(i)
       lats_send(pos) = lats0(i)
       zs_send(pos)   = zs0(i)
       IDs_send(pos)  = IDs0(i)
       tmp_position(rank+1) = tmp_position(rank+1) + 1
    enddo

    call MPI_ScatterV(lons_send, counts_send, disp_send, MPI_REAL,    lons, counts_recv, MPI_REAL,   0, comm, ierror)
    call MPI_ScatterV(lats_send, counts_send, disp_send, MPI_REAL,    lats, counts_recv, MPI_REAL,   0, comm, ierror)
    call MPI_ScatterV(zs_send,   counts_send, disp_send, MPI_REAL,    zs,   counts_recv, MPI_REAL,   0, comm, ierror)
    call MPI_ScatterV(ids_send,  counts_send, disp_send, MPI_INTEGER, IDs,  counts_recv, MPI_INTEGER,0, comm, ierror)

  end subroutine scatter_parcels

  ! gather parcels to root for writing
  subroutine gather_parcels(num_parcels0, lons0, lats0, zs0, IDs0, comm, lons, lats, zs, IDs, num_parcels)
    integer, intent(out) :: num_parcels0
    real, dimension(:), allocatable, intent(out) :: lons0, lats0,zs0
    integer, dimension(:), allocatable, intent(out) :: IDs0
    integer, intent(in) :: comm
    real, dimension(:), intent(in) :: lons, lats, zs
    integer, dimension(:), intent(in) :: IDs
    integer, intent(in) :: num_parcels

    integer :: i, npes, ierror, my_rank
    integer, allocatable :: nums_all(:), displ(:)

    call MPI_Comm_size(comm, npes, ierror)
    call MPI_Comm_rank(comm, my_rank, ierror)

    allocate(nums_all(npes), source = 0)
    call MPI_Gather(num_parcels, 1, MPI_INTEGER, nums_all, 1, MPI_INTEGER, 0, comm, ierror)

    num_parcels0 = sum(nums_all)

    allocate(lons0(num_parcels0))
    allocate(lats0(num_parcels0))
    allocate(  zs0(num_parcels0))
    allocate( IDS0(num_parcels0))
    allocate(  displ(npes), source =0)
    do i =2, npes
       displ(i) = displ(i-1)+nums_all(i-1)
    enddo

    call MPI_GatherV(lons, num_parcels, MPI_REAL, lons0,  nums_all, displ, MPI_REAL, 0, comm,ierror)
    call MPI_GatherV(lats, num_parcels, MPI_REAL, lats0,  nums_all, displ, MPI_REAL, 0, comm,ierror)
    call MPI_GatherV(zs,   num_parcels, MPI_REAL, zs0,    nums_all, displ, MPI_REAL, 0, comm,ierror)
    call MPI_GatherV(IDS,  num_parcels, MPI_INTEGER, IDs0,nums_all, displ, MPI_INTEGER, 0, comm,ierror)

  end subroutine gather_parcels

  subroutine gather_onefield(num_parcels0, field0, comm, field, num_parcels)
    integer, intent(out) :: num_parcels0
    real, dimension(:), allocatable, intent(out) :: field0
    integer, intent(in) :: comm
    real, dimension(:), intent(in) :: field
    integer, intent(in) :: num_parcels

    integer :: i, npes, ierror, my_rank
    integer, allocatable :: nums_all(:), displ(:)

    call MPI_Comm_size(comm, npes, ierror)
    call MPI_Comm_rank(comm, my_rank, ierror)

    allocate(nums_all(npes), source = 0)
    call MPI_Gather(num_parcels, 1, MPI_INTEGER, nums_all, 1, MPI_INTEGER, 0, comm, ierror)

    num_parcels0 = sum(nums_all)

    allocate(field0(num_parcels0))
    allocate(  displ(npes), source =0)
    do i =2, npes
       displ(i) = displ(i-1)+nums_all(i-1)
    enddo

    call MPI_GatherV(field, num_parcels, MPI_REAL, field0,  nums_all, displ, MPI_REAL, 0, comm,ierror)

  end subroutine gather_onefield

  subroutine write_parcels(GC, state, CLOCK, currentTime, rc)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state  ! Import state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    type(ESMF_TIME), intent(in) :: currentTime 
    integer, optional, intent(out) :: rc
  
    character(len=:), allocatable :: Iam
    type (ESMF_VM)  :: vm
    type(Netcdf4_fileformatter) :: formatter
    integer :: comm, my_rank, total_num, i, k,status, last_time, ierror, count3
    real, allocatable :: lats0(:), lons0(:), zs0(:), ids0_r(:), values(:), values0(:)
    real,target, allocatable ::  values_2d(:,:)
    real,pointer ::  field(:,:,:)
    integer, allocatable :: ids0(:)
    type(ESMF_Alarm)  :: GigaTrajOutAlarm
    type(FileMetadata) :: meta
    real(ESMF_KIND_R8) :: tint_d
    type(ESMF_TimeInterval) :: tint
    type(MAPL_MetaComp),pointer  :: MPL
    character(len=ESMF_MAXSTR)   :: parcels_file, other_fields
    character(len=ESMF_MAXSTR), allocatable   :: varnames(:)
    character(len=:), allocatable:: var_name, var_
    
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)      :: wrap
    real :: lon_start
    character(len=20), target :: ctime

    Iam = "write_parcels"
    call ESMF_ClockGetAlarm(clock, 'GigatrajOut', GigaTrajOutAlarm, _RC)

    if ( .not. ESMF_AlarmIsRinging(GigaTrajOutAlarm)) then
       RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_GetObjectFromGC ( GC, MPL, _RC)
    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC)
    call MAPL_GetResource(MPL, lon_start   , "GIGATRAJ_LON_START:",    default= 0. , _RC)

    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)

    call gather_parcels(total_num, lons0, lats0, zs0, IDs0, &
                          comm,                             & 
                          GigaTrajInternalPtr%parcels%lons, &
                          GigaTrajInternalPtr%parcels%lats, &
                          GigaTrajInternalPtr%parcels%zs,   &
                          GigaTrajInternalPtr%parcels%IDS,  &
                          GigaTrajInternalPtr%parcels%num_parcels )

    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)
    if (my_rank ==0) then
       ! reorder 
       lats0 = lats0(ids0(:)+1) ! id is zero-bases, plus 1 Fortran
       lons0 = lons0(ids0(:)+1)
       !lons0 = lons0 + lon_start 
       if (lon_start < 0.) then
         where (lons0 > 180.0) lons0 = lons0 - 360.
       endif
       zs0   = zs0(ids0(:)+1)
       call formatter%open(trim(parcels_file), pFIO_WRITE, _RC)
       meta = formatter%read(_RC)
       last_time = meta%get_dimension('time', _RC)
       tint = CurrentTime - GigaTrajInternalPtr%startTime
       call ESMF_TimeIntervalGet(tint,d_r8=tint_d,rc=status)

       call formatter%put_var('lat', lats0, start=[1, last_time+1], _RC)
       call formatter%put_var('lon', lons0, start=[1, last_time+1], _RC)
       call formatter%put_var('P',   zs0,   start=[1, last_time+1], _RC)
       call formatter%put_var('time',  [tint_d],  start=[last_time+1], _RC)
     endif

     ! extra fields
     if (allocated(GigaTrajInternalPtr%ExtraFieldNames)) then
       call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime)
       ctime(20:20) = c_null_char
       do k = 1, size(GigaTrajInternalPtr%ExtraFieldNames)
         select case (trim(GigaTrajInternalPtr%ExtraFieldNames(k)))
         case('TH')
            var_name = 'theta'
         case('T')
            var_name = 't'
         case('PALT')
            var_name = 'palt'
         case default
            var_name = trim(GigaTrajInternalPtr%ExtraFieldNames(k))
         end select
         
         if ( index(var_name, 'bcDP') /= 0 .or.   &
              index(var_name, 'ocDP') /= 0 .or.   &
              index(var_name, 'bcWT') /= 0 .or.   &
              index(var_name, 'ocWT') /= 0 .or.   &
              index(var_name, 'bcSD') /= 0 .or.   &
              index(var_name, 'ocSD') /= 0 .or.   &
              index(var_name, 'bcSV') /= 0 .or.   &
              index(var_name, 'ocSV') /= 0 ) then

            call MAPL_GetPointer(state, field, var_name, _RC)
            count3 = size(field,3)
            allocate(values_2d(GigaTrajInternalPtr%parcels%num_parcels, count3))
            call get_metsrc_data2d (GC, state, ctime, var_name, values_2d, RC )
            do i = 1, count3
               call gather_onefield(total_num, values0, comm, values_2d(:,i), GigaTrajInternalPtr%parcels%num_parcels)

               if (my_rank == 0) then
                 values0 = values0(ids0(:)+1)
                 var_ = var_name //'00'//i_to_string(i)
                 if ( meta%has_variable(var_)) then
                    call formatter%put_var( var_,   values0,   start=[1, last_time+1], _RC)
                 else
                    print*, "Please provide "//var_ // " in the file "//trim(parcels_file)
                 endif
               endif
            enddo
            deallocate(values_2d)
         else if ( index(var_name, 'TRI') /= 0) then
            allocate(values(GigaTrajInternalPtr%parcels%num_parcels))
            allocate(varnames(4))
            varnames(1) = "CA.bc::CA.bcphilicIT"
            varnames(2) = "CA.bc::CA.bcphobicIT"
            varnames(3) = "CA.oc::CA.ocphilicIT"
            varnames(4) = "CA.oc::CA.ocphobicIT"
            do i = 1, 4
              var_ = varnames(i)(8:)
              call get_metsrc_data (GC, state, ctime, varnames(i), values, _RC)
              call gather_onefield(total_num, values0, comm, values, GigaTrajInternalPtr%parcels%num_parcels)

              if (my_rank == 0) then
                values0 = values0(ids0(:)+1)
                if ( meta%has_variable(var_)) then
                  call formatter%put_var( var_,   values0,   start=[1, last_time+1], _RC)
                else
                  print*, "Please provide "//var_ // " in the file "//trim(parcels_file)
                endif
              endif
            enddo
            deallocate(values, varnames)
         else
            allocate(values(GigaTrajInternalPtr%parcels%num_parcels))
            call get_metsrc_data (GC, state, ctime, var_name, values, _RC )
            call gather_onefield(total_num, values0, comm, values, GigaTrajInternalPtr%parcels%num_parcels)

            if (my_rank == 0) then
              values0 = values0(ids0(:)+1)
              if ( meta%has_variable(var_name)) then
                call formatter%put_var( var_name,   values0,   start=[1, last_time+1], _RC)
              else
                print*, "Please provide "//var_name // " in the file "//trim(parcels_file)
              endif
            endif
            deallocate(values)
         endif
       enddo
     endif

     if (my_rank ==0) then 
       call formatter%close(_RC)
     endif
     RETURN_(ESMF_SUCCESS)
     contains

  end subroutine write_parcels

  subroutine read_parcels(GC,internal, rc)
     type(ESMF_GridComp), intent(inout) :: GC      
     type(GigaTrajInternal), intent(inout) :: internal
     integer, optional, intent(out) :: rc

     type(Netcdf4_fileformatter) :: formatter
     type(FileMetadata) :: meta 
     integer :: comm, my_rank, total_num, ierror, last_time, DIMS(3)
     real, allocatable :: lats(:), lons(:), zs(:), lats0(:), lons0(:), zs0(:)
     real(kind=ESMF_KIND_R8), allocatable :: ids0_r(:)
     integer, allocatable :: ids0(:)
     integer :: status
     character(len=ESMF_MAXSTR) :: parcels_file
     character(len=ESMF_MAXSTR) :: regrid_to_latlon
     type(MAPL_MetaComp),pointer  :: MPL
     type (ESMF_VM)   :: vm
     class(Variable), pointer :: t
     type(Attribute), pointer :: attr
     real :: lon_start
     class(*), pointer :: units
     character(len=ESMF_MAXSTR) :: Iam ="read_parcels"

     call ESMF_VMGetCurrent(vm, _RC)
     call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
     
     call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)
     call MAPL_GetObjectFromGC ( GC, MPL, _RC)
     call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC)
     call MAPL_GetResource(MPL, lon_start   , "GIGATRAJ_LON_START:",  default= 0. , _RC)
     total_num = 0
     if (my_rank ==0) then
        call formatter%open(parcels_file, pFIO_READ, _RC)
        meta = formatter%read(_RC)
        total_num = meta%get_dimension('id', _RC)
        last_time = meta%get_dimension('time', _RC)
        t => meta%get_variable('time', _RC)
        attr => t%get_attribute('long_name')
        units => attr%get_value()
        select type(units)
        type is (character(*))
           internal%startTime = parse_time_string(units, _RC)
        class default
           _FAIL('unsupported subclass for units')
        end select
     endif

     allocate(lats0(total_num), lons0(total_num), zs0(total_num), ids0_r(total_num))

     if (my_rank ==0) then
        call formatter%get_var('lat', lats0, start = [1,last_time], _RC)
        call formatter%get_var('lon', lons0, start = [1,last_time], _RC)
        call formatter%get_var('P',   zs0,   start = [1,last_time], _RC)
        call formatter%get_var('id',  ids0_r,start = [1,last_time], _RC)
        call formatter%close(_RC)
        ids0 = int(ids0_r)
        !lons0 = lons0 - lon_start
        where (lons0<0) lons0 = lons0 + 360.
     endif
     call MAPL_GetResource(MPL, regrid_to_latlon, "GIGATRAJ_REGRID_TO_LATLON:",  default= 'YES' , _RC)
     if (trim (regrid_to_latlon) == 'YES') then
        call scatter_parcels(total_num, lons0, lats0, zs0, IDs0, internal%CellToRank, internal%LatLonGrid, comm, &
                              Internal%parcels%lons, &
                              Internal%parcels%lats, &
                              Internal%parcels%zs,   &
                              Internal%parcels%IDS,  & 
                              Internal%parcels%num_parcels) 
     else
        call scatter_parcels(total_num, lons0, lats0, zs0, IDs0, internal%CellToRank, internal%CubedGrid, comm, &
                              Internal%parcels%lons, &
                              Internal%parcels%lats, &
                              Internal%parcels%zs,   &
                              Internal%parcels%IDS,  & 
                              Internal%parcels%num_parcels) 
     endif

     deallocate(lats0, lons0, zs0, ids0_r)
     RETURN_(ESMF_SUCCESS)
     contains
          ! a copy from MAPL_TimeMod
          function parse_time_string(timeUnits,rc) result(time)
             character(len=*), intent(inout) :: timeUnits
             integer, optional, intent(out) :: rc
         
             type(ESMF_Time) :: time
             integer :: status
         
             integer        year               ! 4-digit year
             integer        month              ! month
             integer        day                ! day
             integer        hour               ! hour
             integer        min                ! minute
             integer        sec                ! second
         
             integer ypos(2), mpos(2), dpos(2), hpos(2), spos(2)
             integer strlen
             integer firstdash, lastdash
             integer firstcolon, lastcolon
             integer lastspace
             strlen = LEN_TRIM (TimeUnits)
         
             firstdash = index(TimeUnits, '-')
             lastdash  = index(TimeUnits, '-', BACK=.TRUE.)
             if (firstdash .LE. 0 .OR. lastdash .LE. 0) then
                _FAIL('time string is not a valid format')
             endif
             ypos(2) = firstdash - 1
             mpos(1) = firstdash + 1
             ypos(1) = ypos(2) - 3
         
             mpos(2) = lastdash - 1
             dpos(1) = lastdash + 1
             dpos(2) = dpos(1) + 1
         
             read ( TimeUnits(ypos(1):ypos(2)), * ) year
             read ( TimeUnits(mpos(1):mpos(2)), * ) month
             read ( TimeUnits(dpos(1):dpos(2)), * ) day
         
             firstcolon = index(TimeUnits, ':')
             if (firstcolon .LE. 0) then
         
                ! If no colons, check for hour.
         
                ! Logic below assumes a null character or something else is after the hour
                ! if we do not find a null character add one so that it correctly parses time
                if (TimeUnits(strlen:strlen) /= char(0)) then
                   TimeUnits = trim(TimeUnits)//char(0)
                   strlen=len_trim(TimeUnits)
                endif
                lastspace = index(TRIM(TimeUnits), ' ', BACK=.TRUE.)
                if ((strlen-lastspace).eq.2 .or. (strlen-lastspace).eq.3) then
                   hpos(1) = lastspace+1
                   hpos(2) = strlen-1
                   read (TimeUnits(hpos(1):hpos(2)), * ) hour
                   min  = 0
                   sec  = 0
                else
                   hour = 0
                   min  = 0
                   sec  = 0
                endif
             else
                hpos(1) = firstcolon - 2
                hpos(2) = firstcolon - 1
                lastcolon =  index(TimeUnits, ':', BACK=.TRUE.)
                if ( lastcolon .EQ. firstcolon ) then
                   mpos(1) = firstcolon + 1
                   mpos(2) = firstcolon + 2
                   read (TimeUnits(hpos(1):hpos(2)), * ) hour
                   read (TimeUnits(mpos(1):mpos(2)), * ) min
                   sec = 0
                else
                   mpos(1) = firstcolon + 1
                   mpos(2) = lastcolon - 1
                   spos(1) = lastcolon + 1
                   spos(2) = lastcolon + 2
                   read (TimeUnits(hpos(1):hpos(2)), * ) hour
                   read (TimeUnits(mpos(1):mpos(2)), * ) min
                   read (TimeUnits(spos(1):spos(2)), * ) sec
                endif
             endif
         
             call ESMF_TimeSet(time,yy=year,mm=month,dd=day,h=hour,m=min,s=sec,rc=status)
             _VERIFY(status)
             RETURN_(ESMF_SUCCESS)
          end function parse_time_string
  end subroutine read_parcels

  subroutine get_metsrc_data (GC, state, ctime, fieldname, values, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state
    character(*), target,   intent(in) :: ctime
    character(*), intent(in) :: fieldname
    real, target, intent(inout) :: values(:)
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)     :: wrap
    type (ESMF_FieldBundle)         :: TRI
    type (ESMF_Field)               :: field
    real, dimension(:,:,:), pointer     :: ptr3d
    real, dimension(:,:,:), allocatable :: field_latlon
    real, dimension(:,:,:), allocatable, target  :: haloField
    integer :: counts(3), dims(3), d1,d2,km, count3
    character(len=:), target, allocatable :: field_

    Iam = "get_metsrc_data"

    if (index(fieldname,'philicIT') /=0 .or. index(fieldname,'phobicIT') /=0) then
       call ESMF_StateGet(state, 'PHYSICS_Exports/TURBULENCE_Exports/TRI', TRI, _RC)
       call ESMF_FieldBundleGet(TRI, fieldname, field=field, _RC)
       call ESMF_FieldGet(field,farrayPtr=ptr3d, _RC)
    else
       call MAPL_GetPointer(state, ptr3d, fieldname, _RC)
    endif

    count3 = size(ptr3d,3)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr
    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

    allocate(field_latlon(counts(1),counts(2),count3))
    allocate(haloField(counts(1)+2, counts(2)+2,count3), source = 0.0)

    call GigaTrajInternalPtr%cube2latlon%regrid(ptr3d, Field_latlon, _RC)

    call esmf_halo(GigaTrajInternalPtr%LatLonGrid, field_Latlon, haloField, _RC)

    field_ = trim(fieldname)//c_null_char
    call setData( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(field_), c_loc(haloField))
    call getData(GigaTrajInternalPtr%metSrc,  c_loc(ctime),  c_loc(field_),  &
              GigaTrajInternalPtr%parcels%num_parcels,  &
              c_loc(GigaTrajInternalPtr%parcels%lons),         &
              c_loc(GigaTrajInternalPtr%parcels%lats),         &
              c_loc(GigaTrajInternalPtr%parcels%zs),           &
              c_loc(values))

    deallocate(field_latlon, haloField)
    RETURN_(ESMF_SUCCESS)

  end subroutine get_metsrc_data

  subroutine get_metsrc_data2d (GC, state, ctime, fieldname, values, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state
    character(*), target,   intent(in) :: ctime
    character(*), intent(in) :: fieldname
    real, target, intent(inout) :: values(:,:)
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)     :: wrap
    real, dimension(:,:,:), pointer     :: field
    real, dimension(:,:,:), allocatable :: field_latlon
    real, dimension(:,:,:), allocatable, target  :: haloField
    integer :: counts(3), i, count3
    character(len=:), target, allocatable :: field_

    Iam = "get_metsrc_data"

    call MAPL_GetPointer(state, field, fieldname, _RC)
    count3 = size(field,3)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr
    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts, _RC)

    allocate(field_latlon(counts(1),counts(2),count3))
    allocate(haloField(counts(1)+2, counts(2)+2,count3), source = 0.0)

    call GigaTrajInternalPtr%cube2latlon%regrid(field, Field_latlon, _RC)

    call esmf_halo(GigaTrajInternalPtr%LatLonGrid, field_Latlon, haloField, _RC)

    field_ = trim(fieldname)//'_2D'//c_null_char

    do i = 1,count3

      call setData( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(field_), c_loc(haloField(1,1,i)))
      call getData2d(GigaTrajInternalPtr%metSrc,  c_loc(ctime),  c_loc(field_),  &
                 GigaTrajInternalPtr%parcels%num_parcels,  &
                 c_loc(GigaTrajInternalPtr%parcels%lons),         &
                 c_loc(GigaTrajInternalPtr%parcels%lats),         &
                 c_loc(values(1,i)))

    enddo
    deallocate(field_latlon, haloField)
    RETURN_(ESMF_SUCCESS)

  end subroutine get_metsrc_data2d
 
end module GEOS_GigatrajGridCompMod 
