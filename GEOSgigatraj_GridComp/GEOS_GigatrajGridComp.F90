#include "MAPL_Generic.h"

module GEOS_GigatrajGridCompMod
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated, c_null_char
  use, intrinsic :: iso_c_binding, only : c_loc
  use ESMF
  use MAPL
  use MAPL_VerticalDataMod
  use mpi
  use GEOS_Giga_interOpMod
  use Gigatraj_UtilsMod
  implicit none

  public :: SetServices

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
         SHORT_NAME = 'ZL',                                        &
         LONG_NAME  = 'mid_layer_heights',                         &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter, _RC)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                     & 
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter, _RC)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,  _RC)

    call MAPL_AddInternalSpec ( gc,                                  &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K s-1',                                       &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

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
    integer :: DIMS(3)
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
    type (GigatrajInternalWrap)   :: wrap
    type (ESMF_TIME) :: CurrentTime
    type(ESMF_Alarm)  :: GigaTrajOutAlarm, GigaTrajRebalanceAlarm, GigaTrajIntegrateAlarm
    type(ESMF_TimeInterval)      :: parcelsOut_DT, Rebalance_DT, Integrate_DT
    type(ESMF_TimeInterval)      :: ModelTimeStep
    integer :: HH, MM, SS
    integer :: integrate_time, r_time, o_time
    character(len=ESMF_MAXSTR) :: parcels_file
    character(len=ESMF_MAXSTR) :: grid_name, vCoord
    character(len=ESMF_MAXSTR) :: regrid_to_latlon
    character(len=ESMF_MAXSTR), allocatable  :: cName(:), bName(:), fName(:), aName(:)
    type(ESMF_Grid) :: grid_

    call ESMF_GridCompGet ( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // "Initialize"

    call MAPL_GetObjectFromGC ( GC, MPL, _RC)

    call MAPL_TimerOn(MPL,"TOTAL")
    call MAPL_TimerOn(MPL,"INITIALIZE")

    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, _RC)

    call ESMF_TimeIntervalGet(ModelTimeStep, h = hh, m = mm, s = ss, _RC)
    
    call MAPL_GetResource(MPL, integrate_time, "GIGATRAJ_INTEGRATE_DT:", default = hh*10000+mm*100+ss, _RC) 
    hh = integrate_time/10000
    mm = mod(integrate_time, 10000)/100
    ss = mod(integrate_time, 100)
    call ESMF_TimeIntervalSet(Integrate_DT,  h = hh, m = mm, s = ss, _RC)

    call MAPL_GetResource(MPL, r_time, "GIGATRAJ_REBALANCE_DT:", default = integrate_time, _RC) 
    hh = r_time/10000
    mm = mod(r_time, 10000)/100
    ss = mod(r_time, 100)
    call ESMF_TimeIntervalSet(Rebalance_DT, h = hh, m = mm, s = ss, _RC)

    call MAPL_GetResource(MPL, o_time, "GIGATRAJ_OUTPUT_DT:", default = integrate_time, _RC)
    hh = o_time/10000
    mm = mod(o_time, 10000)/100
    ss = mod(o_time, 100)
    call ESMF_TimeIntervalSet(parcelsOut_DT, h = hh, m = mm, s = ss, _RC)

    GigaTrajOutAlarm = ESMF_AlarmCreate(                   &
         clock,                                            &
         name='GigatrajOut',                               &
         ringTime= CurrentTime + parcelsOut_DT-ModelTimeStep, &
         ringInterval=parcelsOut_DT,                          &
         ringTimeStepCount=1,                              &
         sticky=.false., _RC)

    GigaTrajRebalanceAlarm = ESMF_AlarmCreate(             &
         clock,                                            &
         name='GigatrajRebalance',                         &
         ringTime= CurrentTime + Rebalance_DT-ModelTimeStep, &
         ringInterval=Rebalance_DT,                          &
         ringTimeStepCount=1,                              &
         sticky=.false., _RC)

    GigaTrajIntegrateAlarm = ESMF_AlarmCreate(             &
         clock,                                            &
         name='GigatrajIntegrate',                         &
         ringTime= CurrentTime + integrate_DT-ModelTimeStep, &
         ringInterval=integrate_DT,                          &
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

    call MAPL_GetResource(MPL, vCoord, "GIGATRAJ_VERTICAL_COORD:", default='DYN%%PL|P', rc=status)
    call parseCompsAndFieldsName(vCoord, cName, bName, fName, aName)
    GigaTrajInternalPtr%vCoord = trim(fName(1))
    GigaTrajInternalPtr%vAlias = trim(aName(1))
    select case(GigaTrajInternalPtr%vCoord)
    case ('PL')
       GigaTrajInternalPtr%vTendency = 'OMEGA'
    case('TH')
       GigaTrajInternalPtr%vTendency = 'DTDTDYN'
    case('ZL')
       GigaTrajInternalPtr%vTendency = 'W'
    case default
       _ASSERT(.false., "vertical coordinate is needed")
    end select    

    npz = Dims(3)
    GigaTrajInternalPtr%npz = npz
    GigaTrajInternalPtr%Integrate_DT = Integrate_DT

    call MAPL_GetResource(MPL, NX, "NX:", _RC)
    call MAPL_GetResource(MPL, grid_name, "AGCM_GRIDNAME:", _RC)

    ! the level is differtent from the original grid
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

        grid_ = GigaTrajInternalPtr%LatLonGrid
        call MAPL_GridGet(grid_, globalCellCountPerDim=DIMS, _RC)
    else
        grid_  = CubedGrid
    endif
 
    call MAPL_Grid_interior(grid_, i1,i2,j1,j2)

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

   integer :: i, status, k
   character(len=ESMF_MAXSTR) :: IAm 
   type (ESMF_State)          :: INTERNAL, leaf_export
   type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
   type (GigatrajInternalWrap)   :: wrap
   character(len=ESMF_MAXSTR)    :: GigaRstFile
   character(len=ESMF_MAXSTR)    :: other_fields
   type(ESMF_Field) :: tmp_field
   type (ESMF_FieldBundle)       :: bdle
   type (ESMF_TIME) :: CurrentTime
   character(len=20), target :: ctime 
   type (MAPL_MetaComp),  pointer  :: MPL
   logical, save :: init = .false.
   real, dimension(:,:,:) , pointer :: ptr3d
   type(Netcdf4_fileformatter) :: formatter
   type(FileMetadata) :: meta
   character(len=ESMF_MAXSTR) :: parcels_file
   character(len=:), allocatable :: fieldname, tmp_name
   character(len=20) :: diffusions(4)
   character (len=ESMF_MAXSTR), allocatable  :: itemNameList(:)
   character (len=ESMF_MAXSTR) :: LONG_NAME, UNITS
   character (len=ESMF_MAXSTR), allocatable  :: fieldnames(:)
   character (len=ESMF_MAXSTR), allocatable  :: bundlenames(:)
   character (len=ESMF_MAXSTR), allocatable  :: compnames(:)
   character (len=ESMF_MAXSTR), allocatable  :: aliasnames(:)
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
      call init_metsrc_field0(GC,  IMPORT,  ctime,  _RC)
   else
      INQUIRE(FILE= GigaRstFile, EXIST=file_exists) 
      _ASSERT( file_exists, " GIGATRAJ_INTERNAL_RESTART_FILE does not exist")
      call init_metsrc_field0(GC, INTERNAL, ctime,  _RC)
   endif

   call MAPL_GetResource(MPL, other_fields, "GIGATRAJ_EXTRA_FIELDS:", default='NONE', _RC)

   if (other_fields /= 'NONE') then
      call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC)
      if (MAPL_AM_I_ROOT()) then
        call formatter%open(trim(parcels_file), pFIO_WRITE, _RC)
        meta = formatter%read(_RC)
      endif

      call parseCompsAndFieldsName(other_fields, compnames, bundlenames, fieldnames, aliasnames)
      GigaTrajInternalPtr%ExtraCompNames   = compnames
      GigaTrajInternalPtr%ExtraBundleNames = bundlenames
      GigaTrajInternalPtr%ExtraFieldNames  = fieldnames
      GigaTrajInternalPtr%ExtraAliasNames  = aliasnames

      do i = 1, size(FieldNames)
         call MAPL_ExportStateGet([import], trim(compnames(i)), leaf_export, _RC)
         if ( trim(bundlenames(i)) == 'NONE') then
            call ESMF_StateGet(leaf_export, trim(FieldNames(i)), tmp_field, _RC)
         else
            call ESMF_StateGet(leaf_export, trim(bundlenames(i)), bdle, _RC)
            call ESMFL_BundleGetPointerToData(bdle   , trim(FieldNames(i)) , ptr3d, _RC)
            if ( .not. associated(ptr3d)) then
               _ASSERT(.false., trim(FieldNames(i)) // " in bundle "//trim(bundlenames(i)) // " is not allocated, gigatraj cannot output this field")
            endif
            call ESMF_FieldBundleGet(bdle, trim(FieldNames(i)), field=tmp_field, _RC)
         endif
         call MAPL_AllocateCoupling(tmp_field, _RC)
         call ESMF_AttributeGet(tmp_field, NAME='LONG_NAME', VALUE=LONG_NAME, _RC)
         call ESMF_AttributeGet(tmp_field, NAME='UNITS', VALUE=UNITS, _RC)
    
         call create_new_vars( meta, formatter, trim(long_name), trim(aliasnames(i)), trim(units))

     enddo
     if (MAPL_AM_I_Root()) then
        call formatter%close()
     endif
   endif

   init = .true.

   RETURN_(ESMF_SUCCESS)

 end subroutine GetInitVars

 subroutine Init_metsrc_field0 (GC, state, ctime, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state
    character(*), target,   intent(in) :: ctime
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm 
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap
    real, dimension(:,:,:), pointer     :: U, V, W, P, PL0, PLE, TH
    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable, target  :: haloU, haloV, haloW, haloP
    integer :: counts(3), dims(3), d1,d2,km,lm, i1,i2,j1,j2,i
    real, allocatable, target :: lats_center(:), lons_center(:), levs_center(:)
    real, allocatable, target :: cube_lats_center(:, :), cube_lons_center(:,:)
    integer :: comm
    real :: delt, High, low
    type(ESMF_VM) :: vm
    type(ESMF_Grid) :: grid_

    Iam = "init_metsrc_field0"

    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    call MAPL_GetPointer(state, U, "U", _RC)
    call MAPL_GetPointer(state, V, "V", _RC)
    call MAPL_GetPointer(state, W, trim(GigaTrajInternalPtr%vTendency), _RC)
    call MAPL_GetPointer(state, P, trim(GigaTrajInternalPtr%vCoord), _RC)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       grid_ = GigaTrajInternalPtr%LatLonGrid
    else
       grid_ = GigaTrajInternalPtr%CubedGrid
    endif

    call MAPL_GridGet( grid_, localCellCountPerDim=counts, globalCellCountPerDim=dims, _RC)

    select case ( trim(GigaTrajInternalPtr%vCoord))
    case ("PL")
       High = 100000.
       Low  = 2.
    case ("TH")
       High = 5000.
       Low  = 200.
    case ("ZL")
       High = 78000.
       Low  = 1.
    end select

    delt = (log(High)-log(low))/dims(3)
    levs_center=[(exp(log(High)-(i-1)*delt), i=1, dims(3))]

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       call get_latlon_centers(gc, lons_center, lats_center, _RC)
       GigaTrajInternalPtr%metSrc = initMetGEOSDistributedLatLonData(comm, c_loc(GigaTrajInternalPtr%CellToRank), DIMS(1), DIMS(2),   &
                                   dims(3), counts(1)+2, counts(2)+2, dims(3), &
                                   c_loc(lons_center), c_loc(lats_center), c_loc(levs_center), c_loc(ctime))
       deallocate(lons_center, lats_center)
    else
       call MAPL_Grid_interior(grid_, i1, i2, j1, j2)
       call get_cube_centers(gc, cube_lons_center, cube_lats_center, _RC)
       GigaTrajInternalPtr%metSrc = initMetGEOSDistributedCubedData(comm, c_loc(GigaTrajInternalPtr%CellToRank), DIMS(1),   &
                                   dims(3), i1, i2, j1, j2, dims(3), &
                                   c_loc(cube_lons_center), c_loc(cube_lats_center), c_loc(levs_center), c_loc(ctime))

       deallocate(cube_lons_center, cube_lats_center)

    endif
    deallocate(levs_center)

    lm = dims(3)
    d1 = counts(1)
    d2 = counts(2)

    allocate(haloU(d1+2, d2+2,lm), source = 0.0)
    allocate(haloV(d1+2, d2+2,lm), source = 0.0)
    allocate(haloW(d1+2, d2+2,lm), source = 0.0)
    allocate(haloP(d1+2, d2+2,lm), source = 0.0)

    if ( GigaTrajInternalPtr%regrid_to_latlon) then
       allocate(U_latlon(d1,d2,lm))
       allocate(V_latlon(d1,d2,lm))
       allocate(W_latlon(d1,d2,lm))
       allocate(P_latlon(d1,d2,lm))
       call GigaTrajInternalPtr%cube2latlon%regrid(U, V, U_latlon, V_latlon, _RC)
       call GigaTrajInternalPtr%cube2latlon%regrid(W, W_latlon, _RC)
       call GigaTrajInternalPtr%cube2latlon%regrid(P, P_latlon, _RC)

       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, U_Latlon, haloU, _RC)
       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, V_Latlon, haloV, _RC)
       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, W_Latlon, haloW, _RC)
       call esmf_halo(GigaTrajInternalPtr%LatLonGrid, P_Latlon, haloP, _RC)
       deallocate(U_latlon, V_latlon, W_latlon, P_latlon)
    else
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, U, haloU, _RC)
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, V, haloV, _RC)
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, W, haloW, _RC)
       call esmf_halo(GigaTrajInternalPtr%CubedGrid, P, haloP, _RC)
    endif

    call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(haloU), c_loc(haloV), c_loc(haloW), c_loc(haloP))

    if(associated(PL0)) deallocate(PL0)
    deallocate(haloU, haloV, haloW, haloP)
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
    integer        :: CSTAT, ESTAT, YY, DD
    character(512) :: CMSG
    character(256) :: command_line
    character(19)  :: begdate, enddate
    character(64)  :: format_string
    type(ESMF_TimeInterval) :: ModelTimeStep
    type(ESMF_Time)         :: CurrentTime, preTime
    type(ESMF_Grid)         :: grid_
  
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap
  
    integer :: lm, d1, d2, k
    integer ::counts(3), DIMS(3), comm, ierror
    type (ESMF_VM)   :: vm
  
    real, dimension(:,:,:), pointer     :: U_cube, V_cube, W_cube, P_cube, PLE_Cube, with_halo
    real, dimension(:,:,:), pointer     :: internal_field, model_field
    real, dimension(:,:,:), pointer     :: tmp_ptr

    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable :: U_inter, V_inter, W_inter, P_inter
    real, dimension(:,:,:), allocatable, target  :: haloU, haloV, haloW, haloP

    real(ESMF_KIND_R8) :: DT

    character(len=20), target :: ctime, ctime0
    type(ESMF_State)            :: INTERNAL
    type(MAPL_MetaComp),pointer :: MPL  
    type(ESMF_Alarm)  :: GigaTrajIntegrateAlarm
    type(MAPL_VarSpec ), pointer:: internal_specs(:)
    character(len=ESMF_MAXSTR)  :: SHORT_NAME

!---------------
!  Update internal 
!---------------
    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr
    call MAPL_GetObjectFromGC ( GC, MPL, _RC)
    call MAPL_Get (MPL, INTERNAL_ESMF_STATE=INTERNAL, _RC)

    call MAPL_GridCompGetVarSpecs(GC, INTERNAL=internal_specs, _RC)
    do K=1,size(internal_specs)
       call MAPL_VarSpecGet(internal_specs(k), SHORT_NAME=SHORT_NAME, _RC)
       call MAPL_GetPointer(Import, model_field, trim(short_name), _RC)
       call MAPL_GetPointer(INTERNAL, internal_field, trim(short_name), _RC)
       internal_field(:,:,:) = model_field(:,:,:)
    enddo

    call ESMF_ClockGetAlarm(clock, 'GigatrajIntegrate', GigaTrajIntegrateAlarm, _RC)

    if ( .not. ESMF_AlarmIsRinging(GigaTrajIntegrateAlarm)) then
       RETURN_(ESMF_SUCCESS)
    endif

    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)


    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, _RC)
 
    ! W.J note: this run is after agcm's run. The clock is not yet ticked
    !           So the values we are using are at (CurrentTime + ModelTimeStep)

    CurrentTime = CurrentTime + ModelTimeStep
    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime) 
    ctime(20:20) = c_null_char

    preTime = CurrentTime - GigaTrajInternalPtr%Integrate_DT 
    call ESMF_TimeGet(preTime, timeStringISOFrac=ctime0) 
    ctime0(20:20) = c_null_char

    call ESMF_TimeIntervalGet(GigaTrajInternalPtr%Integrate_DT, d_r8=DT, _RC)

!---------------
! Step 1) Regrid the metData field from cubed to lat-lon
!---------------
    call MAPL_GetPointer(Import, U_cube, "U", _RC)
    call MAPL_GetPointer(Import, V_cube, "V", _RC)
    call MAPL_GetPointer(Import, W_cube, GigaTrajInternalPtr%vTendency, _RC)
    call MAPL_GetPointer(Import, P_cube, GigaTrajInternalPtr%vCoord, _RC)

    lm = size(u_cube,3)
    d1 = size(u_cube,1)
    d2 = size(u_cube,2)

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

!---------------
! Step 2) Get halo
!---------------
    if (GigaTrajInternalPtr%regrid_to_latlon) then

      call GigaTrajInternalPtr%cube2latlon%regrid(U_cube,V_cube, U_latlon, V_latlon, _RC)
      call GigaTrajInternalPtr%cube2latlon%regrid(W_cube, W_latlon, _RC)
      call GigaTrajInternalPtr%cube2latlon%regrid(P_cube, P_latlon, _RC)

      call esmf_halo(grid_, U_Latlon, haloU, _RC)
      call esmf_halo(grid_, V_Latlon, haloV, _RC)
      call esmf_halo(grid_, W_Latlon, haloW, _RC)
      call esmf_halo(grid_, P_Latlon, haloP, _RC)

      deallocate( U_Latlon, V_latlon, W_latlon, P_latlon) 
    else
      call esmf_halo(grid_, U_cube, haloU, _RC)
      call esmf_halo(grid_, V_cube, haloV, _RC)
      call esmf_halo(grid_, W_cube, haloW, _RC)
      call esmf_halo(grid_, P_cube, haloP, _RC)
    endif

!---------------
! Step 3) Update
!---------------

    call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(haloU), c_loc(haloV), c_loc(haloW), c_loc(haloP))

!---------------
! Step 4) Time advance 
!---------------
    call RK4_advance( GigaTrajInternalPtr%metSrc, c_loc(ctime0), DT, GigaTrajInternalPtr%parcels%num_parcels,  &
                       c_loc(GigaTrajInternalPtr%parcels%lons), &
                       c_loc(GigaTrajInternalPtr%parcels%lats), &
                       c_loc(GigaTrajInternalPtr%parcels%zs))
   
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
    real, dimension(:,:), pointer  :: with_halo

    Iam = "Gigatraj ESMF Halo"
    call MAPL_GridGet(grid, localCellCountPerDim=counts, _RC)

    count3 = size(field,3) ! may be nbins

    halo_field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R4, name='halo_field',  &
                                totalLWidth=[1,1],totalUWidth=[1,1])
    call ESMF_FieldGet(halo_field, farrayPtr=with_halo, _RC)
    with_halo = 0.0
    call ESMF_FieldHaloStore(halo_field, rh, _RC)
    !
    ! W.Y note, the pointer with_halo's lbound is 0
    !
    do k = 1, count3
       with_halo(1:counts(1), 1:counts(2)) = Field(:,:,k)
       call ESMF_FieldHalo(halo_field, rh, _RC)
       haloField(:,:,k) = with_halo
    enddo

    call ESMF_FieldDestroy(halo_field)
    call ESMF_FieldHaloRelease(rh, _RC)

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

    where (lons0 < -180.0)  lons0 =lons0 + 360.0
    where (lons0 > 180.0)   lons0 =lons0 - 360.0

    allocate(II(num_parcels0), JJ(num_parcels0))
    call MAPL_GridGet(Grid, globalCellCountPerDim=DIMS)
    if (DIMS(2) == 6*DIMS(1)) then
       call MAPL_GetGlobalHorzIJIndex(num_parcels0, II, JJ, lons0/180.0*MAPL_PI, lats0/180.0*MAPL_PI, Grid=Grid, rc=status)
    else

       dlon = 360.0 / DIMS(1)
       dlat = 180.0 / (DIMS(2)-1)
       ! DC
       II = min( max(ceiling ((lons0+dlon/2.+180.0)/dlon),1), DIMS(1))
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
    where (lons0 < -180.0)  lons0 =lons0 + 360.0
    where (lons0 > 180.0 )  lons0 =lons0 - 360.0

    if (my_rank == 0) then
       if (DIMS(2) == 6*DIMS(1)) then
         call MAPL_GetGlobalHorzIJIndex(num_parcels0, II, JJ, lons0/180.0*MAPL_PI, lats0/180.0*MAPL_PI, Grid=Grid, rc=status)
       else
          dlon = 360.0 / DIMS(1)
          dlat = 180.0 / (DIMS(2)-1) !PC

          II = min( max(ceiling ((lons0+dlon/2.0 + 180.0)/dlon),1), DIMS(1))
          JJ = min( max(ceiling ((lats0+dlat/2.0 + 90.0 )/dlat),1), DIMS(2))
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

    if (my_rank == 0) then
      do i = 1, num_parcels0
         rank   = ranks(i)
         pos = tmp_position(rank+1) +1
         lons_send(pos) = lons0(i)
         lats_send(pos) = lats0(i)
         zs_send(pos)   = zs0(i)
         IDs_send(pos)  = IDs0(i)
         tmp_position(rank+1) = tmp_position(rank+1) + 1
      enddo
    endif

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
    real, allocatable :: lats0(:), lons0(:), zs0(:), values(:), values0(:)
    real,target, allocatable ::  values_2d(:,:)
    real,pointer ::  field(:,:,:)
    integer, allocatable :: ids0(:), ids0_in(:)
    type(ESMF_Alarm)  :: GigaTrajOutAlarm
    type(FileMetadata) :: meta
    real(ESMF_KIND_R8) :: tint_d
    type(ESMF_TimeInterval) :: tint
    type(MAPL_MetaComp),pointer  :: MPL
    character(len=ESMF_MAXSTR)   :: parcels_file, other_fields
    character(len=ESMF_MAXSTR), allocatable   :: varnames(:)
    character(len=:), allocatable:: var_name, var_, comp_name, var_alias, bdlename
    
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)      :: wrap
    character(len=20), target :: ctime
    character(len=:), allocatable :: vAlias

    Iam = "write_parcels"
    call ESMF_ClockGetAlarm(clock, 'GigatrajOut', GigaTrajOutAlarm, _RC)

    if ( .not. ESMF_AlarmIsRinging(GigaTrajOutAlarm)) then
       RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_GetObjectFromGC ( GC, MPL, _RC)
    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC)

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
       if (GigaTrajInternalPtr%vCoord == 'PL') then
         zs0 = zs0 / 100.0 ! hard coded, conert Pa back to hPa
       endif
       ! reorder
       ids0       = ids0 + 1 ! element start 0, make it to 1 for ordering
       ids0(ids0) = [(k, k=1,size(ids0))]

      ! test if ordering is right
      ! ids0_in = ids0
      ! ids0_in = ids0_in(ids0)
      ! do k = 1, size(ids0)
      !    if (k /= ids0_in(k)) then
      !      RETURN_(-1)
      !    endif 
      ! enddo

       lats0 = lats0(ids0(:)) ! id is zero-bases, plus 1 Fortran
       lons0 = lons0(ids0(:))
       where (lons0 >  180.0) lons0 = lons0 - 360.
       where (lons0 < -180.0) lons0 = lons0 + 360.
       zs0   = zs0(ids0(:))
       call formatter%open(trim(parcels_file), pFIO_WRITE, _RC)
       meta = formatter%read(_RC)
       last_time = meta%get_dimension('time', _RC)
       tint = CurrentTime - GigaTrajInternalPtr%startTime
       call ESMF_TimeIntervalGet(tint,d_r8=tint_d,rc=status)

       call formatter%put_var('lat', lats0, start=[1, last_time+1], _RC)
       call formatter%put_var('lon', lons0, start=[1, last_time+1], _RC)
       call formatter%put_var(GigaTrajInternalPtr%vAlias,   zs0,   start=[1, last_time+1], _RC)
       call formatter%put_var('time',  [tint_d],  start=[last_time+1], _RC)
     endif

     ! extra fields
     if (allocated(GigaTrajInternalPtr%ExtraFieldNames)) then
       call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime)
       ctime(20:20) = c_null_char
       do k = 1, size(GigaTrajInternalPtr%ExtraFieldNames)
         comp_name = trim(GigaTrajInternalPtr%ExtraCompNames(k))
         var_name  = trim(GigaTrajInternalPtr%ExtraFieldNames(k))
         bdlename  = trim(GigaTrajInternalPtr%ExtraBundleNames(k))
         var_alias = trim(GigaTrajInternalPtr%ExtraAliasNames(k))
         
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
                 values0 = values0(ids0(:))
                 var_ = var_alias //'00'//i_to_string(i)
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
              call get_metsrc_data (GC, state, ctime, comp_name, bdlename, varnames(i), values, _RC)
              call gather_onefield(total_num, values0, comm, values, GigaTrajInternalPtr%parcels%num_parcels)

              if (my_rank == 0) then
                values0 = values0(ids0(:))
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
            call get_metsrc_data (GC, state, ctime, comp_name, bdlename, var_name, values, _RC )
            call gather_onefield(total_num, values0, comm, values, GigaTrajInternalPtr%parcels%num_parcels)
            if (my_rank == 0) then
              values0 = values0(ids0(:))
              if ( meta%has_variable(var_alias)) then
                if(var_alias == 'P') values0 = values0/100.0 ! hard coded to hPa
                call formatter%put_var( var_alias,   values0,   start=[1, last_time+1], _RC)
              else
                print*, "Please provide "//var_alias // " in the file "//trim(parcels_file)
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
     integer :: comm, my_rank, total_num, ierror, last_time
     real, allocatable :: lats0(:), lons0(:), zs0(:)
     !real(kind=ESMF_KIND_R8), allocatable :: ids0_r(:)
     integer, allocatable :: ids0(:)
     integer :: status, k
     character(len=ESMF_MAXSTR) :: parcels_file
     character(len=ESMF_MAXSTR) :: regrid_to_latlon
     type(MAPL_MetaComp),pointer  :: MPL
     type (ESMF_VM)   :: vm
     type (ESMF_GRID) :: grid_
     class(Variable), pointer :: t
     type(Attribute), pointer :: attr
     class(*), pointer :: units
     character(len=ESMF_MAXSTR) :: Iam ="read_parcels"

     call ESMF_VMGetCurrent(vm, _RC)
     call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
     
     call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)
     call MAPL_GetObjectFromGC ( GC, MPL, _RC)
     call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC)
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

     allocate(lats0(total_num), lons0(total_num), zs0(total_num),ids0(total_num))
        
     if (my_rank ==0) then
        call formatter%get_var('lat', lats0, start = [1,last_time], _RC)
        call formatter%get_var('lon', lons0, start = [1,last_time], _RC)
        call formatter%get_var(internal%vAlias,zs0,   start = [1,last_time], _RC)
        if (internal%vCoord == 'PL') zs0 = zs0*100.0 ! hard coded from hPa to Pa
        call formatter%close(_RC)
        ids0 = [(k, k=0,total_num-1)]
     endif
     call MAPL_GetResource(MPL, regrid_to_latlon, "GIGATRAJ_REGRID_TO_LATLON:",  default= 'YES' , _RC)
     if (trim (regrid_to_latlon) == 'YES') then
        grid_ = internal%LatLonGrid
     else
        grid_ = internal%CubedGrid
     endif
     call scatter_parcels(total_num, lons0, lats0, zs0, IDs0, internal%CellToRank, grid_, comm, &
                              Internal%parcels%lons, &
                              Internal%parcels%lats, &
                              Internal%parcels%zs,   &
                              Internal%parcels%IDS,  & 
                              Internal%parcels%num_parcels) 

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

  subroutine get_metsrc_data (GC, state, ctime, compname, bundlename, fieldname, values, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state
    character(*), target,   intent(in) :: ctime
    character(*), intent(in) :: compname
    character(*), intent(in) :: bundlename
    character(*), intent(in) :: fieldname
    real, target, intent(inout) :: values(:)
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)     :: wrap
    type (ESMF_FieldBundle)         :: bdle
    type (ESMF_GRID)                :: grid_
    type(ESMF_State)  :: leaf_export
    real, dimension(:,:,:), pointer     :: ptr3d
    real, dimension(:,:,:), allocatable :: field_latlon

    real, dimension(:,:,:), allocatable, target  :: haloField
    integer :: counts(3), dims(3), d1, d2, lm, count3
    character(len=:), target, allocatable :: field_

    Iam = "get_metsrc_data"

    !if (index(fieldname,'philicIT') /=0 .or. index(fieldname,'phobicIT') /=0) then
    !   call ESMF_StateGet(state, 'PHYSICS_Exports/TURBULENCE_Exports/TRI', TRI, _RC)
    !   call ESMF_FieldBundleGet(TRI, fieldname, field=field, _RC)
    !   call ESMF_FieldGet(field,farrayPtr=ptr3d, _RC)
    !else
    !   call MAPL_GetPointer(state, ptr3d, fieldname, _RC)
    !endif

    call MAPL_ExportStateGet([state], trim(compname), leaf_export, _RC)
    if (trim(bundlename) /= 'NONE') then
       call ESMF_StateGet(leaf_export,  trim(bundlename), bdle,  _RC)
       call ESMFL_BundleGetPointerToData(bdle, fieldname, ptr3d, _RC)
    else
       call MAPL_GetPointer(leaf_export, ptr3d, fieldname, _RC)
    endif

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    lm = size(ptr3d,3)
    d1 = size(ptr3d,1)
    d2 = size(ptr3d,2)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       grid_ = GigaTrajInternalPtr%LatLonGrid
       call MAPL_GridGet(grid_, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)
       allocate(field_latlon(counts(1),counts(2), lm))
    else
       grid_ = GigaTrajInternalPtr%CubedGrid
       call MAPL_GridGet(grid_, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)
    endif

    allocate(haloField(counts(1)+2, counts(2)+2, lm), source = 0.0)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       call GigaTrajInternalPtr%cube2latlon%regrid(ptr3d, Field_latlon, _RC)
       call esmf_halo(grid_, Field_latlon, haloField, _RC)
       deallocate(Field_latlon)
    else
       call esmf_halo(grid_, ptr3d, haloField, _RC)
    endif

    field_ = trim(fieldname)//c_null_char
    call setData( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(field_), c_loc(haloField))
    call getData(GigaTrajInternalPtr%metSrc,  c_loc(ctime),  c_loc(field_),  &
              GigaTrajInternalPtr%parcels%num_parcels,  &
              c_loc(GigaTrajInternalPtr%parcels%lons),         &
              c_loc(GigaTrajInternalPtr%parcels%lats),         &
              c_loc(GigaTrajInternalPtr%parcels%zs),           &
              c_loc(values))

    deallocate(haloField)
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
    type (ESMF_GRID) :: grid_
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
    if (GigaTrajInternalPtr%regrid_to_latlon) then
       grid_ = GigaTrajInternalPtr%LatLonGrid
    else
       grid_ = GigaTrajInternalPtr%CubedGrid
    endif

    call MAPL_GridGet(grid_, localCellCountPerDim=counts, _RC)

    allocate(field_latlon(counts(1),counts(2),count3))
    allocate(haloField(counts(1)+2, counts(2)+2,count3), source = 0.0)

    if (GigaTrajInternalPtr%regrid_to_latlon) then
       call GigaTrajInternalPtr%cube2latlon%regrid(field, Field_latlon, _RC)
       call esmf_halo(grid_, Field_latlon, haloField, _RC)
    else
       call esmf_halo(grid_, field, haloField, _RC)
    endif

    field_ = trim(fieldname)//'_2D'//c_null_char

    do i = 1,count3

      call setData( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(field_), c_loc(haloField(1,1,i)))
      call getData2d(GigaTrajInternalPtr%metSrc,  c_loc(ctime),  c_loc(field_),  &
                 GigaTrajInternalPtr%parcels%num_parcels,  &
                 c_loc(GigaTrajInternalPtr%parcels%lons),         &
                 c_loc(GigaTrajInternalPtr%parcels%lats),         &
                 c_loc(values(1,i)))

    enddo
    deallocate(field_latlon, haloField, field_)
    RETURN_(ESMF_SUCCESS)

  end subroutine get_metsrc_data2d
 
end module GEOS_GigatrajGridCompMod 
