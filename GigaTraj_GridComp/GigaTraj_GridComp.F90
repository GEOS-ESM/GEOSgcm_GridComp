#include "MAPL_Generic.h"

module GigaTraj_GridCompMod
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated
  use, intrinsic :: iso_c_binding, only : c_loc
  use ESMF
  use MAPL
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
    class (AbstractRegridder), pointer :: cube2latlon => null()
    integer, allocatable :: CellToRank(:,:)
    type(horde) :: parcels
    type(c_ptr) :: metSrc
    type(ESMF_Time) :: startTime
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
    integer :: I1, I2, J1, J2, comm, npes, rank, ierror, NX, NY, NPZ
    type(ESMF_Grid) :: CubedGrid
    integer, allocatable :: I1s(:), J1s(:), I2s(:),J2s(:)
    integer :: DIMS(3), counts(3), i,j,k
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
    type (GigatrajInternalWrap)   :: wrap
    real :: dlat, dlon
    real, allocatable, target :: lats_center(:), lons_center(:), levs_center(:)
    type (ESMF_TIME) :: CurrentTime
    character(len=20), target :: ctime 
    character(len=:), allocatable, target :: name_, unit_
    type(ESMF_Alarm)  :: GigaTrajOutAlarm, GigaTrajRebalanceAlarm
    type(ESMF_TimeInterval)      :: parcels_DT, Rebalance_DT
    type(ESMF_TimeInterval)      :: ModelTimeStep
    type(ESMF_State)             :: INTERNAL
    integer :: minutes_
    character(len=ESMF_MAXSTR) :: parcels_file

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

    call ESMF_GridCompGet(GC, grid=CubedGrid, _RC)
    call MAPL_GridGet(CubedGrid, globalCellCountPerDim=DIMS, _RC)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, status); _VERIFY(STATUS)
    GigaTrajInternalPtr => wrap%ptr

    GigaTrajInternalPtr%npes = npes
    call MAPL_MakeDecomposition(NX,NY,_RC)
    call MAPL_GridGet(CubedGrid, globalCellCountPerDim=DIMS, _RC)

    GigaTrajInternalPtr%LatLonGrid = grid_manager%make_grid( &
                 LatLonGridFactory(im_world=DIMS(1)*4, jm_world=DIMS(1)*2, lm=DIMS(3),  &
                 nx=NX, ny=NY, pole='PE', dateline= 'DE', rc=status) ); _VERIFY(status) 
   
    GigaTrajInternalPtr%cube2latlon => new_regridder_manager%make_regridder(CubedGrid,  GigaTrajInternalPtr%LatLonGrid, REGRID_METHOD_CONSERVE, _RC)

    call MAPL_Grid_interior(GigaTrajInternalPtr%LatLonGrid ,i1,i2,j1,j2)
    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

    ! lat and lon centers need to hold the halo
    dlon = 360.0/dims(1) 
    lons_center = [(-dlon/2+ dlon*i, i= i1-1, i2+1)]
    where(lons_center < 0. ) lons_center  = 0.
    where(lons_center >360.) lons_center = 360.

    dlat = 180.0/dims(2) 
    lats_center = [(-dlat/2. + dlat*j-90.0, j= j1-1, j2+1)] 
    where(lats_center <-90.) lats_center = -90.
    where(lats_center >90. ) lats_center =  90.

    levs_center = [1000.,  975.,  950.,  925.,  900., 875., 850., 825., 800., 775., 750., 725.,700., 650., 600., 550., 500., &
                   450., 400., 350., 300., 250., 200., 150., 100.,  70.,  50.,  40.,  30., 20., 10., 7., 5., 4., 3., 2., &
                   1.  , 0.6,  0.3, 0.1, 0.07, 0.04, 0.02, 0.01]

    npz = size(levs_center, 1)
    GigaTrajInternalPtr%npz = npz

    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime)
    ctime(20:20) = c_null_char

    allocate(I1s(npes),J1s(npes))
    allocate(I2s(npes),J2s(npes))

    call MPI_Allgather(i1, 1, MPI_INTEGER, I1s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)
    call MPI_Allgather(i2, 1, MPI_INTEGER, I2s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)
    call MPI_Allgather(j1, 1, MPI_INTEGER, J1s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)
    call MPI_Allgather(j2, 1, MPI_INTEGER, J2s, 1, MPI_INTEGER, comm, ierror); _VERIFY(ierror)

    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, globalCellCountPerDim=DIMS, _RC)

    allocate(GigaTrajInternalPtr%CellToRank(DIMS(1),DIMS(2))) 
    do rank = 0, npes -1
       I1 = I1s(rank+1)
       I2 = I2s(rank+1)      
       J1 = J1s(rank+1)
       J2 = J2s(rank+1)
       GigaTrajInternalPtr%CellToRank(I1:I2,J1:J2) = rank
    enddo

    GigaTrajInternalPtr%metSrc = initMetGEOSDistributedData(comm, c_loc(GigaTrajInternalPtr%CellToRank), DIMS(1), DIMS(2),   &
                                   dims(3), counts(1)+2, counts(2)+2, npz, &
                                   c_loc(lons_center), c_loc(lats_center), c_loc(levs_center), c_loc(ctime))

    deallocate(lons_center, lats_center,levs_center)

    call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC) 

    call read_parcels(parcels_file, GigaTrajInternalPtr, _RC) 

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

   integer :: status
   character(len=ESMF_MAXSTR) :: IAm 
   type (ESMF_State)          :: INTERNAL
   type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
   type (GigatrajInternalWrap)   :: wrap
   character(len=ESMF_MAXSTR)    :: GigaRstFile
   type (ESMF_TIME) :: CurrentTime
   character(len=20), target :: ctime 
   type (MAPL_MetaComp),  pointer  :: MAPL
   logical, save :: init = .false.

   Iam = "getInitVars"

   if (init) then
     RETURN_(ESMF_SUCCESS)
   endif

   call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

   call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
   call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime) 
   ctime(20:20) = c_null_char

   call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, _RC)

   call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
   GigaTrajInternalPtr => wrap%ptr

   call MAPL_GetResource(MAPL, GigaRstFile, 'GIGATRAJ_INTERNAL_RESTART_FILE:', default="NONE", RC=STATUS )

   if (trim(GigaRstFile) == 'NONE') then
      ! without restart file, get value from import
      call init_metsrc_field0(GC,  IMPORT,  ctime, 'PLE', _RC)
   else
      call init_metsrc_field0(GC, INTERNAL, ctime, 'PL',  _RC)
   endif

   init = .true.

   RETURN_(ESMF_SUCCESS)
 end subroutine GetInitVars

 subroutine Init_metsrc_field0 (GC, state, ctime, PL, RC )
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: state
    character(*), target,   intent(in) :: ctime
    character(*),           intent(in) :: PL
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm 
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaInternalPtr
    type (GigatrajInternalWrap)   :: wrap
    real, dimension(:,:,:), pointer     :: U, V, W, P, PL0, PLE
    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable, target  :: haloU, haloV, haloW, haloP
    integer :: counts(3), dims(3), d1,d2,km

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
    GigaInternalPtr => wrap%ptr
    call MAPL_GridGet(GigaInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

    allocate(U_latlon(counts(1),counts(2),counts(3)))
    allocate(V_latlon(counts(1),counts(2),counts(3)))
    allocate(W_latlon(counts(1),counts(2),counts(3)))
    allocate(P_latlon(counts(1),counts(2),counts(3)))

    allocate(haloU(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(haloV(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(haloW(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(haloP(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)

    call GigaInternalPtr%cube2latlon%regrid(U, V, U_latlon, V_latlon, _RC)
    call GigaInternalPtr%cube2latlon%regrid(W, W_latlon, _RC)
    call GigaInternalPtr%cube2latlon%regrid(P, P_latlon, _RC)

    call esmf_halo(GigaInternalPtr%LatLonGrid, U_Latlon, haloU, _RC)
    call esmf_halo(GigaInternalPtr%LatLonGrid, V_Latlon, haloV, _RC)
    call esmf_halo(GigaInternalPtr%LatLonGrid, W_Latlon, haloW, _RC)
    call esmf_halo(GigaInternalPtr%LatLonGrid, P_Latlon, haloP, _RC)

    call updateFields( GigaInternalPtr%metSrc, c_loc(ctime), c_loc(haloU), c_loc(haloV), c_loc(haloW), c_loc(haloP))

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
  
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap
    type(ESMF_Grid) :: CubedGrid
  
    integer :: num_parcels, my_rank
    real, allocatable, target :: lats(:), lons(:), zs(:)
    real, allocatable, target :: U(:), V(:), W(:)
    integer ::counts(3), DIMS(3), rank, comm, ierror
    type (ESMF_VM)   :: vm
  
    real, dimension(:,:,:), pointer     :: U_cube, V_cube, W_cube, p_cube, with_halo
    real, dimension(:,:,:), pointer     :: U_internal, V_internal, W_internal, P_internal

    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable, target  :: haloU, haloV, haloW, haloP

    real(ESMF_KIND_R8) :: DT

    type(ESMF_Field)   :: halo_field
    type(ESMF_RouteHandle) :: rh
    character(len=20), target :: ctime
    type(ESMF_State)            :: INTERNAL
    type(MAPL_MetaComp),pointer :: MPL  
    character(len=ESMF_MAXSTR) :: parcels_file


    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)

    call MAPL_GetObjectFromGC ( GC, MPL, _RC)
    call MAPL_Get (MPL, INTERNAL_ESMF_STATE=INTERNAL, _RC)

    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, _RC)
 
    ! W.J note: this run is after agcm's run. The clock is not yet ticked
    !           So the values we are using are at (CurrentTime + ModelTimeStep)
    CurrentTime = CurrentTime + ModelTimeStep
    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime) 
    ctime(20:20) = c_null_char

    call ESMF_TimeIntervalGet(ModelTimeStep,d_r8=DT, _RC)

    call ESMF_GridCompGet(GC, grid=CubedGrid, _RC)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts, &
                      globalCellCountPerDim=DIMS, _RC)

    allocate(U_latlon(counts(1), counts(2),counts(3)), source = 0.0)
    allocate(V_latlon(counts(1), counts(2),counts(3)), source = 0.0)
    allocate(W_latlon(counts(1), counts(2),counts(3)), source = 0.0)
    allocate(P_latlon(counts(1), counts(2),counts(3)), source = 0.0)

    allocate(haloU(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(haloV(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(haloW(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(haloP(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)

!---------------
! Step 1) Regrid the metData field from cubed to lat-lon
!---------------

    call MAPL_GetPointer(Import, U_cube, "U", _RC)
    call MAPL_GetPointer(Import, V_cube, "V", _RC)
    call MAPL_GetPointer(Import, W_cube, "OMEGA", _RC)
    call MAPL_GetPointer(Import, P_cube, "PL", _RC)

    call GigaTrajInternalPtr%cube2latlon%regrid(U_cube,V_cube, U_latlon, V_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(W_cube, W_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(P_cube, P_latlon, _RC)

!---------------
! Step 2) Get halo
!---------------

    call esmf_halo(GigaTrajInternalPtr%LatLonGrid, U_Latlon, haloU, _RC)
    call esmf_halo(GigaTrajInternalPtr%LatLonGrid, V_Latlon, haloV, _RC)
    call esmf_halo(GigaTrajInternalPtr%LatLonGrid, W_Latlon, haloW, _RC)
    call esmf_halo(GigaTrajInternalPtr%LatLonGrid, P_Latlon, haloP, _RC)


!---------------
! Step 3) Update
!---------------
    call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(haloU), c_loc(haloV), c_loc(haloW), c_loc(haloP))

!---------------
! Step 3) Time advance 
!---------------
    call rk4a_advance( GigaTrajInternalPtr%metSrc, c_loc(ctime), DT, GigaTrajInternalPtr%parcels%num_parcels,  &
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
 
    U_internal = U_cube
    V_internal = V_cube
    W_internal = W_cube
    P_internal = P_cube
    
    deallocate( U_Latlon, V_latlon, W_latlon, P_latlon, haloU, haloV, haloW, haloP)
  
!---------------
! Step 5) rebalance parcels among processors ( configurable with alarm) 
!---------------
    call rebalance_parcels(clock, GigaTrajInternalPtr%parcels, GigaTrajInternalPtr%CellToRank, comm, DIMS, _RC)

!---------------
! Step 6) write out parcel positions ( configurable with alarm) 
!---------------
    call MAPL_GetResource(MPL, parcels_file, "GIGATRAJ_PARCELS_FILE:", default='parcels.nc4', _RC) 
    call write_parcels(clock, parcels_file, GigaTrajInternalPtr%parcels,currentTime, GigaTrajInternalPtr%startTime, _RC)


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
    integer :: counts(3), dims(3), k
    integer :: status
    type(ESMF_Field)   :: halo_field
    type(ESMF_RouteHandle) :: rh
    real, dimension(:,:,:), pointer  :: with_halo

    Iam = "Gigatraj ESMF Halo"
    call MAPL_GridGet(grid, localCellCountPerDim=counts, &
                      globalCellCountPerDim=DIMS, _RC)

    halo_field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R4, name='halo_field',  &
                                ungriddedLBound=[1],ungriddedUBound=[counts(3)], &
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
  subroutine rebalance_parcels(clock, parcels, CellToRank, comm, DIMS, rc)
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    type(horde), intent(inout) :: parcels
    integer, dimension(:,:), intent(in) :: CellToRank
    integer :: comm, DIMS(3)
    integer, optional, intent(out) :: rc

    integer :: status
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

    dlon = 360.0 / DIMS(1)
    dlat = 180.0 / DIMS(2)

    II = min( max(ceiling (lons0         /dlon),1), DIMS(1))
    JJ = min( max(ceiling ((lats0 + 90.0)/dlat),1), DIMS(2))

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
  subroutine scatter_parcels(num_parcels0, lons0, lats0, zs0, IDs0, CellToRank, DIMS, comm, lons, lats, zs, IDs, num_parcels)
    integer :: num_parcels0
    real, dimension(:), intent(inout) :: lons0
    real, dimension(:), intent(in)    :: lats0, zs0
    integer, dimension(:), intent(in) :: IDs0
    integer, dimension(:,:), intent(in) :: CellToRank
    integer :: comm, DIMS(3)
    real, dimension(:), allocatable, intent(out) :: lons, lats, zs
    integer, dimension(:), allocatable, intent(out) :: IDs
    integer, intent(out) :: num_parcels

    integer :: i, npes, ierror, rank, my_rank, counts_recv, pos
    real :: dlon, dlat

    real, allocatable :: lons_send(:), lats_send(:), zs_send(:)
    integer, allocatable :: ids_send(:)

    integer, allocatable :: counts_send(:), II(:), JJ(:), ranks(:)
    integer, allocatable :: disp_send(:), tmp_position(:)

    call MPI_Comm_size(comm, npes, ierror)
    call MPI_Comm_rank(comm, my_rank, ierror)

    allocate(counts_send(npes), source = 0)
    allocate(disp_send(npes), source = 0)
    if (my_rank == 0) then
       dlon = 360.0 / DIMS(1)
       dlat = 180.0 / DIMS(2)

       where (lons0 < 0) lons0 =lons0 + 360.0
       II = min( max(ceiling (lons0         /dlon),1), DIMS(1))
       JJ = min( max(ceiling ((lats0 + 90.0)/dlat),1), DIMS(2))

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

  subroutine write_parcels(CLOCK, fname, parcels, currentTime, startTime, rc)
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    character(*), intent(in):: fname
    type(horde), intent(in) :: parcels
    type(ESMF_TIME), intent(in) :: currentTime 
    type(ESMF_TIME), intent(in) :: startTime 
    integer, optional, intent(out) :: rc
  
    character(len=:), allocatable :: Iam
    type (ESMF_VM)   :: vm
    type(Netcdf4_fileformatter) :: formatter
    integer :: comm, my_rank, total_num, status, last_time, ierror
    real, allocatable :: lats0(:), lons0(:), zs0(:), ids0_r(:)
    integer, allocatable :: ids0(:)
    type(ESMF_Alarm)  :: GigaTrajOutAlarm
    type(FileMetadata) :: meta
    real(ESMF_KIND_R8) :: tint_d
    type(ESMF_TimeInterval) :: tint

    call ESMF_ClockGetAlarm(clock, 'GigatrajOut', GigaTrajOutAlarm, _RC)

    if ( .not. ESMF_AlarmIsRinging(GigaTrajOutAlarm)) then
       RETURN_(ESMF_SUCCESS)
    endif

    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)

    call gather_parcels(total_num, lons0, lats0, zs0, IDs0, &
                          comm,                             & 
                          parcels%lons, &
                          parcels%lats, &
                          parcels%zs,   &
                          parcels%IDS,  &
                          parcels%num_parcels )

    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)
    if (my_rank ==0) then
       ! reorder 
       lats0 = lats0(ids0(:)+1) ! id is zero-bases, plus 1 Fortran
       lons0 = lons0(ids0(:)+1)
       zs0   = zs0(ids0(:)+1)
       call formatter%open(fname, pFIO_WRITE, _RC)
       meta = formatter%read(_RC)
       last_time = meta%get_dimension('time', _RC)
       tint = CurrentTime - startTime
       call ESMF_TimeIntervalGet(tint,d_r8=tint_d,rc=status)

       call formatter%put_var('lat', lats0, start=[1, last_time+1], _RC)
       call formatter%put_var('lon', lons0, start=[1, last_time+1], _RC)
       call formatter%put_var('pressure',   zs0,   start=[1, last_time+1], _RC)
       call formatter%put_var('time',  [tint_d],  start=[last_time+1], _RC)
       call formatter%close(_RC)
     endif

     RETURN_(ESMF_SUCCESS)
  end subroutine write_parcels

  subroutine read_parcels(fname, internal, rc)
     character(*), intent(in) :: fname
     type(GigaTrajInternal), intent(inout) :: internal
     integer, optional, intent(out) :: rc

     type(Netcdf4_fileformatter) :: formatter
     type(FileMetadata) :: meta 
     integer :: comm, my_rank, total_num, ierror, last_time, DIMS(3)
     real, allocatable :: lats(:), lons(:), zs(:), lats0(:), lons0(:), zs0(:)
     real(kind=ESMF_KIND_R8), allocatable :: ids0_r(:)
     integer, allocatable :: ids0(:)
     integer :: status
     type (ESMF_VM)   :: vm
     class(Variable), pointer :: v
     type(Attribute), pointer :: attr
     class(*), pointer :: units
     character(len=ESMF_MAXSTR) :: Iam ="read_parcels"

     call ESMF_VMGetCurrent(vm, _RC)
     call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
     
     call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)

      total_num = 0
      if (my_rank ==0) then
        call formatter%open(fname, pFIO_READ, _RC)
        meta = formatter%read(_RC)
        total_num = meta%get_dimension('id', _RC)
        last_time = meta%get_dimension('time', _RC)
        v => meta%get_variable('time', _RC)
        attr => v%get_attribute('long_name')
        units => attr%get_value()
        select type(units)
        type is (character(*))
           internal%startTime = parse_time_string(units, _RC)
        class default
           _FAIL('unsupported subclass for units')
         end select
      endif

      allocate(lats0(total_num), lons0(total_num), zs0(total_num), ids0_r(total_num))

      if  (my_rank ==0) then
        call formatter%get_var('lat', lats0, start = [1,last_time], _RC)
        call formatter%get_var('lon', lons0, start = [1,last_time], _RC)
        call formatter%get_var('pressure',   zs0,   start = [1,last_time], _RC)
        call formatter%get_var('id',  ids0_r,start = [1,last_time], _RC)
        call formatter%close(_RC)
        ids0 = int(ids0_r)
      endif
      call MAPL_GridGet(internal%LatLonGrid, globalCellCountPerDim=DIMS, _RC)
      call scatter_parcels(total_num, lons0, lats0, zs0, IDs0, internal%CellToRank, DIMS, comm, &
                              Internal%parcels%lons, &
                              Internal%parcels%lats, &
                              Internal%parcels%zs,   &
                              Internal%parcels%IDS,  & 
                              Internal%parcels%num_parcels) 

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
    character(*), target,   intent(in) :: fieldname
    real, target, intent(inout) :: values(:)
    integer, optional,     intent(out) :: RC     ! Error code

    character(len=ESMF_MAXSTR)    :: IAm
    integer                       :: STATUS

    type(GigaTrajInternal), pointer :: GigaInternalPtr
    type (GigatrajInternalWrap)     :: wrap
    real, dimension(:,:,:), pointer     :: field
    real, dimension(:,:,:), allocatable :: field_latlon
    real, dimension(:,:,:), allocatable, target  :: haloField
    integer :: counts(3), dims(3), d1,d2,km

    Iam = "get_metsrc_data"

    call MAPL_GetPointer(state, field, fieldname, _RC)
    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaInternalPtr => wrap%ptr
    call MAPL_GridGet(GigaInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

    allocate(field_latlon(counts(1),counts(2),counts(3)))
    allocate(haloField(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)

    call GigaInternalPtr%cube2latlon%regrid(field, Field_latlon, _RC)

    call esmf_halo(GigaInternalPtr%LatLonGrid, field_Latlon, haloField, _RC)

    call setData( GigaInternalPtr%metSrc, c_loc(ctime), c_loc(fieldname), c_loc(haloField))

    call getData(GigaInternalPtr%metSrc,  c_loc(ctime),  c_loc(fieldname),  &
                 GigaInternalPtr%parcels%num_parcels,  &
                 c_loc(GigaInternalPtr%parcels%lons),         &
                 c_loc(GigaInternalPtr%parcels%lats),         &
                 c_loc(GigaInternalPtr%parcels%zs),           &
                 c_loc(values))

    RETURN_(ESMF_SUCCESS)

  end subroutine get_metsrc_data

end module GigaTraj_GridCompMod
