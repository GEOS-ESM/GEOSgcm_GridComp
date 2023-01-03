#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GigaTraj_GridCompMod -- A Module to run gigatraj

! !INTERFACE:

module GigaTraj_GridCompMod

! !USES:
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated
  use, intrinsic :: iso_c_binding, only : c_loc
  use ESMF
  use MAPL
  use mpi
  use GEOS_Giga_interOpMod
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public :: SetServices

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
  end type

  type GigatrajInternalWrap
    type (GigaTrajInternal), pointer :: PTR
  end type

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    type (ESMF_Config)                  :: CF
    type (MAPL_MetaComp),  pointer      :: MAPL

    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)      :: wrap

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, _RC )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  GetInitVars , _RC )
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run         , _RC )

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, _RC )

! Set the state variable specs.
! -----------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_bias_tendency',              &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             _RC  )


    allocate(GigaTrajInternalPtr)
    wrap%ptr => GigaTrajInternalPtr
    call ESMF_UserCompSetInternalState ( GC, 'GigaTrajInternal', wrap, status ); _VERIFY(STATUS) 

    call MAPL_GenericSetServices    ( GC, _RC )

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,_RC)
    call MAPL_TimerAdd(GC, name="RUN"           ,_RC)

! All done
!---------

    RETURN_(ESMF_SUCCESS)  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
 
! !IROUTINE: Initialize -- Initialize method for the composite GigaTraj Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)           :: IAm 
  integer                              :: STATUS
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases
  type (MAPL_MetaComp),  pointer  :: STATE 
  type (ESMF_VM)                  :: vm
  integer :: I1, I2, J1, J2, comm, npes, rank, ierror, NX, NY, NPZ
  type(ESMF_Grid) :: CubedGrid
  integer, allocatable :: I1s(:), J1s(:), I2s(:),J2s(:)
  integer :: DIMS(3), counts(3), i,j,k
  type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
  type (GigatrajInternalWrap)   :: wrap
  real :: dlat, dlon
  real, pointer :: lats_center(:), lons_center(:), levs_center(:)
  type (ESMF_TIME) :: CurrentTime
  character(len=20), target :: ctime 
  character(len=:), allocatable, target :: name_, unit_

! =============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, _RC)


    call MAPL_TimerOn(STATE,"TOTAL")
    call MAPL_TimerOn(STATE,"INITIALIZE")

! Call Initialize for every Child

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
    GigaTrajInternalPtr%LatLonGrid = grid_manager%make_grid(                           &
                 LatLonGridFactory(im_world=DIMS(1)*4, jm_world=DIMS(1)*2+1, lm=DIMS(3),  &
                 nx=NX, ny=NY, pole='PC', dateline= 'DE', rc=status) ); _VERIFY(status) 


    GigaTrajInternalPtr%cube2latlon => new_regridder_manager%make_regridder(CubedGrid,  GigaTrajInternalPtr%LatLonGrid, REGRID_METHOD_CONSERVE, _RC)

    call MAPL_Grid_interior(GigaTrajInternalPtr%LatLonGrid ,i1,i2,j1,j2)
    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

    ! lat and lon centers need to hold the halo
    allocate(lons_center(counts(1)+2))
    allocate(lats_center(counts(2)+2))
    allocate(levs_center(counts(3)  ))
    dlon = 360.0/dims(1) 
    lons_center = [(-dlon/2+ dlon*i, i= i1-1, i2+1)]
    ! for PC grid
    dlat = 180.0/(dims(2)-1) 
    lats_center = [(-dlat + dlat*j-90.0, j= j1-1, j2+1)] 
    npz = 42
    GigaTrajInternalPtr%npz = npz
    allocate(levs_center(GigaTrajInternalPtr%npz))   
    levs_center = [1000.,  975.,  950.,  925.,  900., 875., 850., 825., 800., 775., 750., 725.,700., 650., 600., 550., 500., &
                   450., 400., 350., 300., 250., 200., 150., 100.,  70.,  50.,  40.,  30., 20., 10., 7., 5., 4., 3., 2., &
                   1.  , 0.7,  0.5, 0.4,  0.3, 0.1]

    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
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


! initialize partical positions. It will be read and distribted across processors.
! for now, it is genrated randomly
    block
      real :: rparcels
      integer :: num_parcels, my_rank, i
      real, allocatable :: lats(:), lons(:), zs(:)
      integer, allocatable :: nums_all(:)
 

      call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)

      !call random_seed()
      !call random_number(rparcels)

      !num_parcels = nint(rparcels*10)
      num_parcels = 1

      allocate(lats(num_parcels), lons(num_parcels), zs(num_parcels))

      !call random_number(lats)
      !call random_number(lons)
      !call random_number(zs)

      !lons = lons*360.0 
      !lats = lats*180.0  - 90.0
      !lats = lats*20.0  - 10.0
      !zs   = (999.9*zs) + 0.1 

      call MAPL_Grid_interior(GigaTrajInternalPtr%LatLonGrid ,i1,i2,j1,j2)
      dlon  = 360.0/DIMS(1)
      dlat = 180.0/(DIMS(2)-1)
      lons = [(i1+i2)*dlon/2.0]
      lats = [(j1+j2)*dlat/2.0 - 90.0]
      zs  = [my_rank*10.0+10.0]

      GigaTrajInternalPtr%parcels%num_parcels = num_parcels
      GigaTrajInternalPtr%parcels%lats = lats
      GigaTrajInternalPtr%parcels%lons = lons
      GigaTrajInternalPtr%parcels%zs   = zs
      GigaTrajInternalPtr%parcels%IDS = [(i, i=1, num_parcels)]

      allocate(nums_all(npes))
      call MPI_AllGather(num_parcels, 1, MPI_INTEGER, nums_all, 1, MPI_INTEGER, comm, ierror)

      if (my_rank > 0) then
        GigaTrajInternalPtr%parcels%IDS = GigaTrajInternalPtr%parcels%IDS + sum(nums_all(1:my_rank))
      endif

    end block

    deallocate(lons_center, lats_center,levs_center)

    call MAPL_TimerOff(STATE,"INITIALIZE")
    call MAPL_TimerOff(STATE,"TOTAL")


    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: GetInitVars -- GetInitVars method for Gigatraj GridComp to get initial state from AGCM's export

! !INTERFACE:

  subroutine GetInitVars ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 
!EOP

    character(len=ESMF_MAXSTR)       :: IAm 
    type (ESMF_VM)                :: VM
    integer                       :: STATUS, Comm
    logical, save :: first_time_run = .true.
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap

    real, dimension(:,:,:), pointer     :: U, V, W, with_halo, P
    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable, target :: preU, preV, preW, preP
    integer :: counts(3), dims(3)
    type(ESMF_Field)   :: field
    type(ESMF_RouteHandle) :: rh
    real, pointer :: lats_center(:), lons_center(:), levs_center(:)
    integer :: i1, i2, j1, j2, i,j,k
    type (ESMF_TIME) :: CurrentTime
    character(len=20), target :: ctime 

    if (.not. first_time_run) then
       RETURN_(ESMF_SUCCESS)
    endif

    call ESMF_VmGetCurrent(VM, _RC)
    call ESMF_VMGet(VM, mpiCommunicator=Comm, _RC)
    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
    call ESMF_TimeGet(CurrentTime, timeStringISOFrac=ctime)
    ctime(20:20) = c_null_char

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr
    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=counts,globalCellCountPerDim=DIMS,  _RC)

    call MAPL_GetPointer(Import, U, "U", _RC)
    call MAPL_GetPointer(Import, V, "V", _RC)
    call MAPL_GetPointer(Import, W, "OMEGA", _RC)
    call MAPL_GetPointer(Import, P, "PL", _RC)

    allocate(U_latlon(counts(1),counts(2),counts(3)))
    allocate(V_latlon(counts(1),counts(2),counts(3)))
    allocate(W_latlon(counts(1),counts(2),counts(3)))
    allocate(P_latlon(counts(1),counts(2),counts(3)))

    call GigaTrajInternalPtr%cube2latlon%regrid(U,V, U_latlon, V_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(W, W_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(P, P_latlon, _RC)

    field = ESMF_FieldCreate(GigaTrajInternalPtr%LatLonGrid, ESMF_TYPEKIND_R4, name='halo_field',  &
                               ungriddedLBound=[1],ungriddedUBound=[counts(3)], &
                               totalLWidth=[1,1],totalUWidth=[1,1])

    call ESMF_FieldHaloStore(field,rh,_RC)

    call ESMF_FieldGet(field, farrayPtr=with_halo, _RC)
    with_halo(2:counts(1)+1, 2:counts(2)+1, :) = U_latlon
    call ESMF_FieldHalo(field,rh,_RC)
    preU = with_halo

    with_halo(2:counts(1)+1, 2:counts(2)+1, :) = V_latlon
    call ESMF_FieldHalo(field,rh,_RC)
    preV = with_halo

    with_halo(2:counts(1)+1, 2:counts(2)+1, :) = W_latlon
    call ESMF_FieldHalo(field,rh,_RC)
    preW = with_halo

    with_halo(2:counts(1)+1, 2:counts(2)+1, :) = p_latlon
    call ESMF_FieldHalo(field,rh,_RC)
    preP = with_halo

    call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime), c_loc(preU), c_loc(preV), c_loc(preW), c_loc(preP))

    call ESMF_FieldDestroy(field)
 
    first_time_run = .false.
    
    RETURN_(ESMF_SUCCESS)
    
  end subroutine GetInitVars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Run -- Run method for Gigatraj GridComp

! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 
!EOP

! ErrLog Variables

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
    real, dimension(:,:,:), allocatable :: U_latlon, V_latlon, W_latlon, P_latlon
    real, dimension(:,:,:), allocatable, target  :: U_latlon_halo, V_latlon_halo, W_latlon_halo, p_latlon_halo

    real(ESMF_KIND_R8) :: DT

    type(ESMF_Field)   :: halo_field
    type(ESMF_RouteHandle) :: rh
    character(len=20), target :: ctime


    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)

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

    allocate(U_latlon_halo(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(V_latlon_halo(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(W_latlon_halo(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)
    allocate(P_latlon_halo(counts(1)+2, counts(2)+2,counts(3)), source = 0.0)

!---------------
! Step 1) Regrid the metData field from cubed to lat-lon
!---------------

    call MAPL_GetPointer(Import, U_cube, "U", _RC)
    call MAPL_GetPointer(Import, V_cube, "V", _RC)
    call MAPL_GetPointer(Import, W_cube, "OMEGA", _RC)
    call MAPL_GetPointer(Import, P_cube, "PL", _RC)

    call GigaTrajInternalPtr%cube2latlon%regrid(U_cube, U_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(V_cube, V_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(W_cube, W_latlon, _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(P_cube, P_latlon, _RC)

!---------------
! Step 2) Get halo of latlon metData field
!         After this step, the local field has distributed horizonal + halo
!---------------
  
     halo_field = ESMF_FieldCreate(GigaTrajInternalPtr%LatLonGrid, ESMF_TYPEKIND_R4, name='halo_field',  &
                                ungriddedLBound=[1],ungriddedUBound=[counts(3)], &
                                totalLWidth=[1,1],totalUWidth=[1,1])
     call ESMF_FieldHaloStore(halo_field, rh, _RC)

     call ESMF_FieldGet(halo_field, farrayPtr=with_halo, _RC)

     ! get U + halo
     with_halo(2:counts(1)+1, 2:counts(2)+1, :) = U_latlon
     call ESMF_FieldHalo(halo_field, rh, _RC)
     U_latlon_halo = with_halo
     
     ! get V + halo
     with_halo(2:counts(1)+1, 2:counts(2)+1, :) = V_latlon
     call ESMF_FieldHalo(halo_field, rh, _RC)
     V_latlon_halo = with_halo

     ! get W + halo
     with_halo(2:counts(1)+1, 2:counts(2)+1, :) = W_latlon
     call ESMF_FieldHalo(halo_field, rh, _RC)
     W_latlon_halo = with_halo

     ! get W + halo
     with_halo(2:counts(1)+1, 2:counts(2)+1, :) = P_latlon
     call ESMF_FieldHalo(halo_field, rh, _RC)
     P_latlon_halo = with_halo

     call updateFields( GigaTrajInternalPtr%metSrc, c_loc(ctime),c_loc(U_latlon_halo),  &
                       c_loc(V_latlon_halo), c_loc(W_latlon_halo), c_loc(P_latlon_halo))


     num_parcels = GigaTrajInternalPtr%parcels%num_parcels

    !lons = GigaTrajInternalPtr%parcels%lons
    !lats = GigaTrajInternalPtr%parcels%lats
    !zs   = GigaTrajInternalPtr%parcels%zs

    !call test_dataflow(num_parcels, lons, lats, zs, GigaTrajInternalPtr%CellToRank, DIMS, comm)

    !allocate(U(num_parcels))
    !allocate(V(num_parcels))
    !allocate(W(num_parcels)) 

    !call test_metData( GigaTrajInternalPtr%metSrc, 0.01d0, num_parcels, c_loc(lons), c_loc(lats), c_loc(zs), c_loc(U), c_loc(V), c_loc(W))

    call rk4a_advance( GigaTrajInternalPtr%metSrc, c_loc(ctime), DT, num_parcels,  &
                       c_loc(GigaTrajInternalPtr%parcels%lons), &
                       c_loc(GigaTrajInternalPtr%parcels%lats), &
                       c_loc(GigaTrajInternalPtr%parcels%zs))
 
  
    call ESMF_FieldDestroy( halo_field)

    if (MAPL_AM_I_ROOT())  print*," Great, the end of the GigaTraj_GridCompMod run"

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

end module GigaTraj_GridCompMod
