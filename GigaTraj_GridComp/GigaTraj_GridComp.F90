#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GigaTraj_GridCompMod -- A Module to run gigatraj

! !INTERFACE:

module GigaTraj_GridCompMod

! !USES:

  use ESMF
  use MAPL
  use mpi
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public :: SetServices

  type GigaTrajInternal
    integer :: npes
    type (ESMF_Grid) :: LatLonGrid
    class (AbstractRegridder), pointer :: cube2latlon => null()
    integer, allocatable :: CellToRank(:,:)
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

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run       , RC=STATUS )
    VERIFY_(STATUS)

! Get the configuration from the component
!-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_bias_tendency',              &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         FIELD_TYPE = MAPL_VectorField,                            &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)


    allocate(GigaTrajInternalPtr)
    wrap%ptr => GigaTrajInternalPtr
    call ESMF_UserCompSetInternalState ( GC, 'GigaTrajInternal', wrap, status )
    VERIFY_(STATUS) 

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)
! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

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
  type (ESMF_VM)                      :: vm
  integer :: I1, I2, J1, J2, comm, npes, rank, ierror, NX, NY
  type(ESMF_Grid) :: CubedGrid
  integer, allocatable :: I1s(:), J1s(:), I2s(:),J2s(:)
  integer :: DIMS(3)
  type (GigaTrajInternal), pointer :: GigaTrajInternalPtr 
  type (GigatrajInternalWrap)   :: wrap
  

! =============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_TimerOn(STATE,"TOTAL")
    call MAPL_TimerOn(STATE,"INITIALIZE")

! Call Initialize for every Child

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGetCurrent(vm, rc=status)
    call ESMF_VMGet(vm, mpiCommunicator=comm, __RC__)
    call MPI_Comm_size(comm, npes, ierror); _VERIFY(ierror)

    call ESMF_GridCompGet(GC, grid=CubedGrid, rc=status)
    call MAPL_GridGet(CubedGrid, globalCellCountPerDim=DIMS, RC=status)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, status)
    VERIFY_(STATUS)
    GigaTrajInternalPtr => wrap%ptr

    GigaTrajInternalPtr%npes = npes
    call MAPL_MakeDecomposition(NX,NY,rc=status)
    VERIFY_(status)
    GigaTrajInternalPtr%LatLonGrid = grid_manager%make_grid(                           &
                 LatLonGridFactory(im_world=DIMS(1)*4, jm_world=DIMS(1)*2, lm=DIMS(3),  &
                 nx=NX, ny=NY, pole='PC', dateline= 'DC', rc=status) ) 


    GigaTrajInternalPtr%cube2latlon => new_regridder_manager%make_regridder(CubedGrid,  GigaTrajInternalPtr%LatLonGrid, REGRID_METHOD_CONSERVE, rc=status)
    call MAPL_Grid_interior(GigaTrajInternalPtr%LatLonGrid ,i1,i2,j1,j2)

    allocate(I1s(npes),J1s(npes))
    allocate(I2s(npes),J2s(npes))

    call MPI_Allgather(i1, 1, MPI_INTEGER, I1s, 1, MPI_INTEGER, comm, ierror)
    _VERIFY(ierror)
    call MPI_Allgather(i2, 1, MPI_INTEGER, I2s, 1, MPI_INTEGER, comm, ierror)
    _VERIFY(ierror)
    call MPI_Allgather(j1, 1, MPI_INTEGER, J1s, 1, MPI_INTEGER, comm, ierror)
    _VERIFY(ierror)
    call MPI_Allgather(j2, 1, MPI_INTEGER, J2s, 1, MPI_INTEGER, comm, ierror)
    _VERIFY(ierror)

    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, globalCellCountPerDim=DIMS, RC=status)

    allocate(GigaTrajInternalPtr%CellToRank(DIMS(1),DIMS(2))) 
    
    do rank = 0, npes -1
       I1 = I1s(rank+1)
       I2 = I2s(rank+1)      
       J1 = J1s(rank+1)
       J2 = J2s(rank+1)
       GigaTrajInternalPtr%CellToRank(I1:I2,J1:J2) = rank
    enddo
   
    call MAPL_TimerOff(STATE,"INITIALIZE")
    call MAPL_TimerOff(STATE,"TOTAL")


    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

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
    integer        :: CSTAT, ESTAT, YY, MM, HH, DD, H, M,S, model_dtstep_m
    character(512) :: CMSG
    character(256) :: command_line
    character(19)  :: begdate, enddate
    character(64)  :: format_string
    type(ESMF_TimeInterval) :: ModelTimeStep
    type(ESMF_Time)         :: CurrentTime
  
    type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
    type (GigatrajInternalWrap)   :: wrap
    type(ESMF_Grid) :: CubedGrid
  
    integer :: num_parcels, i
    integer,parameter :: seed = 86456
    real, allocatable :: lats(:), lons(:), lons_send(:), lats_send(:), U(:), pos(:)
    real, allocatable :: lons_recv(:), lats_recv(:), U_recv(:), U_send(:)
    real :: rparcels, dlat, dlon
    integer, allocatable :: counts_send(:),counts_recv(:), II(:), JJ(:), ranks(:)
    integer, allocatable :: disp_send(:), disp_recv(:), tmp_position(:)
    integer :: DIMS(3), rank, comm, ierror
    type (ESMF_VM)   :: vm
  
    real, dimension(:,:,:), pointer     :: U_cube
    real, dimension(:,:,:), allocatable :: U_latlon
    real, dimension(:,:,:), pointer     :: U_latlon_halo
    integer :: halowidth(3)
    type(ESMF_Field)   :: U_field
    type(ESMF_RouteHandle) :: rh

    call ESMF_VMGetCurrent(vm, _RC)
    call ESMF_VMGet(vm, mpiCommunicator=comm, _RC)
    call MPI_Comm_rank(comm, my_rank, ierror); _VERIFY(ierror)

    call ESMF_ClockGet(clock, currTime=CurrentTime, _RC)
  
    call ESMF_ClockGet(clock, timeStep=ModelTimeStep, _RC)

    call ESMF_GridCompGet(GC, grid=CubedGrid, _RC)

    call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, _RC)
    GigaTrajInternalPtr => wrap%ptr

    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, localCellCountPerDim=DIMS, _RC)

    allocate(U_latlon(DIMS(1), DIMS(2),DIMS(3)), source = 0.0)

!---------------
! Step 1) Regrid the metData field from cubed to lat-lon
!---------------

    call MAPL_GetPointer(Import, U_cube, "U", _RC)
    call GigaTrajInternalPtr%cube2latlon%regrid(U_cube, U_latlon, _RC)

!---------------
! Step 2) Get halo of latlon metData field
!         After this step, the local field has distributed horizonal + halo
!---------------
  
     U_field = ESMF_FieldCreate(GigaTrajInternalPtr%LatLonGrid, ESMF_TYPEKIND_R4, name='U',  &
                                ungriddedLBound=[1],ungriddedUBound=[DIMS(3)], &
                                totalLWidth=[1,1],totalUWidth=[1,1])
     call ESMF_FieldHaloStore(U_field,rh,rc=status)
     _VERIFY(status)

     call ESMF_FieldGet(U_field, farrayPtr=U_latlon_halo, _RC)
     U_latlon_halo(2:DIMS(1)+1, 2:DIMS(2)+1, :) = U_latlon
     call ESMF_FieldHalo(U_field,rh,rc=status)
      _VERIFY(status)
 
!-----------------
! Step 3) Given partical position (lat lon), find out which processor it should go ?
!-----------------

    call random_seed()
    call random_number(rparcels)

    num_parcels = nint(rparcels*10)

    allocate(lats(num_parcels), lons(num_parcels))

    call random_number(lats)
    call random_number(lons)

    lons = lons*2*MAPL_PI
    lats = lats*MAPL_PI

    allocate(II(num_parcels),JJ(num_parcels), ranks(num_parcels))
    allocate(counts_send(GigaTrajInternalPtr%npes))
    allocate(counts_recv(GigaTrajInternalPtr%npes))
    allocate(disp_send(GigaTrajInternalPtr%npes))
    allocate(disp_recv(GigaTrajInternalPtr%npes))

    call MAPL_GridGet(GigaTrajInternalPtr%LatLonGrid, globalCellCountPerDim=DIMS, RC=status)

    dlon = 2*MAPL_PI / DIMS(1)
    dlat = MAPL_PI / DIMS(2)

    II = ceiling (lons/dlon)
    JJ = ceiling (lats/dlat)

    if (any(II > DIMS(1)) .or. any(II<0)) stop ("wrong II")
    if (any(JJ > DIMS(2)) .or. any(JJ<0)) stop ("wrong JJ")
    do i = 1, num_parcels
       ranks(i) = GigaTrajInternalPtr%CellToRank(II(i), JJ(i))
    enddo

!---------------------
!step 4) Pack the location data and send them to where the metData sit
!---------------------

    do rank = 0, GigaTrajInternalPtr%npes-1
       counts_send(rank+1) = count(ranks == rank)
    enddo

    call ESMF_VMGetCurrent(vm, rc=status)
    call ESMF_VMGet(vm, mpiCommunicator=comm, __RC__)
    call MPI_AllToALL(counts_send, 1, MPI_INTEGER, counts_recv, 1, MPI_INTEGER, comm, ierror)
 
    disp_send = 0
    do rank = 1, GigaTrajInternalPtr%npes-1
       disp_send(rank+1) = disp_send(rank)+ counts_send(rank)
    enddo
    disp_recv = 0
    do rank = 1, GigaTrajInternalPtr%npes-1
       disp_recv(rank+1) = disp_recv(rank)+ counts_recv(rank)
    enddo

    ! re-arranged lats lons, and ids
    tmp_position = disp_send
    allocate(lons_send(num_parcels))
    allocate(lons_recv(sum(counts_recv)))
    allocate(pos(num_parcels))
    do i = 1, num_parcels
       rank   = ranks(i)
       pos(i) = tmp_position(rank+1) +1
       lons_send(pos(i)) = lons(i)
       tmp_position(rank+1) = tmp_position(rank+1) + 1
    enddo

    call MPI_AllToALLv(lons_send, counts_send, disp_send, MPI_REAL, lons_recv, counts_recv, disp_recv, MPI_REAL, comm, ierror)

!---------------------
!step 5) Interpolate the data ( horiontally and vertically) and send back where they are from
!---------------------

    allocate(U_recv(sum(counts_recv)), source = rank*1.0)
    allocate(U_send(num_parcels), source = -1.0)
    !
    ! Horizontal and vertical interpolator here
    !
    call MPI_AllToALLv(U_recv, counts_recv, disp_recv, MPI_REAL, U_send, counts_send, disp_send, MPI_REAL, comm, ierror)


!---------------------
!step 6) Rearrange data ( not necessary if ids was rearranged ins step 4)
    allocate(U(num_parcels))
    U(:) = U_send(pos(:))
  
    !deallocate(U_latlon_halo)
    call ESMF_FieldDestroy(U_field)

    print*," Great, I am still alive" 
    RETURN_(ESMF_SUCCESS)
  end subroutine Run

end module GigaTraj_GridCompMod
