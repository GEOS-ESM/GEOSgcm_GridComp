#include "MAPL_Generic.h"

module GEOS_RouteGridCompMod

!BOP
! !MODULE: GEOS_Route -- child to the "Land" gridded component.  

! !DESCRIPTION:
!   {\tt GEOS\_Route} is a gridded component to route total runoff produced in
!   {\tt GEOS\_Catch} (RUNOFF in {\tt GEOScatch\_GridComp} or {\tt GEOScatchCN\_GridComp}) through a 290,188   
!   watersheds globally (excluding submerged catchments (watersheds) from the global list of 291,284. 
!   All of its calculations are done on Pfafstetter watershed space. {\tt GEOS\_Route} has no children. \\
!
!   IMPORTS   : RUNOFF \\

! !USES: 

  use ESMF
  use MAPL_Mod
  use MAPL_ConstantsMod
  use ROUTING_MODEL,          ONLY: river_routing_lin, river_routing_hyd, ROUTE_DT
  use reservoirMod,           ONLY: RES_STATE, Reservoir
  use catch_constants,        ONLY: N_pfaf_g => CATCH_N_PFAFS

  use, intrinsic :: iso_c_binding
  
  implicit none

  logical, parameter :: use_res  = .True.

  private

  type T_RROUTE_STATE !routing related variables
     private
     type (ESMF_RouteHandle) :: routeHandle
     type (ESMF_Field)       :: field_src
     type (ESMF_Field)       :: field
     type (RES_STATE)        :: reservoir

     integer :: n_pfaf_local
     integer :: nt_global
     integer :: nt_local
     integer :: comm
     integer :: nDes
     integer :: myPe
     integer :: minCatch
     integer :: maxCatch

     real,    allocatable :: areacat(:)  ! m2
     real,    allocatable :: lengsc(:)   ! m
     integer, allocatable :: downid(:) 

     real,    allocatable :: runoff_acc(:) 
     real,    allocatable :: wriver_acc(:)    
     real,    allocatable :: wstream_acc(:)        
     real,    allocatable :: qoutflow_acc(:)
     real,    allocatable :: qsflow_acc(:)  

     real,    allocatable :: lstr(:)         ! m
     real,    allocatable :: qri_clmt(:)     ! m3/s
     real,    allocatable :: qin_clmt(:)     ! m3/s
     real,    allocatable :: qstr_clmt(:)    ! m3/s
     real,    allocatable :: K(:)            
     real,    allocatable :: Kstr(:)         

     integer, allocatable :: re_order(:), to_down(:)
     integer, allocatable :: send_count(:), displ_send(:)
     integer, allocatable :: recv_count(:), displ_recv(:)
     integer              :: total_send, total_recv
  end type T_RROUTE_STATE

  ! Wrapper for extracting internal state
  ! -------------------------------------
  type RROUTE_WRAP
     type (T_RROUTE_STATE), pointer :: PTR => null()
  end type RROUTE_WRAP

  include "mpif.h"

  
! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (ESMF_Config          )            :: CF
    
    type (T_RROUTE_STATE), pointer          :: route_internal_state => null()
    type (RROUTE_wrap)                      :: wrap

    integer      :: RUN_DT
    real         :: DT

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                                     ,&
                          NAME=COMP_NAME			 ,&
                          RC=STATUS )

    VERIFY_(STATUS)

    Iam = trim(COMP_NAME) // 'SetServices'

! -----------------------------------------------------------
! Set the Initialize and Run entry points
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run,        RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,   Finalize,   RC=status)
    VERIFY_(status)
    
! -----------------------------------------------------------
! Get the configuration
! -----------------------------------------------------------
! 
    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)
!
! -----------------------------------------------------------
! Get the intervals
! -----------------------------------------------------------
!
    call ESMF_ConfigGetAttribute (CF, DT                         ,&
                                  Label="RUN_DT:"                ,&
                                  RC=STATUS)

    VERIFY_(STATUS)

    RUN_DT = nint(DT)

! -----------------------------------------------------------
! At the moment, this will refresh when the land parent 
! needs to refresh.

    call ESMF_ConfigGetAttribute ( CF, DT, Label=trim(COMP_NAME)//"_DT:", &
         default=DT, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Set the state variable specs.
! -----------------------------------------------------------

!BOS

! -----------------------------------------------------------
!   Import States
! -----------------------------------------------------------
! WY note: Here TileOnly is on tile space

    call MAPL_AddImportSpec(GC,                            &
         LONG_NAME          = 'runoff_total_flux'         ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RUNOFF'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC ) 

!!!!!!!!!!!!!!!!
! Internal
!!!!!!!!!!!!!!
! WY note: Here TileOnly is on *Pfafstetter* catchment space

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_local_stream',&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WSTREAM'                  ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         _RC )

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_river' ,&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WRIVER'                   ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         _RC  )


    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'reservoir_storage' ,&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WRES'                   ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         _RC  )

!!!!!!!!!!!!!!!!
! Export
!!!!!!!!!!!!!!!
! WY note: Here TileOnly is on *Pfafstetter* catchment space

    call MAPL_AddExportSpec(GC,                        &
         LONG_NAME          = 'transfer_of_moisture_from_stream_variable_to_river_variable' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QSFLOW'                   ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'transfer_of_river_water_from_upstream_catchments' ,&
         UNITS              = 'm+3 s-1'                   ,&
         SHORT_NAME         = 'QINFLOW'                   ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         _RC )

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'transfer_of_river_water_to_downstream_catchments' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QOUTFLOW'                 ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )

    call MAPL_AddExportSpec(GC,                           &
         LONG_NAME          = 'reservoir_discharge'      ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QRES'                     ,&
         DIMS               = MAPL_DimsTileOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         _RC )
!EOS

    call MAPL_TimerAdd(GC,    name="-RRM" ,RC=STATUS)
    VERIFY_(STATUS)



! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( route_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => route_internal_state
    
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)

! All done
!---------
    
    call MAPL_GenericSetServices(GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

! -----------------------------------------------------------
! INITIALIZE -- Initialize method for the route component
! -----------------------------------------------------------

  subroutine INITIALIZE (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:
! -----------------------------------------------------------

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

! -----------------------------------------------------------
! ErrLog Variables
! -----------------------------------------------------------

    character(len=ESMF_MAXSTR)          :: IAm="Initialize"
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------
    type (ESMF_VM) :: VM
    integer        :: comm
    integer        :: nDEs
    integer        :: myPE
    integer        :: beforeMe, minCatch, maxCatch, i
    integer        :: n_pfaf_local, nt_global

    type(ESMF_Grid)     :: tileGrid
    type(ESMF_Grid)     :: newTileGrid, catch_grid

    type(MAPL_MetaComp), pointer   :: MAPL
    type(MAPL_LocStream)           :: locstream, catch_LocStream

    character(len=ESMF_MAXSTR)     :: RIVER_INPUT_FILE    
    character(len=ESMF_MAXSTR)     :: TILE_PFAF_FILE
    
    type(ESMF_Grid)     :: agrid 
    type(ESMF_DELayout) :: layout
    type(ESMF_DistGrid) :: dist_grid 

    integer, pointer :: ims(:) => NULL()
    integer, pointer :: local_id(:)  => NULL()

    type (T_RROUTE_STATE), pointer :: route => null()
    type (RROUTE_wrap)             :: wrap
    real, allocatable              :: tmp_real(:)
    integer, allocatable           :: tmp_int(:)

    type(ESMF_Time)  :: CurrentTime
    type(ESMF_Alarm) :: CollectWaterAlarm
    type(ESMF_TimeInterval) :: CollectWater_DT, ModelTimeStep
    character(len=3) :: resname
    type(Netcdf4_Fileformatter)  :: formatter 
    integer          :: j,nt_local, mpierr

    ! ------------------
    ! begin

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)
    ! get LocStream
    call MAPL_Get(MAPL, LocStream = locstream, RC=status)
    VERIFY_(STATUS) 
    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)

    route => wrap%ptr
    ! get vm
    ! extract comm
    call ESMF_VMGetCurrent(VM,                                RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM,       mpiCommunicator =comm,   RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_VMGet       (VM, localpet=MYPE, petcount=nDEs,  RC=STATUS)
    VERIFY_(STATUS)


    route%comm = comm
    route%ndes = ndes
    route%mype = mype
 
    allocate(ims(1:ndes))
    ! define catchment space for this processor
    call MAPL_DecomposeDim ( N_pfaf_g,ims,ndes ) ! ims(mype+1) gives the size of my partition
    ! myPE is 0-based!
    beforeMe = sum(ims(1:mype))
    minCatch = beforeMe + 1
    maxCatch = beforeMe + ims(myPe+1)
    
    ! Get grid info from the gridcomp
    call ESMF_GridCompGet(gc, grid=agrid, rc=status)
    VERIFY_(status)
    call ESMF_GridGet(agrid,   distGrid=dist_grid, _RC)
    call ESMF_DistGridGet(dist_grid, delayout=layout, _RC) 

    call MAPL_LocStreamGet(locstream, &
         tileGrid=tilegrid, nt_local=nt_local, nt_global=nt_global, &
         LOCAL_ID=local_id, RC=status)
    VERIFY_(STATUS)

    route%nt_global   = nt_global
    route%nt_local    = nt_local       
    n_pfaf_local      = maxCatch-minCatch+1
    route%n_pfaf_local= n_pfaf_local  
    route%minCatch    = minCatch
    route%maxCatch    = maxCatch 

    call MAPL_GetResource (MAPL, RIVER_INPUT_FILE,  label = 'RIVER_INPUT_FILE:',    default = '../input/river_input.nc', RC=STATUS )
    if (MAPL_AM_I_Root()) then
       call formatter%open(RIVER_INPUT_FILE, PFIO_READ, _RC)
    endif
    allocate(tmp_real(n_pfaf_g))
    allocate(tmp_int(n_pfaf_g))

    ! read areacat
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('area_catch', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%areacat(n_pfaf_local), source = tmp_real(minCatch:maxCatch))
    route%areacat=route%areacat*1.e6

   !read lengsc
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('lengsc', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%lengsc(n_pfaf_local), source = tmp_real(minCatch:maxCatch))
    route%lengsc=route%lengsc*1.e3

   !read downid
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('downid', tmp_int(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_int,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%downid(n_pfaf_local), source = tmp_int(minCatch:maxCatch))

   !read lstr
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('lstr', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%lstr(n_pfaf_local), source = tmp_real(minCatch:maxCatch))
    route%lstr = route%lstr*1.e3

    !read K
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('K', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%K(n_pfaf_local), source = tmp_real(minCatch:maxCatch))

    !read Kstr
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('Kstr', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%Kstr(n_pfaf_local), source = tmp_real(minCatch:maxCatch))

    !read qri_clmt
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('qri_clmt', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%qri_clmt(n_pfaf_local), source = tmp_real(minCatch:maxCatch))

    !read qin_clmt
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('qin_clmt', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%qin_clmt(n_pfaf_local), source = tmp_real(minCatch:maxCatch))

  !read qstr_clmt
    if (MAPL_AM_I_Root()) then
       call formatter%get_var('qstr_clmt', tmp_real(:), _RC)
    endif
    call MAPL_CommsBcast(layout, tmp_real,   n_pfaf_g,  MAPL_Root, status)
    allocate(route%qstr_clmt(n_pfaf_local), source = tmp_real(minCatch:maxCatch))

    deallocate(tmp_real, tmp_int)

    if (MAPL_AM_I_Root()) then
       call formatter%close()
    endif

    allocate(route%runoff_acc(nt_local), source = 0.)

    ! accumulated variables for output 
    allocate(route%wriver_acc(n_pfaf_local),route%wstream_acc(n_pfaf_local),route%qoutflow_acc(n_pfaf_local),route%qsflow_acc(n_pfaf_local))
    route%wriver_acc  =0.
    route%wstream_acc =0.
    route%qoutflow_acc=0.
    route%qsflow_acc  =0.

    !Initial reservoir module
    route%reservoir = Reservoir(layout, trim(RIVER_INPUT_FILE), N_pfaf_g,n_pfaf_local, minCatch,maxCatch,use_res, _RC)
    allocate(route%reservoir%qres_acc(n_pfaf_local), source =0.0)

    if(mapl_am_I_root()) print *,"reservoir init success" 

    call create_catchment_grid(catch_grid, catch_locstream, _RC)

    call MAPL_GetResource (MAPL, TILE_PFAF_File, label = 'TILE_PFAF_FILE:',  default = '../input/tile_pfaf.nc4', RC=STATUS )
    call create_mapping_handler(trim(TILE_PFAF_File), _RC)

    call MAPL%grid%set(catch_grid, _RC)

    call ESMF_GridCompSet(gc, grid=catch_grid, RC=status)
    VERIFY_(STATUS)

    call MAPL_set(MAPL, locstream = catch_locstream, rc=status)
    VERIFY_(STATUS)

    call setup_exchange_water(_RC)
    
    call ESMF_TimeIntervalSet(CollectWater_DT, s=ROUTE_DT, rc=status)
    VERIFY_(status)
    call ESMF_ClockGet(clock, currTime=CurrentTime, timeStep=ModelTimeStep, rc=status) 
    CollectWaterAlarm = ESMF_AlarmCreate(                                       &
         clock,                                                                 &
         name='CollectWater',                                                   &
         ringTime=CurrentTime - ModelTimeStep+CollectWater_DT,                  &
         ringInterval=CollectWater_DT,                                          &
         ringTimeStepCount=1,                                                   &
         sticky=.false.,                                                        &
         rc=status                                                              &
         )
    VERIFY_(status)
 

    deallocate(ims)
    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)
    RETURN_(ESMF_SUCCESS)

  contains

    subroutine create_catchment_grid(catch_Grid, catch_Locstream, rc)
      type (ESMF_Grid), intent(out)  :: catch_grid
      type (MAPL_LocStream), intent(out) :: catch_Locstream
      integer, optional, intent(out) :: rc
      integer :: status
      real(kind=8), pointer :: centers(:,:)
      ! create catchment grid and it is tile space
      catch_Grid = ESMF_GridCreate(       &
           name='CATCHMENT_GRID',         &
           countsPerDEDim1=IMs,           &
           countsPerDEDim2=[1],           &
           indexFlag=ESMF_INDEX_DELOCAL,  &
           coordDep1 = (/1,2/),           &
           coordDep2 = (/1,2/),           &
           gridEdgeLWidth = (/0,0/),      &
           gridEdgeUWidth = (/0,0/),      &
           _RC)
      ! coord and centers are required for a valid grid,
      ! even if their values don't make sense;
      ! later on, the coord will be set to catchment's lat lon.
      call ESMF_GridAddCoord(catch_Grid, _RC)
      _VERIFY(status)

      call ESMF_GridGetCoord(catch_Grid, coordDim=1, localDE=0, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=centers, _RC)
      centers = 0 ! ?? just assign
      call ESMF_GridGetCoord(catch_Grid, coordDim=2, localDE=0, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=centers, _RC)
      centers = 0 ! 

      call MAPL_LocstreamCreate(catch_Locstream, catch_Grid, rc=status)
      _VERIFY(STATUS)
      _RETURN(_SUCCESS)
    end subroutine create_catchment_grid

    ! ----------------------------------------------------

    subroutine create_mapping_handler(tile_pfaf_file, rc)

      character(*),      intent(in)  :: tile_pfaf_file
      integer, optional, intent(out) :: rc

      integer :: status, nWeights, nLocal_weights
      integer, allocatable :: global_id(:)
      logical, allocatable :: mask(:)
      integer, allocatable :: srcIndices(:), positions(:), factorIndexList(:,:)
      real,    allocatable :: weights(:), global_frac(:)
      integer, allocatable :: local_src(:), local_dst(:), global_src(:), global_dst(:)
      integer :: unit
      type(Netcdf4_Fileformatter) :: formatter
      type(Filemetadata)          :: meta


      ! create source for orignal tile space
      route%field_src = ESMF_FieldCreate(grid=tilegrid, typekind=ESMF_TYPEKIND_R4, _RC)

      call MAPL_LocStreamGet(catch_LocStream, TILEGRID=newtilegrid, RC=status)
      route%field     = ESMF_FieldCreate(grid=newtilegrid, typekind=ESMF_TYPEKIND_R4,_RC)

      if (MAPL_AM_I_ROOT()) then
         call formatter%open(tile_pfaf_file, PFIO_READ, _RC)
         meta     = formatter%read(rc=status)
         nweights = meta%get_dimension('tile')
      endif
      call MAPL_CommsBcast(layout, nWeights, 1, MAPL_Root, status)
      allocate(global_src(nWeights), global_dst(nWeights), global_frac(nWeights))
      if (MAPL_AM_I_ROOT()) then
         call formatter%get_var('tile_id',    global_src,  _RC)
         call formatter%get_var('pfaf_index', global_dst,  _RC)
         call formatter%get_var('pfaf_frac',  global_frac, _RC)
      endif
      call MAPL_CommsBcast(layout, global_src, nWeights, MAPL_Root, status)
      call MAPL_CommsBcast(layout, global_dst, nWeights, MAPL_Root, status)
      call MAPL_CommsBcast(layout, global_frac,nWeights, MAPL_Root, status)
      allocate(mask(nWeights))
      mask = minCatch <= global_dst(:) .and. global_dst(:) <=maxCatch
      local_src = pack(global_src(:),  mask)
      local_dst = pack(global_dst(:),  mask)
      weights   = pack(global_frac(:), mask)
      deallocate(global_src, global_dst, global_frac)

      !UNIT = GETFILE(route_file, form='FORMATTED', _RC)
      !call READ_PARALLEL(layout, nWeights, UNIT=UNIT, _RC)
      !allocate(mapping(3, nWeights))
      !call READ_PARALLEL(layout, mapping(:,:), unit=UNIT, _RC)
      !call FREE_FILE(unit)

      ! get local  number of weight         
      !allocate(mask(nWeights))
      !mask = minCatch <= mapping(2,:) .and. mapping(2,:) <=maxCatch
      !local_src = nint(pack(mapping(1,:), mask))
      !local_dst = nint(pack(mapping(2,:), mask))
      !weights   = pack(mapping(3,:), mask)
      !deallocate(mapping)

      ! ESMF use global indices increasing with mpi_rank, no mask here for tile grid 
      allocate(global_id(nt_global))
      call ESMFL_Fcollect(tilegrid, global_id, local_id, _RC)

      ! mapping form local to global index
      nLocal_weights = count(mask)
      allocate(srcIndices(nLocal_weights))
      do i =1, nLocal_weights
         positions = pack([(j, j=1, nt_global)], global_id == local_src(i))
         srcIndices(i) = positions(1)
      enddo

      allocate(factorIndexList(2, nlocal_weights))
      factorIndexList(1,:) = srcIndices
      factorIndexList(2,:) = local_dst
      call ESMF_FieldSMMStore(route%field_src, route%field, &
           routeHandle=route%routeHandle, &
           factorList=weights, &
           factorIndexList= factorIndexList, &
           _RC)

      ! testing
      !call ESMF_FieldGet(route%field_src, farrayPtr=ptr, rc=status)
      !VERIFY_(STATUS)
      !ptr(:) = 1.

      !call ESMF_FieldSMM(route%field_src, route%Field, &
      !             route%routeHandle, rc=status)
      !VERIFY_(STATUS)

      !call ESMF_FieldGet(route%Field, farrayPtr=ptr, rc=status)
      !VERIFY_(STATUS)
      !! after remapping, all values should be 1
      !if (route%mype == 20) then
      !  do i = 1, n_pfaf_local
      !    print* , 'ptr(i): ',i,  ptr(i)
      !  enddo
      !endif
      _RETURN(_SUCCESS)
      
    end subroutine create_mapping_handler

    ! -------------------------------------------------------------------------------------------
    
    subroutine setup_exchange_water(rc)
      integer, optional, intent(out) :: rc
      integer :: pf, down_id, rank
      integer, allocatable :: cat_to_ranks_global(:)
      integer, allocatable :: cat_to_ranks_local(:)
      integer :: status, mype,mpierr, k
      integer, allocatable :: to_downstream(:), positions(:)

      mype = route%mype
      allocate(cat_to_ranks_local(route%n_pfaf_local))
      allocate(cat_to_ranks_global(n_pfaf_g))
      cat_to_ranks_local = mype

      call ESMFL_FCollect(newtilegrid, cat_to_ranks_global, cat_to_ranks_local, _RC)

      allocate(route%send_count(route%ndes), source=0)
      allocate(route%recv_count(route%ndes), source=0)

      do pf = 1, route%n_pfaf_local
         down_id = route%downid(pf)
         if (down_id == -1) cycle
         ! down stream is not in the process
         if (down_id < minCatch .or. maxCatch < down_id) then
            rank = cat_to_ranks_global(down_id)
            route%send_count(rank+1) = route%send_count(rank+1) + 1
         endif
      enddo

      call MPI_AlltoALL(route%send_count, 1, MPI_INTEGER,  &
           route%recv_count, 1, MPI_INTEGER, route%comm, status)
      VERIFY_(STATUS)

      allocate(route%displ_send(ndes), source = 0)
      allocate(route%displ_recv(ndes), source = 0)
      do rank = 1, ndes-1
         route%displ_send(rank+1) = route%displ_send(rank) + route%send_count(rank)
         route%displ_recv(rank+1) = route%displ_recv(rank) + route%recv_count(rank)
      enddo
      k = sum(route%send_count)
      route%total_send = k
      allocate(to_downstream(k))
      allocate(route%re_order(k))
      allocate(positions(ndes), source = route%displ_send )
      positions = positions + 1
      k = 0
      do pf = 1, route%n_pfaf_local
         down_id = route%downid(pf)
         if (down_id == -1) cycle
         if (down_id < minCatch .or. maxCatch < down_id) then
            rank = cat_to_ranks_global(down_id)
            k    = k + 1
            to_downstream(k) = down_id
            route%re_order(k)= positions(rank+1)
            positions(rank+1)= positions(rank+1) + 1
         endif
      enddo
      if (route%total_send > 0) then
         to_downstream = to_downstream(route%re_order)
      endif
      k = sum(route%recv_count)
      route%total_recv = k
      allocate(route%to_down(k))
      call MPI_AllToAllv(to_downstream, route%send_count, route%displ_send, MPI_Integer,  &
           route%to_down, route%recv_count, route%displ_recv, MPI_INTEGER, &
           route%comm, status)

      ! testing
      !do i = 1, route%total_recv
      !   down_id = route%to_down(i)
      !   if (down_id < minCatch .or. maxCatch < down_id) then
      !     _ASSERT(.false., "Got the down_id that does not belong to me")
      !   endif
      !enddo

      _VERIFY(STATUS)
      _RETURN(_SUCCESS)
    end subroutine setup_exchange_water

  end subroutine INITIALIZE

  ! --------------------------------------------------------------------------------  
  
  subroutine RUN (GC,IMPORT, EXPORT, CLOCK, RC )
    
! -----------------------------------------------------------
! !ARGUMENTS:
! -----------------------------------------------------------

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

! -----------------------------------------------------------
! ErrLog Variables
! -----------------------------------------------------------

    character(len=ESMF_MAXSTR)          :: IAm='RUN'
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------
    type (ESMF_State       )            :: INTERNAL
    type (MAPL_MetaComp),     pointer   :: MAPL
!    type(ESMF_Alarm)                    :: ALARM
    type (ESMF_Config )                 :: CF

! -----------------------------------------------------
! IMPORT pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer :: RUNOFF_SRC0   

! -----------------------------------------------------
! INTERNAL pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer :: AREACAT
    real, dimension(:), pointer :: LENGSC
    real, dimension(:), pointer :: DNSTR
    real, dimension(:), pointer :: WSTREAM
    real, dimension(:), pointer :: WRIVER
    real, dimension(:), pointer :: WRES
    real, dimension(:), pointer :: LRIVERMOUTH
    real, dimension(:), pointer :: ORIVERMOUTH

! -----------------------------------------------------
! EXPORT pointers 
! -----------------------------------------------------

    real, dimension(:), pointer :: QSFLOW
    real, dimension(:), pointer :: QINFLOW
    real, dimension(:), pointer :: QOUTFLOW
    real, dimension(:), pointer :: QRES
  
! Time attributes and placeholders

!    type(ESMF_Time) :: CURRENT_TIME

! Others

    type(ESMF_Grid)                    :: TILEGRID
    type (MAPL_LocStream)              :: LOCSTREAM
 
    integer                            :: n_pfaf_local, N_CYC

    integer                                  :: I
    REAL                                     :: HEARTBEAT 
    REAL, ALLOCATABLE, DIMENSION(:)          :: RUNOFF_ACT,AREACAT_ACT,& 
         LENGSC_ACT, QSFLOW_ACT,QOUTFLOW_ACT,QRES_ACT,QOUT_CAT
    type(ESMF_Field) :: runoff_src

    integer                                :: ndes, mype
    type (T_RROUTE_STATE), pointer         :: route 
    type (RROUTE_wrap)                     :: wrap

    integer :: nt_global,nt_local
    integer,save :: nstep_per_day

    type(ESMF_Time) :: CurrentTime, nextTime
    integer :: YY,MM,DD,HH,MMM,SS,YY_next,MM_next,DD_next
    character(len=4) :: yr_s
    character(len=2) :: mon_s,day_s

    real,             pointer     :: arrayPtr(:)
    type (RES_STATE), pointer     :: res 

    real,             allocatable :: WTOT_BEFORE(:),QINFLOW_LOCAL(:)
   
    type(ESMF_Alarm) :: CollectWaterAlarm 
    
    ! ------------------
    ! begin    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)
    route => wrap%ptr

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS) 
    Iam = trim(COMP_NAME) // "::RUN"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)
    call MAPL_Get(MAPL, HEARTBEAT = HEARTBEAT, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ClockGetAlarm(clock, 'CollectWater', CollectWaterAlarm, _RC)
    !if (mapl_am_I_root()) print *, "HEARTBEAT=",HEARTBEAT 

! internal
    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL,  _RC)
    call MAPL_GetPointer(INTERNAL, WRIVER, 'WRIVER',   _RC )
    call MAPL_GetPointer(INTERNAL, WSTREAM,'WSTREAM',  _RC)
    call MAPL_GetPointer(INTERNAL, WRES,    'WRES',    _RC)
 
! export
    call MAPL_GetPointer(EXPORT, QSFLOW,   'QSFLOW',   _RC)
    call MAPL_GetPointer(EXPORT, QINFLOW,  'QINFLOW' , _RC)
    call MAPL_GetPointer(EXPORT, QOUTFLOW, 'QOUTFLOW', _RC)
    call MAPL_GetPointer(EXPORT, QRES,     'QRES',     _RC)

! Start timers
! ------------
    call MAPL_TimerOn(MAPL,"RUN")
! Get parameters from generic state
! ---------------------------------

    !   call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    !   VERIFY_(STATUS)
    
! get pointers to inputs variables
! ----------------------------------

    ndes          =  route%ndes
    mype          =  route%mype  
    n_pfaf_local  =  route%n_pfaf_local  
    nt_global     =  route%nt_global  
    nt_local      =  route%nt_local
    res           => route%reservoir

    ! get the field from IMPORT
    call ESMF_StateGet(IMPORT, 'RUNOFF', field=runoff_src, RC=STATUS)
    VERIFY_(STATUS)    
    call ESMF_FieldGet(runoff_src, farrayPtr=RUNOFF_SRC0, rc=status)   
    VERIFY_(STATUS) 
    call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)
    VERIFY_(STATUS)   
    call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
    VERIFY_(STATUS)    
    call MAPL_TimerOn  ( MAPL, "-RRM" )

    ! For efficiency, the time step to call the river routing model is set at ROUTE_DT 
    N_CYC = ROUTE_DT/HEARTBEAT    
    if (ESMF_AlarmIsRinging(CollectWaterAlarm)) then

       !accumulates runoff
       route%runoff_acc = (route%runoff_acc + RUNOFF_SRC0)/real (N_CYC)

       !Gets time used for output and restart 
       call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
       call ESMF_TimeGet(CurrentTime, yy=YY, mm=MM, dd=DD, h=HH, m=MMM, s=SS, rc=status)  
       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=status)
       call ESMF_TimeGet(nextTime, yy=YY_next, mm=MM_next, dd=DD_next, rc=status) 
       write(yr_s, '(I4.4)')YY
       write(mon_s,'(I2.2)')MM
       write(day_s,'(I2.2)')DD

       ! redistribute runoff from tile space to catchment space
       call ESMF_FieldGet(route%field_src, farrayPtr=arrayPtr, rc=status)
       VERIFY_(STATUS)
       ArrayPtr = route%runoff_acc(:)
       call ESMF_FieldSMM(srcField=route%field_src, dstField=route%Field, &
            routeHandle=route%routeHandle, rc=rc)
       call ESMF_FieldGet(route%field, farrayPtr=arrayPtr, rc=status)
       VERIFY_(STATUS)
       RUNOFF_ACT = arrayPtr * route%areacat/1000.


       ! Prepares to conduct routing model
       allocate (AREACAT_ACT (n_pfaf_local))       
       allocate (LENGSC_ACT  (n_pfaf_local))
       allocate (QSFLOW_ACT  (n_pfaf_local))
       allocate (QOUTFLOW_ACT(n_pfaf_local),QRES_ACT(n_pfaf_local),QOUT_CAT(n_pfaf_local))  

       QRES_ACT=0.
       LENGSC_ACT =route%lengsc/1.e3 !m->km
       AREACAT_ACT=route%areacat/1.e6 !m2->km2

       allocate(WTOT_BEFORE(n_pfaf_local))
       WTOT_BEFORE=WSTREAM + WRIVER + WRES

       ! Call river_routing_model
       ! ------------------------     
       !CALL RIVER_ROUTING_LIN  (n_pfaf_local, RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,  &
       !     WSTREAM_ACT,WRIVER_ACT, QSFLOW_ACT,QOUTFLOW_ACT) 

       CALL RIVER_ROUTING_HYD  (n_pfaf_local, &
            RUNOFF_ACT, route%lengsc, route%lstr, &
            route%qstr_clmt, route%qri_clmt, route%qin_clmt, &
            route%K, route%Kstr, &
            WSTREAM,WRIVER, &
            QSFLOW_ACT,QOUTFLOW_ACT)  
       ! Call reservoir module        

       call res%calc( QOUTFLOW_ACT, QRES_ACT, WRES, real(route_dt), _RC)
       QOUT_CAT = QOUTFLOW_ACT              
       where(res%active_res==1) QOUT_CAT=QRES_ACT
       allocate(QINFLOW_LOCAL(n_pfaf_local))
       call exchange_water(QOUT_CAT, QINFLOW_LOCAL, _RC)
       WRIVER = WRIVER + QINFLOW_LOCAL*real(route_dt)

       ! Check balance if needed
       !call check_balance(route,n_pfaf_local,nt_local,runoff_acc,WRIVER_ACT,WSTREAM_ACT,WTOT_BEFORE,RUNOFF_ACT,QINFLOW_LOCAL,QOUT_CAT,FirstTime,yr_s,mon_s)

       ! Update accumulated variables for output
       nstep_per_day = 86400/route_dt

       route%wriver_acc   = route%wriver_acc   + WRIVER      /real(nstep_per_day)
       route%wstream_acc  = route%wstream_acc  + WSTREAM     /real(nstep_per_day)
       route%qoutflow_acc = route%qoutflow_acc + QOUTFLOW_ACT/real(nstep_per_day)
       route%qsflow_acc   = route%qsflow_acc   + QSFLOW_ACT  /real(nstep_per_day)

       res%qres_acc       = res%qres_acc       + QRES_ACT    /real(nstep_per_day)       

       if (associated(QSFLOW))    QSFLOW  = QSFLOW_ACT
       if (associated(QOUTFLOW)) QOUTFLOW = QOUTFLOW_ACT
       if (associated(QRES)) then
          _ASSERT(use_res, "Set use_res be true to get QRES export")
          QRES = QRES_ACT
       endif
       deallocate(RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,QOUTFLOW_ACT,QINFLOW_LOCAL,QSFLOW_ACT,WTOT_BEFORE,QRES_ACT,QOUT_CAT)

       route%runoff_acc = 0.
       if (HH == 23) then
          route%wriver_acc   = 0.
          route%wstream_acc  = 0.
          route%qoutflow_acc = 0.
          route%qsflow_acc   = 0.
          res%qres_acc       = 0.
       endif
    else

       route%runoff_acc = route%runoff_acc + RUNOFF_SRC0

    endif ! alarm

    ! All done
    ! --------
    call MAPL_TimerOff ( MAPL, "-RRM" ) 
    call MAPL_TimerOff(MAPL,"RUN")
    !call MPI_Barrier(route%comm, mpierr)

    RETURN_(ESMF_SUCCESS)
    
  contains

    subroutine exchange_water(cat_out, cat_in, rc)
      real, intent(in)  :: cat_out(:)
      real, intent(out) :: cat_in(:)
      integer, optional, intent(out) :: rc

      integer :: status
      real, allocatable :: send_data(:), recv_data(:)
      integer :: pf, down_id, k, pos, total_send, total_recv, mpierr

      cat_in = 0.0
      total_send = route%total_send
      allocate(send_data(total_send))
      pos = 0
      do pf = 1, route%n_pfaf_local
         down_id = route%downid(pf)
         if (down_id == -1) cycle
         if (route%minCatch <= down_id .and. down_id <= route%maxCatch) then
            ! from local
            k = down_id - route%minCatch + 1
            cat_in(k) = cat_in(k) + cat_out(pf)
         else ! from the other process
            pos = pos + 1
            send_data(pos) = cat_out(pf)
         endif
      enddo
      if (total_send > 0) then
         send_data = send_data(route%re_order)
      endif
      total_recv = route%total_recv
      allocate(recv_data(total_recv))
      call MPI_AllToAllv(send_data, route%send_count, route%displ_send, MPI_REAL, &
           recv_data, route%recv_count, route%displ_recv, MPI_REAL, &
           route%comm, mpierr)
      do i = 1, total_recv
         down_id = route%to_down(i)
         k = down_id - route%minCatch + 1
         cat_in(k) = cat_in(k) + recv_data(i)
      enddo
      RETURN_(ESMF_SUCCESS)
      
    end subroutine exchange_water
    
  end subroutine RUN

  ! ----------------------------------------------------------------------------------------

  subroutine Finalize(gc, import, export, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Clean-up.
    !EOP
    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap

    ! Begin...
    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"

    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap, _RC)
    route => wrap%ptr    

    CALL ESMF_FieldSMMRelease(routeHandle=route%routeHandle, _RC)

    ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, _RC)
    ! End
    RETURN_(ESMF_SUCCESS)

  end subroutine Finalize

end module GEOS_RouteGridCompMod

! ======================= EOF =========================================================

