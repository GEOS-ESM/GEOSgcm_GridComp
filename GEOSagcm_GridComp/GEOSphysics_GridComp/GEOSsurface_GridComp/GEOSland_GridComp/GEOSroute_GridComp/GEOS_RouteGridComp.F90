
#include "MAPL_Generic.h"

#ifndef RUN_FOR_REAL
#define MAPL_DimsCatchOnly MAPL_DimsTileOnly
#endif

!=============================================================================
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
!   INTERNALS : AREACAT, LENGSC2, DNSTR, WSTREAM, WRIVER, LRIVERMOUTH, ORIVERMOUTH \\
!   EXPORTS   : QSFLOW, QINFLOW, QOUTFLOW \\

! !USES: 

  use ESMF
  use MAPL_Mod
  use MAPL_ConstantsMod
  use ROUTING_MODEL,          ONLY:     &
       river_routing, ROUTE_DT
#if 0
  USE catch_constants, ONLY:          &
       N_CatG => N_Pfaf_Catchs  
#endif
  
  implicit none
  integer, parameter :: N_CatG = 291284
  private

  type T_RROUTE_STATE
     private
     type (ESMF_RouteHandle) :: routeHandle
     type (ESMF_Field)       :: field
     integer :: nTiles
     integer :: comm
     integer :: nDes
     integer :: myPe
     integer :: minCatch
     integer :: maxCatch
     integer, pointer :: pfaf(:) => NULL()
     real,    pointer :: tile_area(:) => NULL()
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
    
    type (T_RROUTE_STATE), pointer         :: route_internal_state => null()
    type (RROUTE_wrap)                     :: wrap

    integer      :: RUN_DT
    real         :: DT

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                                 ,&
                          NAME=COMP_NAME			   ,&
                          RC=STATUS )

    VERIFY_(STATUS)

    Iam = trim(COMP_NAME) // 'SetServices'

! -----------------------------------------------------------
! Set the Initialize and Run entry points
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, RC=STATUS)
    VERIFY_(STATUS)

!------------------------------------------------------------
! Set generic final method 
!------------------------------------------------------------

    
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

    call MAPL_AddImportSpec(GC,                          &
         LONG_NAME          = 'runoff_flux'               ,&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'RUNOFF'                    ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)

! -----------------------------------------------------------
!   INTERNAL STATE
! -----------------------------------------------------------

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'area_of_catchment'        ,&
         UNITS              = 'km+2'                     ,&
         SHORT_NAME         = 'AREACAT'                  ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 
 
    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'length_of_channel_segment',&
         UNITS              = 'km+2'                     ,&
         SHORT_NAME         = 'LENGSC'                   ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'index_of_downtream_catchment',&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'DNSTR'                    ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_local_stream',&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WSTREAM'                  ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  )
    
    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_river' ,&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WRIVER'                   ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'TileID_of_the_lake_tile_at_the_river_mouth' ,&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'LRIVERMOUTH'              ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'TileID_of_the_ocean_tile_at_the_river_mouth' ,&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'ORIVERMOUTH'              ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   
! -----------------------------------------------------------
!  EXPORT STATE:
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC,                        &
         LONG_NAME          = 'transfer_of_moisture_from_stream_variable_to_river_variable' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QSFLOW'                   ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
                                         RC=STATUS  ) 

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'transfer_of_river_water_from_upstream_catchments' ,&
         UNITS              = 'm+3 s-1'                   ,&
         SHORT_NAME         = 'QINFLOW'                   ,&
         DIMS               = MAPL_DimsCatchOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
                                          RC=STATUS  ) 

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'transfer_of_river_water_to_downstream_catchments' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'QOUTFLOW'                 ,&
         DIMS               = MAPL_DimsCatchOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
                                         RC=STATUS  ) 

!EOS

    call MAPL_TimerAdd(GC,    name="RUN"  ,RC=STATUS)
    VERIFY_(STATUS)
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
    integer        :: beforeMe, minCatch, maxCatch, pf, i
    integer        :: ntiles, nt_global

    type(ESMF_Grid) :: tileGrid
    type(ESMF_Grid) :: newTileGrid
    type(ESMF_Grid) :: catchGrid
    type(ESMF_DistGrid) :: distGrid
    type(ESMF_Field) :: field, field0

    type(MAPL_MetaComp), pointer   :: MAPL
    type(MAPL_LocStream) :: locstream
    
    integer, pointer :: ims(:) => NULL()
    integer, pointer :: pfaf(:) => NULL()
    integer, pointer :: arbSeq(:) => NULL()
    integer, allocatable :: arbIndex(:,:)
    real, pointer :: tile_area_src(:) => NULL()
    real, pointer :: tile_area(:) => NULL()
    real, pointer :: ptr2(:) => NULL()
    
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap

    ! ------------------
    ! begin
    
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

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    route%comm = comm
    route%ndes = ndes
    route%mype = mype
    
    ! define minCatch, maxCatch
    call MAPL_DecomposeDim ( n_catg,ims,ndes ) ! ims(mype+1) gives the size of my partition
    ! myPE is 0-based!
    beforeMe = sum(ims(1:mype))
    minCatch = beforeMe + 1
    maxCatch = beforeMe + ims(myPe+1)
    
    ! get LocStream
    call MAPL_Get(MAPL, LocStream = locstream, RC=status)
    VERIFY_(STATUS)
    ! extract Pfaf (TILEI on the "other" grid)
    call MAPL_LocStreamGet(locstream, tilei=pfaf, OnAttachedGrid=.false., &
         tileGrid=tilegrid, nt_global=nt_global, RC=status)
    VERIFY_(STATUS)
    
    ! exchange Pfaf across PEs

    ntiles = 0
    !loop over total_n_tiles
    do i = 1, nt_global
       pf = pfaf(i)
       if (pf >= minCatch .and. pf <= maxCatch) then ! I want this!
          ntiles = ntiles+1
          !realloc if needed
          arbSeq(ntiles) = i
       end if
    end do ! global tile loop

    distgrid = ESMF_DistGridCreate(arbSeqIndexList=arbSeq, rc=status)
    VERIFY_(STATUS)

    newTileGRID = ESMF_GridEmptyCreate(rc=status)
    VERIFY_(STATUS)
         
    allocate(arbIndex(nTiles,1), stat=status)
    VERIFY_(STATUS)

    arbIndex(:,1) = arbSeq

    call ESMF_GridSet(newTileGrid,  &
         name='redist_tile_grid_for_'//trim(COMP_NAME),    &
         distgrid=distgrid, & 
         gridMemLBound=(/1/), &
         indexFlag=ESMF_INDEX_USER, &
         distDim = (/1/), &
         localArbIndexCount=ntiles, &
         localArbIndex=arbIndex, &
         minIndex=(/1/), &
         maxIndex=(/NT_GLOBAL/), &
         rc=status)
    VERIFY_(STATUS)

    deallocate(arbIndex)

    call ESMF_GridCommit(newTileGrid, rc=status)
    VERIFY_(STATUS)


    ! now create a "catch" grid to be the "native" grid for this component
    distgrid = ESMF_DistGridCreate(arbSeqIndexList=(/minCatch:maxCatch/), &
         rc=status)
    VERIFY_(STATUS)

    catchGRID = ESMF_GridEmptyCreate(rc=status)
    VERIFY_(STATUS)

    allocate(arbIndex(ims(myPE+1),1), stat=status)
    VERIFY_(STATUS)

    arbIndex(:,1) = (/minCatch:maxCatch/)
         
    call ESMF_GridSet(catchGrid,  &
         name='catch_grid_for_'//trim(COMP_NAME),    &
         distgrid=distgrid, & 
         gridMemLBound=(/1/), &
         indexFlag=ESMF_INDEX_USER, &
         distDim = (/1/), &
         localArbIndexCount=ims(myPE+1), &
         localArbIndex=arbIndex, &
         minIndex=(/1/), &
         maxIndex=(/N_CatG/), &
         rc=status)
    VERIFY_(STATUS)

    deallocate(arbIndex)

    call ESMF_GridCommit(catchGrid, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSet(gc, grid=catchGrid, RC=status)
    VERIFY_(STATUS)

    call MAPL_LocStreamGet(locstream, TILEAREA = tile_area_src, RC=status)
    VERIFY_(STATUS)

    field0 = ESMF_FieldCreate(grid=tilegrid, datacopyflag=ESMF_DATACOPY_VALUE, &
         farrayPtr=tile_area_src, name='TILE_AREA_SRC', RC=STATUS)
    VERIFY_(STATUS)
    ! create field on the "new" tile grid
    allocate(tile_area(ntiles), stat=status)
    VERIFY_(STATUS)
    field = ESMF_FieldCreate(grid=newtilegrid, datacopyflag=ESMF_DATACOPY_VALUE, &
         farrayPtr=tile_area, name='TILE_AREA', RC=STATUS)
    VERIFY_(STATUS)

    ! create routehandle
    call ESMF_FieldRedistStore(srcField=field0, dstField=field, &
                routehandle=route%routehandle, rc=status)
    VERIFY_(STATUS)
 
    ! redist tile_area
    call ESMF_FieldRedist(srcField=FIELD0, dstField=FIELD, &
         routehandle=route%routehandle, rc=status)
    VERIFY_(STATUS)

    call ESMF_FieldDestroy(field, rc=status)
    VERIFY_(STATUS)
    call ESMF_FieldDestroy(field0, rc=status)
    VERIFY_(STATUS)

#ifdef __GFORTRAN__
    deallocate(tile_area_src)
#endif
    
    ! redist pfaf (NOTE: me might need a second routehandle for integers)

    route%pfaf => arbSeq
    route%ntiles = ntiles
    route%minCatch = minCatch
    route%maxCatch = maxCatch

    allocate(ptr2(ntiles), stat=status)
    VERIFY_(STATUS)
    route%field = ESMF_FieldCreate(grid=newtilegrid, datacopyflag=ESMF_DATACOPY_VALUE, &
         farrayPtr=ptr2, name='RUNOFF', RC=STATUS)
    VERIFY_(STATUS)
    
    deallocate(ims)
    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine INITIALIZE
  
! -----------------------------------------------------------
! RUN -- Run method for the route component
! -----------------------------------------------------------

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

    character(len=ESMF_MAXSTR)          :: IAm="Run"
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------

    type (MAPL_MetaComp),     pointer   :: MAPL
    type (ESMF_State       )            :: INTERNAL
!    type(ESMF_Alarm)                    :: ALARM
    type (ESMF_Config )                 :: CF
    type(ESMF_VM)                       :: VM

! -----------------------------------------------------
! IMPORT pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer :: RUNOFF 

! -----------------------------------------------------
! INTERNAL pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer :: AREACAT
    real, dimension(:), pointer :: LENGSC
    real, dimension(:), pointer :: DNSTR
    real, dimension(:), pointer :: WSTREAM
    real, dimension(:), pointer :: WRIVER
    real, dimension(:), pointer :: LRIVERMOUTH
    real, dimension(:), pointer :: ORIVERMOUTH

! -----------------------------------------------------
! EXPORT pointers 
! -----------------------------------------------------

    real, dimension(:), pointer :: QSFLOW
    real, dimension(:), pointer :: QINFLOW
    real, dimension(:), pointer :: QOUTFLOW
  
! Time attributes and placeholders

!    type(ESMF_Time) :: CURRENT_TIME

! Others

    type(ESMF_Grid)                    :: TILEGRID
    type (MAPL_LocStream)              :: LOCSTREAM
    integer                            :: NTILES, N_CatL, N_CYC
    logical, save                      :: FirstTime=.true.
    real, pointer, dimension(:)    :: tile_area
    integer, pointer, dimension(:) :: pfaf_code

    INTEGER, DIMENSION(:,:), POINTER, SAVE   :: AllActive,DstCatchID 
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: srcProcsID, LocDstCatchID  
    integer, dimension (:),allocatable, SAVE :: GlbActive
    INTEGER, SAVE                            :: N_Active, ThisCycle   
    INTEGER                                  :: Local_Min, Local_Max
    integer                                  :: K, N, I, req
    REAL                                     :: mm2m3, rbuff, HEARTBEAT 
    REAL, ALLOCATABLE, DIMENSION(:)          :: RUNOFF_CATCH, RUNOFF_ACT,AREACAT_ACT,& 
         LENGSC_ACT, WSTREAM_ACT,WRIVER_ACT, QSFLOW_ACT,QOUTFLOW_ACT, runoff_save
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_index
    type(ESMF_Field) :: runoff_src

    integer                                :: ndes, mype
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap

    ! ------------------
    ! begin
    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)

    route => wrap%ptr

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
  
    Iam = trim(COMP_NAME) // "RUN"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, HEARTBEAT = HEARTBEAT, RC=STATUS)
    VERIFY_(STATUS)

! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"RUN")

! Get parameters from generic state
! ---------------------------------

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS) 
    
! get pointers to inputs variables
! ----------------------------------

    ! get the field from IMPORT
    call ESMF_StateGet(IMPORT, 'RUNOFF', field=runoff_src, RC=STATUS)
    VERIFY_(STATUS)

    ! redist RunOff
    call ESMF_FieldRedist(srcField=runoff_src, dstField=route%field, &
                routehandle=route%routehandle, rc=status)
    VERIFY_(STATUS)

    call ESMF_FieldGet(route%field, farrayPtr=RUNOFF, rc=status)
    VERIFY_(STATUS)

    pfaf_code => route%pfaf
    tile_area => route%tile_area
    
! get pointers to internal variables
! ----------------------------------
  
    call MAPL_GetPointer(INTERNAL, AREACAT , 'AREACAT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, LENGSC  , 'LENGSC',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DNSTR   , 'DNSTR'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, WSTREAM , 'WSTREAM', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, WRIVER  , 'WRIVER' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, LRIVERMOUTH, 'LRIVERMOUTH' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ORIVERMOUTH, 'ORIVERMOUTH' , RC=STATUS)
    VERIFY_(STATUS)

! get pointers to EXPORTS
! -----------------------

    call MAPL_GetPointer(EXPORT, QSFLOW,   'QSFLOW'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QINFLOW,  'QINFLOW' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QOUTFLOW, 'QOUTFLOW', RC=STATUS)
    VERIFY_(STATUS)
 
    call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
    VERIFY_(STATUS)    

    call MAPL_TimerOn  ( MAPL, "-RRM" )

    call MAPL_LocStreamGet(LocStream, NT_LOCAL=NTILES, RC=STATUS )
    N_CatL  = size(AREACAT)

!@@    ALLOCATE (pfaf_code (1:NTILES)) ! 9th_coulumn_in_TILFILE

    ! NOTES    :    
    !Need below area and pfaf_index from the .til file (Maybe, they are already in LocStream)
    !
    !	 TILFILE: /discover/nobackup/smahanam/bcs/Heracles-4_3/Heracles-4_3_MERRA-3/CF0090x6C_DE1440xPE0720/CF0090x6C_DE1440xPE0720-Pfafstetter.til
    !	 The 8-line header is followed by 1061481 number of rows.
    !	 do n = 1,475330
    !	        read (10,*)type,area, longitude, latitude, ig, jg, cell_frac, integer,   & 
    !                      pfaf_code, pfaf_index, pfaf_frac
    !	 end do
    !	 
    !	where for each tile:
    !	 (1)    type      [-]      tile type (100-land; 19-lakes; 20-ice)
    !	 (2)    area      [x EarthRadius^2 km2]  tile area
    !	 (3)    longitude [degree] longitude at the centroid of the tile
    !	 (4)    latitude  [degree] latitude at the centroid of the tile
    !	 (5)    ig        [-]      i-index of the AGCM grid cell where the tile is located
    !	 (6)    jg        [-]      j-index of the AGCM grid cell where the tile is located
    !	 (7)    cell_frac [-]      fraction of the AGCM grid cell  
    !    (8)    integer            some integer that matters only for OGCM tiles, I suppose.  
    !	 (9)    pfaf_code [-]      catchment index (1-291284) after sorting Pfafstetter codes in ascending order 
    !	 (10)    pfaf_index[-]      catchment index (1-290188) after sorting Pfafstetter codes 
    !                                  and removing submerged in ascending order 
    !    (11)   pfaf_frac [-]      fraction of the pfafstetter catchment

    !call MAPL_LocStreamGet(LocStream, 9th_coulumn_in_TILFILE=pfaf_code, RC=STATUS )

    Local_Min = route%minCatch
    Local_Max = route%maxCatch
 
    FIRST_TIME : IF (FirstTime) THEN

       ! Pfafstetter catchment Domain Decomposition :         
       ! --------------------------------------------

       ! AllActive      : Processor(s) where the catchment is active  (identical in any processor). 
       ! srcProcsID     : For all active catchments anywhere tells which processor is the principal owner of the catchment (identical in any processor).
       ! DstCatchID  : 2-D array contains downstream catchID and downstream processor (identical in any processor)
       ! LocDstCatchID  : Downstream catchID when for catchments that are local to the processor.

       ndes = route%ndes
       mype = route%mype
       allocate (AllActive    (1:N_CatG, 1: nDEs))
       allocate (DstCatchID(1:N_CatG, 1: nDEs)) 
       allocate (srcProcsID   (1:N_CatG ))
       allocate (LocDstCatchID(1:N_CatG ))

       AllActive       = -9999
       srcProcsID      = -9999
       DstCatchID      = -9999
       LocDstCatchID   = NINT(DNSTR)

       call InitializeRiverRouting(MYPE, nDEs, MAPL_am_I_root(vm),pfaf_code, & 
            AllActive, DstCatchID, srcProcsID, LocDstCatchID, rc=STATUS)

       VERIFY_(STATUS)

       N_Active = count (srcProcsID == MYPE)

       allocate (GlbActive(1 : N_Active))
       allocate (tmp_index(1 : N_CatG  ))

       forall (N=1:N_CatG) tmp_index(N) = N

       GlbActive = pack (tmp_index, mask = (srcProcsID == MYPE))

       ! Initialize the cycle counter and sum (runoff) 

       allocate (runoff_save (1:NTILES))

       runoff_save = 0.
       ThisCycle   = 1

       FirstTime = .false.

       deallocate (tmp_index)
       
    ENDIF FIRST_TIME

    ! For efficiency, the time step to call the river routing model is set at ROUTE_DT 

    N_CYC = ROUTE_DT/HEARTBEAT

    RUN_MODEL : if (ThisCycle == N_CYC) then  

       runoff_save = runoff_save + runoff/real (N_CYC)

       ! Here we aggreagate GEOS_Catch/GEOS_CatchCN produced RUNOFF from TILES to CATCHMENTS
       ! Everything is local to the parallel block. Units: RUNOFF [kg m-2 s-1], 
       !        RUNOFF_CATCH [m3 s-1]
       ! -----------------------------------------------------------------------------------
       
       ! Unit conversion 
       
       mm2m3 = MAPL_RADIUS * MAPL_RADIUS / 1000.
       
       ALLOCATE (RUNOFF_CATCH(1:N_CatG))
    
       RUNOFF_CATCH = 0.
       
       DO N = 1, NTILES 
          RUNOFF_CATCH (pfaf_code(n)) = RUNOFF_CATCH (pfaf_code(n)) + mm2m3 * RUNOFF_SAVE (N) * TILE_AREA (N)
       END DO
       
       ! Inter-processor communication 1
       ! For catchment-tiles that contribute to the main catchment in some other processor, 
       ! send runoff to the corresponding srcProcsID(N)    
       ! -----------------------------------------------------------------------------------
       
       do N = Local_Min, Local_Max
          
          if ((AllActive (N,MYPE+1) > 0).and.(srcProcsID(N) /= MYPE)) then
             
             rbuff = RUNOFF_CATCH (N)
             
             call MPI_ISend(rbuff,1,MPI_real,srcProcsID(N),999,MPI_COMM_WORLD,req,status)
             call MPI_WAIT (req  ,MPI_STATUS_IGNORE,status)
             
             RUNOFF_CATCH (N) = 0.
             
          else
             
             if(srcProcsID(N) == MYPE) then 
                
                do i = 1,nDEs                
                   if((i-1 /= MYPE).and.(AllActive (N,i) > 0))  then                   
                      
                      call MPI_RECV(rbuff,1,MPI_real,i-1,999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)
                      RUNOFF_CATCH (N) = RUNOFF_CATCH (N) + rbuff
                      
                   endif
                end do
             endif
          endif
       end do
       
       ! Now compress and create subsets of arrays that only contain active catchments 
       !    in the local processor
       ! -----------------------------------------------------------------------------
       
       if(allocated (LENGSC_ACT ) .eqv. .false.) allocate (LENGSC_ACT  (1:N_Active))
       if(allocated (AREACAT_ACT ) .eqv. .false.) allocate (AREACAT_ACT (1:N_Active))
       if(allocated (WSTREAM_ACT ) .eqv. .false.) allocate (WSTREAM_ACT (1:N_Active))
       if(allocated (WRIVER_ACT  ) .eqv. .false.) allocate (WRIVER_ACT  (1:N_Active))
       if(allocated (QSFLOW_ACT  ) .eqv. .false.) allocate (QSFLOW_ACT  (1:N_Active))
       if(allocated (QOUTFLOW_ACT) .eqv. .false.) allocate (QOUTFLOW_ACT(1:N_Active))  
       if(allocated (RUNOFF_ACT  ) .eqv. .false.) allocate (RUNOFF_ACT  (1:N_Active))  
       
       DO N = 1, size (GlbActive)
          
          I = GlbActive (N)
          RUNOFF_ACT  (N) = RUNOFF_CATCH (I)
          
          I = GlbActive (N) - Local_Min + 1
          WSTREAM_ACT (N) = WSTREAM (I)
          WRIVER_ACT  (N) = WRIVER  (I)
          LENGSC_ACT  (N) = LENGSC  (I)
          AREACAT_ACT (N) = AREACAT (I)
          
       END DO
       
       QSFLOW_ACT   = 0.
       QOUTFLOW_ACT = 0.
       QSFLOW       = 0.
       QOUTFLOW     = 0.
       QINFLOW      = 0.
       
       ! Call river_routing_model
       ! ------------------------
       
       CALL RIVER_ROUTING  (N_Active, RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,  &
            WSTREAM_ACT,WRIVER_ACT, QSFLOW_ACT,QOUTFLOW_ACT) 
       
       DO N = 1, size (GlbActive)
          
          I = GlbActive (N) - Local_Min + 1
          
          WSTREAM (I) = WSTREAM_ACT (N) 
          WRIVER  (I) = WRIVER_ACT  (N)  
          QSFLOW  (I) = QSFLOW_ACT  (N)
          QOUTFLOW(I) = QOUTFLOW_ACT(N)

          if (LocDstCatchID (GlbActive (N)) ==  GlbActive (N)) then

             ! This catchment drains to the ocean, lake or a sink 
             ! if(ORIVERMOUTH(... ) > 0) send QOUTFLOW(I) [m3/s] to ORIVERMOUTH(N) th ocean tile
             ! if(LRIVERMOUTH(... ) > 0) send QOUTFLOW(I) [m3/s] to LRIVERMOUTH(N) th lake tile

          endif
       END DO
       
       ! Inter-processor communication-2
       ! Update down stream catchments
       ! -------------------------------
       
       do N = 1,N_CatG 
          
          if ((srcProcsID (N) == MYPE).and.(srcProcsID (LocDstCatchID (N)) == MYPE)) then ! destination is local
             
             I = LocDstCatchID (N) - Local_Min + 1 ! Downstream index in the local processor
             K = N - Local_Min + 1                 ! Source index in the local processor  
             
             if(LocDstCatchID (N) /= N) then ! ensure not to refill the reservoir by itself
                
                QINFLOW(I) = QINFLOW(I) + QOUTFLOW (K)
                WRIVER (I) = WRIVER (I) + QOUTFLOW (K) * real(route_dt)

             endif
                             
          elseif ((srcProcsID (N) == MYPE).and.(srcProcsID (LocDstCatchID (N)) /= MYPE)) then 
             
             if(srcProcsID (LocDstCatchID (N)) >= 0) then
                
                ! Send to downstream processor
                
                K = N - Local_Min + 1                 ! Source index in the local processor  
                
                call MPI_ISend(QOUTFLOW(K),1,MPI_real,srcProcsID (LocDstCatchID (N)),999,MPI_COMM_WORLD,req,status)
                call MPI_WAIT(req,MPI_STATUS_IGNORE,status) 
                
             endif
             
          elseif ((srcProcsID (N) /= MYPE).and.(srcProcsID (N) >= 0)) then
             
             K = srcProcsID (dstCatchID(N,srcProcsID (N)+1))
             
             if (k == MYPE) then
                
                do i = 1,nDEs
                   
                   if(MYPE /= i-1) then 
                      
                      if((srcProcsID  (n) == i-1).and.(srcProcsID (dstCatchID(N, i)) == MYPE))then                       
                         call MPI_RECV(rbuff,1,MPI_real, srcProcsID (N),999,MPI_COMM_WORLD,MPI_STATUS_IGNORE,status)
                         K = dstCatchID(N,i) - Local_Min + 1
                         QINFLOW (K) = QINFLOW (K) + rbuff
                         WRIVER  (K) = WRIVER (K)  + rbuff * real(route_dt)
                         
                      endif
                   endif
                end do
             endif
             
          endif
          
       end do

       ! initialize the cycle counter and sum (runoff_tile) 

       runoff_save = 0.
       ThisCycle   = 1         


    else
       
       runoff_save = runoff_save + runoff/real (N_CYC)
       
       ThisCycle = ThisCycle + 1
       
    endif RUN_MODEL

 call MAPL_TimerOff ( MAPL, "-RRM" )

! All done
! --------

    call MAPL_TimerOff(MAPL,"RUN")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

! ---------------------------------------------------------------------------

  subroutine InitializeRiverRouting(MYPE, numprocs, master_proc,           &
       pfaf_code, AllActive, AlldstCatchID, srcProcsID, LocDstCatchID, rc)
    
    implicit none
    INTEGER, INTENT (IN)                             :: MYPE, numprocs
    LOGICAL, INTENT (IN)                             :: master_proc
    INTEGER, DIMENSION (:),  INTENT (IN)             :: pfaf_code
    INTEGER, DIMENSION (N_CatG),          INTENT (INOUT) :: srcProcsID, LocDstCatchID
    INTEGER, DIMENSION (N_CatG,numprocs), INTENT (INOUT) :: Allactive,  AlldstCatchID

    INTEGER, DIMENSION(:)  ,ALLOCATABLE  :: global_buff, scounts, rdispls, rcounts, LocalActive
    INTEGER                              :: N_active, I,J,K,N,i1,i2,NProcs, Local_Min, Local_Max

    integer, optional,    intent(OUT):: rc

    integer :: mpierr
    character(len=ESMF_MAXSTR), parameter :: Iam='InitializeRiverRouting'

    ! STEP 1: Identify active catchments within the local processor. If the catchment is active in 
    !         more than 1 processor, choose an owner.
    ! --------------------------------------------------------------------------------------------

    allocate (LocalActive (1:N_CatG))
    LocalActive = -9999
 
    Local_Min = minval (pfaf_code)
    Local_Max = maxval (pfaf_code)
 
    do N = 1, size (pfaf_code)               
       LocalActive(pfaf_code(n)) = pfaf_code(n) 
    end do

    allocate (global_buff (N_CatG * numprocs))
    allocate (scounts(numprocs),rdispls(numprocs),rcounts(numprocs))  

    scounts = N_CatG
    rcounts = N_CatG
    
    rdispls(1) = 0
    global_buff= 0
    
    do i=2,numprocs
       rdispls(i)=rdispls(i-1)+rcounts(i-1)
    enddo
    
    call MPI_allgatherv  (                          &
         LocalActive, scounts         ,MPI_INTEGER, &
         global_buff, rcounts, rdispls,MPI_INTEGER, &
         MPI_COMM_WORLD, mpierr)
    
    do i=1,numprocs
       Allactive (:,i) = global_buff((i-1)*N_CatG+1:i*N_CatG)
    enddo

    if (master_proc) then

       DO N = 1, N_CatG
          NPROCS = count(Allactive(N,:) >= 1)
          if(NPROCS > 0)then
             if (NPROCS == 1) then
                srcProcsID (N) = maxloc(Allactive(N,:),dim=1) - 1
             else
                i1 = MAX(N - 5,1)
                i2 = MIN(N + 5, N_CatG)
                N_active = 0
                do I = 1,numprocs
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

    call MPI_BCAST (srcProcsID, N_CatG, MPI_INTEGER, 0,MPI_COMM_WORLD,mpierr)

    ! STEP 2: reset downstream catchment indeces (from -1 OR 1:291284) of catchments that are
    !            in the local processor to full domain indeces.
    ! ------------------------------------------------------------------------------------------

    do N = Local_Min, Local_Max
       
       if(LocalActive (N) >=1) then 
          
          if (LocDstCatchID (N) == -1) then
             ! (a) DNST Catch is a sink hole, ocean or lake so water drains to self 
             LocDstCatchID (N)    = N 
             
          endif
          
       else
          
          LocDstCatchID (N) = -9999 ! is inactive
          
       endif
    end do

    global_buff= 0
    
    call MPI_allgatherv  (                          &
         LocDstCatchID,  scounts      ,MPI_INTEGER, &
         global_buff, rcounts, rdispls,MPI_INTEGER, &
         MPI_COMM_WORLD, mpierr)
    
    do i=1,numprocs
       AlldstCatchID (:,i) = global_buff((i-1)*N_CatG+1:i*N_CatG)
    enddo
    
    deallocate (global_buff, scounts, rdispls, rcounts, LocalActive)

    RETURN_(ESMF_SUCCESS)
  end subroutine InitializeRiverRouting
end module GEOS_RouteGridCompMod
