#include "MAPL_Generic.h"

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
!   EXPORTS   : SFLOW, INFLOW, DISCHARGE \\

! !USES: 

  use ESMF
  use MAPL_Mod
  use MAPL_ConstantsMod
  use ROUTING_MODEL,          ONLY:     &
       river_routing, RRM_TIMESTEP
  USE catch_constants, ONLY: N_CatG => N_Pfaf_Catchs
  USE MAPL_MinMaxMod
  implicit none

  private

  type T_RROUTE_STATE
     private
     integer :: nTiles
     integer :: comm
     integer :: nDes
     integer :: myPe
     integer :: minCatch
     integer :: maxCatch
     integer, pointer :: pfaf(:) => NULL()
     real,    pointer :: tile_area(:) => NULL()
     type(MAPL_LocStreamXForm) :: xform
     type(MAPL_LocStream) :: route_ls
     logical :: ease_tiles
     type(MAPL_LocStream) :: application_ls
     type(ESMF_Grid) :: application_grid
     character(len=:), allocatable :: routing_tile_file
  end type T_RROUTE_STATE

  ! Wrapper for extracting internal state
  ! -------------------------------------
  type RROUTE_WRAP
     type (T_RROUTE_STATE), pointer :: PTR => null()
  end type RROUTE_WRAP

  character(len=ESMF_MAXSTR)        :: TILFILE
  integer                           :: RUN_DT, RRM_DT
  REAL                              :: DT
  
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

    type (ESMF_Config          )           :: CF 
    type (T_RROUTE_STATE), pointer         :: route_internal_state => null()
    type (RROUTE_wrap)                     :: wrap
    
!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                                 ,&
                          NAME=COMP_NAME	             ,&
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

! -----------------------------------------------------------
! Get the intervals
! -----------------------------------------------------------

    call ESMF_ConfigGetAttribute (CF, RRM_DT, Label="RRM_DT:", &
         default=RRM_TIMESTEP, RC=STATUS)
    VERIFY_(STATUS)
    
! -----------------------------------------------------------
! At the moment, this will refresh when the land parent 
! needs to refresh.

    call ESMF_ConfigGetAttribute ( CF, DT, Label=trim(COMP_NAME)//"_DT:", &
         default=DT, RC=STATUS)
    VERIFY_(STATUS)
    
! -----------------------------------------------------------
! get the TILFILE    
    
    call ESMF_ConfigGetAttribute ( CF, TILFILE, Label="TILING_FILE:", &
         default='tile.bin', RC=STATUS)
    VERIFY_(STATUS)    
    
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
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'length_of_channel_segment',&
         UNITS              = 'km+2'                     ,&
         SHORT_NAME         = 'LENGSC'                   ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'index_of_downtream_catchment',&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'DNSTR'                    ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  ) 

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_local_stream',&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WSTREAM'                  ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
         RC=STATUS  )
    
    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'volume_of_water_in_river' ,&
         UNITS              = 'm+3'                      ,&
         SHORT_NAME         = 'WRIVER'                   ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'TileID_of_the_lake_tile_at_the_river_mouth' ,&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'LRIVERMOUTH'              ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   

    call MAPL_AddInternalSpec(GC                     ,&
         LONG_NAME          = 'TileID_of_the_ocean_tile_at_the_river_mouth' ,&
         UNITS              = '1'                        ,&
         SHORT_NAME         = 'ORIVERMOUTH'              ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
         RESTART            = MAPL_RestartRequired       ,&
                                         RC=STATUS  )   
! -----------------------------------------------------------
!  EXPORT STATE:
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC,                        &
         LONG_NAME          = 'RUNOFF_IMPORT_from_tiles_to_catchments' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'RUNOFFIN'                 ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
                                         RC=STATUS  ) 
    
    call MAPL_AddExportSpec(GC,                        &
         LONG_NAME          = 'water_transfer_rate_from_streams_to_local_river_channel' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'SFLOW'                    ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
         VLOCATION          = MAPL_VLocationNone         ,&
                                         RC=STATUS  ) 

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'inflow_of_river_water_from_upstream_catchments' ,&
         UNITS              = 'm+3 s-1'                   ,&
         SHORT_NAME         = 'INFLOW'                    ,&
         DIMS               = MAPL_DimsHorzOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
                                          RC=STATUS  ) 

    call MAPL_AddExportSpec(GC,                    &
         LONG_NAME          = 'discharge_to_the_downstream_catchment' ,&
         UNITS              = 'm+3 s-1'                  ,&
         SHORT_NAME         = 'DISCHARGE'                ,&
         DIMS               = MAPL_DimsHorzOnly          ,&
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
    
!! Save pointer to the wrapped internal state in the GC
!! ----------------------------------------------------

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

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------

    type (ESMF_VM) :: VM
    integer        :: comm
    integer        :: nDEs
    integer        :: myPE
    integer        :: beforeMe, minCatch, maxCatch,  i

    type(ESMF_Grid)     :: catchGrid
    type(ESMF_DistGrid) :: distGrid
    type(ESMF_DELayout) :: layout_ease,layout_catch
    type (MAPL_MetaComp), pointer  :: MAPL
    type(MAPL_LocStream)           :: locstream,route_ls,application_ls
    type (T_RROUTE_STATE), pointer :: route => null()
    type (RROUTE_wrap)             :: wrap
    real(ESMF_KIND_R8), pointer    :: coords(:,:)
    integer, allocatable :: ims(:)
    integer              :: jms(1), ii,counts(3)
    character(len=ESMF_MAXSTR) :: master_gridname
    type(ESMF_Config) :: cf
    logical :: ease_grid_found,alternate_tile_found
     
    ! ------------------
    ! begin
    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)
    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
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

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

    route%comm = comm
    route%ndes = ndes
    route%mype = mype
    
    ! define minCatch, maxCatch
    allocate(ims(ndes))
    call MAPL_DecomposeDim ( n_catg,ims,ndes ) ! ims(mype+1) gives the size of my partition
    jms=1
    ! myPE is 0-based!
    beforeMe = sum(ims(1:mype))
    minCatch = beforeMe + 1
    maxCatch = beforeMe + ims(myPe+1)
    
    ! get LocStream
    call MAPL_Get(MAPL, LocStream = locstream, RC=status)
    VERIFY_(STATUS)
    block
     type(ESMF_Grid) :: mygrid
     integer :: mydims(3)
     call ESMF_GridCompGet(GC,grid=mygrid)
     call MAPL_GridGet(mygrid,globalcellcountperdim=mydims)
     if (mapl_am_I_root()) write(*,*)'bmaa mydims: ',mydims(1),mydims(2)
    end block

    !    CatchGrid = ESMF_GridCreate(name="ROUTE_GRID",countsPerDEDim1=ims,countsPerDeDim2=jms,indexFlag=ESMF_INDEX_DELOCAL,rc=status)
    catchGrid = grid_manager%make_grid(LatLonGridFactory(grid_name="ROUTE_GRID",im_world=n_catg,jm_world=1,lm=1, &
         nx=ndes,ny=1,pole='XY',dateline='XY',lon_range=RealMinMax(-1.0,1.0), &
         lat_range=RealMinMax(-1.0,10.),ims=ims,jms=jms,force_decomposition=.true.,rc=status))
    _VERIFY(status)
    call ESMF_GridAddCoord(CatchGrid,rc=status)
    _VERIFY(status)
    call ESMF_GridGetCoord(CatchGrid,coordDim=1,farrayptr=coords,rc=status)
    _VERIFY(status)
    coords=0.0d0
    call ESMF_GridGetCoord(CatchGrid,coordDim=2,farrayptr=coords,rc=status)
    _VERIFY(status)
    coords=0.0d0
   
    call ESMF_ConfigFindLabel ( CF, "GEOSldas.GRIDNAME:", isPresent=ease_grid_found,rc=status)
    if (ease_grid_found) then
       call ESMF_ConfigGetAttribute ( CF, master_gridname, Label="GEOSldas.GRIDNAME:", RC=STATUS)
       _VERIFY(status)
    else
       master_gridname =""
    end if
    route%ease_tiles = (index(master_gridname,"EASE")/=0)

    alternate_tile_found=.false.
    if (route%ease_tiles) then
       call ESMF_ConfigGetAttribute ( CF, route%routing_tile_file, Label="EASE_PFAFFSETTER_FILE:", RC=STATUS)
       _VERIFY(status)
       alternate_tile_found=.true.
    else 
       call ESMF_ConfigFindLabel( CF, "routing_tile:", isPresent=alternate_tile_found,rc=status)
       if (alternate_tile_found) then
          call ESMF_ConfigGetAttribute(CF, route%routing_tile_file, label="routing_tile:", rc=status)
          _VERIFY(status)
       end if
    end if

    _ASSERT(.not.(alternate_tile_found .and. ease_grid_found),"can not use alternative tile file for route when running ease grid")

    if (alternate_tile_found) then

       if (route%ease_tiles) then
          route%application_grid = grid_manager%make_grid(cf,prefix="GEOSldas"//".",rc=status)
          _VERIFY(status)
       else
          call ESMF_GridCompGet(GC,grid=route%application_grid,rc=status)
          _VERIFY(status)
       end if

       call ESMF_GridGet(route%application_grid,distgrid=distgrid,rc=status)
       _VERIFY(status)
       call ESMF_DistGridGet(distgrid,deLayout=layout_ease,rc=status)
       _VERIFY(status)
       call ESMF_GridGet(CatchGrid,distgrid=distgrid,rc=status)
       _VERIFY(status)
       call ESMF_DistGridGet(distgrid,deLayout=layout_catch,rc=status)
       _VERIFY(status)

       call MAPL_LocStreamCreate(application_ls,layout_ease,trim(route%routing_tile_file),'application_ls',mask=[MAPL_LAND], &
            grid=route%application_grid,rc=status)
       _VERIFY(status)
       route%application_ls=application_ls
       call MAPL_LocstreamCreate(route_ls,layout_catch,trim(route%routing_tile_file),'route_ls', mask=[MAPL_LAND], &
            grid=CatchGrid,use_pfaf=.true.,rc=status)
       _VERIFY(status)
       route%route_ls=route_ls
       call MAPL_LocStreamCreateXform(xform=route%xform,locstreamOut=route_ls, &
            locstreamIn=application_ls,name="application_to_route",rc=status)
       _VERIFY(status)
        
    else
       call ESMF_GridGet(CatchGrid,distgrid=distgrid,rc=status)
       _VERIFY(status)
       call ESMF_DistGridGet(distgrid,deLayout=layout_catch,rc=status)
       _VERIFY(status)
       call MAPL_LocstreamCreate(route_ls,layout_catch,trim(TILFILE),'catch_ls', mask=[MAPL_LAND], &
            grid=CatchGrid,use_pfaf=.true.,rc=status)
       _VERIFY(status)
       route%route_ls=route_ls

       call MAPL_LocStreamCreateXForm(xform=route%xform, locstreamOut=route_ls, &
             locStreamIn=locstream,name="tile_to_route",rc=status)
    end if


    call MAPL_LocStreamGet(route%route_ls,nt_local=route%ntiles,rc=status)
    _VERIFY(status)

    call MAPL_GridGet(catchGrid,localCellCountPerDim=counts,rc=status)
    allocate (route%pfaf(counts(1)))
    ii=1

    do i=MinCatch, MaxCatch
       route%pfaf(ii)=i
       ii=ii+1
    enddo
    
    route%minCatch = minCatch
    route%maxCatch = maxCatch

    if (allocated(ims)) deallocate(ims)
    call ESMF_GridCompSet(gc,grid=CatchGrid,rc=status)
    _VERIFY(status)
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

    character(len=ESMF_MAXSTR)         :: IAm="Run"
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! -----------------------------------------------------------
! Locals
! -----------------------------------------------------------

    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State   )              :: INTERNAL
    type (ESMF_Config  )              :: CF
    type(ESMF_VM)                     :: VM

! -----------------------------------------------------
! IMPORT pointers
! ----------------------------------------------------- 

    real, dimension(:), pointer   :: RUNOFF 

! -----------------------------------------------------
! INTERNAL pointers
! ----------------------------------------------------- 

    real, dimension(:,:), pointer :: AREACAT_2D
    real, dimension(:,:), pointer :: LENGSC_2D
    real, dimension(:,:), pointer :: DNSTR_2D
    real, dimension(:,:), pointer :: WSTREAM_2D
    real, dimension(:,:), pointer :: WRIVER_2D
    real, dimension(:,:), pointer :: LRIVERMOUTH_2D
    real, dimension(:,:), pointer :: ORIVERMOUTH_2D

    real, dimension(:),   pointer :: AREACAT
    real, dimension(:),   pointer :: LENGSC
    real, dimension(:),   pointer :: DNSTR
    real, dimension(:),   pointer :: WSTREAM
    real, dimension(:),   pointer :: WRIVER
    real, dimension(:),   pointer :: LRIVERMOUTH
    real, dimension(:),   pointer :: ORIVERMOUTH

! -----------------------------------------------------
! EXPORT pointers 
! -----------------------------------------------------
    
    real, dimension(:,:), pointer :: RUNOFFIN_2D
    real, dimension(:,:), pointer :: SFLOW_2D
    real, dimension(:,:), pointer :: INFLOW_2D
    real, dimension(:,:), pointer :: DISCHARGE_2D
    real, dimension(:),   pointer :: RUNOFFIN
    real, dimension(:),   pointer :: SFLOW
    real, dimension(:),   pointer :: INFLOW
    real, dimension(:),   pointer :: DISCHARGE
  
! Others

    type(ESMF_Grid)                          :: ESMFGRID
    type (MAPL_LocStream)                    :: LOCSTREAM
    integer                                  :: N_CatL, N_CYC
    logical, save                            :: FirstTime=.true.
    INTEGER, DIMENSION(:,:), POINTER,   SAVE :: AllActive,DstCatchID 
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: srcProcsID, LocDstCatchID  
    INTEGER,                            SAVE :: ThisCycle   
    INTEGER,                            SAVE :: Local_Min, Local_Max, Pfaf_Min, Pfaf_Max
    REAL, ALLOCATABLE, DIMENSION(:),    SAVE :: runoff_save   
    integer                                  :: K, N, I, req
    REAL                                     :: rbuff 
    real, allocatable, DIMENSION(:)          :: runoff_catch_dist
    real, allocatable, DIMENSION(:,:)        :: runoff_on_catch
    integer                                  :: ndes, mype, nt_local, counts(3)
    type (T_RROUTE_STATE), pointer           :: route => null()
    type (RROUTE_wrap)                       :: wrap
    real, allocatable :: ease_array(:,:), tile_in(:), tile_out(:)
    integer :: ntiles_in,ntiles_out
    integer, pointer :: local_i(:), local_j(:)

    ! ------------------
    ! begin

    if (mapl_am_i_root()) write(*,*)'bmaa running route'
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)    
    route => wrap%ptr
    call ESMF_VMGetCurrent(VM, RC=STATUS) ; VERIFY_(STATUS)

!! Get the target components name and set-up traceback handle.
!! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
  
    Iam = trim(COMP_NAME) // "RUN"

!! Get my internal MAPL_Generic state
!! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)
    call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    RUN_DT = NINT (DT)
    
! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"RUN")
    
! Get parameters from generic state
! ---------------------------------

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS) 

! Get runoff from import
    call MAPL_GetPointer(import,runoff,"RUNOFF",rc=status)
    _VERIFY(status)
    
!! get pointers to internal variables (note they are 2D in the internal state)
!! ----------------------------------

    call MAPL_GetPointer(INTERNAL, AREACAT_2D , 'AREACAT', RC=STATUS)
    VERIFY_(STATUS)
    AREACAT => AREACAT_2D(:,1)
    call MAPL_GetPointer(INTERNAL, LENGSC_2D  , 'LENGSC',  RC=STATUS)
    VERIFY_(STATUS)
    LENGSC => LENGSC_2D(:,1)
    call MAPL_GetPointer(INTERNAL, DNSTR_2D   , 'DNSTR'  , RC=STATUS)
    VERIFY_(STATUS)
    DNSTR => DNSTR_2D(:,1)
    call MAPL_GetPointer(INTERNAL, WSTREAM_2D , 'WSTREAM', RC=STATUS)
    VERIFY_(STATUS)
    WSTREAM => WSTREAM_2D(:,1)
    call MAPL_GetPointer(INTERNAL, WRIVER_2D  , 'WRIVER' , RC=STATUS)
    VERIFY_(STATUS)
    WRIVER => WRIVER_2D(:,1)
    call MAPL_GetPointer(INTERNAL, LRIVERMOUTH_2D, 'LRIVERMOUTH' , RC=STATUS)
    VERIFY_(STATUS)
    LRIVERMOUTH => LRIVERMOUTH_2D(:,1)
    call MAPL_GetPointer(INTERNAL, ORIVERMOUTH_2D, 'ORIVERMOUTH' , RC=STATUS)
    VERIFY_(STATUS)
    ORIVERMOUTH => ORIVERMOUTH_2D(:,1)

!! get pointers to EXPORTS
!! -----------------------
    call MAPL_GetPointer(EXPORT, RUNOFFIN_2D, 'RUNOFFIN'  , alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    RUNOFFIN => RUNOFFIN_2D(:,1)
    call MAPL_GetPointer(EXPORT, SFLOW_2D, 'SFLOW'  , alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    SFLOW => SFLOW_2D(:,1)
    call MAPL_GetPointer(EXPORT, INFLOW_2D, 'INFLOW' , alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    INFLOW => INFLOW_2D(:,1)
    call MAPL_GetPointer(EXPORT, DISCHARGE_2D, 'DISCHARGE', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    DISCHARGE => DISCHARGE_2D(:,1)

! Begin water transfer across the river network
! ---------------------------------------------
    
    call MAPL_TimerOn  ( MAPL, "-RRM" )
    N_CatL  = size(AREACAT)
    ndes = route%ndes
    mype = route%mype

    FIRST_TIME : IF (FirstTime) THEN

       Local_min  = MINVAL (route%pfaf)
       Local_max  = MAXVAL (route%pfaf)

       call MPI_Reduce(Local_Min,Pfaf_Min,1,MPI_INTEGER,MPI_MIN,0,route%comm,STATUS) ; VERIFY_(STATUS)
       call MPI_Reduce(Local_Max,Pfaf_Max,1,MPI_INTEGER,MPI_MAX,0,route%comm,STATUS) ; VERIFY_(STATUS)    
       call MPI_BARRIER( route%comm, STATUS ) ; VERIFY_(STATUS)    
       call MPI_BCAST (Pfaf_Min , 1, MPI_INTEGER, 0,route%comm,STATUS) ; VERIFY_(STATUS)
       call MPI_BCAST (Pfaf_Max , 1, MPI_INTEGER, 0,route%comm,STATUS) ; VERIFY_(STATUS)
       
       !! Pfafstetter catchment Domain Decomposition :         
       !! --------------------------------------------

       !! AllActive      : Processor(s) where the catchment is active  (identical in any processor). 
       !! srcProcsID     : For all active catchments anywhere tells which processor is the principal owner of the catchment (identical in any processor).
       !! DstCatchID  : 2-D array contains downstream catchID and downstream processor (identical in any processor)
       !! LocDstCatchID  : Downstream catchID when for catchments that are local to the processor.

       allocate (AllActive    (1:N_CatG, 1: nDEs))
       allocate (DstCatchID   (1:N_CatG, 1: nDEs)) 
       allocate (srcProcsID   (1:N_CatG ))
       allocate (LocDstCatchID(1:N_CatG ))

       AllActive       = -9999
       srcProcsID      = -9999
       DstCatchID      = -9999
       LocDstCatchID   = -9999
       LocDstCatchID (Local_Min - Pfaf_Min + 1 :Local_Max - Pfaf_Min + 1) = NINT(DNSTR)

       call InitializeRiverRouting(MYPE, nDEs, route%comm, MAPL_am_I_root(vm),route%pfaf, & 
            Pfaf_Min, Pfaf_Max, AllActive, DstCatchID, srcProcsID, LocDstCatchID, rc=STATUS)

       VERIFY_(STATUS)
       ASSERT_ ((count (srcProcsID == MYPE) - N_CatL) == 0)

       ! Initialize the cycle counter and sum (runoff) 

       allocate (runoff_save (1:size (runoff)))

       runoff_save = 0.
       ThisCycle   = 0

       FirstTime = .false.
       
    ENDIF FIRST_TIME

    N_CYC = RRM_DT/RUN_DT

    ! save mean tile-space runoff during RRM_DT in [kg/m2/s] 
    runoff_save (:) = runoff_save (:) + runoff (:)/real (N_CYC)
    ThisCycle = ThisCycle + 1
    
    RUN_MODEL : if (ThisCycle == N_CYC) then

       call ESMF_GridCompGet(gc,grid=esmfgrid,rc=status)
       _VERIFY(status)
       call MAPL_GridGet(esmfgrid,localCellCountPerDim=counts,rc=status)
       _VERIFY(status)
       allocate(runoff_on_catch(counts(1),counts(2)))    
       if (allocated(route%routing_tile_file)) then
          call MAPL_GridGet(route%application_grid,localCellCountPerDim=counts,rc=status)
          _VERIFY(status)
          allocate(ease_array(counts(1),counts(2)))
          ease_array=0.
          call MAPL_Get(MAPL, LocStream = locstream, RC=status)
          _VERIFY(status)
          allocate(ease_array(counts(1),counts(2)))
          if (route%ease_tiles) then
             call MAPL_LocStreamGet(locstream,local_i=local_i,local_j=local_j,rc=status)
             _VERIFY(status)
             do i=1,size(runoff_save)
                ease_array(local_i(i),local_j(i))=runoff_save(i)
             enddo
          else
             call MAPL_LocStreamTransform(locstream,ease_array,runoff_save,rc=status)
             _VERIFY(status)
          end if

          ! now that we are on cube, do cube to tile, tile-tile, tile-catchment grid
          call mapl_locstreamget(route%application_ls,nt_local=ntiles_in,rc=status)
          _VERIFY(status)
          call mapl_locstreamget(route%route_ls,nt_local=ntiles_out,rc=status)
          _VERIFY(status)
          allocate(tile_in(ntiles_in))
          allocate(tile_out(ntiles_out))
          
          call MAPL_LocStreamTransform(route%application_ls,tile_in,ease_array,rc=status)
          _VERIFY(status)
          call MAPL_LocStreamTransform(tile_in,route%xform,tile_out,rc=status)
          _VERIFY(status)
          call MAPL_LocStreamTransform(route%route_ls,runoff_on_catch,tile_out,rc=status)
          _VERIFY(status)
          deallocate(ease_array,tile_in,tile_out)
       else 
          ! Call Locstream transform, takes land tiles that are on the distribution used by
          ! catchgridcomp (based on atmosphere) and shuffles them to be on catchment distribution
          ! in other words the tile is on the processor contains the catchment making up the tile
          call MAPL_LocStreamGet(route%route_ls,nt_local=nt_local,rc=status)
          _VERIFY(status)
          allocate(runoff_catch_dist(nt_local))
          call MAPL_LocStreamTransform(runoff_catch_dist,route%xform,runoff_save,rc=status)
          _VERIFY(status)
          ! Finall we can take the tiles and do the conservative locstream transform to 
          ! aggregate them on the actual catchments. I am allocataing a temporary array
          ! get the result of this computation
          call MAPL_LOcStreamTransform(route%route_ls,runoff_on_catch,runoff_catch_dist,rc=status)
          _VERIFY(status)
       end if

       ! unit conversion of spatially averaged catchment runoff to [m3/s]
       runoff_on_catch(:,1) = runoff_on_catch (:,1) * AREACAT(:) * 1000.
       RUNOFFIN             = runoff_on_catch (:,1) 
       
       SFLOW       = 0.
       DISCHARGE   = 0.
       INFLOW      = 0.
       
       ! Call river_routing_model
       ! ------------------------
       
       CALL RIVER_ROUTING  (N_catL, real(RRM_DT), runoff_on_catch(:,1),AREACAT,LENGSC,  &
            WSTREAM, WRIVER, SFLOW, DISCHARGE) 

       ! Inter-processor communication: Update downstream catchments
       ! -----------------------------------------------------------
       
       do N = 1,N_CatG 
          
          if ((srcProcsID (N) == MYPE).and.(srcProcsID (LocDstCatchID (N)) == MYPE)) then ! destination is local
             
             I = LocDstCatchID (N) - Local_Min + 1 ! Downstream index in the local processor
             K = N - Local_Min + 1                 ! Source index in the local processor  
             
             if(LocDstCatchID (N) /= N) then ! ensure not to refill the reservoir by itself
                INFLOW(I)  = INFLOW(I) + DISCHARGE (K)
                WRIVER (I) = WRIVER (I) + DISCHARGE (K) * real(RRM_DT)                
             endif
             
          elseif ((srcProcsID (N) == MYPE).and.(srcProcsID (LocDstCatchID (N)) /= MYPE)) then 
             
             if(srcProcsID (LocDstCatchID (N)) >= 0) then
                
                ! Send to downstream processor                
                K = N - Local_Min + 1         ! Source index in the local processor                  
                call MPI_ISend(DISCHARGE(K),1,MPI_real,srcProcsID (LocDstCatchID (N)),999,route%comm,req,status)
                call MPI_WAIT(req,MPI_STATUS_IGNORE,status)                 
             endif
             
          elseif ((srcProcsID (N) /= MYPE).and.(srcProcsID (N) >= 0)) then
             
             K = srcProcsID (dstCatchID(N,srcProcsID (N)+1))
             
             if (k == MYPE) then                
                do i = 1,nDEs                   
                   if(MYPE /= i-1) then                       
                      if((srcProcsID  (n) == i-1).and.(srcProcsID (dstCatchID(N, i)) == MYPE))then                       
                         call MPI_RECV(rbuff,1,MPI_real, srcProcsID (N),999,route%comm,MPI_STATUS_IGNORE,status)
                         K = dstCatchID(N,i) - Local_Min + 1
                         INFLOW (K)  = INFLOW (K) + rbuff
                         WRIVER  (K) = WRIVER (K) + rbuff * real(RRM_DT)                    
                      endif
                   endif
                end do
             endif             
          endif          
       end do

       ! initialize the cycle counter and sum (runoff_tile) 

       runoff_save = 0.
       ThisCycle   = 0         
       deallocate (runoff_catch_dist, runoff_on_catch)
       
    end if RUN_MODEL

    call MAPL_TimerOff ( MAPL, "-RRM" )

    ! All done
    ! --------

    call MAPL_TimerOff(MAPL,"RUN")

    RETURN_(ESMF_SUCCESS)

  contains

    ! ---------------------------------------------------------------------------

    subroutine InitializeRiverRouting(MYPE, numprocs, comm, root_proc,           &
         pfaf_code, Pfaf_Min, Pfaf_Max, AllActive, AlldstCatchID, srcProcsID, LocDstCatchID, rc)
      
      implicit none
      INTEGER, INTENT (IN)                             :: MYPE, numprocs, comm, Pfaf_Min, Pfaf_Max
      LOGICAL, INTENT (IN)                             :: root_proc
      INTEGER, DIMENSION (:),  INTENT (IN)             :: pfaf_code
      INTEGER, DIMENSION (:),           INTENT (INOUT) :: srcProcsID, LocDstCatchID
      INTEGER, DIMENSION (:,:),         INTENT (INOUT) :: Allactive,  AlldstCatchID
      
      INTEGER, DIMENSION(:)  ,ALLOCATABLE  :: global_buff, scounts, rdispls, rcounts, LocalActive
      INTEGER                              :: N_active, I,J,K,N,i1,i2,NProcs, Local_Min, Local_Max      
      integer, optional,    intent(OUT):: rc      
      integer :: status
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
           COMM, status)
      
      do i=1,numprocs
         Allactive (:,i) = global_buff((i-1)*N_CatG+1:i*N_CatG)
      enddo
      
      if (root_proc) then
         
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
      
      call MPI_BCAST (srcProcsID, N_CatG, MPI_INTEGER, 0,COMM, STATUS)
      
      ! STEP 2: reset downstream catchment indeces (from -1 OR 1:291284) of catchments that are
      !            in the local processor to full domain indeces.
      ! ------------------------------------------------------------------------------------------
      
      do N = Local_Min, Local_Max
         
         if(LocalActive (N) >=1) then 
            
            if ((LocDstCatchID (N) == -1).OR.(LocDstCatchID (N) < Pfaf_Min) &
                 .OR. (LocDstCatchID (N) > Pfaf_Max)) then
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
           COMM, STATUS)
      
      do i=1,numprocs
         AlldstCatchID (:,i) = global_buff((i-1)*N_CatG+1:i*N_CatG)
      enddo
      
      deallocate (global_buff, scounts, rdispls, rcounts, LocalActive)
      
      RETURN_(ESMF_SUCCESS)

    end subroutine InitializeRiverRouting
        
  end subroutine RUN

  
end module GEOS_RouteGridCompMod
