
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

! !USES: 

  use ESMF
  use MAPL_Mod
  use MAPL_ConstantsMod
  use ROUTING_MODEL,          ONLY:     &
       river_routing_lin, river_routing_hyd, ROUTE_DT
#if 0
  USE catch_constants, ONLY:          &
       N_CatG => N_Pfaf_Catchs  
#endif
  use, intrinsic :: iso_c_binding
  
  implicit none
  integer, parameter :: N_CatG = 291284
  integer,parameter :: upmax=34
  character(len=500) :: inputdir="/discover/nobackup/yzeng3/data/river_input/"
  integer,save :: nmax 

  private

  type T_RROUTE_STATE
     private
     type (ESMF_RouteHandle) :: routeHandle
     type (ESMF_Field)       :: field
     integer :: nTiles
     integer :: nt_global
     integer :: nt_local
     integer :: comm
     integer :: nDes
     integer :: myPe
     integer :: minCatch
     integer :: maxCatch
     integer, pointer :: pfaf(:) => NULL()
     real,    pointer :: tile_area(:) => NULL() !m2
     integer, pointer :: nsub(:) => NULL()
     integer, pointer :: subi(:,:) => NULL()
     real,    pointer :: subarea(:,:) => NULL() !m2

     integer, pointer :: scounts_global(:) => NULL()
     integer, pointer :: rdispls_global(:) => NULL()
     integer, pointer :: scounts_cat(:) => NULL()
     integer, pointer :: rdispls_cat(:) => NULL()
 
     real,    pointer :: runoff_save(:) => NULL()
     real,    pointer :: areacat(:) => NULL() !m2
     real,    pointer :: lengsc(:) => NULL() !m

     real,    pointer :: wstream(:) => NULL() !m3
     real,    pointer :: wriver(:)  => NULL() !m3
     integer, pointer :: downid(:) => NULL()
     integer, pointer :: upid(:,:) => NULL()

     real,    pointer :: wriver_acc(:)  => NULL()
     real,    pointer :: wstream_acc(:) => NULL()     
     real,    pointer :: qoutflow_acc(:) => NULL()
     real,    pointer :: qsflow_acc(:)  => NULL()

     real,    pointer :: lstr(:) => NULL() !m
     real,    pointer :: qri_clmt(:) => NULL() !m3/s
     real,    pointer :: qin_clmt(:) => NULL() !m3/s
     real,    pointer :: qstr_clmt(:) =>NULL() !m3/s
     real,    pointer :: K(:) => NULL()
     real,    pointer :: Kstr(:) => NULL()

  end type T_RROUTE_STATE


  interface
    function mkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: mkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t), value :: mode
    end function mkdir
  end interface

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
!    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, RC=STATUS)
!    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN2, RC=STATUS )
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
    call MAPL_TimerAdd(GC, name="RUN1"           ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN2"           ,RC=STATUS)
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
    integer, pointer :: arbSeq_pf(:) => NULL()    
    integer, pointer :: arbSeq_ori(:) => NULL()    
    integer, allocatable :: arbIndex(:,:)
    real, pointer :: tile_area_src(:) => NULL()
    integer,pointer :: local_id(:)  => NULL()
    real, pointer :: tile_area_local(:) => NULL(), tile_area_global(:) => NULL()
    real, pointer :: tile_area(:) => NULL()    
    real, pointer :: ptr2(:) => NULL()

    real,pointer :: subarea_global(:,:)=> NULL(),subarea(:,:)=> NULL() ! Arrays for sub-area and fractions
    integer,pointer :: subi_global(:,:)=> NULL(),subi(:,:)=> NULL()
    integer,pointer :: nsub_global(:)=> NULL(),nsub(:)=> NULL()
    real,pointer :: area_cat_global(:)=> NULL(),area_cat(:)=> NULL()
    integer,pointer :: scounts(:)=>NULL()
    integer,pointer :: scounts_global(:)=>NULL(),rdispls_global(:)=>NULL()
    integer,pointer :: scounts_cat(:)=>NULL(),rdispls_cat(:)=>NULL()    

    real,pointer :: runoff_save(:)=>NULL(), areacat(:)=>NULL()
    real,pointer :: lengsc_global(:)=>NULL(), lengsc(:)=>NULL(), buff_global(:)=>NULL()
    integer,pointer :: downid_global(:)=>NULL(), downid(:)=>NULL()
    integer,pointer :: upid_global(:,:)=>NULL(), upid(:,:)=>NULL()    

    real,pointer :: wstream(:)=>NULL(),wriver(:)=>NULL()
    real,pointer :: wstream_global(:)=>NULL(),wriver_global(:)=>NULL()    
    
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap

    type(ESMF_Time) :: CurrentTime
    integer :: YY,MM,DD,HH,MMM,SS
    character(len=4) :: yr_s
    character(len=2) :: mon_s,day_s    
    character(len=3) :: resname

    real, pointer :: dataPtr(:)
    integer :: j,nt_local,mpierr,it   
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
 
    allocate(ims(1:ndes))
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
    call MAPL_LocStreamGet(locstream, &
         tileGrid=tilegrid, nt_global=nt_global, RC=status)
    VERIFY_(STATUS)     
    route%nt_global = nt_global

    if(nt_global==112573)then
      resname="M36"
      nmax=150
    else if(nt_global==1684725)then
      resname="M09"
      nmax=458
    else
      if(mapl_am_I_root())then
        print *,"unknown grid for routing model"
        stop
      endif
    endif
    ! exchange Pfaf across PEs

    call MAPL_LocStreamGet(locstream, TILEAREA = tile_area_src, LOCAL_ID=local_id, RC=status)
    VERIFY_(STATUS) 
    nt_local=size(tile_area_src,1) 
    route%nt_local=nt_local       

    ntiles = maxCatch-minCatch+1
    allocate(arbSeq_pf(maxCatch-minCatch+1))
    arbSeq_pf = [(i, i = minCatch, maxCatch)]
   
    ! redist pfaf (NOTE: me might need a second routehandle for integers)

    route%pfaf => arbSeq_pf
    route%ntiles = ntiles  
    route%minCatch = minCatch
    route%maxCatch = maxCatch 
  ! Read sub-area data from text files
    allocate(nsub_global(N_CatG),subarea_global(nmax,N_CatG))
    open(77,file=trim(inputdir)//"/Pfaf_nsub_"//trim(resname)//".txt",status="old",action="read"); read(77,*)nsub_global; close(77)
    open(77,file=trim(inputdir)//"/Pfaf_asub_"//trim(resname)//".txt",status="old",action="read"); read(77,*)subarea_global; close(77)       
    allocate(nsub(ntiles),subarea(nmax,ntiles))
    nsub=nsub_global(minCatch:maxCatch)
    subarea=subarea_global(:,minCatch:maxCatch)
    subarea=subarea*1.e6 !km2->m2
    deallocate(nsub_global,subarea_global)

    route%nsub => nsub
    route%subarea => subarea
    
    allocate(subi_global(nmax,N_CatG),subi(nmax,ntiles))
    open(77,file=trim(inputdir)//"/Pfaf_isub_"//trim(resname)//".txt",status="old",action="read");read(77,*)subi_global;close(77)
    subi=subi_global(:,minCatch:maxCatch)
    route%subi => subi
    deallocate(subi_global)


    allocate(scounts(ndes),scounts_global(ndes),rdispls_global(ndes))
    scounts=0
    scounts(mype+1)=nt_local  
    call MPI_Allgather(scounts(mype+1), 1, MPI_INTEGER, scounts_global, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr) 
    rdispls_global(1)=0
    do i=2,nDes
      rdispls_global(i)=rdispls_global(i-1)+scounts_global(i-1)
    enddo
    deallocate(scounts)
    route%scounts_global=>scounts_global
    route%rdispls_global=>rdispls_global

    allocate(scounts(ndes),scounts_cat(ndes),rdispls_cat(ndes))
    scounts=0
    scounts(mype+1)=ntiles  
    call MPI_Allgather(scounts(mype+1), 1, MPI_INTEGER, scounts_cat, 1, MPI_INTEGER, MPI_COMM_WORLD, mpierr) 
    rdispls_cat(1)=0
    do i=2,nDes
      rdispls_cat(i)=rdispls_cat(i-1)+scounts_cat(i-1)
    enddo
    deallocate(scounts)
    route%scounts_cat=>scounts_cat
    route%rdispls_cat=>rdispls_cat

    allocate(runoff_save(1:nt_local))
    route%runoff_save => runoff_save
    route%runoff_save=0.

    allocate(tile_area_local(nt_local),tile_area_global(nt_global))  
    open(77,file=trim(inputdir)//"/area_"//trim(resname)//"_1d.txt",status="old",action="read");read(77,*)tile_area_global;close(77)
    tile_area_local=tile_area_global(rdispls_global(mype+1)+1:rdispls_global(mype+1)+nt_local)*1.e6 !km2->m2
    route%tile_area => tile_area_local
    deallocate(tile_area_global)

    allocate(areacat(1:ntiles))
    areacat=0. 
    do i=1,ntiles
      do j=1,nmax
        it=route%subi(j,i) 
        if(it>0)then 
          areacat(i)=areacat(i)+route%subarea(j,i)
        endif
        if(it==0)exit
      enddo
    enddo  
    route%areacat=>areacat

    allocate(lengsc_global(n_catg),lengsc(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_lriv_PR.txt",status="old",action="read");read(77,*)lengsc_global;close(77)
    lengsc=lengsc_global(minCatch:maxCatch)*1.e3 !km->m
    route%lengsc=>lengsc
    deallocate(lengsc_global)

    allocate(downid_global(n_catg),downid(ntiles))
    open(77,file=trim(inputdir)//"/downstream_1D_new_noadj.txt",status="old",action="read");read(77,*)downid_global;close(77)    
    downid=downid_global(minCatch:maxCatch)
    route%downid=>downid
    deallocate(downid_global)

    allocate(upid_global(upmax,n_catg),upid(upmax,ntiles))   
    open(77,file=trim(inputdir)//"/upstream_1D.txt",status="old",action="read");read(77,*)upid_global;close(77)  
    upid=upid_global(:,minCatch:maxCatch)   
    route%upid=>upid
    deallocate(upid_global)

    call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
    call ESMF_TimeGet(CurrentTime, yy=YY, mm=MM, dd=DD, h=HH, m=MMM, s=SS, rc=status) 
    write(yr_s,'(I4.4)')YY
    write(mon_s,'(I2.2)')MM
    write(day_s,'(I2.2)')DD    
    if(mapl_am_I_root())print *, "init time is ", YY, "/", MM, "/", DD, " ", HH, ":", MMM, ":", SS    
    allocate(wriver(ntiles),wstream(ntiles))
    allocate(wriver_global(n_catg),wstream_global(n_catg))
    open(77,file="../input/restart/river_storage_rs_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",status="old",action="read",iostat=status)
    if(status==0)then
      read(77,*)wriver_global;close(77)
    else
      close(77)
      open(78,file=trim(inputdir)//"/river_storage_rs_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",status="old",action="read",iostat=status)
      if(status==0)then   
        read(78,*)wriver_global;close(78)  
      else
        close(78)      
        open(79,file=trim(inputdir)//"/river_storage_rs.txt",status="old",action="read",iostat=status)      
        if(status==0)then 
          read(79,*)wriver_global;close(79) 
        else
          close(79)    
          wriver_global=0.
        endif
      endif
    endif
    open(77,file="../input/restart/stream_storage_rs_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",status="old",action="read",iostat=status)
    if(status==0)then
      read(77,*)wstream_global;close(77)
    else
      close(77)
      open(78,file=trim(inputdir)//"/stream_storage_rs_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",status="old",action="read",iostat=status)
      if(status==0)then   
        read(78,*)wstream_global;close(78)  
      else
        close(78)      
        open(79,file=trim(inputdir)//"/stream_storage_rs.txt",status="old",action="read",iostat=status)      
        if(status==0)then 
          read(79,*)wstream_global;close(79) 
        else
          close(79)    
          wstream_global=0.
        endif
      endif
    endif
    if(mapl_am_I_root())print *, "init river storage is: ",sum(wriver_global)/1.e9
    if(mapl_am_I_root())print *, "init stream storage is: ",sum(wstream_global)/1.e9        
    wriver=wriver_global(minCatch:maxCatch)
    wstream=wstream_global(minCatch:maxCatch)
    deallocate(wriver_global,wstream_global)
    route%wstream=>wstream
    route%wriver=>wriver

    allocate(route%wriver_acc(ntiles),route%wstream_acc(ntiles),route%qoutflow_acc(ntiles),route%qsflow_acc(ntiles))
    route%wriver_acc=0.
    route%wstream_acc=0.
    route%qoutflow_acc=0.
    route%qsflow_acc=0.

   !input for geometry hydraulic
    allocate(buff_global(n_catg),route%lstr(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_lstr_PR.txt",status="old",action="read");read(77,*)buff_global;close(77)
    route%lstr=buff_global(minCatch:maxCatch)*1.e3 !km->m
    deallocate(buff_global)   

    allocate(buff_global(n_catg),route%K(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt",status="old",action="read");read(77,*)buff_global;close(77)
    route%K=buff_global(minCatch:maxCatch) 
    deallocate(buff_global)  

    allocate(buff_global(n_catg),route%Kstr(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt",status="old",action="read");read(77,*)buff_global;close(77)
    route%Kstr=buff_global(minCatch:maxCatch)
    deallocate(buff_global)     

    allocate(buff_global(n_catg),route%qri_clmt(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_qri.txt",status="old",action="read");read(77,*)buff_global;close(77)
    route%qri_clmt=buff_global(minCatch:maxCatch) !m3/s
    deallocate(buff_global)      

    allocate(buff_global(n_catg),route%qin_clmt(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_qin.txt",status="old",action="read");read(77,*)buff_global;close(77)
    route%qin_clmt=buff_global(minCatch:maxCatch) !m3/s
    deallocate(buff_global)  

    allocate(buff_global(n_catg),route%qstr_clmt(ntiles))   
    open(77,file=trim(inputdir)//"/Pfaf_qstr.txt",status="old",action="read");read(77,*)buff_global;close(77)
    route%qstr_clmt=buff_global(minCatch:maxCatch) !m3/s
    deallocate(buff_global) 

    !if (mapl_am_I_root())then
    !  open(88,file="nsub.txt",action="write")
    !  open(89,file="subarea.txt",action="write")
    !  open(90,file="subi.txt",action="write")
    !  open(91,file="tile_area.txt",action="write")
    !  do i=1,nTiles
    !    write(88,*)route%nsub(i)
    !    write(89,'(150(1x,f10.4))')route%subarea(:,i)
    !    write(90,'(150(i7))')route%subi(:,i)
    !    write(91,*)route%tile_area(i)
    !  enddo
    !  stop
    !endif
   
    deallocate(ims)
    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)      
  end subroutine INITIALIZE
  
! -----------------------------------------------------------
! RUN -- Run method for the route component
! -----------------------------------------------------------
  subroutine RUN1 (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:
! -----------------------------------------------------------

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC
  end subroutine RUN1


  subroutine RUN2 (GC,IMPORT, EXPORT, CLOCK, RC )

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

    character(len=ESMF_MAXSTR)          :: IAm="Run2"
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
    real, dimension(:), pointer :: RUNOFF_SRC0   

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
    INTEGER, SAVE                            :: N_Active, ThisCycle=1  
    INTEGER                                  :: Local_Min, Local_Max
    integer                                  :: K, N, I, req
    REAL                                     :: mm2m3, rbuff, HEARTBEAT 
    REAL, ALLOCATABLE, DIMENSION(:)          :: RUNOFF_CATCH, RUNOFF_ACT,AREACAT_ACT,& 
         LENGSC_ACT, QSFLOW_ACT,QOUTFLOW_ACT
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_index
    type(ESMF_Field) :: runoff_src

    integer                                :: ndes, mype
    type (T_RROUTE_STATE), pointer         :: route => null()
    type (RROUTE_wrap)                     :: wrap
    INTEGER, DIMENSION(:)  ,ALLOCATABLE  :: scounts, scounts_global,rdispls, rcounts  
    real, dimension(:), pointer :: runoff_global,runoff_local,area_local,runoff_cat_global    

    integer :: mpierr, nt_global,nt_local, it, j, upid,cid,temp(1),tid,istat
    integer,save :: nstep_per_day

    type(ESMF_Time) :: CurrentTime, nextTime
    integer :: YY,MM,DD,HH,MMM,SS,YY_next,MM_next,DD_next
    character(len=4) :: yr_s
    character(len=2) :: mon_s,day_s

    real,pointer :: runoff_save(:)=>NULL()
    real,pointer :: WSTREAM_ACT(:)=>NULL()
    real,pointer :: WRIVER_ACT(:)=>NULL()
    real,allocatable :: runoff_save_m3(:),runoff_global_m3(:),QOUTFLOW_GLOBAL(:)
    real,allocatable :: WTOT_BEFORE(:),WTOT_AFTER(:),QINFLOW_LOCAL(:),UNBALANCE(:),UNBALANCE_GLOBAL(:),ERROR(:),ERROR_GLOBAL(:)
    real,allocatable :: QFLOW_SINK(:),QFLOW_SINK_GLOBAL(:),WTOT_BEFORE_GLOBAL(:),WTOT_AFTER_GLOBAL(:)
    real,allocatable :: wriver_global(:),wstream_global(:),qsflow_global(:)
    
    ! ------------------
    ! begin    
    call ESMF_UserCompGetInternalState ( GC, 'RiverRoute_state',wrap,status )
    VERIFY_(STATUS)
    route => wrap%ptr

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS) 
    Iam = trim(COMP_NAME) // "RUN2"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)
    call MAPL_Get(MAPL, HEARTBEAT = HEARTBEAT, RC=STATUS)
    VERIFY_(STATUS)
    !if (mapl_am_I_root()) print *, "HEARTBEAT=",HEARTBEAT 
! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"RUN2")
! Get parameters from generic state
! ---------------------------------

 !   call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
 !   VERIFY_(STATUS) 
! get pointers to inputs variables
! ----------------------------------

    ndes = route%ndes
    mype = route%mype  
    ntiles = route%ntiles  
    nt_global = route%nt_global  
    runoff_save => route%runoff_save
    nt_local = route%nt_local

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
    RUN_MODEL : if (ThisCycle == N_CYC) then   

       runoff_save = runoff_save + RUNOFF_SRC0/real (N_CYC)

       call ESMF_ClockGet(clock, currTime=CurrentTime, rc=status)
       call ESMF_TimeGet(CurrentTime, yy=YY, mm=MM, dd=DD, h=HH, m=MMM, s=SS, rc=status)  
       call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=status)
       call ESMF_TimeGet(nextTime, yy=YY_next, mm=MM_next, dd=DD_next, rc=status) 
       write(yr_s,'(I4.4)')YY
       write(mon_s,'(I2.2)')MM
       write(day_s,'(I2.2)')DD

       allocate(runoff_global(nt_global))
       call MPI_allgatherv  (                          &
          runoff_save,  route%scounts_global(mype+1)      ,MPI_REAL, &
          runoff_global, route%scounts_global, route%rdispls_global,MPI_REAL, &
          MPI_COMM_WORLD, mpierr) 

       if(FirstTime.and.mapl_am_I_root()) print *,"nmax=",nmax
       allocate(RUNOFF_ACT(ntiles))
       RUNOFF_ACT=0.
       do i=1,ntiles
         do j=1,nmax
           it=route%subi(j,i) 
           if(it>0)then
             RUNOFF_ACT(i)=RUNOFF_ACT(i)+route%subarea(j,i)*runoff_global(it)/1000.   
           endif
           if(it==0)exit
         enddo
       enddo  

       deallocate(runoff_global) 


       allocate (AREACAT_ACT (1:ntiles))       
       allocate (LENGSC_ACT  (1:ntiles))
       allocate (QSFLOW_ACT  (1:ntiles))
       allocate (QOUTFLOW_ACT(1:ntiles))     

       LENGSC_ACT=route%lengsc/1.e3 !m->km
       AREACAT_ACT=route%areacat/1.e6 !m2->km2

       WSTREAM_ACT => route%wstream
       WRIVER_ACT => route%wriver

      
       allocate(WTOT_BEFORE(ntiles))
       WTOT_BEFORE=WSTREAM_ACT+WRIVER_ACT

       ! Call river_routing_model
       ! ------------------------     
       !CALL RIVER_ROUTING_LIN  (ntiles, RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,  &
       !     WSTREAM_ACT,WRIVER_ACT, QSFLOW_ACT,QOUTFLOW_ACT) 

       CALL RIVER_ROUTING_HYD  (ntiles, &
            RUNOFF_ACT, route%lengsc, route%lstr, &
            route%qstr_clmt, route%qri_clmt, route%qin_clmt, &
            route%K, route%Kstr, &
            WSTREAM_ACT,WRIVER_ACT, &
            QSFLOW_ACT,QOUTFLOW_ACT)        

       allocate(QOUTFLOW_GLOBAL(n_catg))
       call MPI_allgatherv  (                          &
            QOUTFLOW_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
            QOUTFLOW_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
            MPI_COMM_WORLD, mpierr) 

       allocate(QINFLOW_LOCAL(ntiles))
       QINFLOW_LOCAL=0.
       do i=1,nTiles
         do j=1,upmax
           if(route%upid(j,i)>0)then
             upid=route%upid(j,i)
             WRIVER_ACT(i)=WRIVER_ACT(i)+QOUTFLOW_GLOBAL(upid)*real(route_dt)
             QINFLOW_LOCAL(i)=QINFLOW_LOCAL(i)+QOUTFLOW_GLOBAL(upid)
           else
             exit
           endif
         enddo
       enddo

       !call check_balance(route,ntiles,nt_local,runoff_save,WRIVER_ACT,WSTREAM_ACT,WTOT_BEFORE,RUNOFF_ACT,QINFLOW_LOCAL,QOUTFLOW_ACT,FirstTime,yr_s,mon_s)

       if(FirstTime) nstep_per_day = 86400/route_dt
       route%wriver_acc = route%wriver_acc + WRIVER_ACT/real(nstep_per_day)
       route%wstream_acc = route%wstream_acc + WSTREAM_ACT/real(nstep_per_day)
       route%qoutflow_acc = route%qoutflow_acc + QOUTFLOW_ACT/real(nstep_per_day)
       route%qsflow_acc = route%qsflow_acc + QSFLOW_ACT/real(nstep_per_day)

       deallocate(RUNOFF_ACT,AREACAT_ACT,LENGSC_ACT,QOUTFLOW_ACT,QINFLOW_LOCAL,QOUTFLOW_GLOBAL,QSFLOW_ACT,WTOT_BEFORE)
      !initialize the cycle counter and sum (runoff_tile)       
       WSTREAM_ACT=>NULL()
       WRIVER_ACT=>NULL()      

       runoff_save = 0.
       ThisCycle   = 1           

      ! output
       !if(mapl_am_I_root())print *, "nstep_per_day=",nstep_per_day
       if(mapl_am_I_root())print *, "Current time is ", YY, "/", MM, "/", DD, " ", HH, ":", MMM, ":", SS, ", next MM_next:",MM_next
       if(FirstTime)then
         if(mapl_am_I_root()) istat = mkdir("../river", int(o'755',c_int16_t))  
       endif
       if(HH==23)then
         allocate(wriver_global(n_catg),wstream_global(n_catg),qoutflow_global(n_catg),qsflow_global(n_catg))       
         !call MPI_allgatherv  (                          &
         !     route%wriver_acc,  route%scounts_cat(mype+1)      ,MPI_REAL, &
         !     wriver_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
         !     MPI_COMM_WORLD, mpierr)    
         !call MPI_allgatherv  (                          &
         !     route%wstream_acc,  route%scounts_cat(mype+1)      ,MPI_REAL, &
         !     wstream_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
         !     MPI_COMM_WORLD, mpierr)    
         call MPI_allgatherv  (                          &
              route%qoutflow_acc,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              qoutflow_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)        
         !call MPI_allgatherv  (                          &
         !     route%qsflow_acc,  route%scounts_cat(mype+1)      ,MPI_REAL, &
         !     qsflow_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
         !     MPI_COMM_WORLD, mpierr)                              
         if(mapl_am_I_root())then   
              !open(88,file="../river/river_storage_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              !open(89,file="../river/stream_storage_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              open(90,file="../river/river_flow_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")              
              !open(91,file="../river/stream_flow_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              do i=1,n_catg
                !write(88,*)wriver_global(i)
                !write(89,*)wstream_global(i)
                write(90,*)qoutflow_global(i)
                !write(91,*)qsflow_global(i)
              enddo
              close(90)
              !close(88);close(89);close(90)!;close(91)
              !print *, "output river storage is: ",sum(wriver_global)/1.e9
              !print *, "output stream storage is: ",sum(wstream_global)/1.e9                
         endif           
         deallocate(wriver_global,wstream_global,qoutflow_global,qsflow_global)
         route%wriver_acc = 0.
         route%wstream_acc = 0.
         route%qoutflow_acc = 0.
         route%qsflow_acc = 0.
       endif
 
       !restart
       if(MM_next/=MM)then
         allocate(wriver_global(n_catg),wstream_global(n_catg))
         call MPI_allgatherv  (                          &
              route%wstream,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              wstream_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)
         call MPI_allgatherv  (                          &
              route%wriver,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              wriver_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)    
         if(mapl_am_I_root())then
              write(yr_s,'(I4.4)')YY_next
              write(mon_s,'(I2.2)')MM_next
              write(day_s,'(I2.2)')DD_next
              open(88,file="../input/restart/river_storage_rs_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")
              open(89,file="../input/restart/stream_storage_rs_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt",action="write")   
              do i=1,n_catg
                write(88,*)wriver_global(i)
                write(89,*)wstream_global(i)
              enddo   
              close(88);close(89) 
              print *, "saved river storage is: ",sum(wriver_global)/1.e9
              print *, "saved stream storage is: ",sum(wstream_global)/1.e9                                  
         endif
         deallocate(wriver_global,wstream_global)
       endif

       if(FirstTime) FirstTime=.False.

    else
       
       runoff_save = runoff_save + RUNOFF_SRC0/real (N_CYC)
       
       ThisCycle = ThisCycle + 1

    endif RUN_MODEL 

    runoff_save => NULL()

! All done
! --------
    call MAPL_TimerOff ( MAPL, "-RRM" ) 
    call MAPL_TimerOff(MAPL,"RUN2")
    !call MPI_Barrier(MPI_COMM_WORLD, mpierr)


    RETURN_(ESMF_SUCCESS)
  end subroutine RUN2

! --------------------------------------------------------


  subroutine check_balance(route,ntiles,nt_local,runoff_save,WRIVER_ACT,WSTREAM_ACT,WTOT_BEFORE,RUNOFF_ACT,QINFLOW_LOCAL,QOUTFLOW_ACT,FirstTime,yr_s,mon_s)
      
      type(T_RROUTE_STATE), intent(in) :: route 
      integer, intent(in) :: ntiles,nt_local
      real,intent(in) :: runoff_save(nt_local),WRIVER_ACT(ntiles),WSTREAM_ACT(ntiles),WTOT_BEFORE(ntiles),RUNOFF_ACT(ntiles)
      real,intent(in) :: QINFLOW_LOCAL(ntiles),QOUTFLOW_ACT(ntiles)
      logical,intent(in) :: FirstTime
      character(len=*), intent(in) :: yr_s,mon_s
   
      real,allocatable :: runoff_cat_global(:)
      real,allocatable :: runoff_save_m3(:),runoff_global_m3(:)
      real,allocatable :: WTOT_AFTER(:),UNBALANCE(:),UNBALANCE_GLOBAL(:),ERROR(:),ERROR_GLOBAL(:)
      real,allocatable :: QFLOW_SINK(:),QFLOW_SINK_GLOBAL(:),WTOT_BEFORE_GLOBAL(:),WTOT_AFTER_GLOBAL(:)

      integer :: i, nt_global,mype,cid,temp(1),tid,mpierr
      real :: wr_error, wr_tot, runf_tot

      nt_global = route%nt_global
      mype = route%mype   

         allocate(WTOT_AFTER(ntiles),UNBALANCE(ntiles),UNBALANCE_GLOBAL(n_catg),runoff_cat_global(n_catg))
         allocate(QFLOW_SINK(ntiles),QFLOW_SINK_GLOBAL(n_catg),WTOT_BEFORE_GLOBAL(n_catg),WTOT_AFTER_GLOBAL(n_catg))
         allocate(runoff_save_m3(nt_local),runoff_global_m3(nt_global),ERROR(ntiles),ERROR_GLOBAL(n_catg))

    

         WTOT_AFTER=WRIVER_ACT+WSTREAM_ACT
         ERROR = WTOT_AFTER - (WTOT_BEFORE + RUNOFF_ACT*route_dt + QINFLOW_LOCAL*route_dt - QOUTFLOW_ACT*route_dt)
         where(QOUTFLOW_ACT>0.) UNBALANCE = abs(ERROR)/(QOUTFLOW_ACT*route_dt)
         where(QOUTFLOW_ACT<=0.) UNBALANCE = 0.
         call MPI_allgatherv  (                          &
              UNBALANCE,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              UNBALANCE_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)           
         QFLOW_SINK=0.
         do i=1,ntiles
           if(route%downid(i)==-1)then
              QFLOW_SINK(i) = QOUTFLOW_ACT(i)
           endif
         enddo
         call MPI_allgatherv  (                          &
              QFLOW_SINK,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              QFLOW_SINK_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)
         call MPI_allgatherv  (                          &
              WTOT_BEFORE,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              WTOT_BEFORE_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)   
         call MPI_allgatherv  (                          &
              WTOT_AFTER,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              WTOT_AFTER_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr) 
         runoff_save_m3=runoff_save*route%tile_area/1000. 
         call MPI_allgatherv  (                          &
              runoff_save_m3,  route%scounts_global(mype+1)      ,MPI_REAL, &
              runoff_global_m3, route%scounts_global, route%rdispls_global,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)     
         call MPI_allgatherv  (                          &
              RUNOFF_ACT,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              runoff_cat_global, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)     
         if(mapl_am_I_root())then 
             open(88,file="../runoff_tile_global_"//trim(yr_s)//"_"//trim(mon_s)//".txt",status="unknown", position="append")
             write(88,*)sum(runoff_global_m3)
             close(88)
             open(88,file="../runoff_cat_global_"//trim(yr_s)//"_"//trim(mon_s)//".txt",status="unknown", position="append")
             write(88,*)sum(runoff_cat_global)
             close(88)  
             !print *,"sum(runoff_global_m3)=",sum(runoff_global_m3)
             !print *,"sum(runoff_cat_global)=",sum(runoff_cat_global)   
         endif                   
         if(mapl_am_I_root())then 
             open(88,file="../WTOT_AFTER_"//trim(yr_s)//"_"//trim(mon_s)//".txt",status="unknown", position="append")
             write(88,*)sum(WTOT_AFTER_GLOBAL)
             close(88)
             open(88,file="../WTOT_BEFORE_RUNOFF_QSINK_"//trim(yr_s)//"_"//trim(mon_s)//".txt",status="unknown", position="append")
             write(88,*) sum(WTOT_BEFORE_GLOBAL)+sum(runoff_global_m3)*route_dt-sum(QFLOW_SINK_GLOBAL)*route_dt
             close(88)  
             wr_error=sum(WTOT_AFTER_GLOBAL)-(sum(WTOT_BEFORE_GLOBAL)+sum(runoff_global_m3)*route_dt-sum(QFLOW_SINK_GLOBAL)*route_dt)
             runf_tot=sum(runoff_global_m3)*route_dt
             wr_tot=sum(WTOT_AFTER_GLOBAL)
             open(88,file="../WTOT_ERROR_2_RUNOFF_"//trim(yr_s)//"_"//trim(mon_s)//".txt",status="unknown", position="append")
             write(88,*) wr_error/runf_tot
             close(88)
             open(88,file="../WTOT_ERROR_2_WTOT_"//trim(yr_s)//"_"//trim(mon_s)//".txt",status="unknown", position="append")
             write(88,*) wr_error/wr_tot
             close(88)                 
             !print *,"WTOT_ERROR_2_RUNOFF:",(sum(WTOT_AFTER_GLOBAL)-(sum(WTOT_BEFORE_GLOBAL)+sum(runoff_global_m3)*route_dt-sum(QFLOW_SINK_GLOBAL)*route_dt))/(sum(runoff_global_m3)*route_dt)          
         endif                     

         call MPI_allgatherv  (                          &
              ERROR,  route%scounts_cat(mype+1)      ,MPI_REAL, &
              ERROR_GLOBAL, route%scounts_cat, route%rdispls_cat,MPI_REAL, &
              MPI_COMM_WORLD, mpierr)
         temp = maxloc(abs(ERROR_GLOBAL))
         cid = temp(1)
         if(cid>=route%minCatch.and.cid<=route%maxCatch)then
           tid=cid-route%minCatch+1
           print *,"my PE is:",mype,", max abs value of ERROR=", ERROR(tid)," at pfafid: ",route%minCatch+tid-1,", W_BEFORE=",WTOT_BEFORE(tid),", RUNOFF=",RUNOFF_ACT(tid)*route_dt,", QINFLOW=",QINFLOW_LOCAL(tid)*route_dt,", QOUTFLOW=",QOUTFLOW_ACT(tid)*route_dt,", W_AFTER=",WTOT_AFTER(tid)
         endif  
         !if(FirstTime)then     
         !  if(mapl_am_I_root())then  
         !    open(88,file="ERROR_TOTAL.txt",action="write")
         !    do i=1,n_catg
         !       write(88,*)ERROR_GLOBAL(i)
         !    enddo
         !  endif
         !endif

         deallocate(WTOT_AFTER,UNBALANCE,UNBALANCE_GLOBAL,ERROR,QFLOW_SINK,QFLOW_SINK_GLOBAL,WTOT_BEFORE_GLOBAL,WTOT_AFTER_GLOBAL)
         deallocate(runoff_save_m3,runoff_global_m3,ERROR_GLOBAL,runoff_cat_global)


  end subroutine check_balance


end module GEOS_RouteGridCompMod
