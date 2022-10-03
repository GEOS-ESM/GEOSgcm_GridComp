#include "MAPL_Generic.h"

!=============================================================================
module GEOS_LandGridCompMod

!BOP
! !MODULE: GEOS_LandGridCompMod -- A Module to combine VegDyn and Catch Gridded Components

! !DESCRIPTION: This gridded component operates on the land tiles as
! as child of GEOS\_SurfaceGridComp.  The core functionality is the
! calculation of energy and water fluxes impacting the lowest layer
! of the atmospheric grid.  In order to operate on the tilespace
! specified by its parent, GEOS\_LandGridComp runs its child VegdynGridComp
! to determine relevant time-dependent land-surface characteristics.
! All parameters calculated in VegdynGridComp are required by CatchGridComp.
! Furthermore, several exports of the Vegdyn routines are also exports
! from the Land composite, for use in other modules, such as the case
! for lai and grn needed in radiation.  Vegdyn will be updated first.
! Then the catchment call will be issued.  The composite exports
! consist of the union of the catchment exports with a subset of the 
! vegdyn exports.  All imports and exports are on the prescribed tile
! grid in the (IM, JM)=(NTILES, 1) convention.  

!
! !USES:

  use ESMF
  use MAPL

  use GEOS_VegdynGridCompMod,  only : VegdynSetServices   => SetServices
  use GEOS_CatchGridCompMod,   only : CatchSetServices    => SetServices
  use GEOS_CatchCNGridCompMod, only : CatchCNSetServices  => SetServices
!  use GEOS_RouteGridCompMod,   only : RouteSetServices    => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices


!EOP


  integer                                 :: VEGDYN
  integer, allocatable                    :: CATCH(:), ROUTE (:), CATCHCN (:)
  INTEGER                                 :: LSM_CHOICE, RUN_ROUTE, DO_GOSWIM

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Physics GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (VegDyn, Catch, CatchCN, Route) and runs their respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals
    
    character(len=ESMF_MAXSTR)              :: GCName
    type(ESMF_Config)                       :: CF, SCF
    integer                                 :: NUM_CATCH
    integer                                 :: I
    character(len=ESMF_MAXSTR)              :: TMP
    type(MAPL_MetaComp),pointer             :: MAPL=>null()
    integer                                 :: NUM_LDAS_ENSEMBLE, ens_id_width
    character(len=ESMF_MAXSTR)              :: SURFRC

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                                 ,&
                          NAME=COMP_NAME	             ,&
                          CONFIG=CF                          ,&
                                                    RC=STATUS )
    VERIFY_(STATUS)

    Iam = trim(COMP_NAME) // 'SetServices'

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, NUM_LDAS_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, ens_id_width, Label="ENS_ID_WIDTH:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    tmp = ''
    if (NUM_LDAS_ENSEMBLE >1) then
        ! land_exxxx
        ! here ENS_ID_WIDTH is the width of digits. add 2 for '_e'
        tmp(1:ens_id_width+2)=COMP_NAME(5:5+ens_id_width-1+2)
    endif

!------------------------------------------------------------
! Register services for this component
!------------------------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run2, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute ( CF, NUM_CATCH, Label="NUM_CATCH_ENSEMBLES:", default=1, RC=STATUS)
    VERIFY_(STATUS)

!------------------------------------------------------------
! Create children's gridded components and invoke their 
! SetServices
!------------------------------------------------------------

    VEGDYN  = MAPL_AddChild(GC, NAME='VEGDYN'//trim(tmp), SS=VegdynSetServices, RC=STATUS)
    VERIFY_(STATUS)

! Get CHOICE OF  Land Surface Model (1:Catch, 2:Catch-CN)
! and Runoff Routing Model (0: OFF, 1: ON)
! -------------------------------------------------------

    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call MAPL_GetResource (SCF, RUN_ROUTE, label='RUN_ROUTE:',           DEFAULT=0, __RC__ )
    call MAPL_GetResource (SCF, DO_GOSWIM, label='N_CONST_LAND4SNWALB:', DEFAULT=0, __RC__ )
    call ESMF_ConfigDestroy      (SCF, __RC__)

    SELECT CASE (LSM_CHOICE)

    CASE (1) 
    
       allocate (CATCH(NUM_CATCH), stat=status)
       VERIFY_(STATUS)
       if (NUM_CATCH == 1) then
          CATCH(1) = MAPL_AddChild(GC, NAME='CATCH'//trim(tmp), SS=CatchSetServices, RC=STATUS)
          VERIFY_(STATUS)
       else
          do I = 1, NUM_CATCH
             WRITE(TMP,'(I3.3)') I
             GCName  = 'ens' // trim(TMP) // ':CATCH'
             CATCH(I) = MAPL_AddChild(GC, NAME=GCName, SS=CatchSetServices, RC=STATUS)
             VERIFY_(STATUS)
          end do
       end if
       
    CASE (2,3) 
       
       allocate (CATCHCN(NUM_CATCH), stat=status)
       VERIFY_(STATUS)
       if (NUM_CATCH == 1) then
          CATCHCN(1) = MAPL_AddChild(GC, NAME='CATCHCN'//trim(tmp), SS=CatchCNSetServices, RC=STATUS)
          VERIFY_(STATUS)
       else
          do I = 1, NUM_CATCH
             WRITE(TMP,'(I3.3)') I
             GCName  = 'ens' // trim(TMP) // ':CATCHCN'
             CATCHCN(I) = MAPL_AddChild(GC, NAME=GCName, SS=CatchCNSetServices, RC=STATUS)
             VERIFY_(STATUS)
          end do
       end if
       
    END SELECT

!    IF(RUN_ROUTE == 1) THEN
!       if (NUM_CATCH == 1) then
!          ROUTE(1) = MAPL_AddChild(GC, NAME='ROUTE', SS=RouteSetServices, RC=STATUS)
!          VERIFY_(STATUS)
!       else
!          do I = 1, NUM_CATCH
!             WRITE(TMP,'(I3.3)') I
!             GCName  = 'ens' // trim(TMP) // ':ROUTE'
!             ROUTE(I) = MAPL_AddChild(GC, NAME=GCName, SS=RouteSetServices, RC=STATUS)
!             VERIFY_(STATUS)
!          end do
!       end if
!    ENDIF
    
!BOS

!------------------------------------------------------------
! Set UNIQUE IMPORT or EXPORT specs for the Composite GC
! There are none of these at the moment...
!------------------------------------------------------------

! !IMPORT STATE:

! There are no explicit imports

! !EXPORT STATE:

! Selected exports of the children are made explicit exports of land

! These are from RUN2 of the first catchment instance
    SELECT CASE (LSM_CHOICE)

    CASE (1)
       call MAPL_AddExportSpec(GC, CATCH(1),_RC) 
    CASE (2,3) 
       call MAPL_AddExportSpec(GC, CATCHCN(1),_RC) 
    END SELECT

! These are from RUN1 of vegdyn and the first catchment instance
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'LAI', &
                              CHILD_ID = VEGDYN, &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GRN', &
                              CHILD_ID = VEGDYN, &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'ROOTL', &
                              CHILD_ID = VEGDYN, &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'NDVI', &
                              CHILD_ID = VEGDYN,&
                              RC=STATUS  )
    VERIFY_(STATUS) 
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'Z2CH', &
                              CHILD_ID = VEGDYN,&
                              RC=STATUS  )
    VERIFY_(STATUS) 
!    IF(RUN_ROUTE == 1) THEN
!       call MAPL_AddExportSpec ( GC, &
!            SHORT_NAME = 'QOUTFLOW', &
!            CHILD_ID = ROUTE(1),     &
!            RC=STATUS  )
!       VERIFY_(STATUS)       
!    ENDIF

!EOS
    
!------------------------------------------------------------
! Set internal connections between the children's IMPORTS 
! and EXPORTS... this will be a one-way street between the
! vegdyn routine and the catchment model
!------------------------------------------------------------

! !CONNECTIONS:

    DO I = 1, NUM_CATCH

       SELECT CASE (LSM_CHOICE)

       CASE (1) 
          call MAPL_AddConnectivity (                                    &
            GC                                                 ,         &
            SHORT_NAME  = (/'LAI    ', 'GRN    ', 'ROOTL  ', 'Z2CH   ',  &
                            'ITY    ', 'ASCATZ0', 'NDVI   '/) ,          &
            DST_ID =  CATCH(I)                                 ,         &
            SRC_ID =  VEGDYN                                   ,         &
                                                      RC=STATUS )
          VERIFY_(STATUS)

!          IF(RUN_ROUTE == 1) THEN
!             call MAPL_AddConnectivity (                              &
!                  GC                                                 ,&
!                  SHORT_NAME  = (/'RUNOFF  '/)                       ,&
!                  SRC_ID =  CATCH(I)                                 ,&
!                  DST_ID =  ROUTE(I)                                 ,&
!                  
!                  RC=STATUS )
!             VERIFY_(STATUS)            
!          ENDIF

       CASE (2,3)
          call MAPL_AddConnectivity (                                    & 
            GC                                                 ,         &
            SHORT_NAME  = (/'LAI    ', 'GRN    ', 'ROOTL  ', 'Z2CH   ',  &
                            'ITY    ', 'ASCATZ0', 'NDVI   ' /) ,         &
            DST_ID =  CATCHCN(I)                               ,         &
            SRC_ID =  VEGDYN                                   ,         &
                                                      RC=STATUS ) 

!          IF(RUN_ROUTE == 1) THEN
!             call MAPL_AddConnectivity (                              &
!                  GC                                                 ,&
!                  SHORT_NAME  = (/'RUNOFF  '/)                       ,&
!                  SRC_ID =  CATCHCN(I)                               ,&
!                  DST_ID =  ROUTE(I)                                 ,&
!                  
!                  RC=STATUS )
!             VERIFY_(STATUS)            
!          ENDIF
       END SELECT
    END DO


    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN1"          ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN2"          ,RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_GenericSetServices(GC, RC=STATUS )
    VERIFY_(STATUS)

    if (allocated(CATCH)) deallocate(CATCH)
    if (allocated(CATCHCN)) deallocate(CATCHCN)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! !IROUTINE: Initialize -- Initialize method for the composite Surface Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Land Composite Gridded Component.
 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)           :: IAm 
    integer                              :: STATUS
    character(len=ESMF_MAXSTR)           :: COMP_NAME
    
! Local derived type aliases

    type (MAPL_MetaComp    ), pointer       :: MAPL
    type (MAPL_MetaComp    ), pointer       :: CHILD_MAPL 
    type (MAPL_LocStream       )            :: LOCSTREAM
    type (ESMF_DELayout        )            :: LAYOUT
    type (ESMF_Config          )            :: CF
    type (ESMF_GridComp        ), pointer   :: GCS(:)
  
    integer                                 :: I

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOn(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)

! Get the land tilegrid and the child components
!----------------------------------------------- 

    call MAPL_Get (MAPL, LOCSTREAM=LOCSTREAM, GCS=GCS, RC=STATUS )
    VERIFY_(STATUS)

! Place the land tilegrid in the generic state of each child component
!---------------------------------------------------------------------

    do I = 1, SIZE(GCS)
       call MAPL_GetObjectFromGC( GCS(I), CHILD_MAPL, RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_Set (CHILD_MAPL, LOCSTREAM=LOCSTREAM, RC=STATUS )
       VERIFY_(STATUS)
    end do

    call MAPL_TimerOff(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)

! Call Initialize for every Child
!--------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"INITIALIZE", RC=STATUS ); VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


!=============================================================================

!BOP
! !IROUTINE: Run1 -- First Run method for the composite Land Gridded Component
! !INTERFACE:
  subroutine Run1(GC, IMPORT, EXPORT, CLOCK, RC )
! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: This first run method calls the children's
!   first run methods. VEGDYN has only one, and it is called here.
!EOP

! ErrLog Variables

   character(len=ESMF_MAXSTR)          :: IAm 
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),      pointer  :: MAPL
   type (ESMF_GridComp),      pointer  :: GCS(:)
   type (ESMF_State),         pointer  :: GIM(:)
   type (ESMF_State),         pointer  :: GEX(:)
   character(len=ESMF_MAXSTR),pointer  :: GCnames(:)
   integer                             :: I 

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run1"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOn(MAPL,"RUN1", RC=STATUS ); VERIFY_(STATUS)

    call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCnames=GCnames,rc=STATUS)
    VERIFY_(STATUS)

! Call the children's RUN methods
!--------------------------------

    DO I = 1, size(GCS)
       call MAPL_TimerOn(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
       call ESMF_GridCompRun(GCS(I), importState=GIM(I), exportState=GEX(I), &
                             CLOCK=CLOCK, PHASE=1, userRC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
    END DO

    call MAPL_TimerOff(MAPL,"RUN1", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine Run1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: Run2 -- Second Run method for the composite Land Gridded Component
! !INTERFACE:
  subroutine Run2(GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: This second run method call only the catchments second method.
!EOP

! ErrLog Variables

   character(len=ESMF_MAXSTR)          :: IAm 
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),  pointer  :: MAPL
   type (ESMF_GridComp),  pointer  :: GCS(:)
   type (ESMF_State),     pointer  :: GIM(:)
   type (ESMF_State),     pointer  :: GEX(:)
   character(len=ESMF_MAXSTR),pointer  :: GCnames(:)

   integer :: I
  
!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run2"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOn(MAPL,"RUN2", RC=STATUS ); VERIFY_(STATUS)

    call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCnames=GCnames,rc=STATUS)
    VERIFY_(STATUS)

! Call the children's RUN methods
!--------------------------------
    DO I=1,size(GCS)
       if (I == VEGDYN) cycle
       call MAPL_TimerOn(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
       call ESMF_GridCompRun(GCS(I), importState=GIM(I), exportState=GEX(I), &
                             CLOCK=CLOCK, PHASE=2, userRC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
    END DO

    call MAPL_TimerOff(MAPL,"RUN2", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)
    
    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_LandGridCompMod
