#include "MAPL_Generic.h"

!=============================================================================
module GEOS_LandGridCompMod

!BOP
! !MODULE: GEOS_LandGridCompMod -- A Module to combine VegDyn and Catch, and Igni Gridded Components

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
  use GEOS_IgniGridCompMod,    only : IgniSetServices     => SetServices
!  use GEOS_RouteGridCompMod,   only : RouteSetServices    => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices


!EOP


  integer                                 :: VEGDYN
  integer, allocatable                    :: CATCH(:), ROUTE (:), CATCHCN (:)
  integer                                 :: LSM_CHOICE, RUN_ROUTE, DO_GOSWIM
  integer                                 :: IGNI
  logical                                 :: DO_FIRE_DANGER

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
        tmp(1:ens_id_width)=COMP_NAME(5:5+ens_id_width-1)
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
    call MAPL_GetResource (SCF, DO_FIRE_DANGER, label='FIRE_DANGER:',    DEFAULT=.false., __RC__ )
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
   
    if (DO_FIRE_DANGER) then
        IGNI = MAPL_AddChild(GC, NAME='IGNI'//trim(tmp), SS=IgniSetServices, RC=STATUS)
        VERIFY_(STATUS)
    else
        IGNI = -1
    end if

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

    CASE (1)                                    ! Catchment model
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LST', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TST', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'QST', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DELTS', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DELQS', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ALBVR', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ALBVF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ALBNR', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ALBNF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EMIS', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SNOWMASS', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SNOWDP', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'GHFLX', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPUNST', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPSURF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPWLT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPSAT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ASNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SHSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'AVETSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRSAT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRUST', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRWLT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WET1', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WET2', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WET3', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WCSF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WCRZ', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WCPR', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TP1', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TP2', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TP3', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TP4', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TP5', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TP6', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVAPOUT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SUBLIM', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SHOUT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'RUNOFF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVPINT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVPSOI', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVPVEG',  &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVPICE', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WAT10CM', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WATSOI', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ICESOI', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVPSNO', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'BASEFLOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'RUNSURF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SMELT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'HLWUP',&
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SWNDSRF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LWNDSRF', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'HLATN', &
            CHILD_ID = CATCH(1),&
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'QINFIL', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ACCUM', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'EVLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'PRLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SNOLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DRPARLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DFPARLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LHSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SWNETSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LWUPSNOW', &
            CHILD_ID = CATCH(1) ,&
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LWDNSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TCSORIG', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPSN1IN', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TPSN1OUT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LHLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SHLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SWLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SWDOWNLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LWLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'GHLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'GHSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'GHTSKIN', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SMLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TWLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TELAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TSLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DWLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DHLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SPLAND', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SPLH', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SPWATR', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SPSNOW', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WESNN1', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WESNN2', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WESNN3', &
            CHILD_ID = CATCH(1),   &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'CAPAC',&
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)

    ! the following constants are needed by GEOSldas (to assemble the catparam structure)   

    call MAPL_AddExportSpec(GC, SHORT_NAME = 'COND' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'PSIS' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BEE'  , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'GNU'  , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'VGWMAX',CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF1'  , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF2'  , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF3'  , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'CDCR1', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS1' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS2' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS3' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA1' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA2' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA3' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA4' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW1' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW2' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW3' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW4' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSA1' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSA2' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSB1' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSB2' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ATAU' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BTAU' , CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)

    ! the following constants are needed by GEOSldas and for the "land constants" output collection

    call MAPL_AddExportSpec(GC, SHORT_NAME = 'WPWET', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'CDCR2', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'POROS', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)

    ! the following constants are needed for the "land constants" output collection

    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZGT1', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZGT2', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZGT3', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZGT4', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZGT5', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZGT6', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZPR',  CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZRZ',  CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZSF',  CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'DZTS',  CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)  
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'WPEMW', CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'WPMC',  CHILD_ID = CATCH(1), RC=STATUS); VERIFY_(STATUS)


!   From catment grid internal to be perturbed by land_pert grid
!   WESNN1-3 are originally exported
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'TC', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'QC', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'CATDEF', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'RZEXC', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'SRFEXC', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'HTSNNN1', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'HTSNNN2', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'HTSNNN3', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'SNDZN1', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'SNDZN2', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'SNDZN3', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GHTCNT1', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GHTCNT2', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GHTCNT3', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GHTCNT4', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GHTCNT5', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, &
                              SHORT_NAME = 'GHTCNT6', &
                              CHILD_ID = CATCH(1), &
                              RC=STATUS  )
    VERIFY_(STATUS)
!   Ther are from the first catchment instance
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TH', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'QH', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'CHT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'CQT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'CMT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'CNT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'RIT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'Z0', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'D0', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'Z0H', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'VENT', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'GUST', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOU50M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOV50M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOT10M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOQ10M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOU10M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOV10M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOT2M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOQ2M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOU2M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'MOV2M', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'ITY', &
            CHILD_ID = CATCH(1), &
            RC=STATUS  )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU001', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU002', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU003', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)      
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU004', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU005', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTBC001', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTBC002', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTOC001', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTOC002', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)  
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PEATCLSM_WATERLEVEL',CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PEATCLSM_FSWCHANGE', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)

       if (DO_GOSWIM /= 0) then
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU001', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU002', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU003', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)      
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU004', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU005', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RBC001', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RBC002', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ROC001', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ROC002', CHILD_ID = CATCH(1), RC=STATUS) ; VERIFY_(STATUS)     
       end if

    CASE (2,3)           ! CatchmentCN model
       
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LST',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TST',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QST',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DELTS',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DELQS',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBVR',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBVF',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBNR',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBNF',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EMIS',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNOWMASS', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNOWDP',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHFLX',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPUNST',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSURF',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSNOW',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPWLT',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSAT',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ASNOW',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SHSNOW',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'AVETSNOW', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FRSAT',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FRUST',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FRWLT',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WET1',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WET2',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WET3',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WCSF',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WCRZ',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WCPR',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP1',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP2',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP3',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP4',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP5',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP6',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVAPOUT',  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SUBLIM' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SHOUT'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RUNOFF' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPINT' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPSOI' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPVEG' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPICE' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WAT10CM',  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WATSOI' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ICESOI' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPSNO' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BASEFLOW', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RUNSURF',  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SMELT'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HLWUP'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWNDSRF',  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWNDSRF',  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)                                                                         
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HLATN'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QINFIL' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ACCUM'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PRLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNOLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DRPARLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DFPARLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LHSNOW' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWNETSNOW' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWUPSNOW' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWDNSNOW' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TCSORIG' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSN1IN' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSN1OUT' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LHLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SHLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWDOWNLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHSNOW' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTSKIN' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SMLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TWLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TELAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TSLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DWLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DHLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPLAND' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )              ! a.k.a. SPSHLAND
       VERIFY_(STATUS)
! will need later for CatchCN:
!       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPLH'   ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
!       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPWATR' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )              ! a.k.a. SPEVLAND
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPSNOW' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WESNN1' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WESNN2' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WESNN3' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CAPAC'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)

    ! the following constants are needed by GEOSldas (to assemble the catparam structure)   

    call MAPL_AddExportSpec(GC, SHORT_NAME = 'COND' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'PSIS' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BEE'  , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'GNU'  , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'VGWMAX',CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF1'  , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF2'  , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF3'  , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'CDCR1', CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS1' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS2' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS3' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA1' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA2' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA3' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA4' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW1' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW2' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW3' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW4' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSA1' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSA2' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSB1' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSB2' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ATAU' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BTAU' , CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, SHORT_NAME = 'WPWET', CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'CDCR2', CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'POROS', CHILD_ID = CATCHCN(1), RC=STATUS); VERIFY_(STATUS)

!   From catmentcn grid internal to be perturbed by land_pert grid
!   WESNN1-3 are originally exported
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TC'     , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TG'     , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QC'     , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CATDEF' , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RZEXC'  , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SRFEXC' , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HTSNNN1', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HTSNNN2', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HTSNNN3', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNDZN1' , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNDZN2' , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNDZN3' , CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT1', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT2', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT3', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT4', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT5', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT6', CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       
       ! Unified CN from RUN2 of the first catchment instance
       
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNLAI'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNTLAI' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNSAI'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNTOTC' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNVEGC' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNROOT' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       if (LSM_CHOICE == 3) then
         call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNFROOTC' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
         VERIFY_(STATUS)
       endif
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNNPP'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNGPP'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNSR'   ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNNEE'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNXSMR' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNADD'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNLOSS' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNBURN' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PARABS' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PARINC' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SCSAT'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SCUNS'  ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BTRANT' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SIF'    ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNFSEL' ,  CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TH',       CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QH',       CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CHT',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CQT',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CMT',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNT',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RIT',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'Z0',       CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'D0',       CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'Z0H',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'VENT',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GUST',     CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOU50M',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOV50M',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOT10M',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOQ10M',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOU10M',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOV10M',   CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOT2M',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOQ2M',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOU2M',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOV2M',    CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ITY',      CHILD_ID = CATCHCN(1), RC=STATUS  )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU001', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU002', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU003', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)      
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU004', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU005', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTBC001', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTBC002', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTOC001', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTOC002', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)  
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PEATCLSM_WATERLEVEL',CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PEATCLSM_FSWCHANGE', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)

       if (DO_GOSWIM /= 0) then
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU001', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU002', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU003', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)      
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU004', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU005', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RBC001', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RBC002', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ROC001', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
          call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ROC002', CHILD_ID = CATCHCN(1), RC=STATUS) ; VERIFY_(STATUS)     
       endif

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


    if (DO_FIRE_DANGER) then
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FFMC',        CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GFMC',        CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DMC',         CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DC',          CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ISI',         CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BUI',         CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FWI',         CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DSR',         CHILD_ID = IGNI,  __RC__ )

       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FFMC_DAILY',  CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DMC_DAILY',   CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DC_DAILY',    CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ISI_DAILY',   CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BUI_DAILY',   CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FWI_DAILY',   CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DSR_DAILY',   CHILD_ID = IGNI,  __RC__ )

       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FFMC_DAILY_', CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DMC_DAILY_',  CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DC_DAILY_',   CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ISI_DAILY_',  CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BUI_DAILY_',  CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FWI_DAILY_',  CHILD_ID = IGNI,  __RC__ )
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DSR_DAILY_',  CHILD_ID = IGNI,  __RC__ )

       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'VPD',         CHILD_ID = IGNI,  __RC__ )
    end if


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

          if (DO_FIRE_DANGER) then
              call MAPL_AddConnectivity (                      &
                GC,                                            &
                SHORT_NAME = (/ 'MOT2M     ', 'MOQ2M     ',    &
                                'MOU10M    ', 'MOV10M    ',    &
                                'PRLAND    ', 'SWDOWNLAND',    &
                                'ASNOW     ', 'SNOWDP    ' /), &
                DST_ID = IGNI,                                 &
                SRC_ID = CATCH(I),                             &
                RC = STATUS )
              VERIFY_(STATUS)
          end if

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

          if (DO_FIRE_DANGER) then
              call MAPL_AddConnectivity (                      &
                GC,                                            &
                SHORT_NAME = (/ 'MOT2M     ', 'MOQ2M     ',    &
                                'MOU10M    ', 'MOV10M    ',    &
                                'PRLAND    ', 'SWDOWNLAND',    &
                                'ASNOW     ', 'SNOWDP    ' /), &
                DST_ID = IGNI,                                 &
                SRC_ID = CATCHCN(I),                           &
                RC = STATUS )
              VERIFY_(STATUS)
          end if

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
