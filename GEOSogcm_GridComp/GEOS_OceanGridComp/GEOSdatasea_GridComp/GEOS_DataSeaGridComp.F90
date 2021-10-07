!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_DataSea -- A fake ocean surface

! !INTERFACE:

module GEOS_DataSeaGridCompMod

! !USES: 

  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_DataSea} is a gridded component that reads the 
!   ocean\_bcs file 
!   This module interpolates the SST and SSS data from 
!   either daily or monthly values to the correct time of the simulation.
!   Data are read only if the simulation time is not in the save interval.
!   It also sets surface currents US and VS to constant value (=0.) in 
!   Data Mode.
!

!EOP

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

!  !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.
!
!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = "SetServices"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=STATUS)
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

!BOS

! !Import state:
! None. Its only job is to simply export: SST, SSS, US, VS

!  !Export state:

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'UW',                                &
    LONG_NAME          = 'zonal_velocity_of_surface_water',   &
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'VW',                                &
    LONG_NAME          = 'meridional_velocity_of_surface_water',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'TW',                          &
    LONG_NAME          = 'foundation_temperature_for_interface_layer',&
    UNITS              = 'K',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'SW',                          &
    LONG_NAME          = 'foundation_salinity_for_interface_layer',&
    UNITS              = 'PSU',                               &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

! import

  call MAPL_AddImportSpec(GC,                            &
        SHORT_NAME         = 'FRACICE',                           &
        LONG_NAME          = 'fractional_cover_of_seaice',        &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

!EOS

    call MAPL_TimerAdd(GC,    name="RUN"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-UPDATE" ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run stage for the DataSea component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )


! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the SST and Ice information.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_Time)                    :: CurrentTime
  character(len=ESMF_MAXSTR)          :: DATASeaFILE
! character(len=ESMF_MAXSTR)          :: DATASeaSalFILE
  integer                             :: IFCST
  logical                             :: FCST
  integer                             :: adjSST
  real, pointer, dimension(:,:)       :: SST
  integer                             :: IM
  integer                             :: JM
  real                                :: TICE
  real                                :: CTB  ! Ocean-ice turbulent mixing coefficient (m/sec)
  real                                :: DT
  real                                :: RUN_DT

  real, pointer, dimension(:,:)       :: TNEW   => null()
  real, pointer, dimension(:,:)       :: F1     => null()

! pointers to export

   real, pointer, dimension(:,:)  :: UW
   real, pointer, dimension(:,:)  :: VW
   real, pointer, dimension(:,:)  :: TW
   real, pointer, dimension(:,:)  :: SW

! pointers to import

   real, pointer, dimension(:,:)  :: FI

!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN" )

! Pointers to Imports
!--------------------
    call MAPL_GetPointer(IMPORT,      FI  , 'FRACICE'  , RC=STATUS)
    VERIFY_(STATUS)

!  Pointers to Exports
!---------------------

    call MAPL_GetPointer(EXPORT,      UW  , 'UW'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,      VW  , 'VW'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,      TW  , 'TW'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,      SW  , 'SW'       , RC=STATUS)
    VERIFY_(STATUS)

! Set current time and calendar
!------------------------------

    call ESMF_ClockGet(CLOCK, currTime=CurrentTime, RC=STATUS)
    VERIFY_(STATUS)

! Get the SST bcs file name from the resource file
!-------------------------------------------------

    call MAPL_GetResource(MAPL,DATASeaFILE,LABEL="DATA_SST_FILE:", RC=STATUS)
    VERIFY_(STATUS)

! Get the SSS bcs file name from the resource file
!-------------------------------------------------

!   call MAPL_GetResource(MAPL,DATASeaSalFILE,LABEL="DATA_SSS_FILE:",  RC=STATUS)
!   VERIFY_(STATUS)

! In atmospheric forecast mode we do not have future SST and SSS
!--------------------------------------------------------------

    call MAPL_GetResource(MAPL,IFCST,LABEL="IS_FCST:",default=0,    RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,adjSST,LABEL="SST_ADJ_UND_ICE:",default=0,    RC=STATUS)
    VERIFY_(STATUS)

    FCST = IFCST==1

! SST is usually Reynolds/OSTIA SST or bulk SST
!------------------------------------------------

   call MAPL_Get(MAPL, IM=IM, JM=JM, RC=STATUS)
   VERIFY_(STATUS)

   allocate(SST(IM,JM), stat=STATUS)
   VERIFY_(STATUS)

! SSS is usually bulk SSS
!--------------------------

!  allocate(SSS(IM, JM), stat=STATUS)
!  VERIFY_(STATUS)

!  Update the friendly skin values
!---------------------------------

   call MAPL_TimerOn(MAPL,"-UPDATE" )

!  Read bulk SST from retrieval
!------------------------------

   call MAPL_ReadForcing(MAPL,'SST',DATASeaFILE, CURRENTTIME, SST, INIT_ONLY=FCST, RC=STATUS)
   VERIFY_(STATUS)

!  Read bulk SSS from retrieval
!------------------------------

!  call MAPL_ReadForcing(MAPL,'SSS',DATASeaSalFILE, CURRENTTIME, SSS, INIT_ONLY=FCST, RC=STATUS)
!  VERIFY_(STATUS)

   call MAPL_TimerOff(MAPL,"-UPDATE" )

!  Update the exports
!--------------------

   if(associated(UW)) UW = 0.0
   if(associated(VW)) VW = 0.0

   TICE   = MAPL_TICE-1.8
   if (adjSST == 1) then
      SST = max(SST, TICE)
      SST = (1.-FI)*SST+FI*TICE
   endif

   if (adjSST == 2) then

      call MAPL_GetResource(MAPL,CTB    , LABEL="CTB:"   , default=1.0e-4, RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetResource(MAPL,RUN_DT , LABEL="RUN_DT:" ,                 RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetResource(MAPL,DT     , LABEL="DT:"    , default=RUN_DT, RC=STATUS)
      VERIFY_(STATUS)

      allocate(TNEW(size(TW,1),size(TW,2)), stat=STATUS)
      VERIFY_(STATUS)
      allocate(F1  (size(TW,1),size(TW,2)), stat=STATUS)
      VERIFY_(STATUS)

      TNEW=0.0
      F1  =0.0

!     ! SST below freezing point is set to freezing temperature
      TNEW   = max( SST,TICE)

      where(FI == 1.0)
!     ! if fraction of ice is 1, set SST to freezing temperature        
        TNEW   =  TICE
      elsewhere
        F1=FI*CTB/(2.0*(1.0-FI))
        TNEW=(TNEW+TICE*F1*DT)/(1.0+F1*DT)
      end where

      SST = TNEW

      deallocate( TNEW)
      deallocate( F1)
   endif

   if(associated(TW)) then
        TW = SST        ! SA: SST is in deg Kelvin, hence no need for abs(SST)
   end if

   if(associated(SW)) then
      SW = 30.0       ! SA: for now
!     SW = SSS        ! SA: every SST data point must have SSS (in PSU) as well
   end if

! Clean-up
!---------

   deallocate(SST,     STAT=STATUS); VERIFY_(STATUS)
!  deallocate(SSS,     STAT=STATUS); VERIFY_(STATUS)

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN"  )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)
end subroutine RUN

end module GEOS_DataSeaGridCompMod
