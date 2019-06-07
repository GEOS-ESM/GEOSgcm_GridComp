
#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: DFI -- A Module to handle Digital Filter Initialization

! !INTERFACE:

module DFI_GridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!

!EOP

! Common block (for now)
  integer,  parameter :: r8 = 8

  real,allocatable,dimension(:) :: dfi(:) ! _RT temporary common block
  integer :: idfi=0 ! _RT temporary counter indexing steps in DFI procedure
  integer :: NSTEPS ! _RT temporary too

contains

!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

  subroutine SetServices ( GC, RC )

! ! ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    type (MAPL_MetaComp),         pointer   :: MAPL

    integer DFI

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run,  &
                                      RC=STATUS)
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

! !IMPORT STATE:

!   call MAPL_AddImportSpec ( gc,                                  &
!        SHORT_NAME = 'AK',                                        &
!        LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
!        UNITS      = 'Pa',                                        &
!        DIMS       = MAPL_DimsVertOnly,                           &
!        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
!    VERIFY_(STATUS)

!   call MAPL_AddImportSpec ( gc,                                  &
!        SHORT_NAME = 'BK',                                        &
!        LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
!        UNITS      = '1',                                         &
!        DIMS       = MAPL_DimsVertOnly,                           &
!        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
!    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'U_DGRID',                                   &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'V_DGRID',                                   &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         & !_RT: check name/unit
         UNITS      = 'kg kg^${-1}$',                              &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'OX',                                        &
         LONG_NAME  = 'molecular_oxigen',                          & !_RT: check name/unit
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

! !INTERNAL STATE:

!   call MAPL_AddInternalSpec ( gc,                                &
!        SHORT_NAME = 'AK',                                        &
!        LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
!        UNITS      = 'Pa',                                        &
!        PRECISION  = ESMF_KIND_R8,                                &
!        DIMS       = MAPL_DimsVertOnly,                           &
!        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
!    VERIFY_(STATUS)

!   call MAPL_AddInternalSpec ( gc,                                &
!        SHORT_NAME = 'BK',                                        &
!        LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
!        UNITS      = '1',                                         &
!        PRECISION  = ESMF_KIND_R8,                                &
!        DIMS       = MAPL_DimsVertOnly,                           &
!        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
!    VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'U_DGRID',                                   &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'V_DGRID',                                   &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         & !_RT: check name/unit
         UNITS      = 'kg kg^${-1}$',                              &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'OX',                                        &
         LONG_NAME  = 'molecular_oxigen',                          & !_RT: check name/unit
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, name="INITIALIZE" ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the DFI

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

   type (MAPL_MetaComp),  pointer      :: MAPL
   type (ESMF_Time)                    :: CurrTime, RingTime, BegTime
   type (ESMF_TimeInterval)            :: TIMEINT
   type (ESMF_Alarm)                   :: ALARM
   type (ESMF_Config)                  :: cf
   integer                             :: I, NQ
   real                                :: DT             
   integer                             :: DFI_BEG_DATE
   integer                             :: DFI_BEG_TIME
   integer                             :: BEG_TIME(6)
   character(len=ESMF_MAXSTR)          :: FilterType
   character(len=ESMF_MAXSTR)          :: FilterVerb
   real                                :: TAUANL
   real                                :: POFFSET

! =============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Call Initialize for every Child

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

!  Initialize digital filter/iau coefficients
!  ------------------------------------------
   call MAPL_GetResource( MAPL, FilterType, Label="FILTER_TYPE:"   , default="NULL" , RC=STATUS)
   VERIFY_(STATUS)
   if (trim(FilterType) /= 'DFI' ) then ! nothing to do otherwise
       RETURN_(ESMF_SUCCESS)
   endif 

   call MAPL_GetResource( MAPL, DT, Label="RUN_DT:", RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetResource( MAPL, TAUANL, label="TAUANL:", default=21600., rc=status)
   VERIFY_(STATUS)
   call MAPL_GetResource( MAPL, POFFSET, label="PREDICTOR_OFFSET:", default=21600., rc=status)
   VERIFY_(STATUS)
   call MAPL_GetResource( MAPL, FilterVerb, Label="FILTER_VERBOSE:", default="NO"  , RC=STATUS)
   VERIFY_(STATUS)
   nsteps = POFFSET/DT+1
   allocate(dfi(nsteps))
   call dfi_coeffs (FilterType,DT,POFFSET,nsteps,dfi)
   if (MAPL_am_I_root().and.(trim(FilterVerb)=="YES".or.trim(FilterVerb)=="yes")) then
      print*, 'DFI initialized for this many steps: ', nsteps 
      do i=1,nsteps
         print*, 'i,dfi-coeff=', i, dfi(i)
      enddo
   endif

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"INITIALIZE")

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif


    RETURN_(ESMF_SUCCESS)
 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: RUN -- Run method for DFI component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: Run main component of digital filter. Simply accummulate
!                state applying the coefficients of the filter as it goes
!                along.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State)                   :: INTERNAL

  integer                             :: IM, JM, LM
  real, pointer, dimension(:,:,:)     :: u
  real, pointer, dimension(:,:,:)     :: v
  real, pointer, dimension(:,:,:)     :: pt
  real, pointer, dimension(:,:,:)     :: pe
  real, pointer, dimension(:,:,:)     :: qv
  real, pointer, dimension(:,:,:)     :: ox

  real(r8), pointer, dimension(:,:,:)     :: dfu
  real(r8), pointer, dimension(:,:,:)     :: dfv
  real(r8), pointer, dimension(:,:,:)     :: dfpt
  real(r8), pointer, dimension(:,:,:)     :: dfpe
  real    , pointer, dimension(:,:,:)     :: dfqv
  real    , pointer, dimension(:,:,:)     :: dfox

  type(ESMF_Grid)                         :: grid
  type(ESMF_Alarm)                        :: Alarm
  type(ESMF_Time)                         :: currTime
  character(len=ESMF_MAXSTR)              :: FilterVerb
  character(len=ESMF_MAXSTR)              :: FilterType
  integer :: it,L
  logical :: do_dfi


!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run"
   call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

   call MAPL_GetResource( MAPL, FilterType, Label="FILTER_TYPE:", default="NULL" , RC=STATUS)
   VERIFY_(STATUS)
   if (trim(FilterType) /= 'DFI' ) then ! nothing to do otherwise
       RETURN_(ESMF_SUCCESS)
   endif 

   call MAPL_GetResource( MAPL, FilterVerb, Label="FILTER_VERBOSE:", default="NO"  , RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_ClockGet(Clock, CurrTime=currTIME, rc=status)
   VERIFY_(status)

   call check_dfi_time_(do_dfi)

   if(.not. do_dfi) then
      RETURN_(ESMF_SUCCESS)
   endif

! If it's time to accumulate state start by setting counter
! ---------------------------------------------------------
   call set_dfi()

!  call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,    &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)


! **********************************************************************
! ****               Get Pointers to BKG Import Data                ****
! **********************************************************************
#if 0
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
    end if
#endif

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(import,    u, 'U_DGRID',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,    v, 'V_DGRID',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   pt, 'PT'     ,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   pe, 'PE'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   qv, 'Q'      , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,   ox, 'OX'     , RC=STATUS)
    VERIFY_(STATUS)
    
!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(internal,    dfu, 'U_DGRID',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,    dfv, 'V_DGRID',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dfpt, 'PT'     ,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dfpe, 'PE'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dfqv, 'Q'      , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   dfox, 'OX'     , RC=STATUS)
    VERIFY_(STATUS)

!   Apply filter coeffient at this time (indx)
!   ------------------------------------------
    it = idfi
    if (it==0) then
        call WRITE_PARALLEL ( trim(Iam)//": troubled DFI counter" )
        STATUS=99
        VERIFY_(STATUS)
    endif
    if ( trim(FilterVerb)=='yes' .or. trim(FilterVerb)=='YES') then
        call WRITE_PARALLEL ( trim(Iam)//": accumulating DFI now" )
    endif
    dfu    = dfu   + dfi(it) * u
    dfv    = dfv   + dfi(it) * v
    dfpt   = dfpt  + dfi(it) * pt
    dfpe   = dfpe  + dfi(it) * pe
    dfqv   = dfqv  + dfi(it) * qv
    dfox   = dfox  + dfi(it) * ox

    RETURN_(ESMF_SUCCESS)

  contains

  subroutine check_dfi_time_(do_dfi)

     logical, intent(out) :: do_dfi
     type(ESMF_Time) :: BegTime, EndTime
     type(ESMF_TimeInterval) :: RefTGap, RefDT
     integer         :: DFI_BEG_DATE
     integer         :: DFI_BEG_TIME
     integer         :: BEG_TIME(6), REF_TGAP(6), REF_DT(6)
     real            :: DT, POFFSET

     call MAPL_GetResource( MAPL, DT, Label="RUN_DT:", RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetResource( MAPL, POFFSET     , label="PREDICTOR_OFFSET:", default=21600., rc=status)
     VERIFY_(STATUS)
     CALL MAPL_GetResource( MAPL, DFI_BEG_DATE, label='DFI_BEG_DATE:', default=-1, rc=status )
     VERIFY_(STATUS)
     CALL MAPL_GetResource( MAPL, DFI_BEG_TIME, label='DFI_BEG_TIME:', default=-1, rc=status )
     VERIFY_(STATUS)

     BEG_TIME(1) =     DFI_BEG_DATE/10000
     BEG_TIME(2) = mod(DFI_BEG_DATE,10000)/100
     BEG_TIME(3) = mod(DFI_BEG_DATE,100)
     BEG_TIME(4) =     DFI_BEG_TIME/10000
     BEG_TIME(5) = mod(DFI_BEG_TIME,10000)/100
     BEG_TIME(6) = mod(DFI_BEG_TIME,100)

!   set replay time
!   ---------------
     call ESMF_TimeSet(  BegTime, YY =  BEG_TIME(1), &
                                  MM =  BEG_TIME(2), &
                                  DD =  BEG_TIME(3), &
                                  H  =  BEG_TIME(4), &
                                  M  =  BEG_TIME(5), &
                                  S  =  BEG_TIME(6), rc=STATUS ); VERIFY_(STATUS)

!    Offset due to T0
!    ----------------
     REF_TGAP    = 0
     REF_TGAP(4) = nint(POFFSET)/3600
     REF_TGAP(5) = nint(POFFSET-REF_TGAP(4)*3600.)/60
     REF_TGAP(6) = nint(POFFSET-REF_TGAP(4)*3600.-REF_TGAP(5)*60.)
     call ESMF_TimeIntervalSet(  RefTGap, YY = REF_TGAP(1), &
                                          MM = REF_TGAP(2), &
                                           D = REF_TGAP(3), &
                                           H = REF_TGAP(4), &
                                           M = REF_TGAP(5), &
                                           S = REF_TGAP(6), &
                                   startTime = BegTime,    &
                                               rc = STATUS  ); VERIFY_(STATUS)
!    Offset due to DT
!    ----------------
     REF_DT   = 0
     REF_DT(4) = nint(DT)/3600
     REF_DT(5) = nint(DT-REF_DT(4)*3600.)/60
     REF_DT(6) = nint(DT-REF_DT(4)*3600.-REF_DT(5)*60.)
     call ESMF_TimeIntervalSet(  RefDT, YY = REF_DT(1), &
                                        MM = REF_DT(2), &
                                         D = REF_DT(3), &
                                         H = REF_DT(4), &
                                         M = REF_DT(5), &
                                         S = REF_DT(6), &
                                 startTime = BegTime,   &
                                             rc = STATUS  ); VERIFY_(STATUS)

    BegTime = BegTime - RefDT
    EndTime = BegTime + RefTGap

    do_dfi=.false.
    if(currTime>=BegTime.and.currTime<=EndTime) do_dfi=.true.
    end subroutine check_dfi_time_

  end subroutine RUN

  subroutine set_dfi
    implicit none
    idfi=idfi+1
    if(idfi>nsteps) idfi=0
  end subroutine set_dfi

  subroutine dfi_coeffs (FilterOpt,DT,TC,nsteps,dfi)
! This subroutine belongs to GEOS_Shared, but for now it lives here
   implicit none

   character(len=*),intent(in):: FilterOpt ! IAU or DFI
   real,   intent(in)  :: DT     ! model time step
   real,   intent(in)  :: TC     ! cutoff time (typically 21600 sec)
   integer,intent(in)  :: nsteps ! number of steps TC/DT+1
   real,   intent(out) :: dfi(nsteps)

   real pi,arg,wc,thetac
   integer n,i,k,nhlf,np1

!  If regular IAU, simply set to constant and return
!  -------------------------------------------------
   if (trim(FilterOpt)=='IAU' .or. trim(FilterOpt)=='iau' ) then
      dfi=1.0
      return
   endif

!  Calculate DFI coefficients
!  --------------------------
   pi = 4.0*atan(1.)
   thetac = TC ! cutoff
   nhlf = (nsteps+1)/2
   do k = 1, nhlf-1
      n   = k-nhlf
      arg = n*pi/nhlf            
      wc  = sin(arg)/arg ! Lanczos window
      dfi(k) = wc*sin(n*2.0*pi*DT/tc)/(n*pi)
   end do
   dfi(nhlf) = 2*DT/TC
   do i = nhlf+1, nsteps
      dfi(i) = dfi(nsteps-i+1)
   end do

!  Normalize coefficients
!  ----------------------
   dfi = dfi/sum(dfi)

  end subroutine dfi_coeffs

end module DFI_GridCompMod
