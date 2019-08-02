
#include "MAPL_Generic.h"

!-----------------------------------------------------------------------
!              ESMA - Earth System Modeling Applications
!-----------------------------------------------------------------------
   Module CreateInterpWeights_GridCompMod

!BOP
!
! !MODULE: CreateInterpWeights_GridCompMod
!
! !USES:

   use ESMF                ! ESMF base class
   use MAPL_Mod            ! GEOS base class

! !PUBLIC MEMBER FUNCTIONS:

  implicit none
  private

  public  SetServices      ! Register component methods
!EOP

contains

!----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices

! !DESCRIPTION:  SetServices registers Initialize, Run, and Finalize
!   methods for FV. Two stages of the FV run method are registered. The
!   first one does the dynamics calculations, and the second adds 
!   increments from external sources that appear in the Import state.
!   SetServices also creates a private internal state in which FV
!   keeps invariant or auxilliary state variables, as well as pointers to
!   the true state variables. The MAPL internal state contains the
!   true state variables and is managed by MAPL.
!
! !INTERFACE:

   Subroutine SetServices ( gc, rc )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gc     ! gridded component
   integer, intent(out), optional     :: rc     ! return code
    

!EOP         
!----------------------------------------------------------------------
  
    integer                          :: status
    character(len=ESMF_MAXSTR)       :: IAm
    character(len=ESMF_MAXSTR)       :: COMP_NAME

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    _VERIFY(STATUS)
    Iam = trim(COMP_NAME) // "SetServices"
 
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize, rc=status)
    _VERIFY(STATUS)
 
! Generic SetServices
!--------------------

    call MAPL_GenericSetServices( GC, RC=STATUS )
    _VERIFY(STATUS)

    _RETURN(ESMF_SUCCESS)

  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component 
  type(ESMF_State),    intent(inout) :: import   ! import state
  type(ESMF_State),    intent(inout) :: export   ! export state
  type(ESMF_Clock),    intent(inout) :: clock    ! the clock
  
  integer, intent(out), OPTIONAL     :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error

  integer                            :: status
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME

  type (MAPL_MetaComp), pointer      :: MAPL
  type (ESMF_Config)                 :: CF
  type (ESMF_VM)                     :: VM
  character (len=ESMF_MAXSTR)        :: strTxt
  character (len=ESMF_MAXSTR)        :: layout_file
  integer :: npx, npy
  integer :: nlon, nlat
  real, allocatable :: data_ll(:,:), data_cs(:,:)

  integer :: numtasks
  integer :: comm, gid, pe0, pe1, pe2, pe3, pe4, pe5, pe6, pe7

! Begin
!------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
    _VERIFY(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Call Generic Initialize
!------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    _VERIFY(STATUS)

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    _VERIFY(STATUS)

    call ESMF_VMGetCurrent(VM, rc=status)

    call ESMF_VMGet(VM,mpiCommunicator=comm,rc=status)

  write(strTxt,'(A,i5.5)') trim(Iam), __LINE__
  call MAPL_MemUtilsWrite(VM, strTxt, RC=STATUS )
  _VERIFY(STATUS)

    call MAPL_GetResource ( MAPL, layout_file, 'LAYOUT:', default='weights.rc', rc=status )
    _VERIFY(STATUS)
    call ESMF_ConfigLoadFile( cf, LAYOUT_FILE, rc = rc )
    call ESMF_ConfigGetAttribute   ( cf, npx, label = 'npx:', default=180, rc = rc )
    npy = npx*6
    call ESMF_ConfigGetAttribute   ( cf, nlon, label = 'nlon:', default=360, rc = rc )
    call ESMF_ConfigGetAttribute   ( cf, nlat, label = 'nlat:', default=181, rc = rc )

    call write_parallel(npx)
    call write_parallel(npy)
    call write_parallel(nlon)
    call write_parallel(nlat)

    call ESMF_VMGet(vm, petCount=numtasks, rc=status)
    call ESMF_VMGet(vm, localPet=gid, rc=status)

   !call GetWeights_init (6,1,npx,npy,1,&
   !     1,6,.false.,.true.,comm)

    call write_parallel('Generating Weights')
    allocate ( data_ll(nlon,nlat) )
    allocate ( data_cs(npx ,npy ) )
    data_cs(:,:) = 1.0
    data_ll(:,:) = 1.0
    call cube2latlon(npx, npy, nlon, nlat, data_cs, data_ll)
    deallocate( data_ll, data_cs )

    _RETURN(ESMF_SUCCESS)
  end subroutine Initialize
 
end module CreateInterpWeights_GridCompMod
 
