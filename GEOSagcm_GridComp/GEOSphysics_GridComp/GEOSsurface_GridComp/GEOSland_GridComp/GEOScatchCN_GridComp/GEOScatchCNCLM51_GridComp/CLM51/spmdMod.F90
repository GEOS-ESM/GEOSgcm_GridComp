
module spmdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdMod
!
! !DESCRIPTION:
! SPMD initialization
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  use ESMF
  use MAPL
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: iulog
  implicit none

  private

#include <mpif.h>

  save

  ! Default settings valid even if there is no spmd 

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: iam             ! processor number
  integer, public :: npes            ! number of processors for clm
  integer, public :: mpicom          ! communicator group for clm
  integer, public :: comp_id         ! component id

  !
  ! Public methods
  !
  public :: spmd_init                ! Initialization

  !
  ! Values from mpif.h that can be used
  !
  public :: MPI_INTEGER
  public :: MPI_REAL8
  public :: MPI_LOGICAL
  public :: MPI_SUM
  public :: MPI_MIN
  public :: MPI_MAX
  public :: MPI_LOR
  public :: MPI_STATUS_SIZE
  public :: MPI_ANY_SOURCE
  public :: MPI_CHARACTER
  public :: MPI_COMM_WORLD
  public :: MPI_MAX_PROCESSOR_NAME

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init( clm_mpicom )
!
! !INTERFACE:
  subroutine spmd_init()
!
! !DESCRIPTION:
! MPI initialization (number of cpus, processes, tids, etc)
!
! !USES
!
! !ARGUMENTS:
    implicit none
     type(ESMF_VM) :: vm
     integer :: status     ! Error code
!    integer, intent(in) :: clm_mpicom
!    integer, intent(in) :: LNDID
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i     ! indices
    integer  :: npes  ! MPI size
    integer  :: MYID  ! MPI Rank
!-----------------------------------------------------------------------

    call ESMF_VmGetCurrent(VM, rc=status)

    ! Get MPI communicator

    call ESMF_VmGet(VM, mpicommunicator=mpicom, RC=status) 

    ! Get my processor id and number of processors

    call ESMF_VmGet(VM, localPet=MYID, petCount=npes, RC=status)

    ! determine master process
    if (MAPL_Am_I_Root(vm)) then
       masterproc = .true.
    else
       masterproc = .false.
    end if

    if (masterproc) then
       write(iulog,100)npes
       write(iulog,200)
       write(iulog,220)
       do i=0,npes-1
          write(iulog,250)i,MYID
       end do
    endif


100 format(//,i3," pes participating in computation for CLM")
200 format(/,35('-'))
220 format(/,"NODE#",2x,"NAME")
250 format("(",i5,")",2x,100a1,//)

  end subroutine spmd_init

end module spmdMod
