#include "MAPL_Generic.h"

module ncdio_pio

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: ncdio_pioMod
  !
  ! !DESCRIPTION:
  ! Generic interfaces to write fields to netcdf files for CLM
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8, i4=>shr_kind_i4, shr_kind_cl, r4 => shr_kind_r4
  !use shr_infnan_mod , only : nan => shr_infnan_nan,  isnan => shr_infnan_isnan
  use nanMod         , only : nan
  use shr_sys_mod    , only : shr_sys_abort
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use MAPL           , only : file_desc_t =>  NetCDF4_FileFormatter, pFIO_READ
  use MAPL_ExceptionHandling

  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public :: ncd_pio_openfile   ! open a file
  public :: ncd_pio_closefile  ! close a file
  public :: ncd_io             ! write local data

  public file_desc_t
  !

  interface ncd_io

    module procedure ncd_io_char_0d
    module procedure ncd_io_char_1d
   ! module procedure ncd_io_log_1d
    module procedure ncd_io_r4_0d
    module procedure ncd_io_r4_1d
    module procedure ncd_io_r4_2d
    module procedure ncd_io_r4_3d
    module procedure ncd_io_r4_4d
    module procedure ncd_io_r8_0d
    module procedure ncd_io_r8_1d
    module procedure ncd_io_r8_2d
    module procedure ncd_io_r8_3d
    module procedure ncd_io_r8_4d
    module procedure ncd_io_i4_0d
    module procedure ncd_io_i4_1d
    module procedure ncd_io_i4_2d
    module procedure ncd_io_i4_3d
    module procedure ncd_io_i4_4d

  end interface

 contains 

 subroutine ncd_io_char_0d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  character(len=*),  intent(inout) :: data
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer,optional,  intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_char_0d

 subroutine ncd_io_char_1d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  character(len=*),  intent(inout) :: data(:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer,optional,  intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_char_1d

! subroutine ncd_io_log_1d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)
!
! ! ARGUMENTS:
! !-------------
!  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
!  logical,           intent(inout) :: data(:)
!  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
!  character(len=*),  intent(in)    :: varname      ! variable name
!  logical,           intent(out)   :: readvar
!  integer,optional,  intent(out)   :: rc
!  integer, optional  , intent(in)    :: nt        ! time sample index
!  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file
!
!  ! LOCAL:
!
!  integer :: status
!
!  !-------------------------------------
!
!   if (flag == 'read') then
!      readvar = .false.
!     ! call ncid%get_var(varname, data, rc=status)
!     ! call MAPL_VarRead(ncid,varname,data,status)
!      if (status ==0) readvar = .true.
!   endif
!
! end subroutine ncd_io_log_1d

!----------------------------------------------------
 subroutine ncd_io_r4_0d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r4),          intent(inout) :: data
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer,optional,  intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r4_0d

!----------------------------------------------------
 subroutine ncd_io_r4_1d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r4),          intent(inout) :: data(:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer,optional,  intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r4_1d

  !-----------------------------------------------------------------------

 subroutine ncd_io_r4_2d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r4),          intent(inout) :: data(:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status                 
  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r4_2d

  !-----------------------------------------------------------------------

 subroutine ncd_io_r4_3d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r4),          intent(inout) :: data(:,:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r4_3d

  !-----------------------------------------------------------------------

 subroutine ncd_io_r4_4d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r4),          intent(inout) :: data(:,:,:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r4_4d

 !-----------------------------------------------------------------------

 subroutine ncd_io_r8_0d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r8),          intent(inout) :: data
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r8_0d

 !-----------------------------------------------------------------------

 subroutine ncd_io_r8_1d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r8),          intent(inout) :: data(:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r8_1d

 !-----------------------------------------------------------------------

 subroutine ncd_io_r8_2d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r8),          intent(inout) :: data(:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
     ! call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r8_2d

  !-----------------------------------------------------------------------


 subroutine ncd_io_r8_3d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r8),          intent(inout) :: data(:,:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r8_3d

  !-----------------------------------------------------------------------


 subroutine ncd_io_r8_4d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  real(r8),          intent(inout) :: data(:,:,:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_r8_4d

 !-----------------------------------------------------------------------
 subroutine ncd_io_i4_0d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  integer(i4),       intent(inout) :: data
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_i4_0d

 !-----------------------------------------------------------------------

 subroutine ncd_io_i4_1d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  integer(i4),       intent(inout) :: data(:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_i4_1d

 !-----------------------------------------------------------------------

 subroutine ncd_io_i4_2d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  integer(i4),       intent(inout) :: data(:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_i4_2d

  !-----------------------------------------------------------------------
 subroutine ncd_io_i4_3d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  integer(i4),       intent(inout) :: data(:,:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_i4_3d

  !-----------------------------------------------------------------------
 subroutine ncd_io_i4_4d ( varname, data, flag, ncid, readvar, rc, nt, posNOTonfile)

 ! ARGUMENTS:
 !-------------
  type(file_desc_t), intent(inout) :: ncid         ! netcdf file id
  integer(i4),       intent(inout) :: data(:,:,:,:)
  character(len=*),  intent(in)    :: flag         ! 'read' or 'write'
  character(len=*),  intent(in)    :: varname      ! variable name
  logical,           intent(out)   :: readvar
  integer, optional, intent(out)   :: rc
  integer, optional  , intent(in)    :: nt        ! time sample index
  logical            , optional, intent(in) :: posNOTonfile ! position is NOT on this file

  ! LOCAL:

  integer :: status

  !-------------------------------------

   if (flag == 'read') then
      readvar = .false.
      call ncid%get_var(varname, data, rc=status)
      !call MAPL_VarRead(ncid,varname,data,status)
      if (status ==0) readvar = .true.
   endif

 end subroutine ncd_io_i4_4d

  !-----------------------------------------------------------------------

  subroutine ncd_pio_openfile(file, fname, mode, rc)
    !
    ! !DESCRIPTION:
    ! Open a NetCDF PIO file
    !
    ! !ARGUMENTS:
    class(file_desc_t) , intent(inout) :: file   ! Output PIO file handle
    character(len=*)   , intent(in)    :: fname  ! Input filename to open
    integer            , intent(in)    :: mode   ! file mode
    integer, optional  , intent(out)   :: rc
    
    ! LOCAL:

    integer :: status

    !
    !-----------------------------------------------------------------------


    if (mode==0) then
       call file%open(trim(fname),pFIO_READ, rc=status)
    else
       _ASSERT(status==0, "Unrecognized netcdf opening mode")
    end if 

  end subroutine ncd_pio_openfile

  !-----------------------------------------------------------------------
  subroutine ncd_pio_closefile(file)
    !
    ! !DESCRIPTION:
    ! Close a NetCDF PIO file
    !
    ! !ARGUMENTS:
    class(file_desc_t), intent(inout) :: file   ! PIO file handle to close
    !-----------------------------------------------------------------------

    call file%close()

  end subroutine ncd_pio_closefile

end module ncdio_pio
