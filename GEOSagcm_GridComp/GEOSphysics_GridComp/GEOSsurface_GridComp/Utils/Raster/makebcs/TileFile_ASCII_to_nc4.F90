!
! Utility program that converts ASCII-formatted *.til file and catchment.def file into a single nc4 file
!
! Usage TileFile_ASCII-to-nc4.x tile_file catchmentdef_file
!
! wjiang, rreichle, 29 Nov 2024

program TileRst_to_nc4
  use, intrinsic :: iso_fortran_env, only: REAL64 
  use LogRectRasterizeMod,           only: MAPL_UNDEF_R8
  use rmTinyCatchParaMod 
  implicit none
  
  character(512)                 :: arg
  integer                        :: nc, nr
  
  character(:),      allocatable :: tile_file
  character(:),      allocatable :: rst_file
  character(:),      allocatable :: gName
  character(:),  allocatable     :: filenameNC4
  integer                        :: n
  logical                        :: file_exists

  !
  ! usage: TileRst_to_nc4.x nx ny gName
  ! ----------------------------------------------------------------------
  !
  ! process command-line arguments

  CALL get_command_argument(1, arg)
  read(arg,'(i6)') nc 
  CALL get_command_argument(2, arg)
  read(arg,'(i6)') nr 
  CALL get_command_argument(3, arg)
  gName = trim(arg)

  tile_file   = 'til/'//gName
  rst_file    = 'rst/'//gName//'.rst'

  filenameNC4 = tile_file //'.nc4'

  inquire(file=filenameNC4, exist=file_exists)
  if (.not. file_exists) then
     call supplemental_tile_attributes(nc,nr, tile_file, rst_file, write_catch=.false.) 
  endif
end program

! ======================= EOF ====================================================
