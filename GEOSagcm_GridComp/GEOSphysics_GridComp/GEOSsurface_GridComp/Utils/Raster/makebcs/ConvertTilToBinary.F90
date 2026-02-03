#define I_AM_MAIN
#include "MAPL_ErrLog.h"

program ConvertTilToBinary

  use MAPL_ExceptionHandling
  use MAPL_Constants

! Convert text .til files to binary .til.bin for faster I/O
! Usage: ConvertTilToBinary input.til output.til.bin

  implicit none

  integer, parameter     :: TILUNIT_IN  = 20
  integer, parameter     :: TILUNIT_OUT = 21
  character*256          :: Iam = "ConvertTilToBinary"

  integer                :: command_argument_count
  integer                :: STATUS
  integer                :: i, j, k, ip, nf, nx, ny, num_grids
  integer                :: argl

  real(REAL64)           :: val_real
  character*256          :: GridName
  character*256          :: InputFile, OutputFile
  character*256          :: arg

  real(REAL64), allocatable :: Table(:,:)

! Get command line arguments

    if (command_argument_count() /= 2) then
       print *, "Usage: ConvertTilToBinary input.til output.til.bin"
       call exit(1)
    end if

    call get_command_argument(1, InputFile)
    call get_command_argument(2, OutputFile)

    print *, "Converting: ", trim(InputFile)
    print *, "        to: ", trim(OutputFile)

! Open input (text) and output (binary) files

    open (TILUNIT_IN, file=trim(InputFile), form='formatted', status='old')
    open (TILUNIT_OUT, file=trim(OutputFile), form='unformatted', &
          convert='little_endian', status='replace')

! Read and write header: ip, nf, nx, ny

    read(TILUNIT_IN,*) ip, nf, nx, ny
    write(TILUNIT_OUT) ip, nf, nx, ny
    print *, "Header: ip=", ip, " nf=", nf, " nx=", nx, " ny=", ny

! Read and write number of grids

    read(TILUNIT_IN,*) num_grids
    write(TILUNIT_OUT) num_grids
    print *, "Number of grids: ", num_grids

! Read and write grid metadata headers

    do j=1,num_grids
       ! Read grid name
       read(TILUNIT_IN,*) GridName
       write(TILUNIT_OUT) trim(GridName)

       ! Read nx, ny for this grid
       read(TILUNIT_IN,*) nx
       read(TILUNIT_IN,*) ny
       write(TILUNIT_OUT) nx, ny

       print *, "Grid ", j, ": ", trim(GridName), " nx=", nx, " ny=", ny
    end do

! Allocate space for one row of tile data

    allocate(Table(1:8, 1:1), stat=STATUS)
    VERIFY_(STATUS)

! Read and write tile data: ip rows of (2 reals, 2 skip, 6 reals)

    print *, "Reading ", ip, " tile records..."

    do k=1,ip
       ! Read: typ(real), tarea(real), lon(real), lat(real, skip), then 6 more reals
       read(TILUNIT_IN,*) Table(1:2,1), val_real, val_real, Table(3:8,1)

       ! Write all 8 values as unformatted binary
       write(TILUNIT_OUT) Table(1:8,1)

       if (mod(k, max(1, ip/10)) == 0) then
          print *, "  Wrote ", k, " records..."
       end if
    end do

    print *, "Conversion complete."
    print *, "Original file: ", trim(InputFile)
    print *, "Binary file:   ", trim(OutputFile)

! Close files

    close(TILUNIT_IN)
    close(TILUNIT_OUT)

    deallocate(Table)

    call exit(0)

end program ConvertTilToBinary

!===================================================================
