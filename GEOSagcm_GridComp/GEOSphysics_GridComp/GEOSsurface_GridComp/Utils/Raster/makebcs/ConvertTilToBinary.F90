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

  real(kind=kind(1.0d0)) :: val1, val2, val3, val4, val5, val6, val7, val8
  character*256          :: GridName
  character*256          :: InputFile, OutputFile
  character*256          :: arg

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

! Read and write grid metadata for first grid only (matching .til file structure)

    ! Read first grid name
    read(TILUNIT_IN,*) GridName
    write(TILUNIT_OUT) trim(GridName)
    
    ! Read nx, ny for first grid
    read(TILUNIT_IN,*) nx
    read(TILUNIT_IN,*) ny
    write(TILUNIT_OUT) nx, ny
    
    print *, "First grid: ", trim(GridName), " nx=", nx, " ny=", ny

    ! Skip remaining grid headers (they're in the text file but we only need first)
    do j=2,num_grids
       read(TILUNIT_IN,*)
       read(TILUNIT_IN,*)
       read(TILUNIT_IN,*)
    end do

! Read and write tile data: ip rows of 8 values

    print *, "Reading ", ip, " tile records..."
    
    do k=1,ip
       ! Read: typ, tarea, lon, lat (skip), then 4 more values
       read(TILUNIT_IN,*) val1, val2, val3, val4, val5, val6, val7, val8
       
       ! Write all 8 values as unformatted binary
       write(TILUNIT_OUT) val1, val2, val3, val4, val5, val6, val7, val8
       
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

    call exit(0)

end program ConvertTilToBinary

!===================================================================
