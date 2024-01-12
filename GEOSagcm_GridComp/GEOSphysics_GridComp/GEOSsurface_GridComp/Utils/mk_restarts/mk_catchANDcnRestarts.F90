#define I_AM_MAIN
#include "MAPL_Generic.h"

PROGRAM mk_catchANDcnRestarts

  use mpi
  use MAPL
  use ESMF
  use CatchmentRstMod
  use CatchmentCNRstMod

  implicit none
  
  character(len=:), allocatable :: out_bcsdir, out_dir, in_tilefile, out_tilefile, YYYYMMDDHHMM
  character(len=:), allocatable :: model, in_rstfile, out_rstfile
  character(len=:), allocatable :: out_File
  real     :: surflay, wemin_in, wemin_out
  integer  :: rc, status
  integer  :: myid, numprocs, mpierr
  class (CatchmentRst), allocatable :: catch
 
  call MPI_INIT(mpierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )

  call ESMF_Initialize(LogKindFlag=ESMF_LOGKIND_NONE)

  call process_cmd()

  if (index(model, 'catchcn') /=0 ) then
     catch = CatchmentCNRst(in_rstfile, model, yyyymmddhhmm, __RC__)
  else
     catch = CatchmentRst(in_rstfile, yyyymmddhhmm, __RC__)
  endif

  call catch%re_tile(in_tilefile, out_bcsdir, out_tilefile, surflay, __RC__)

  if (myid == 0) then
    call catch%add_bcs_to_rst(surflay, out_bcsdir, rc)
    call catch%re_scale(surflay, wemin_in, wemin_out, __RC__)
    call catch%write_nc4(out_file, __RC__)
  endif

  call ESMF_Finalize(endflag=ESMF_END_KEEPMPI)

  call MPI_FINALIZE(mpierr)
 
  contains
     ! process commands
     subroutine process_cmd()
       integer :: nxt
       character(len=256) :: arg
       nxt = 1
       call getarg(nxt,arg)

       do while(trim(arg) /= '')
          select case (trim(arg))
          case ('-h')
              call print_usage()
              call exit(0)
          case ('-out_bcs')
            nxt = nxt + 1
            call getarg(nxt,arg)
            out_bcsdir = trim(arg)
          case ('-time')
            nxt = nxt + 1
            call getarg(nxt,arg)
            YYYYMMDDHHMM = trim(arg)
          case ('-out_dir')
            nxt = nxt + 1
            call getarg(nxt,arg)
            out_dir = trim(arg)
          case ('-model')
            nxt = nxt + 1
            call getarg(nxt,arg)
            model = trim(arg)
          case ('-surflay')
            nxt = nxt + 1
            call getarg(nxt,arg)
            read(arg,*)  surflay
          case ('-in_tilefile')
            nxt = nxt + 1
            call getarg(nxt,arg)
            in_tilefile =  trim(arg)
          case ('-out_tilefile')
            nxt = nxt + 1
            call getarg(nxt,arg)
            out_tilefile =  trim(arg)
          case ('-in_rst')
            nxt = nxt + 1
            call getarg(nxt,arg)
            in_rstfile   =  trim(arg)
          case ('-out_rst')
            nxt = nxt + 1
            call getarg(nxt,arg)
            out_rstfile   =  trim(arg)
          case ('-in_wemin')
            nxt = nxt + 1
            call getarg(nxt,arg)
            read(arg,*)  wemin_in
          case ('-out_wemin')
            nxt = nxt + 1
            call getarg(nxt,arg)
            read(arg,*)  wemin_out
          case default
            print*, trim(arg)
            call print_usage()
            print*, "wong command line"
            call exit(1)
          end select
          nxt = nxt + 1
          call getarg(nxt,arg)
       end do
       out_file = out_dir //'/'//out_rstfile

     end subroutine

     subroutine print_usage()
        print *,'   '
        print *, 'This program can create catchment or catchmentCN restarts'
        print *, 'depending on the command line option "model"             '
        print *,'   '
        print *,'-out_bcs        : BC directory for output restart file'
        print *,'-time           : time for restart, format (yyyymmddhhmm)' 
        print *,'-out_dir        : directory for output restasrt file'
        print *,'-model          : model ( catch, catchcnclm40, catchcnclm45)'
        print *,'-surflay        : surflay value'
        print *,'-in_wemin       : wemin for input restart'
        print *,'-out_wemin      : wemin for output restart'
        print *,'-in_tilefile    : tile_file for input  restart'
        print *,'-out_tilefile   : tile_file for output restart, if none, it will search out_bcs'
        print *,'-in_rst         : input restart file name WITH path'
        print *,'-out_rst        : output restart file name WITHOUT path'
     end subroutine  
end program  
