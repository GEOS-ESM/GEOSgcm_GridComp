#define I_AM_MAIN
#include "MAPL_Generic.h"

PROGRAM mk_CatchmentsRestarts

  use mpi
  use MAPL
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
  type (scale_var) :: old
 
  call MPI_INIT(mpierr)
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )

  call process_cmd()

  if (index(model, 'catchcn') /=0 ) then
     catch = CatchmentCNRst(in_rstfile, model, yyyymmddhhmm, __RC__)
  else
     catch = CatchmentRst(in_rstfile, yyyymmddhhmm, __RC__)
  endif

  call catch%re_tile(in_tilefile, out_bcsdir, out_tilefile, surflay, __RC__)

  call catch%set_scale_var(old)

  call catch%add_bcs_to_rst(surflay, out_bcsdir, rc)

  call catch%re_scale(surflay, wemin_in, wemin_out, old, __RC__)

  call catch%write_nc4(out_file, __RC__)

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
          case ('-bcs_out')
            call getarg(nxt,arg)
            out_bcsdir = trim(arg)
          case ('-time')
            call getarg(nxt,arg)
            YYYYMMDDHHMM = trim(arg)
          case ('-dir_out')
            call getarg(nxt,arg)
            out_dir = trim(arg)
          case ('-model')
            call getarg(nxt,arg)
            model = trim(arg)
          case ('-surflay')
            call getarg(nxt,arg)
            read(arg,*)  surflay
          case ('-tile_in')
            call getarg(nxt,arg)
            in_tilefile =  trim(arg)
          case ('-tile_out')
            call getarg(nxt,arg)
            out_tilefile =  trim(arg)
          case ('-rst_in')
            call getarg(nxt,arg)
            in_rstfile   =  trim(arg)
          case ('-rst_out')
            call getarg(nxt,arg)
            out_rstfile   =  trim(arg)
          case ('-wemin_in')
            call getarg(nxt,arg)
            read(arg,*)  wemin_in
          case ('-wemin_out')
            call getarg(nxt,arg)
            read(arg,*)  wemin_out
          case default
            call print_usage()
            print*, "wong command line"
            call exit(1)
          end select
          nxt = nxt + 1
          call getarg(nxt,arg)
       end do
       if (index(model, 'catchcn') /=0 ) then
         if((INDEX(out_bcsdir, 'NL') == 0).AND.(INDEX(out_bcsdir, 'OutData') == 0)) then
           print *,'Land BCs in : ',trim(out_bcsdir)
           print *,'do not support ',trim (model)
           stop
         endif
       endif
       out_file = out_dir //'/'//out_rstfile

     end subroutine

     subroutine print_usage()
        print *,'   '
        print *,'-bcs_out    : BC directory for output restart file'
        print *,'-time       : time for restart, format (yyyymmddhhmm)' 
        print *,'-dir_out    : directory for output restasrt file'
        print *,'-model      : model ( catch, catchcnclm40, catchcnclm45)'
        print *,'-surflay    : surflay value'
        print *,'-wemin_in   : wemin for input restart'
        print *,'-wemin_out  : wemin for output restart'
        print *,'-tile_in    : tile_file for input  restart'
        print *,'-tile_out   : tile_file for output restart, if none, it will search out_bcs'
        print *,'-rst_in     : input restart file name'
        print *,'-rst_out    : output restart file name'
     end subroutine  
end program  
