program create_example
use netcdf
use, intrinsic :: iso_fortran_env, only: REAL64
implicit none


character(len=512) :: fin,fout,str,fncar
integer :: im_world,jm_world
integer :: varid, lonid, latid
integer :: i, j, nc,xid,yid
integer :: zid,gwdid,trbid
integer :: ncid,rc
integer :: dimids(2)
integer :: status

integer :: nargs

logical :: doNcar,doGEOS

integer :: ntiles
real, allocatable :: z1d(:)
real(REAL64), allocatable :: xdim(:),ydim(:)
real, allocatable :: a(:,:)
logical :: isCube

nargs = command_argument_count()

doNCAR=.false.
doGEOS=.false.
isCube = .true.
do i=1,nargs
     call get_command_argument(i,str)
     select case(trim(str))
     case ('-i','--input')
         call get_command_argument(i+1,fin)
     case ('-o','--output')
         call get_command_argument(i+1,fout)
         doGEOS=.true.
     case ('--im')
          call get_command_argument(i+1,str)
          read(str,'(I10)')im_world
     case ('--jm')
          call get_command_argument(i+1,str)
          read(str,'(I10)')jm_world
          isCube = .false.
     case ('--ncar')
          call get_command_argument(i+1,fncar)
          doNCAR=.true.
     end select
enddo

if (isCube) jm_world = im_world*6

allocate(a(im_world,jm_world))
open(file=fin,unit=21,form='unformatted')
read(21)a
close(21)

if (doGEOS) then

   call check( nf90_create(fout, NF90_NETCDF4,ncid),"error")
   call check( nf90_def_dim(ncid,"Xdim",im_world,lonid),"error")
   call check( nf90_def_var(ncid,"Xdim",NF90_DOUBLE,(/lonid/),xid),"error")
   call check( nf90_put_att(ncid,xid,"units","degrees_east"),"error")
   call check( nf90_def_dim(ncid,"Ydim",jm_world,latid),"error")
   call check( nf90_def_var(ncid,"Ydim",NF90_DOUBLE,(/latid/),yid),"error")
   call check( nf90_put_att(ncid,yid,"units","degrees_north"),"error")
   call check( nf90_def_var(ncid,"z",NF90_FLOAT,(/lonid,latid/),varid),"error")
   call check( nf90_put_att(ncid,varid,"units","m"),"error")
   call check( nf90_put_att(ncid,varid,"long_name","height above sea level"),"error")

   call check( nf90_enddef(ncid),"error")

   allocate(xdim(im_world),ydim(jm_world))
   do i=1,im_world
   xdim(i)=i
   enddo
   do j=1,jm_world
   ydim(j)=j
   enddo

   call check(nf90_put_var(ncid,xid,xdim),"error")
   call check(nf90_put_var(ncid,yid,ydim),"error")
   call check(nf90_put_var(ncid,varid,a),"error")
   call check(nf90_close(ncid),"error")

end if

if (doNCAR) then

   ntiles=im_world*jm_world
   allocate(z1d(ntiles))
   call check( nf90_create(fncar, NF90_NETCDF4,ncid),"error")
   call check( nf90_def_dim(ncid,"ncol",ntiles,xid),"error")
   call check( nf90_def_var(ncid,"PHIS",NF90_DOUBLE,(/xid/),varid),"error")
   call check( nf90_put_att(ncid,varid,"long_name","height"),"error")
   call check( nf90_put_att(ncid,varid,"units","m"),"error")
   call check( nf90_enddef(ncid),"error")

   nc=0
   do j=1,jm_world
   do i=1,im_world
      nc=nc+1
      z1d(nc)=a(i,j)
   enddo
   enddo

   call check(nf90_put_var(ncid,varid,z1d),"error")

end if

contains

subroutine check(status,loc)

      integer, intent ( in) :: status
      character(len=*), intent ( in) :: loc

      if(status /= NF90_noerr) then
         write (*,*) "Error at ", loc
         write (*,*) nf90_strerror(status)
         stop "Stopped"
      end if

end subroutine check

end program create_example

