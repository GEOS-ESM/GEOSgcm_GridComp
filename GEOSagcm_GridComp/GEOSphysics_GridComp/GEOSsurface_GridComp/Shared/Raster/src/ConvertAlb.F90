
integer, parameter :: nx=8640, ny=4320


real, allocatable :: alb_in(:,:), albo(:,:)
integer*1, allocatable :: alb_out(:,:)
character*128 :: ifile, ofile


allocate(alb_in(nx,ny))
allocate(albo(nx,ny))
allocate(alb_out(nx,ny))

call getarg(1,ifile)

ofile = "data/Attic/foo"!trim(ifile)//'int1'

open(10, file="data/Attic/"//ifile, status='old', form='unformatted',convert='big_endian')
open(11, file=ofile, status='unknown', form='unformatted',convert='big_endian')

do m=1,12
   read(10) alb_in
   where(alb_in>1.) alb_in=0.
   alb_out=nint((alb_in-0.50)*255.)
   write(11) alb_out

   albo=(alb_out+128)/255.0_8
   
!!$   do j=1,ny
!!$      do i=1,nx
!!$         if(albo(i,j)>0.9) &
!!$         print *, i, j, albo(i,j),alb_in(i,j),alb_out(i,j)
!!$      end do
!!$   end do

   print *, "month ",m
   print *, minval(alb_out), maxval(alb_out)
   print *, minval(albo)   , maxval(albo)
   print *, minval(alb_in) , maxval(alb_in)
enddo

stop
end
