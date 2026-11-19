PROGRAM read_tile_param

implicit none
character*800 :: cdum
integer :: n,fvalue

open (10,file ='/home/koster/debug3/tile_parameters.dat', form='formatted',&
     status='old',action='read')

do n=1,1000

read(10,'(a)') cdum
read (cdum,*) fvalue

if (fvalue == 1) print *, trim(cdum)


end do

close (10,status='keep')


END PROGRAM read_tile_param
