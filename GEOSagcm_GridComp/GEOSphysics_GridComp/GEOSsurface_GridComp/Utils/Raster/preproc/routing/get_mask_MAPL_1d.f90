program main

use rwncfile
implicit none

integer,parameter :: nt=2592
integer,parameter :: nlon=72
integer,parameter :: nlat=36

real,allocatable,dimension(:,:) :: msk_MAPL
integer,allocatable,dimension(:) :: t2lati,t2loni,msk_tile

integer :: i

allocate(msk_MAPL(nlon,nlat))
allocate(t2lati(nt),t2loni(nt),msk_tile(nt))
call read_ncfile_real2d("inputs/MAPL_Tripolar.nc","mask",msk_MAPL,nlon,nlat)
open(77,file="inputs/TM0072xTM0036_tile_to_MAPL_lati.txt")
read(77,*)t2lati
open(77,file="inputs/TM0072xTM0036_tile_to_MAPL_loni.txt")
read(77,*)t2loni

do i=1,nt
  msk_tile(i)=int(msk_MAPL(t2loni(i),t2lati(i)))
enddo

open(88,file="outputs/mask_MAPL_1d_TM0072xTM0036.txt")
do i=1,nt
  write(88,*)msk_tile(i)
enddo

end