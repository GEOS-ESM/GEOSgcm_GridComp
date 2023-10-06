program main

implicit none

integer,parameter :: nc=291284
integer,parameter :: ns=22612
integer,parameter :: ng=525

integer,allocatable,dimension(:) :: msk,outid,mskall,final,finalall

integer :: k,i,ntot

ntot=nc+ng
allocate(msk(nc),outid(ns),mskall(ntot),final(nc),finalall(ntot))
open(77,file="inputs/Pfaf_msk.txt")
read(77,*)msk
k=0
do i=1,nc
  if(msk(i).eq.2)then
    k=k+1  	
    outid(k)=i
  end if
end do
do i=k+1,ns
  outid(i)=nc+i-k
end do
open(88,file="outputs/outlet_catchindex.txt")
do i=1,ns
  write(88,*)outid(i)
enddo

mskall(1:nc)=msk
mskall(nc+1:)=2
open(88,file="outputs/Pfaf_msk_all.txt")
do i=1,ntot
  write(88,*)mskall(i)
enddo

open(77,file="inputs/Pfaf_finalID.txt")
read(77,*)final
finalall(1:nc)=final
do i=nc+1,ntot
  finalall(i)=i
end do
open(88,file="outputs/Pfaf_finalID_all.txt")
do i=1,ntot
  write(88,*)finalall(i)
enddo



end