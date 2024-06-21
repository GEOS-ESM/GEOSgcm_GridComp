program main
  
  use routing_constant,only : nc,nl,ng
  implicit none
  
  integer,allocatable,dimension(:) :: msk,outid,mskall,final,finalall
  
  integer :: k,i,ntot,ns

!  Get a list of sink catchment IDs (both land and Greenland)
! nc = number of the land catchments (excluding Greenland), 291284 in our case.
! ng = number of the Greenland catchments, 525 in our case.
! nl = number of outlets to ocean in land (not including Greenland outlets) 
! ns = number of the total outlets (including Greenland); note that all the Greenland catchments are sink catchments (ie, no downstream catchment).
! ntot = number of the total catchments (including Greenland)
! outid = a list of indices for sink catchments (to ocean)  

  ntot=nc+ng
  ns=nl+ng
  allocate(msk(nc),outid(ns),mskall(ntot),final(nc),finalall(ntot))
  open(77,file="outputs/Pfaf_msk.txt") !!!
  read(77,*)msk
  k=0
  do i=1,nc
     if(msk(i).eq.2)then ! msk=2 is for a catchment that drains directly into the ocean
        k=k+1  	
        outid(k)=i
     end if
  end do

  ! Add Greenland catchments to list of outlet catchments; write outlet catchments to file
  do i=k+1,ns
     outid(i)=nc+i-k
  end do
  open(88,file="outputs/outlet_catchindex.txt")
  do i=1,ns
     write(88,*)outid(i)
  enddo

! Append msk values for Greenland to original msk array; write out  
  mskall(1:nc)=msk
  mskall(nc+1:)=2
  open(88,file="outputs/Pfaf_msk_all.txt")
  do i=1,ntot
     write(88,*)mskall(i)
  enddo

! Append final values for Greenland to original “final” array; write out  
  open(77,file="outputs/Pfaf_finalID.txt")
  read(77,*)final
  finalall(1:nc)=final
  do i=nc+1,ntot
     finalall(i)=i
  end do
  open(88,file="outputs/Pfaf_finalID_all.txt")
  do i=1,ntot
     write(88,*)finalall(i)
  enddo
  
end program main
