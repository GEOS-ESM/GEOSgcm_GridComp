program main
  
  use routing_constant,only : nc,nl,ng
  implicit none
  
  integer, allocatable, dimension(:) :: id_final,id_outlet,msk 
  integer,allocatable,dimension(:)   :: lati_outlet,loni_outlet
  integer,allocatable,dimension(:)   :: lati_full,loni_full
  
  integer :: i,j,nall,ns

! Assign outlet locations to all upstream catchments to create a 1d list showing the final x and y indexes for each catchment.
!  id_final = index of the final sink catchment for each catchment 
! id_outlet =catchment index for each outlet
! msk = 1 = has downstream catchment; 2 = drains to ocean; 3 = drains to inland lake
! lati_outlet = lat index for each outlet
! loni_outlet =lon index for each outlet
! lati_full = lat index (on a 30s map) of its final sink point for each catchment
! loni_full = lon index (on a 30s map) of its final sink point for each catchment
  
  nall=nc+ng
  ns=nl+ng
  allocate(id_final(nall),id_outlet(ns),msk(nall),&
       lati_outlet(ns),loni_outlet(ns),lati_full(nall),loni_full(nall))
  
  open(77,file="outputs/Pfaf_finalID_all.txt")
  read(77,*)id_final
  open(77,file="outputs/outlet_catchindex.txt")
  read(77,*)id_outlet
  open(77,file="outputs/outlet_sinky.txt")
  read(77,*)lati_outlet
  open(77,file="outputs/outlet_sinkx.txt")
  read(77,*)loni_outlet
  open(77,file="outputs/Pfaf_msk_all.txt")
  read(77,*)msk
  
  lati_full=-999
  loni_full=-999
  
  do i=1,nall
! For each catchment, check to see if its final sink catchment actually drains to the ocean.    
     if(msk(id_final(i)).eq.2)then
    ! For the catchment being tested, loop over the indices of all sink catchments and find the one that matches its own sink catchment.
    ! Assign lat/longs of outlet point to the catchment being tested.  
        do j=1,ns
           if(id_outlet(j).eq.id_final(i))then 
              lati_full(i)=lati_outlet(j)
              loni_full(i)=loni_outlet(j)
           end if
        enddo
     else if(msk(id_final(i)).eq.3)then
        lati_full(i)=-999
        loni_full(i)=-999
     endif
  end do
  
  open(88,file="outputs/outlet_sinky_allcat.txt")
  do i=1,nall
     write(88,*)lati_full(i)
  enddo
  open(88,file="outputs/outlet_sinkx_allcat.txt")
  do i=1,nall
     write(88,*)loni_full(i)
  enddo
  
end program main
