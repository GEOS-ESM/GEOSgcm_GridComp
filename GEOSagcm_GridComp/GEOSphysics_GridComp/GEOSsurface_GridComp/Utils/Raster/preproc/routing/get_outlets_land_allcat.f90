program main
  
  use routing_constant,only : nall,ns
  implicit none
  
  integer, allocatable, dimension(:) :: id_final,id_outlet,msk 
  integer,allocatable,dimension(:)   :: lati_outlet,loni_outlet
  integer,allocatable,dimension(:)   :: lati_full,loni_full
  
  integer :: i,j
  
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
     !if(mod(i,1000)==0) print *,i
     if(msk(id_final(i)).eq.2)then
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
