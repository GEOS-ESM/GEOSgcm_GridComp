program main
  
  use routing_constant,only : nall,nlon,nlat
  implicit none
  
  character(len=100)  :: var1="outlet_sinky_allcat"
  character(len=100)  :: var2="outlet_sinkx_allcat"
  character(len=100)  :: map="Pfafstetter_Greenland_real"
  
  real*8,allocatable  :: lon(:),lat(:)
  integer,allocatable :: catchind(:,:)
  integer,allocatable :: lons(:,:),lats(:,:)
  integer,allocatable :: data_Pfaf(:)
  
  integer :: i,j,xi,yi,id
  
  
  allocate(catchind(nlon,nlat),lons(nlon,nlat),lats(nlon,nlat))
  allocate(lon(nlon),lat(nlat))
  
  open(30,file="outputs/"//trim(map),form="unformatted")
  do j = 1,nlat
     read (30) catchind(:,j)
  end do
  
  allocate(data_Pfaf(nall))
  
  open(77,file="outputs/"//trim(var1)//".txt")
  read(77,*)data_Pfaf
  lats=-999
  do xi=1,nlon
     do yi=1,nlat
        if(catchind(xi,yi)>=1.and.catchind(xi,yi)<=nall)then
           id=catchind(xi,yi)
           lats(xi,yi)=data_Pfaf(id)
        endif
     enddo
  enddo
  
  open(77,file="outputs/"//trim(var2)//".txt")
  read(77,*)data_Pfaf
  lons=-999
  do xi=1,nlon
     do yi=1,nlat
        if(catchind(xi,yi)>=1.and.catchind(xi,yi)<=nall)then
           id=catchind(xi,yi)
           lons(xi,yi)=data_Pfaf(id)
        endif
     enddo
  enddo
  
  open(30,file="Outlet_latlon.43200x21600",form="unformatted")
  do j = 1, nlat
     write (30) lons(:,j)
     write (30) lats(:,j)
  end do

end program main
