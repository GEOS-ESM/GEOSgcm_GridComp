implicit none

integer, parameter :: nx=8640, ny=4320
integer, allocatable  :: lats(:,:), lons(:,:)
integer i,j

allocate(lats(nx,ny), lons(nx,ny))

open(10,file="/land/koster/lonraster.dat",status="old", &
     form="formatted")
open(20,file="/land/koster/latraster.dat",status="old", &
     form="formatted")

do i=1,ny
   read(10,'(100i5)') lons(:,i)
   read(20,'(100i5)') lats(:,i)
end do

close(10)
close(20)

open(30,file="../data/Outlet_latlon.dat",status="unknown", &
     form="unformatted")

write(30) int(lons,kind=2)
write(30) int(lats,kind=2)


print *, maxval(lats), minval(lats)
print *, maxval(lons), minval(lons)


!!$do i=1,4320
!!$print *, i,lats(1,i),lons(1,i)
!!$enddo

close(30)

end
