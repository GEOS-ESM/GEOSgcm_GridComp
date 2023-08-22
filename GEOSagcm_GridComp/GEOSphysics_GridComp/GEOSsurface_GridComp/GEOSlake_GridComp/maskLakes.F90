#define _ASSERT(cond) if(.not. cond) then; print*, 'ERROR at ', __LINE__; stop; endif
program maskLakes
  implicit none

  integer, parameter :: HDR_SIZE=14
  integer :: header(HDR_SIZE)
  real :: hdr(HDR_SIZE)

  integer :: nt, nl
  integer :: ng
  character(len=128), parameter :: filename='mask.bin'
  character(len=128), parameter :: tilefile='tile.bin'
  real, allocatable :: lats(:), lons(:)
  integer, allocatable :: t(:)
  real, allocatable :: buffer(:)
  real, allocatable :: omask(:,:) ! mask on OSTIA grid
  integer :: im, jm
  integer :: i, j, n
  integer :: status
  integer :: unit
  real, allocatable :: mask(:)
  real :: x0, y0, dx, dy

! read OSTIA mask
  unit=10
  open(unit=unit, file=filename, form='unformatted')
  read(unit) hdr
  header = nint(hdr)
  im = header(13)
  jm = header(14)
  print *,'DEBUG: im/jm',im,jm

  allocate(omask(im,jm), stat=status)
  _ASSERT(status == 0)
  read(unit) omask
  close(unit)

! process tile file
  unit=11
  open(unit=unit, file=tilefile, form='unformatted')
  read(unit) nt
  print *,'DEBUG ntiles ',nt
  read(unit) ng
  print *,'DEBUG ngrids ',ng
  do i=1,ng
     ! 3 blank read for im, jm, gridname
     read(unit)
     read(unit)
     read(unit)
  end do
  allocate(buffer(nt), stat=status)
  _ASSERT(status == 0)
! get typetype
  read(unit) buffer
  t = nint(buffer)
  nl = count(t==19) ! number of lake tiles
  print *,'DEBUG lake points ',nl
  allocate(lons(nl), lats(nl), mask(nl), stat=status)
  _ASSERT(status == 0)
  mask = 0.0

! get X
  read(unit) buffer
  lons = pack(buffer, t==19)
! get Y
  read(unit) buffer
  lats = pack(buffer, t==19)
  close(unit)

  ! the OSTIA grid origin is at the pole edge, dateline edge
  x0 = -180.0
  y0 = -90.0
  dx = 360. / im
  dy = 180. / jm

!  convert tile lat/lon to Ostia grid indices
  do n=1,nl
     i = nint((lons(n)-x0)/dx - 0.5)
     j = nint((lats(n)-y0)/dy -0.5)+1
     i = mod(i+im,im)+1
     if (i<=0 .or. j<=0) print *,i,j,n,lats(n),lons(n)
     if (i>im .or. j>jm) print *,i,j,n,lats(n),lons(n)
     _ASSERT(i>0 .and. i<im)
     _ASSERT(j>0 .and. j<jm)

     if(omask(i,j) /= 0.0) then 
        mask(n) = 1.0
        print *,'DEBUG: masked lake at ',n,lons(n),lats(n)
     end if
  end do

! some sanity check/prints
  print *,'done, processed ',nl,' lake points, found ',count(mask/=0.0)

! ultimately, write mask
  unit=12
  open(unit=unit, file='lakemask.bin', form='unformatted')
  write(unit) mask
  close(unit)

! all done
end program maskLakes
