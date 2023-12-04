#define VERIFY_(A) if(A /=0)then;print *,'ERROR code',A,'at',__LINE__;call exit(3);endif

program checkVegDyn
  implicit none

#ifndef __GFORTRAN__
  integer*4 :: iargc
  external  :: iargc
  integer   :: ftell
  external  :: ftell
#endif
  character(256) :: str, f_in, f_out

  integer :: m, n
  integer :: status
  integer :: bpos, epos, nt
  integer, parameter :: unit=10
  real, allocatable :: a(:)
  integer, allocatable :: veg(:)
  integer :: minVegType
  integer :: maxVegType

! Begin

  if (iargc() /= 2) then
     call getarg(0,str)
     write(*,*) "Usage:",trim(str)," <old vegDynData>"," <new vegDynData>"
     call exit(2)
  end if

  call getarg(1,f_in)
  call getarg(2,f_out)

  open(unit=unit, file=trim(f_in),  form='unformatted')

! count the records
  m=0
  do while(.true.)
     read(unit, end=50, err=200) ! skip to next record
     m = m+1
  end do
50 continue
  if (m == 1) then
     print *, 'File ', trim(f_in), 'contains only only record. Exiting ...'
     goto 100
  end if

  rewind(unit)

  open(unit=20, file=trim(f_out), form='unformatted')

! determine number of tiles by the size of the first record

  bpos=0
  read(unit, err=200) ! skip to next record
  epos = ftell(unit)           ! ending position of file pointer
  nt = (epos-bpos)/4-2    ! record size (in 4 byte words; 
  rewind(unit)

  allocate(a(nt), stat=status)
  VERIFY_(status)

! Read and copy first record
  read (unit) a
  write(20) a

  close(20)

! clean up
100 continue
  close(unit)
  stop

! If we are here, something must have gone wrong
200 VERIFY_(200)

end program checkVegDyn

