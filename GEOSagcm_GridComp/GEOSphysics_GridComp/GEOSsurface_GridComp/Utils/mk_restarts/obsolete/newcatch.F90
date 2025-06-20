#define VERIFY_(A) if(A /=0)then;print *,'ERROR code',A,'at',__LINE__;call exit(3);endif

program newcatch
  implicit none

#ifndef __GFORTRAN__
  integer*4 :: iargc
  external  :: iargc
  integer*8 :: ftell
  external  :: ftell
#endif
  character(256) :: str, f_in, f_out

  integer :: m
  integer :: status
  integer*8 :: bpos, epos, rsize
  real, allocatable :: a(:)

! Begin

  if (iargc() /= 2) then
     call getarg(0,str)
     write(*,*) "Usage:",trim(str)," <old catch restart> <new catch restart>"
     call exit(2)
  end if

  call getarg(1,f_in)

  open(unit=10, file=trim(f_in),  form='unformatted')

! Count the records in the files
! ------------------------------
! Valid numbers are:
! 61 - old catch_internal_restart
! 57 - old catch_internal_restart
  m=0
  do while(.true.)
     read(10, end=50, err=200) ! skip to next record
     m = m+1
  end do
50 continue
  rewind(10)

  if (m == 57) then
     print *,'WARNING: this file contains ', m, ' records and appears to have been already convered'
     print *,'Refuse to convert!'
     print *,'Exiting ...'
     call exit(1)
  else if (m /= 61) then
     print *,'ERROR: this file contains ',m, &
          ' records and does not appear to be a valid catchment internal restart'
     print *,'Exiting ...'
     call exit(2)
  end if

! Open the output file
! --------------------
  call getarg(2,f_out)
  open(unit=20, file=trim(f_out), form='unformatted')

  m=0
  bpos=0
  do while(.true.)
     m = m+1
     read(10, end=100, err=200) ! skip to next record
     epos = ftell(10)           ! ending position of file pointer
     backspace(10)
     
     rsize = (epos-bpos)/4-2    ! record size (in 4 byte words; 
     bpos = epos
     allocate(a(rsize), stat=status)
     VERIFY_(status)
     read (10) a
     if (m < 57 .or. m > 60) then
        print *,'Writing record ',m
        write(20) a
     else
        print *,'Skipping record ',m
     end if
     deallocate(a)
  end do
100 continue
  close(10)
  close(20)
  stop

! If we are here something must have gone wrong
200 VERIFY_(200)

end program newcatch

