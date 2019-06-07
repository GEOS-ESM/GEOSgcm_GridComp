program newvegdyn
  implicit none
  
  real, pointer :: var(:)

  integer :: i, bpos, epos, status
  integer :: rsize
  character(256) :: str, f_in, f_out
  integer*4 :: ftell
  external  :: ftell

  integer*4 :: iargc
  external  :: iargc

! Begin

  if (iargc() /= 2) then
     call getarg(0,str)
     write(*,*) "Usage:",trim(str)," <old_style_vegdyn> <new_style_vegdyn>"
     call exit(2)
  end if

  call getarg(1,f_in)
  call getarg(2,f_out)

  open(unit=10, file=trim(f_in),  form='unformatted')
  open(unit=20, file=trim(f_out), form='unformatted')

  print *,'New Restart Format for File: ',trim(f_in)

  bpos=0
  read(10, err=200) ! skip to next record
  epos = ftell(10)          ! ending position of file pointer

  rsize = (epos-bpos)/4-2   ! record size (in 4 byte words; 
                            ! 2 is the number of fortran control words)
  allocate(var(rsize), stat=status)
  if (status  /= 0) then
     print *, 'Error: allocation ', rsize, ' failed!'
     call exit(11)
  end if

  read(10, err=200) ! skip to next record
  read(10, err=200) ! skip to next record
  read(10, err=200) ! skip to next record
! alltogather we  skip 4 record    
  read (10) var
  write(20) var
  deallocate(var)
  close(10)
  close(20)
  stop

200 print *,'Error reading file ',trim(f_in)
    call exit(11)

end
