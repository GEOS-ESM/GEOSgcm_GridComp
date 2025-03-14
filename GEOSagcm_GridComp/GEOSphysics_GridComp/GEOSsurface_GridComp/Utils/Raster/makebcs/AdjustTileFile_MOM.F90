program AdjustTileFile_MOM
  !
  !Adjust MOM tile file
  ! Input  : file1, original tile file without changing raster 
  !          file2, tile file generated after 
  ! Output : file3, copy land tiles appearing in both file1 and file2 from file1 to file3,
  !                 the othe tiles from file2 to file3 
  !

  implicit none


  character(len=64) ::  Usage = "AdjustTileFile_MOM cf_tile mom_tile new_mom_tile"

  character(len=16) :: Iam = "AdjustTileFile"
  character(len=256)::  arg, line1, line2, line3
  character(:), allocatable :: file1, file2, file3
  integer :: unit1, unit2, unit3, ntile, i, i1, j1, i2, j2, pfaf1, pfaf2, tmpint
  logical :: match
  real :: tmpr
  integer :: typ, itile, total
  call get_command_argument(1,arg)
  file1 = 'til/'//trim(arg)//'.til'
  call get_command_argument(2,arg)
  file2 = 'til/'//trim(arg)//'.til'
  call get_command_argument(3,arg)
  file3 = 'til/'//trim(arg)//'.til'

  open (newunit=unit1, file=file1, form='formatted', action='read')
  open (newunit=unit2, file=file2, form='formatted', action='read')
  open (newunit=unit3, file=file3, form='formatted', action='write')

  do i = 1, 8
    read(unit1,'(A)') line1
    read(unit2,'(A)') line2
    if (i == 1) read(line2,*) ntile
    write(unit3, '(A)') trim(line2)
  enddo

  itile = 1
  total = 0
  do 
    read(unit1,'(A)') line1
    read(unit2,'(A)') line2
    read(line1, *) tmpint, tmpr, tmpr, tmpr, i1, j1, tmpr, pfaf1
    read(line2, *) tmpint, tmpr, tmpr, tmpr, i2, j2, tmpr, pfaf2
    if (tmpint/= 100) exit
    match = ( i1 == i2 .and. j1 == j2 .and. pfaf1 == pfaf2 )
    do while ( .not. match)
       read(unit1,'(A)') line1
       read(line1, *) tmpint, tmpr, tmpr, tmpr, i1, j1, tmpr, pfaf1
       match = ( i1 == i2 .and. j1 == j2 .and. pfaf1 == pfaf2 )
    end do
    write(unit3, '(A)') trim(line1)
    itile = itile + 1
    total = total + 1
  enddo
  write(unit3, '(A)') trim(line2)
  total = total + 1
  do i = itile+1, ntile
    read(unit2,'(A)') line2
    write(unit3, '(A)') trim(line2)
    total = total + 1
  enddo
  close(unit1)
  close(unit2)
  close(unit3)
  if ( total /= ntile) then
    print *, "The tile numeber is wrong, cannot pass the sanity check"
    stop
  endif
end program
