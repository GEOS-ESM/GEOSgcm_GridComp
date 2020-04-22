program mk_CiceRestart

! This utility will work with CICE saltwater internal restart from Fortuna-2_4.
! For CICE import restart use mk_LakeLandiceSaltRestarts.

  use MAPL_ConstantsMod,only: MAPL_PI,  MAPL_radius
  use MAPL_HashMod
  use mk_restarts_getidsMod, only: GetIDs,ReadTileFile_IntLatLon

  implicit none

  character*128 :: Usage="mk_LakeLandiceSaltRestarts OutTileFile InTileFile InRestart"
  character*128 :: OutTileFile
  character*128 :: InTileFile
  character*128 :: InRestart
  character*128 :: arg

  integer :: i, iargc, n,j,ntiles,k
  integer, pointer  :: Lono(:), Lato(:), Id(:), Pf(:)
  integer, pointer  :: Loni(:), Lati(:)
  real*4, allocatable :: var4(:)
  real*8, allocatable :: var8(:)

  integer, parameter ::  zoom=2

!---------------------------------------------------------------------------

  I = iargc()

  if(I /= 3) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     call exit(1)
  end if

  call getarg(1,OutTileFile)
  call getarg(2,InTileFile)
  call getarg(3,InRestart)

! Read Output Tile File .til file
! to get the index into the pfafsttater table

  call ReadTileFile_IntLatLon(OutTileFile,Pf,Id,lono,lato,zoom,0)
  deallocate(Pf,Id)

  call ReadTileFile_IntLatLon(InTileFile ,Pf,Id,loni,lati,zoom,0)
  deallocate(Pf,Id)

  nullify(Pf)
  nullify(Id)

  ntiles = size(lono)

  i = index(InRestart,'/',back=.true.)

  open(unit=40,FILE="OutData/"//trim(InRestart(i+1:)),form='unformatted',&
       status='unknown',convert='little_endian')

  open(unit=50,FILE=InRestart,form='unformatted',&
       status='old',convert='little_endian')

  allocate(var4(size(loni)),var8(size(loni)))

  do n=1,124
     read (50)
  end do

  read (50) var4(:)

  rewind 50

  allocate(Id (ntiles))

  call GetIds(loni,lati,lono,lato,zoom,Id)

  do n=1,18
     read (50) var4(:)
     write(40)(var4(id(i)),i=1,ntiles)
  end do

  do n=19,74
     read (50) var8(:)
     write(40)(var8(id(i)),i=1,ntiles)
  end do

  do n=75,125
     read (50) var4(:)
     write(40)(var4(id(i)),i=1,ntiles)
  end do

  deallocate(var4,var8)

  close(40)
  close(50)

end program mk_CiceRestart

