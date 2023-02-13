program Runoff

  use mapl_hashmod
  use mapl_sortmod
  use netcdf

  implicit none

  integer                :: nx, ny, pf
  integer, allocatable   :: lats(:,:), lons(:,:)
  integer, pointer       :: rst(:,:), SortArr(:,:), key(:)
  integer, pointer       :: srctile(:),  srcweight(:), dstweight(:), dsttile(:)
  real, allocatable      :: SrcFraction(:), area(:), in(:), out(:)
  integer                :: i,j,k,l, Hash, HashC, ii,jj,kk
  integer                :: type, np,lnd, is,ie,ww
  integer                :: numtrans,  numclosed
  integer                :: status
  character*100          :: file, fileT, fileR, fileO, fileB, fileBB

  character*400          :: fileLL
  character*400          :: MAKE_BCS_INPUT_DIR


  character*5            :: C_NX, C_NY

  logical                :: adjust_oceanLandSea_mask = .false. ! default is .false.
  integer                :: nxt, command_argument_count
  character*(128)        :: arg, &
                            Usage = "mk_runofftbl.x CF0012x6C_TM0072xTM0036-Pfafstetter", &
                            mapl_tp_file


  call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
  fileLL=trim(MAKE_BCS_INPUT_DIR)//'/land/route/Outlet_latlon.'


! Read inputs -----------------------------------------------------
  I = command_argument_count()
  if (I < 1 .or. I > 3) then
    print *, " "
    print *, "Wrong number of input arguments, got: ", I
    print *, "Example usage with defaults: "
    print *, " "
    print *, trim(Usage)
    print *, " "
    call exit(1)
  end if

  nxt = 1
  call get_command_argument(nxt, file)
  print *, " "
  print*, "Working with input BCs string: ", file
  print *, " "
  if (I > 1) then
    nxt = nxt + 1
    call get_command_argument(nxt, arg)
    if ( trim(arg) .ne. 'yes') then
      print *, "Incorrect optional second argument, should be: yes"
      call exit(2)
    else
      adjust_oceanLandSea_mask = .true.
      nxt = nxt + 1
      call get_command_argument(nxt, mapl_tp_file)
    endif
  endif
! ------------------------------------------------------------------

  fileT = "til/"//trim(file)//".til" ! input
  fileR = "rst/"//trim(file)//".rst" ! input
  fileO = "til/"//trim(file)//".trn" ! output
  fileB = "til/"//trim(file)//".TRN" ! output

! Read I and J indeces of river outlets.
! These should all be ocean pixels
!---------------------------------------

!  print *, "Getting raster size from "//trim(fileT) 

  open(10,file=fileT, form="formatted", status="old")

  read(10,*) np, pf, nx, ny
  close(10)
!  print *, nx, ny

  write (C_NX, '(i5.5)') NX
  write (C_NY, '(i5.5)') NY

  print *, "Reading outlets..."

  allocate(lats(nx,ny), lons(nx,ny),stat=status)
  if(status/=0) then
     print *, "Out of Memory"
     stop __LINE__
  endif

  open (30,file=trim(fileLL)//C_NX//'x'//C_NY,form="unformatted",status="old")
  do j = 1, ny
     read (30) lons(:,j)
     read (30) lats(:,j)
  end do
  close(30)

!  do j=1,ny
!!     if (mod(j,100) == 0) print *,'J=',j
!     do i=1,nx
!        ii = Lons(i,j)
!        jj = lats(i,j)
!
!        if(ii==-999 .or. jj==-999) then
!           !              ii = i
!           !              jj = j
!           cycle
!        endif
!
!        if(ii==i .and. jj==j) then
!           print *, '>>> Inland Ocean Point ', ii, jj, lons(i,j), lats(i,j)
!           stop
!        end if
!
!    end do
! end do
! stop "DONE"
! Count the number of Ocean and land tiles in the tile file
!  All land tiles preceed the ocean tiles.
!----------------------------------------------------------

! If asked for, adjust tiles to be 
! comptabile with ocean model land-sea mask and write ANOTHER output file
!-------------------------------------------------------------------------

  if (adjust_oceanLandSea_mask) then
    fileBB = "til/"//trim(file)//"_oceanMask_adj.TRN" ! output

    print *, " "
    print *, "Accounting for any mismatch between land-sea masks:"
    print *, "- Of GEOS land and external ocean model."
    print *, "- Output file: ", fileB
    print *, " "
    call read_oceanModel_mask( mapl_tp_file)
!   ... some adjustment of following variable: `type` 
!   ... using ocean model land-sea mask should be done here
  endif

!  print *, "Reading til file "//trim(fileT) 

  open(10,file=fileT, form="formatted", status="old")

  read(10,*) np

  allocate(area(np), in(np), out(np))

  out = 0.0
  in  = 0.0

  do i=1,7
     read(10,*)
  enddo

  lnd=-1
  do l=1,np
     read(10,*) type, area(l)
     if(type==0 .and. lnd==-1) lnd = l-1
     if(lnd<0) in(l) = 1.
  end do

  print *, "area of sphere = ", sum(real(area,kind=8))
  print *, "area of land   = ", sum(real(area*in,kind=8))

  close(10)

  print *, "Number of Tiles ", np
  print *, "Ocean Tiles     ", np-lnd
  print *, "Land tiles      ", lnd

! Read the raster file 
!---------------------

  print *, "Reading rst file "//trim(fileR) 

  open(20,file=fileR,form="unformatted",status="old")

  allocate(rst(nx,ny),stat=status)
  if(status/=0) then
     print *, "Out of Memory"
     stop __LINE__
  endif

  do j=1,ny
     read(20) rst(:,j)
  enddo

  close(20)

  allocate(SortArr(2*lnd,3))

  DstTile   => SortArr(:,1)
  SrcTile   => SortArr(:,2)
  SrcWeight => SortArr(:,3)

! Hash the raster 
!----------------

  print *, "Hashing raster... "

  Hash  = MAPL_HashCreate(1024)
  HashC = MAPL_HashCreate(1024)

  NumTrans=0

  do j=1,ny
!     if (mod(j,100) == 0) print *,'J=',j
     do i=1,nx
        if(rst(i,j)<=lnd) then
           ii = Lons(i,j)
           jj = lats(i,j)

           if(ii==-999 .or. jj==-999) then
!              ii = i
!              jj = j
              cycle
           endif

           if(ii==i .and. jj==j) then
              print *, '>>> Inland Ocean Point ', ii, jj, rst(i,j)
              stop
           end if

           k = MAPL_HASHIncrement(HashC,rst(i,j))
           k = MAPL_HASHIncrement(Hash,rst(ii,jj),rst(i,j))

           SrcWeight(k) = SrcWeight(k) + 1

           if(k>NumTrans) then
              if(k/=NumTrans+1) then
                 print *, NumTrans, k
                 stop 666
              end if
              NumTrans=k
              SrcTile(NumTrans) = rst(i ,j )
              DstTile(NumTrans) = rst(ii,jj)
           endif

        end if
     end do
  end do

  DstTile   => SortArr(:NumTrans,1)
  SrcTile   => SortArr(:NumTrans,2)
  SrcWeight => SortArr(:NumTrans,3)

  print *, "Total Transactions ", NumTrans
  print *, MAPL_HashSize(Hash),MAPL_HashSize(HashC)

  call MAPL_HashDestroy(Hash)
  call MAPL_HashDestroy(HashC)

! Allocate space for transanction lists
!--------------------------------------

  allocate(key(numTrans),stat=status)

  if(status/=0) then
     print *, "Out of Memory"
     stop
  endif

! Sort transactions by source tile number to compute source
!  fractions going into each transaction.
!----------------------------------------------------------

  print *, "Sorting transactions by source..."

  Key = SrcTile

  call MAPL_Sort(Key,SortArr(:NumTrans,:),DIM=1)

  print *, "Computing weights..."

  deallocate(key)
  allocate(SrcFraction(numTrans))

! Compute fractions
!------------------

  is = 1
  ie = 1

  do j=2,NumTrans
     if(SrcTile(j)/=SrcTile(is)) then
        SrcFraction(is:ie) = SrcWeight  (is:ie)
        SrcFraction(is:ie) = SrcFraction(is:ie) / float(sum(SrcWeight(is:ie)))
        is = j
        ie = j
     else
        ie = ie + 1
     end if
  end do

  SrcFraction(is:ie) = SrcWeight  (is:ie)
  SrcFraction(is:ie) = SrcFraction(is:ie) / float(sum(SrcWeight(is:ie)))

  print *,"SrcWeight", sum(SrcFraction), lnd

  print *, '<<<', sum(SrcFraction*area(SrcTile))

  SrcFraction = SrcFraction * (area(SrcTile)/area(DstTile))

  print *, '>>>', sum(SrcFraction*area(DstTile))

! Write output files
!-------------------

  print *, "Writing output file..."

  open(10,file=fileO, form="formatted", status="unknown")
  write(10,*) NumTrans
  do k=1,NumTrans
     write(10,"(2I10,f16.8)") SrcTile(k),DstTile(k),SrcFraction(k)
  end do
  close(10)

  call write_route_file( fileB, NumTrans, SrcTile, DstTile, SrcFraction)
  if (adjust_oceanLandSea_mask) &
     call write_route_file( fileBB, NumTrans, SrcTile, DstTile, SrcFraction)

  do j=1,NumTrans
     Out(DstTile(j)) = Out(DstTile(j)) + In(SrcTile(J))*SrcFraction(J)
  enddo
  print *, "area of land   = ",    sum(real(area*out,kind=8))


  print *, "Completed successfully"

  deallocate( SrcFraction)
  deallocate( SortArr)
  deallocate( rst)
  deallocate( area)
  deallocate( lats, lons)

  call exit(0)

! -----------------------------------------------------------------
contains

  subroutine read_oceanModel_mask( mask_file)
  implicit none
  character*128,    intent(in)  :: mask_file
  integer                       :: nx, ny
  real, allocatable             :: wetMask(:,:)

  integer :: ncid, varid

  print *, "Reading ocean model mask from : ", mask_file
  call check( nf90_open(mask_file, nf90_nowrite, ncid)) ! open nc file

  call check( nf90_inq_dimid(ncid, "n_center_x", varid)) ! read dimenstion (x)
  call check( nf90_inquire_dimension(ncid, varid, len=nx))

  call check( nf90_inq_dimid(ncid, "n_center_y", varid)) ! read dimenstion (y)
  call check( nf90_inquire_dimension(ncid, varid, len=ny))

  allocate( wetMask(nx, ny))
  call check( nf90_inq_varid(ncid, "mask", varid))  ! read mask
  call check( nf90_get_var(ncid, varid,  wetMask))

  call check( nf90_close(ncid)) ! close nc file

  deallocate( wetMask)
  end subroutine read_oceanModel_mask 
! ----------------------

  subroutine write_route_file( fileB, NumTrans, SrcTile, DstTile, SrcFraction)
  implicit none
  character*100,    intent(in) :: fileB
  integer,          intent(in) :: NumTrans
  integer, pointer, intent(in) :: srctile(:), dsttile(:)
  real,             intent(in) :: SrcFraction(:)

  open(10,file=fileB, form="unformatted", status="unknown")
  write(10) NumTrans
  write(10) SrcTile
  write(10) DstTile
  write(10) SrcFraction
  close(10)
  end subroutine write_route_file
! ----------------------

  subroutine check(status)
  implicit none
  integer, intent (in) :: status
  if (status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  print *, "Error in reading ocean mask file."
  stop
  end if
  end subroutine check
! -----------------------------------------------------------------

end program Runoff
