program Runoff

  use mapl_hashmod
  use mapl_sortmod

  implicit none

  integer, parameter     :: nx=8640, ny=4320
  integer*2, allocatable :: lats(:,:), lons(:,:)
  integer, pointer       :: rst(:,:), SortArr(:,:), key(:)
  integer, pointer       :: srctile(:),  srcweight(:), dstweight(:), dsttile(:)
  real, allocatable      :: SrcFraction(:), area(:), in(:), out(:)
  integer                :: i,j,k,l, Hash, HashC, ii,jj,kk
  integer                :: type, np,lnd, is,ie,ww
  integer                :: numtrans,  numclosed
  integer                :: status
  character*100          :: file, fileT, fileR, fileO, fileB
  character*100          :: fileLL="data/Outlet_latlon.dat"

  call getarg(1,file)

  fileT = "til/"//trim(file)//".til"
  fileR = "rst/"//trim(file)//".rst"
  fileO = "til/"//trim(file)//".trn"
  fileB = "til/"//trim(file)//".TRN"

! Read I and J indeces of river outlets.
! These should all be ocean pixels
!---------------------------------------

  print *, "Reading outlets..."

  allocate(lats(nx,ny), lons(nx,ny),stat=status)
  if(status/=0) then
     print *, "Out of Memory"
     stop __LINE__
  endif

  open (30,file=fileLL,status="old",form="unformatted")
  read (30) lons
  read (30) lats
  close(30)

! Count the number of Ocean and land tiles in the tile file
!  All land tiles preceed the ocean tiles.
!----------------------------------------------------------

  print *, "Reading til file "//trim(fileT) 

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

  open(20,file=fileR,form="unformatted",status="old",convert="LITTLE_ENDIAN")

  allocate(rst(nx,ny),stat=status)
  if(status/=0) then
     print *, "Out of Memory"
     stop __LINE__
  endif

  do j=1,ny
     read(20) rst(:,j)
  enddo

  close(20)

  allocate(SortArr(1000000,3))

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
     do i=1,nx
        if(rst(i,j)<=lnd) then
           ii = Lons(i,j)
           jj = lats(i,j)

           if(ii==i .and. jj==j) then
              print *, '>>> Inland Ocean Point ', ii, jj, rst(i,j)
              stop
           end if

           if(ii==-999) then
              ii = i
              jj = j
           endif
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

! Allocate space for transanction lists
!--------------------------------------

  allocate(key(numTrans))

  print *, "Total Transactions ", NumTrans
  print *, MAPL_HashSize(Hash),MAPL_HashSize(HashC)


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

! Write output file
!------------------

  print *, "Writing output file..."

  open(10,file=fileO, form="formatted", status="unknown")

  write(10,*) NumTrans

  do k=1,NumTrans
     write(10,"(2I10,f16.8)") SrcTile(k),DstTile(k),SrcFraction(k)
  end do

  close(10)


  open(10,file=fileB, form="unformatted", Convert="LITTLE_ENDIAN", &
       status="unknown")

  write(10) NumTrans
  write(10) SrcTile
  write(10) DstTile
  write(10) SrcFraction

  close(10)

  do j=1,NumTrans
     Out(DstTile(j)) = Out(DstTile(j)) + In(SrcTile(J))*SrcFraction(J)
  enddo

  print *, "area of land   = ",    sum(real(area*out,kind=8))

  print *, "Completed successfully"

  call exit(0)
end program Runoff
