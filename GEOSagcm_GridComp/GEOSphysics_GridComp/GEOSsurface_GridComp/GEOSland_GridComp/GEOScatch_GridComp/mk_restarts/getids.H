

subroutine ReadTileFile(Tf,Pf,Id,lon,lat,mask)
  character*(*), intent(IN) :: Tf
  integer, pointer          :: Pf(:), Id(:), lon(:), lat(:)
  integer, optional, intent(IN) :: mask

  integer, allocatable :: Pf1(:), Id1(:), ln1(:), lt1(:)
  integer :: k, i, nt, pfs, ids,n,msk, umask
  real    :: dum(4),dum1,lnn,ltt
  integer :: de, ce, st
  logical :: old

  de=180*zoom
  ce=360*zoom
  st=2*zoom
  if(present(mask)) then
     umask = mask
  else
     umask = 100
  endif

  print *, "Reading tilefile ",trim(Tf)

  open(unit=20,file=trim(Tf),form='formatted')

  read(20,*,iostat=n) Nt,i,k
  old=n<0
  close(20)

  open(unit=20,file=trim(Tf),form='formatted')

  read(20,*) Nt

  do i=1,7
     read(20,*)
  enddo

  allocate(Pf1(Nt),Id1(Nt),ln1(Nt),lt1(Nt))

  n=0
  do i=1,Nt
     if(old) then
        read(20,*,end=200) msk, Pfs, lnn, ltt
        ids = 0
     else
        read(20,*,end=200) msk, dum1, lnn, ltt, dum, Pfs, Ids
     end if
     if(msk/=umask) cycle
     n = n+1
     pf1(n) = pfs
     Id1(n) = ids
     ln1(n) = nint(lnn*zoom)
     Lt1(n)=max(min(nint(ltt*zoom),90*zoom),-90*zoom)
     if(ln1(n)<-de) ln1(n) = ln1(n) + ce
     if(ln1(n)> de) ln1(n) = ln1(n) - ce
  enddo

200 continue

  close(20)

  Nt=n
  print *, "Found ",nt," land tiles."

  allocate(Pf(Nt),Id(Nt),lon(Nt),lat(Nt))
  Pf  = Pf1(:Nt)
  Id  = Id1(:Nt)
  lon = ln1(:Nt)
  lat = lt1(:Nt)
  deallocate(Pf1,Id1,ln1,lt1)

  return
end subroutine ReadTileFile

subroutine GetStencil(ii,jj,st)
  integer, intent(OUT) :: ii(0:), jj(0:)
  integer, intent( IN) :: st 

  integer :: n, i, j, iz, jz, di, dj

  n=-1
  do i=0,st
     di = 0
     dj = 1
     jz =  0
     iz =  i
     n  = n+1
     ii(n) = iz
     jj(n) = jz

     do k=1,8*i-1
        if    (iz==i.and.jz==-i) then
           di = 0
           dj = 1
        elseif(iz==i.and.jz==i) then
           di = -1
           dj = 0
        elseif(iz==-i.and.jz==i) then
           di = 0
           dj = -1
        elseif(iz==-i.and.jz==-i) then
           di = 1
           dj = 0
        endif

        iz = iz + di
        jz = jz + dj

        if(jz==0 .and. iz == i) exit
        n  = n+1
        ii(n) = iz
        jj(n) = jz
     end do
  end do

!  print *, 'ii = ',ii
!  print *
!  print *, 'jj = ',jj

end subroutine GetStencil

subroutine GetIds(loni,lati,lon,lat,Id)
  integer, dimension(:), intent( IN) :: loni,lati,lon,lat
  integer, dimension(:), intent(OUT) :: Id

  integer, allocatable :: Idx(:)
  integer :: i, k, l, last, iex, lonx, hash
  integer, allocatable :: ii(:)
  integer, allocatable :: jj(:)
  integer :: jx(7) =(/0,1,-1,2,-2,3,-3/)
  integer, allocatable :: ix(:)
  logical :: found
  integer :: de, ce, st

  de=180*zoom
  ce=360*zoom
  st=2*zoom
  allocate(ix(ce),ii(0:(2*st+1)**2-1),jj(0:(2*st+1)**2-1))
  Hash  = MAPL_HashCreate(8*1024)

  n = 1
  do i=1,ce-1,2
     ix(i  ) =  n
     ix(i+1) = -n
     n=n+1
  end do

  call GetStencil(ii,jj,st)

  allocate(Idx(size(loni)))

  do i=1,size(loni)
     k = MAPL_HashIncrement(Hash,loni(i),lati(i))
     idx(k) = i
  end do

  last = MAPL_HashSize(HASH)

  iex = 0

  do i=1,size(lon)
!     k = MAPL_HashIncrement(Hash,lon(i),lat(i))
!     if (k>last) then
        do n=0,size(ii)-1
           lonx=lon(i)+ii(n)
           if(lonx<-de)lonx=lonx+ce
           if(lonx> de)lonx=lonx-ce
           k = MAPL_HashIncrement(Hash,lonx,lat(i)+jj(n))
           if(k<=last) exit
        end do
        if (k>last) then
           iex = iex + 1
           found=.false.
           do l=1,7
              do n=1,ce
                 lonx=lon(i)+ix(n)
                 if(lonx<-de)lonx=lonx+ce
                 if(lonx> de)lonx=lonx-ce
                 lonx=lon(i)+ix(n)
                 k = MAPL_HashIncrement(Hash,lonx,lat(i)+jx(l))
                 if(k<=last) then
                    found=.true.
                    exit
                 end if
              end do
              if(found) exit
           end do
           if(k>last) then
              print *, 'Failed to find valid data for tile ',i, k, InRestart
              print *, 'Thus using last'
	      k = last
           endif
        end if
!     end if
     Id(i) = Idx(k)
  enddo

  deallocate(Idx,ix,ii,jj)

  print *, 'Used extreme measures ', iex, ' times'
  print *

end subroutine GetIds
