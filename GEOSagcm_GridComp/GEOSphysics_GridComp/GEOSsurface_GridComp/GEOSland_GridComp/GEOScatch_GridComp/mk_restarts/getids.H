
! -----------------------------------------------------------------------------------

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

subroutine GetIds_fast_1p (loni,lati,lon,lat,Id)
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

 end subroutine GetIds_fast_1p

 ! ***************************************************************************** 

  subroutine GetIds_accurate_mpi (loni,lati,lono,lato,Id, tid_in)
  
  implicit none

  integer                             :: NT_IN, NT_OUT, n, i, nplus
  real,    dimension (:), intent (in) :: loni,lati,lono,lato
  integer, dimension (:), intent (in) :: tid_in 
  integer, dimension (:), intent (inout) :: id

  logical                             :: tile_found
  logical, allocatable, dimension(:)  :: mask
  integer, allocatable, dimension (:) :: sub_tid
  real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist
  real                                :: dw, dx, dy, min_lon, max_lon, min_lat, max_lat

  NT_IN  = SIZE (loni)
  NT_OUT = SIZE (lono)

  allocate (mask   (1:  NT_IN))

  Id = -9999

  OUT_TILES : do n = 1,  NT_OUT
	 
     dw = 0.5

     ZOOMOUT : do  

        tile_found = .false. 
        
        ! Min/Max lon/lat of the working window
        ! -------------------------------------
        
        min_lon = MAX(lono (n) - dw, -180.)
        max_lon = MIN(lono (n) + dw,  180.)
        min_lat = MAX(lato (n) - dw,  -90.)
        max_lat = MIN(lato (n) + dw,   90.) 

        mask = .false.
        mask =  ((lati >= min_lat .and. lati <= max_lat).and.(loni >= min_lon .and. loni <= max_lon))
        nplus =  count(mask = mask)
        
        if(nplus < 0) then
           dw = dw + 0.5
           CYCLE
        endif
        
        allocate (sub_tid (1:nplus))
        allocate (sub_lon (1:nplus))
        allocate (sub_lat (1:nplus))
        allocate (rev_dist  (1:nplus))
        
        sub_tid = PACK (tid_in , mask= mask) 
        sub_lon = PACK (loni   , mask= mask)
        sub_lat = PACK (lati   , mask= mask)
        
        ! compute distance from the tile
        
        sub_lat = sub_lat * MAPL_PI/180.
        sub_lon = sub_lon * MAPL_PI/180.
        
        SEEK : if(Id (n) < 0) then
           
           rev_dist  = 1.e20
           
           do i = 1,nplus
              
              rev_dist(i) = haversine(to_radian(lato(n)), to_radian(lono(n)), &
                   sub_lat(i), sub_lon(i))
              
           end do
           
           FOUND : if(minval (rev_dist) < 1.e19) then               
              Id (n) = sub_tid(minloc(rev_dist,1)) 
              tile_found = .true.                  
           endif FOUND
                
        endif SEEK
          
        deallocate (sub_tid, sub_lon, sub_lat, rev_dist)
        
        if(tile_found) GO TO 100
        
        ! if not increase the window size
        dw = dw + 0.5
           
     end do ZOOMOUT
 
100  continue
      
      if(mod (n,10000) == 0)  print *, id(n), loni(id(n)), lono(n), lati(id(n)), lato(n)
  END do OUT_TILES

  deallocate (mask)

 end subroutine GetIds_accurate_mpi

  ! *****************************************************************************

   function to_radian(degree) result(rad)

     ! degrees to radians
     real,intent(in) :: degree
     real :: rad

     rad = degree*MAPL_PI/180.

   end function to_radian

   ! *****************************************************************************
   
   real function haversine(deglat1,deglon1,deglat2,deglon2)
     ! great circle distance -- adapted from Matlab 
     real,intent(in) :: deglat1,deglon1,deglat2,deglon2
     real :: a,c, dlat,dlon,lat1,lat2
     real,parameter :: radius = MAPL_radius
     
!     dlat = to_radian(deglat2-deglat1)
!     dlon = to_radian(deglon2-deglon1)
     !     lat1 = to_radian(deglat1)
!     lat2 = to_radian(deglat2)
     dlat = deglat2-deglat1
     dlon = deglon2-deglon1
     lat1 = deglat1
     lat2 = deglat2     
     a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
     if(a>=0. .and. a<=1.) then
        c = 2*atan2(sqrt(a),sqrt(1-a))
        haversine = radius*c / 1000.
     else
        haversine = 1.e20
     endif
   end function
   
   ! *****************************************************************************

   subroutine ReadCNTilFile (InCNTileFile, ntiles, xlon, xlat,mask)

     implicit none
     character(*), intent (in) :: InCNTileFile
     integer , intent (inout)  :: ntiles
     real, pointer, dimension (:)    :: xlon, xlat
     integer, optional, intent(IN) :: mask
     integer :: n,icnt,ityp, nt, umask, i
     real    :: xval,yval, pf
     real,  allocatable :: ln1(:), lt1(:)
   
   if(present(mask)) then
     umask = mask
   else
     umask = 100
   endif

   open(11,file=InCNTileFile, &
        form='formatted',action='read',status='old')
   read (11,*, iostat=n) Nt

   allocate(ln1(Nt),lt1(Nt))
	  
   do n = 1,7 ! skip header
      read(11,*)
   end do
  
   icnt = 0

   do i=1,Nt
      read(11,*) ityp,pf,xval,yval
      if(ityp == umask) then
         icnt = icnt + 1
         ln1(icnt) = xval
         Lt1(icnt) = yval
      endif
   end do

   close(11)
    
   Ntiles = icnt
   allocate(xlon(Ntiles))
   allocate(xlat(Ntiles))
   xlon = ln1(:Ntiles)
   xlat = lt1(:Ntiles) 

   end subroutine ReadCNTilFile


