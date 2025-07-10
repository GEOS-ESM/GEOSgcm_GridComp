module mk_restarts_getidsMod
  use MAPL
  implicit none
! -----------------------------------------------------------------------------------
  private

  public :: GetIds
  public :: ReadTileFile_IntLatLon    ! returns integer lat/lon for fast but inaccurate processing
  public :: ReadTileFile_RealLatLon   ! returns real lat/lon for slow but accurate processing
  public :: to_radian                 ! should really be replaced with "MAPL_DEGREES_TO_RADIANS" 
  public :: haversine

  interface GetIds
     module procedure GetIds_fast_1p
     module procedure GetIds_accurate_mpi
     module procedure GetIds_carbon
  end interface

contains

   subroutine ReadTileFile_IntLatLon(Tf, ntiles, zoom, lon_int, lat_int, mask)
   
     ! Read *.til tile definition file, return integer lat/lon for fast but inaccurate processing.
     ! Can handle "old" format of *.til files, but that is probably obsolete as of March 2020 and
     !   not used in related subroutine ReadTileFile_RealLatLon().
     ! WARNING: Do NOT use returned "Pf" values.  The content of the 2nd column of the *.til file
     !          that is read into "Pf" depends on whether the file is for EASE or cube-sphere grid tiles!
     ! - reichle, 4 Mar 2020
     
     character*(*), intent(IN)  :: Tf
     integer, intent(out)       :: ntiles
     integer, intent(in)        :: zoom
     integer, pointer, optional :: lon_int(:), lat_int(:)
     integer, optional, intent(IN) :: mask
  
     real, pointer :: xlon(:), xlat(:)
 
     real    :: dum(4),dum1,lnn,ltt
     integer :: de, ce, st
   
  
     if (present(lon_int) .and. present(lat_int)) then
        de=180*zoom
        ce=360*zoom
        call ReadTileFile_RealLatLon(Tf, ntiles, xlon=xlon, xlat=xlat, mask=mask)
        allocate(lon_int(ntiles), lat_int(ntiles))
        lon_int = nint(xlon*zoom)
        lat_int = max(min(nint(xlat*zoom),90*zoom),-90*zoom)
        where(lon_int<-de) lon_int = lon_int + ce
        where(lon_int> de) lon_int = lon_int - ce
        deallocate(xlon, xlat)
     else
        call ReadTileFile_RealLatLon(Tf, ntiles, mask=mask)
     endif

   end subroutine ReadTileFile_IntLatLon
   
   subroutine GetStencil(ii,jj,st)
     integer, intent(OUT) :: ii(0:), jj(0:)
     integer, intent( IN) :: st 
   
     integer :: n, i, j,k, iz, jz, di, dj
   
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
   
    ! *****************************************************************************
   
   subroutine GetIds_fast_1p (loni,lati,lon,lat,zoom,Id)
     integer, dimension(:), intent( IN) :: loni,lati,lon,lat
     integer, intent(in) :: zoom
     integer, dimension(:), intent(OUT) :: Id
   
     integer, allocatable :: Idx(:)
     integer :: i, k, l,n, last, iex, lonx, hash
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
                 print *, 'Failed to find valid data for tile ',i, k
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
   
    subroutine GetIds_carbon (loni,lati,lono,lato,Id, tid_in, &
         CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2, &
         fveg_offl, ityp_offl,isCLM51)

    use clm_varpar_shared , only : nveg_40 => NUM_VEG_CN, nveg_51 => NUM_VEG_CN51, &
                                 npft => numpft_CN, npft_51 => numpft_CN51
      implicit none
      integer :: nveg
      real,    parameter :: fmin= 1.e-4 ! ignore vegetation fractions at or below this value
      integer :: iclass_40_45(npft) = (/1,1,2,3,3,4,5,5,6,7,8,9,10,11,12,11,12,11,12/)
      integer :: iclass_51(npft_51) = (/1,1,2,3,3,4,5,5,6,7,9,10,11,11,11/)
      integer, dimension(:), allocatable :: iclass
      integer :: NT_IN, NT_OUT, n, i, nplus,nv, nx, ityp_new 
      integer, dimension (:), intent (in) :: tid_in 
      integer, dimension (:,:), intent (inout) :: id
      real, dimension (:), intent (in) :: loni,lati,lono,lato
      real, dimension (:), intent (in) :: CLMC_pf1, CLMC_pf2, CLMC_sf1, CLMC_sf2, &
           CLMC_pt1, CLMC_pt2,CLMC_st1,CLMC_st2
      real, dimension(:,:), intent (in)   :: fveg_offl, ityp_offl
      logical, intent(in), optional       :: isCLM51
      logical                             :: tile_found
      logical, allocatable, dimension (:) :: mask
      integer, allocatable, dimension (:) :: sub_tid
      real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist, sub_fevg1, sub_fevg2
      integer, allocatable, dimension (:) :: sub_ityp1, sub_ityp2,icl_ityp1
      real                                :: dw, dx, dy, min_lon, max_lon, min_lat, max_lat, fveg_new, sub_dist
      
      NT_IN  = SIZE (loni)
      NT_OUT = SIZE (lono)

      if (isCLM51) then
         allocate(iclass(1:npft_51))
         iclass = iclass_51
         nveg = nveg_51   
      elseif (.not.isCLM51) then
        allocate(iclass(1:npft))
        iclass = iclass_40_45
        nveg = nveg_40
      end if 
      
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
            
            allocate (sub_tid   (1:nplus))
            allocate (sub_lon   (1:nplus))
            allocate (sub_lat   (1:nplus))
            allocate (rev_dist  (1:nplus))
            allocate (sub_ityp1 (1:nplus))
            allocate (sub_fevg1 (1:nplus))
            allocate (sub_ityp2 (1:nplus))
            allocate (sub_fevg2 (1:nplus))
            allocate (icl_ityp1 (1:nplus))
            
            sub_tid = PACK (tid_in , mask= mask) 
            sub_lon = PACK (loni   , mask= mask)
            sub_lat = PACK (lati   , mask= mask)
            
            ! compute distance from the tile
            
            sub_lat = sub_lat * MAPL_PI/180.
            sub_lon = sub_lon * MAPL_PI/180.
   
            NV_LOOP: do nv = 1, nveg
               
               if (isCLM51) then
                 if (nv == 1) ityp_new = CLMC_pt1(n)
                 if (nv == 1) fveg_new = CLMC_pf1(n)
                 if (nv == 2) ityp_new = CLMC_st1(n)
                 if (nv == 2) fveg_new = CLMC_sf1(n)

                else if (.not.isCLM51) then

                 if (nv == 1) ityp_new = CLMC_pt1(n)
                 if (nv == 1) fveg_new = CLMC_pf1(n)
                 if (nv == 2) ityp_new = CLMC_pt2(n)
                 if (nv == 2) fveg_new = CLMC_pf2(n)
                 if (nv == 3) ityp_new = CLMC_st1(n)
                 if (nv == 3) fveg_new = CLMC_sf1(n)
                 if (nv == 4) ityp_new = CLMC_st2(n)
                 if (nv == 4) fveg_new = CLMC_sf2(n)

                 end if


                 SEEK : if((Id (n, nv) < 0).and.(fveg_new > fmin)) then

                 if (isCLM51) then
                    if(nv <= 1) then ! index for secondary PFT index if primary or primary if secondary
                       nx = nv + 1
                    else
                       nx = nv - 1
                    endif

                 else if (.not.isCLM51) then

                    if(nv <= 2) then ! index for secondary PFT index if primary or primary if secondary
                     nx = nv + 2
                    else
                       nx = nv - 2
                    endif
                 endif                  

                  sub_ityp1 = ityp_offl (sub_tid,nv)
                  sub_fevg1 = fveg_offl (sub_tid,nv)
                  sub_ityp2 = ityp_offl (sub_tid,nx)
                  sub_fevg2 = fveg_offl (sub_tid,nx)
                  
                  rev_dist  = 1.e20
                  icl_ityp1 = iclass(sub_ityp1)
                  
                  do i = 1,nplus
                     if(  ( sub_fevg1(i)>fmin .and. ( ityp_new==sub_ityp1(i) .or. iclass(ityp_new)==iclass(sub_ityp1(i)) ) ) &
                          .or.                                                                                               &
                          ( sub_fevg2(i)>fmin .and. ( ityp_new==sub_ityp2(i) .or. iclass(ityp_new)==iclass(sub_ityp2(i)) ) ) &
                          ) then
                        
                        sub_dist = haversine(to_radian(lato(n)), to_radian(lono(n)), &
                             sub_lat(i), sub_lon(i))
                        
                        if(ityp_new == sub_ityp1(i) .and. sub_fevg1(i) >fmin) then
                           rev_dist(i) = 1.*sub_dist     ! give priority to same (primary if primary, secondary if secondary)   
                           ! gkw: these weights are tunable
                        else if(ityp_new ==sub_ityp2(i) .and. sub_fevg2(i)>fmin) then
                           rev_dist(i) = 2.*sub_dist     ! lower priority if not same (secondary if primary, primary if secondary)
                        else if(iclass(ityp_new)==iclass(sub_ityp1(i)) .and. sub_fevg1(i)>fmin) then
                           rev_dist(i) = 3.*sub_dist     ! even lower priority if same of some other PFT in same class
                        else if(sub_fevg2(i)>fmin) then
                           rev_dist(i) = 4.*sub_dist     ! even lower priority if not same of some other PFT in same class
                        else
                           rev_dist(i) = 1.e20
                        endif
                     endif
                     
                  end do
                  
                  FOUND : if(minval (rev_dist) < 1.e19) then               
                     Id (n, nv) = sub_tid(minloc(rev_dist,1)) 
                     
                  endif FOUND
                   
               endif SEEK
            end do NV_LOOP
   
             deallocate (sub_tid, sub_lon, sub_lat, icl_ityp1)
             deallocate (sub_ityp1, sub_fevg1, sub_ityp2, sub_fevg2, rev_dist)  
   
             tile_found = .true.   
             if (isCLM51) then
                if((tile_found).and.((CLMC_pf1(n) > fmin).and.(Id(n,1) < 0))) tile_found = .false.
                if((tile_found).and.((CLMC_sf1(n) > fmin).and.(Id(n,2) < 0))) tile_found = .false.
             else if (.not.isCLM51) then
                if((tile_found).and.((CLMC_pf1(n) > fmin).and.(Id(n,1) < 0))) tile_found = .false.
                if((tile_found).and.((CLMC_pf2(n) > fmin).and.(Id(n,2) < 0))) tile_found = .false.
                if((tile_found).and.((CLMC_sf1(n) > fmin).and.(Id(n,3) < 0))) tile_found = .false.
                if((tile_found).and.((CLMC_sf2(n) > fmin).and.(Id(n,4) < 0))) tile_found = .false.          
             endif

             if(tile_found) GO TO 100
           
             ! if not increase the window size
             dw = dw + 0.5
             
          end do ZOOMOUT
          
   100    continue
         
   !       if(mod (n,10000) == 0)  print *, id(n), loni(id(n)), lono(n), lati(id(n)), lato(n)
       END do OUT_TILES
   
     deallocate (mask)
     
    end subroutine GetIds_carbon
   
     ! *****************************************************************************
   
      function to_radian(degree) result(rad)
   
        ! should really be replaced with "MAPL_DEGREES_TO_RADIANS" (but note single vs double precision)
        ! - reichle, 4 Mar 2020
   
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
   
      subroutine ReadTileFile_RealLatLon (InCNTileFile, ntiles, xlon, xlat, mask)
   
        ! read *.til tile definition file, return *real* lat/lon for slow but accurate processing
      
        implicit none
        character(*), intent (in) :: InCNTileFile
        integer , intent (out)  :: ntiles
        real, pointer, optional, dimension (:) :: xlon, xlat
        integer, optional, intent(IN) :: mask
        integer :: n,icnt,ityp, nt, umask, i, header
        real    :: xval,yval, pf
        real,  allocatable   :: ln1(:), lt1(:)
        real,  pointer       :: AVR(:,:)
        integer              :: filetype, k
        integer, allocatable :: indices(:), indices_tmp(:)
        logical :: isNC4     
 
        if(present(mask)) then
          umask = mask
        else
          umask = 100
        endif
   
        call MAPL_NCIOGetFileType(InCNTileFile, filetype)
        isNC4   = (filetype == MAPL_FILETYPE_NC4)

        if (isNC4) then
          call MAPL_ReadTilingNC4(InCNTileFile, AVR=AVR)
          allocate(indices_tmp(size(AVR,1)))
          k = 0
          do i = 1, size(AVR,1)
            if( int(AVR(i,1)) == umask) then
              k = k+1
              indices_tmp(k) = i
            endif
          enddo
          indices = indices_tmp(1:k)
          Ntiles  = k
          if ( present(xlon) .and. present(xlat)) then
            if(.not.associated (xlon)) allocate(xlon(Ntiles))
            if(.not.associated (xlat)) allocate(xlat(Ntiles))
            xlon = AVR(indices, 3)
            xlat = AVR(indices, 4)
          endif
          deallocate(AVR) 
        else

          open(11,file=InCNTileFile, form='formatted',action='read',status='old')

          ! first read number of lines in the til file header
          ! -------------------------------------------------
          header = 5
          read (11,*, iostat=n) Nt
          do i = 1, header -1
             read (11,*)
          end do
          read (11,*,IOSTAT=n)ityp,pf,xval, yval
          if(n /= 0) header = 8

          rewind (11)

          ! read the tile file
          !-------------------
          read (11,*, iostat=n) Nt
   
          allocate(ln1(Nt),lt1(Nt))

          do n = 1,header-1 ! skip header
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
          if ( present(xlon) .and. present(xlat)) then
             if(.not.associated (xlon)) allocate(xlon(Ntiles))
             if(.not.associated (xlat)) allocate(xlat(Ntiles))
             xlon = ln1(:Ntiles)
             xlat = lt1(:Ntiles)
          endif
        endif !isNC4
 
      end subroutine ReadTileFile_RealLatLon

end module mk_restarts_getidsMod
