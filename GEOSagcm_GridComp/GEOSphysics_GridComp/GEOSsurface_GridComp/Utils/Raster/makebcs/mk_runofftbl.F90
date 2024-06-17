program Runoff
  
  use mapl_hashmod
  use mapl_sortmod
  use netcdf
  
  implicit none
  include 'netcdf.inc'

  integer                :: nx, ny, pf
  integer, allocatable   :: lats(:,:), lons(:,:)
  integer, pointer       :: rst(:,:), SortArr(:,:), key(:)
  integer, pointer       :: srctile(:),  srcweight(:), dstweight(:), dsttile(:)
  real,    allocatable   :: SrcFraction(:), area(:), in(:), out(:)
  integer                :: i,j,k,l, Hash, HashC, ii,jj,kk
  integer                :: type, np,lnd, is,ie,ww
  integer                :: numtrans,  numclosed
  integer                :: status
  character*100          :: Gridname, fileT, fileR, fileO, fileB

  character*400          :: fileLL
  character*400          :: MAKE_BCS_INPUT_DIR

  character*5            :: C_NX, C_NY

  integer                :: nxt, command_argument_count
  character*(128)        :: arg
  character*(128)        :: Usage = "mk_runofftbl.x M6TP0072x0036-Pfafstetter"
  character*(128)        :: mapl_tp_file
  
  ! ------------------------------------------------------------------
  
  call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
  fileLL=trim(MAKE_BCS_INPUT_DIR)//'/land/route/Outlet_latlon.'
  
  ! Read inputs -----------------------------------------------------

  I = command_argument_count()
  if (I /= 1) then
     print *, " "
     print *, "Wrong number of input arguments, got: ", I
     print *, "Example usage with defaults: "
     print *, " "
     print *, trim(Usage)
     print *, " "
     call exit(1)
  end if
  
  nxt = 1
  call get_command_argument(nxt, Gridname)
  print *, " "
  print*, "Working with input BCs string: ", Gridname
  print *, " "

  ! ------------------------------------------------------------------
  
  fileT = "til/"//trim(Gridname)//".til" ! input
  fileR = "rst/"//trim(Gridname)//".rst" ! input
  fileO = "til/"//trim(Gridname)//".trn" ! output
  fileB = "til/"//trim(Gridname)//".TRN" ! output
  
  ! Read I and J indices of river outlets.
  ! These should all be ocean pixels
  ! ---------------------------------------
  
  !  print *, "Getting raster size from "//trim(fileT) 

  open(10,file=fileT, form="formatted", status="old")
  
  read(10,*) np, pf, nx, ny
  close(10)
  
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
  
  print *, " "
  print *, "Determining river outlets to ocean:"
  print *, "- Output file: ", fileB
  print *, " "
  call outlets_to_ocean(Gridname,lons,lats,nx,ny)
  
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
  ! ---------------------
  
  print *, "Reading raster (rst) file "//trim(fileR) 
  
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
  ! ----------------
  
  print *, "Hashing raster... "
  
  Hash  = MAPL_HashCreate(1024)
  HashC = MAPL_HashCreate(1024)
  
  NumTrans=0
  
  do j=1,ny
     do i=1,nx
        if(rst(i,j)<=lnd) then
           ii = Lons(i,j)
           jj = lats(i,j)
           
           if(ii==-999 .or. jj==-999) then
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
  ! --------------------------------------
  
  allocate(key(numTrans),stat=status)
  
  if(status/=0) then
     print *, "Out of Memory"
     stop
  endif
  
  ! Sort transactions by source tile number to compute source
  !  fractions going into each transaction.
  ! ----------------------------------------------------------
  
  print *, "Sorting transactions by source..."
  
  Key = SrcTile

  call MAPL_Sort(Key,SortArr(:NumTrans,:),DIM=1)
  
  print *, "Computing weights..."
  
  deallocate(key)
  allocate(SrcFraction(numTrans))
  
  ! Compute fractions
  ! ------------------
  
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
  ! -------------------
  
  print *, "Writing output file..."
  
  open(10,file=fileO, form="formatted", status="unknown")
  write(10,*) NumTrans
  do k=1,NumTrans
     write(10,"(2I10,f16.8)") SrcTile(k),DstTile(k),SrcFraction(k)
  end do
  close(10)
  
  call write_route_file( fileB, NumTrans, SrcTile, DstTile, SrcFraction)
  
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
  
  subroutine write_route_file( fileB, NumTrans, SrcTile, DstTile, SrcFraction)
    
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
  
  ! ------------------------------------------------------------------------
  ! The subroutine moves outlets of each land&Greenland gridcell (with endpoint to ocean only) 
  ! to the nearest ocean gridcell defiend by "Gridname" and ocean mask if the origional outlets are not in the ocean.
  subroutine outlets_to_ocean(Gridname,lons,lats,nx,ny)
    
    integer,           intent(in)      :: nx,ny !number of lon and lat, eg: 43200 and 21600 in 30s resolution
    character(len=*),  intent(in)      :: Gridname !name of the domain, eg: CF0180x6C_M6TP0072x0036-Pfafstetter
    integer,           intent(inout)   :: lons(nx,ny),lats(nx,ny)
    !lons(nx,ny): lon idx of outlets of each land&Greenland gridcell (with endpoint to ocean only)
    !lats(nx,ny): lat idx of outlets of each land&Greenland gridcell (with endpoint to ocean only)
    !lons(nx,ny) and lats(nx,ny) are got from input binary file, eg: Outlet_latlon.43200x21600 
    ! -----------------------------------------------------------
    
    integer,       allocatable, dimension(:)   :: lati_lnd,loni_lnd  
    integer,       allocatable, dimension(:)   :: msk1d
    integer,       allocatable, dimension(:,:) :: msk2d  
    integer,       allocatable, dimension(:,:) :: mask  
    integer,       allocatable, dimension(:,:) :: boundary  
    real*8,        allocatable, dimension(:)   :: lonsh,latsh
    real*8,        allocatable, dimension(:)   :: lons_adj,lats_adj
    integer,       allocatable, dimension(:)   :: lati_ocn,loni_ocn   
    character*100                              :: Gridname_ocn
    character*100                              :: Gridname_ocn_lnd
    character*100                              :: res_MAPL 
    integer,       allocatable, dimension(:,:) :: rst_ocn,rst_ocn_lnd
    integer                                    :: nt_ocn_lnd,nl_ocn_lnd,nt_ocn,nx_MAPL,ny_MAPL,nsh
    integer,       pointer, dimension(:)       :: t2lati,t2loni
    real*8,        allocatable, dimension(:)   :: lon30s,lat30s
    integer                                    :: ns
    integer,       allocatable, dimension(:,:) :: ns_map
    real*8,        allocatable, dimension(:)   :: lat_lnd,lon_lnd
    integer                                    :: i,j,l,k,status,type,np,flag,flag2
    
    call get_domain_name(Gridname,Gridname_ocn,Gridname_ocn_lnd,res_MAPL,nx_MAPL,ny_MAPL)
    allocate(rst_ocn(nx,ny),rst_ocn_lnd(nx,ny))
    allocate(lon30s(nx),lat30s(ny))
    call read_rst_til_files(nx,ny,Gridname_ocn,Gridname_ocn_lnd,rst_ocn,rst_ocn_lnd,t2loni,t2lati,nt_ocn_lnd,nl_ocn_lnd,nt_ocn,lon30s,lat30s)
    !print *,"running outlets_num() ..."
    call outlets_num(rst_ocn_lnd,nl_ocn_lnd,nt_ocn_lnd,lons,lats,nx,ny,ns)
    !print *,"outlets num is ",ns
    allocate(loni_lnd(ns),lati_lnd(ns))
    allocate(lons_adj(ns),lats_adj(ns))
    allocate(loni_ocn(ns),lati_ocn(ns))  
    allocate(ns_map(nx,ny))
    allocate(lon_lnd(ns),lat_lnd(ns))
    !print *,"running retrieve_outlets() ..."
    call retrieve_outlets(lons,lats,lon30s,lat30s,loni_lnd,lati_lnd,lon_lnd,lat_lnd,ns_map,nx,ny,ns)
    !print *,"running mask_MAPL_1d() ..."
    allocate(msk1d(nt_ocn))  
    call mask_MAPL_1d(msk1d,t2loni,t2lati,nt_ocn,res_MAPL,nx_MAPL,ny_MAPL)
    deallocate(t2loni,t2lati)
    !print *,"running mask_MAPL_2d() ..."
    allocate(msk2d(nx,ny))  
    call mask_MAPL_2d(rst_ocn,msk1d,msk2d,nt_ocn,nx,ny)
    deallocate(rst_ocn,msk1d)
    !print *,"running mask_MAPL_bcs() ..."
    allocate(mask(nx,ny))
    call mask_MAPL_bcs(rst_ocn_lnd,msk2d,mask,nx,ny,nl_ocn_lnd,nt_ocn_lnd)
    deallocate(msk2d,rst_ocn_lnd)
    !print *,"running ocean_boundary() ..."
    allocate(boundary(nx,ny))  
    call ocean_boundary(mask,boundary,nx,ny)
    !print *,"running ocean_boundary_num() ..."
    call ocean_boundary_num(boundary,nx,ny,nsh)  
    !print *,"ocean boundary point num is ",nsh
    allocate(lonsh(nsh),latsh(nsh))  
    !print *,"running ocean_boundary_points() ..."
    call ocean_boundary_points(boundary,lon30s,lat30s,lonsh,latsh,nx,ny,nsh)
    deallocate(boundary)
    !print *,"running move_to_ocean() ..."
    call move_to_ocean(loni_lnd,lati_lnd,lon_lnd,lat_lnd,mask,lonsh,latsh,lons_adj,lats_adj,ns,nx,ny,nsh) 
    deallocate(mask,lonsh,latsh) 
    !print *,"running sinkxy_ocean() ..."
    call sinkxy_ocean(lons_adj,lats_adj,lon30s,lat30s,loni_ocn,lati_ocn,ns,nx,ny)  
    !print *,"running update_outlets() ..."
    call update_outlets(loni_ocn,lati_ocn,ns_map,lons,lats,nx,ny,ns)
    
    deallocate(loni_lnd,lati_lnd,lons_adj,lats_adj,loni_ocn,lati_ocn)
    deallocate(lon30s,lat30s)
    deallocate(ns_map,lon_lnd,lat_lnd)
    
  end subroutine outlets_to_ocean

!-------------------------------------------------------------------------
! This subroutine gets the name of 'Gridname_ocn' and 'file_ocn_lnd' from the input name 'Gridname'.
! It also gets resolution of the ocean domain.
  subroutine get_domain_name(Gridname,Gridname_ocn,Gridname_ocn_lnd,res_MAPL,nx_MAPL,ny_MAPL)

    character(len=*),  intent(in)      :: Gridname !input domain name, eg: CF0180x6C_M6TP0072x0036-Pfafstetter
    character(len=*),  intent(out)     :: Gridname_ocn,Gridname_ocn_lnd,res_MAPL
    !Gridname_ocn: ocean domain name, eg: M6TP0072x0036
    !Gridname_ocn_lnd: ocean-land domain name, eg: M6TP0072x0036-Pfafstetter
    !res_MAPL: ocean resolution name, eg: M6TP0072x0036
    integer,           intent(out)     :: nx_MAPL,ny_MAPL
    !nx_MAPL: number of lon of ocean domain, eg: 72
    !ny_MAPL: number of lat of ocean domain, eg: 36

    character*100   :: nx_str,ny_str,fileT 
    integer         :: px,plats,plate,plons,plone,plonss,pocns,pocne  
    integer         :: nstr1,nstr2
    integer         :: i,length

    Gridname_ocn=""
    Gridname_ocn_lnd=""
    res_MAPL=""
    nx_MAPL=""
    ny_MAPL=""
    fileT = "til/"//trim(Gridname)//".til" 
    open(10,file=fileT, form="formatted", status="old")
    do i=1,5
      read(10,*)      
    enddo
    read(10,*)Gridname_ocn_lnd
    do i=100,1,-1
      if(Gridname_ocn_lnd(i:i).eq."-")then
        Gridname_ocn(1:i-1)=Gridname_ocn_lnd(1:i-1)
        exit
      endif
    enddo
    read(10,*)nx_MAPL
    read(10,*)ny_MAPL
    nx_str=""
    ny_str=""
    write(nx_str,*)nx_MAPL
    write(ny_str,*)ny_MAPL
    res_MAPL=""
    length = len( trim(adjustl(nx_str))//"x"//trim(adjustl(ny_str)) )
    res_MAPL(1:length)=trim(adjustl(nx_str))//"x"//trim(adjustl(ny_str))
  

  end subroutine get_domain_name

!-------------------------------------------------------------------------
! This subroutine reads rst and til files
  subroutine read_rst_til_files(nx,ny,Gridname_ocn,Gridname_ocn_lnd,rst_ocn,rst_ocn_lnd,t2loni,t2lati,nt_ocn_lnd,nl_ocn_lnd,nt_ocn,lon30s,lat30s)

    integer,           intent(in)              :: nx,ny !number of lon and lat, eg: 43200 and 21600 in 30s resolution  
    character(len=*),  intent(in)              :: Gridname_ocn,Gridname_ocn_lnd ! input filename, eg: M6TP0072x0036 and M6TP0072x0036-Pfafstetter

    integer,       intent(out)                 :: rst_ocn(nx,ny),rst_ocn_lnd(nx,ny) !data from rst files
    integer,pointer,intent(out), dimension(:)  :: t2loni,t2lati !relationship between ocn tile idx in M6TP0072x0036.rst and lat/lon idx in MAPL_Tripolar.nc
    integer,       intent(out)                 :: nt_ocn_lnd,nl_ocn_lnd,nt_ocn
    !nt_ocn_lnd: number of total tiles in the M6TP0072x0036-Pfafstetter
    !nl_ocn_lnd: number of land tiles in the M6TP0072x0036-Pfafstetter
    !nt_ocn: number of total tiles in the M6TP0072x0036
    real*8,        intent(out)                 :: lon30s(nx),lat30s(ny)!lon and lat value arrays of the 30s map

    character*100  :: fileT_ocn, fileR_ocn
    character*100  :: fileT_ocn_lnd, fileR_ocn_lnd    
    integer        :: i,j,type,np,k,l
    real           :: num1,num2,num3,num4    
    real*8         :: dx,dy

    fileT_ocn = "til/"//trim(Gridname_ocn)//".til"          ! input
    fileR_ocn = "rst/"//trim(Gridname_ocn)//".rst"          ! input
    fileT_ocn_lnd = "til/"//trim(Gridname_ocn_lnd)//".til"  ! input
    fileR_ocn_lnd = "rst/"//trim(Gridname_ocn_lnd)//".rst"  ! input  
    
    !print *, "Reading rst file "//trim(fileR_ocn) 
    open(20,file=fileR_ocn,form="unformatted",status="old")
    do j=1,ny
       read(20) rst_ocn(:,j)
    enddo
    close(20)
    
    !print *, "Reading rst file "//trim(fileR_ocn_lnd) 
    open(21,file=fileR_ocn_lnd,form="unformatted",status="old")
    do j=1,ny
       read(21) rst_ocn_lnd(:,j)
    enddo
    close(21)  
    
    open(10,file=fileT_ocn, form="formatted", status="old")
    read(10,*) nt_ocn
    allocate(t2lati(nt_ocn),t2loni(nt_ocn))
    do i=1,4
       read(10,*)
    enddo
    do l=1,nt_ocn
       read(10,*)type,num1,num2,num3,t2loni(l),t2lati(l)
    enddo
    close(10)
    
    open(10,file=fileT_ocn_lnd, form="formatted", status="old")
    read(10,*) np
    do i=1,4
       read(10,*)
    enddo
    k=0
    do l=1,np
       read(10,*)type,num1,num2,num3,num4
       if(type/=100)exit
       k=k+1
    enddo
    close(10)
    nt_ocn_lnd=np
    nl_ocn_lnd=k
    
    dx=360.d0/nx
    dy=180.d0/ny
    do i=1,nx
       lon30s(i)=-180.d0+dx/2.d0+dx*(i-1)
    enddo
    do j=1,ny
       lat30s(j)=-90.d0+dy/2.d0+dy*(j-1)
    enddo    

  end subroutine read_rst_til_files

!-------------------------------------------------------------------------
 ! This subroutine counts the number of outlet points from input outlets lon/lat idx map.
  subroutine outlets_num(rst_ocn_lnd,nl,nt,lons,lats,nx,ny,ns)

    integer, intent(in)                  :: nx,ny,nl,nt
    !nx: number of lon of 30s map, 43200
    !ny: number of lat of 30s map, 21600
    !nl: number of land tiles
    !nt: number of total tiles
    integer, intent(inout)               :: lons(nx,ny),lats(nx,ny) !map of lon/lat idx of outlets for each cells on the 30s map
    integer, intent(in)                  :: rst_ocn_lnd(nx,ny) !map of tile idx from ocean-land rst map, eg: from M6TP0072x0036-Pfafstetter.rst
    integer, intent(out)                 :: ns !number of the outlets
    
    integer, allocatable, dimension(:)   :: lonp,latp
    integer, allocatable, dimension(:,:) :: acc,np_map
    integer                              :: i,j,k,l,lonc,latc,flag,maxbak,status,num
    
    allocate(acc(nx,ny))
    
    !It first masks out the outlets map with the ocean-land rst map...
    !If a cell is not defined as land or glacier in the ocean-land rst map, we ignore it.
    do i=1,nx
       do j=1,ny
          if(rst_ocn_lnd(i,j)>nl.and.rst_ocn_lnd(i,j)/=nt)then
             lons(i,j)=-999
             lats(i,j)=-999
          endif
       enddo
    enddo
    
    !Counting the outlets...
    acc=0
    k=0
    do i=1,nx
       do j=1,ny
          if(lons(i,j)/=-999.and.lats(i,j)/=-999)then
             lonc=lons(i,j)
             latc=lats(i,j)
             if(acc(lonc,latc)==0)then
                k=k+1
                acc(lonc,latc)=1
             endif
          endif
       enddo
    enddo
    ns=k
    deallocate(acc)

  end subroutine outlets_num

  !------------------------------------------------------------------------
  ! This subroutine retrives the outlets locations from the outlets map (nx,ny) to a list (ns)
  ! It also stores the outlets idx on the 30s map
  subroutine retrieve_outlets(lons,lats,lon30s,lat30s,lonp,latp,lon_lnd,lat_lnd,ns_map,nx,ny,ns)
    
    integer, intent(in)                 :: nx,ny,ns
    !nx: number of lon of 30s map, 43200
    !ny: number of lat of 30s map, 21600
    !ns: number of the outlets
    integer, intent(in)                 :: lons(nx,ny),lats(nx,ny) !input map of outlets idx
    real*8,  intent(in)                 :: lon30s(nx),lat30s(ny) !lon and lat value arrays of the 30s map 
    integer, intent(out)                :: lonp(ns),latp(ns) !list of lon and lat idx for the ns outlets 
    real*8,  intent(out)                :: lon_lnd(ns),lat_lnd(ns) !list of lon and lat value for the ns outlets 
    integer, intent(out)                :: ns_map(nx,ny) !It stores the outlets idx on the 30s map
    
    integer, allocatable,dimension(:,:) :: acc
    integer                             :: i,j,k,l,lonc,latc
    
    allocate(acc(nx,ny))
    ns_map=-9999
    acc=0
    k=0
    do i=1,nx
       do j=1,ny
          if(lons(i,j)/=-999.and.lats(i,j)/=-999)then
             lonc=lons(i,j)
             latc=lats(i,j)
             if(acc(lonc,latc)==0)then
                k=k+1
                acc(lonc,latc)=1
                lonp(k)=lonc
                latp(k)=latc
                ns_map(lonc,latc)=k
             endif
          endif
       enddo
    enddo
    
    do i=1,ns
       lon_lnd(i)=lon30s(lonp(i))
       lat_lnd(i)=lat30s(latp(i))
    enddo
    
    deallocate(acc)
    
  end subroutine retrieve_outlets

  !------------------------------------------------------------------------
  ! convert the ocean mask in MAPL_Tripolar.nc to 1d list for each ocean tile defined by ocn rst, eg: M6TP0072x0036.rst
  subroutine mask_MAPL_1d(msk_tile,t2loni,t2lati,nt,res_MAPL,nlon,nlat)
    
    integer,           intent(in)     :: nt,nlon,nlat
    !nt: number of ocean tiles
    !nlon: number of lon in MAPL_Tripolar.nc, eg 72
    !nlat: number of lat in MAPL_Tripolar.nc, eg 36   
    integer,           intent(in)     :: t2loni(nt),t2lati(nt) !relationship between ocn tile idx in M6TP0072x0036.rst and lat/lon idx in MAPL_Tripolar.nc
    character(len=*),  intent(in)     :: res_MAPL !name of the ocean resolution, eg M6TP0072x0036
    integer,           intent(out)    :: msk_tile(nt) !1d list for the ocean msk for each ocean tile

    real, allocatable, dimension(:,:) :: msk_MAPL
    integer                           :: i
    
    allocate(msk_MAPL(nlon,nlat))
    call read_oceanModel_mapl(res_MAPL,msk_MAPL,nlon,nlat)
    
    do i=1,nt
       msk_tile(i)=int(msk_MAPL(t2loni(i),t2lati(i)))
    enddo
    
    deallocate(msk_MAPL)

  end subroutine mask_MAPL_1d

  !------------------------------------------------------------------------
  ! read ocean mask from "MAPL_Tripolar.nc"
  subroutine read_oceanModel_mapl(res_MAPL,wetMask,nx,ny)
    

    character(len=*),    intent(in)  :: res_MAPL !ocn resolution name
    integer,             intent(in)  :: nx, ny !ocn resolution numbers
    real,                intent(out) :: wetMask(nx,ny) !ocn mask from MAPL_Tripolar.nc
    
    integer                          :: ncid, varid, ret
    character(len=4)                 :: subname="read"
        
    ! try MOM6 first
    
    ret=nf90_open( trim(MAKE_BCS_INPUT_DIR) // "/ocean/MOM6/" // trim(res_MAPL) // "/MAPL_Tripolar.nc", 0, ncid )
    
    ! if MOM6 did not work, try MOM5, 
    
    if(ret /= NF_NOERR)then
       call check_ret( nf90_open( trim(MAKE_BCS_INPUT_DIR) // "/ocean/MOM5/" // trim(res_MAPL) // "/MAPL_Tripolar.nc", 0, ncid ), subname)
    endif
    
    ! read "mask" from netcdf file into "wetMask"
    
    call check_ret(nf90_inq_varid(ncid,"mask",varid),subname)
    call check_ret(nf90_get_var(ncid,varid,wetMask),subname)
    call check_ret(nf90_close(ncid),subname)    
    
  end subroutine read_oceanModel_mapl

  !------------------------------------------------------------------------

  subroutine check_ret(ret, calling)

    integer, intent(in) :: ret
    character(len=*)    :: calling
    
    if (ret /= NF_NOERR) then 
       write(6,*)'netcdf error from ',trim(calling)
       call endrun(nf_strerror(ret))
    end if
    
  end subroutine check_ret

  !-----------------------------------------------------------------------  

  subroutine endrun(msg,subname)
    
    character(len=*), intent(in), optional :: msg        ! string to be printed
    character(len=*), intent(in), optional :: subname    ! subname
    
    if (present(subname)) then 
       write(6,*) 'ERROR in subroutine :', trim(subname)
    end if
    
    if (present(msg)) then
       write(6,*)'ENDRUN: ', msg
    else
       write(6,*) 'ENDRUN: called without a message string'
    end if
    
    stop 
    
  end subroutine endrun
  
  !------------------------------------------------------------------------
  ! convert 1d list of ocn mask to a 30s map
  subroutine mask_MAPL_2d(ocean,mask1d,msk2d,nt,nlon,nlat)
    
    integer, intent(in)               :: nt,nlon,nlat
    !nt: number of ocean tiles
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600   
    integer, intent(in)               :: ocean(nlon,nlat) !tile idx from ocn rst file, eg: from M6TP0072x0036.rst
    integer, intent(in)               :: mask1d(nt) !1d list of ocn mask got from subroutine mask_MAPL_1d
    integer, intent(out)              :: msk2d(nlon,nlat) !output of the ocn mask on the 30s map
    
    real*8,  allocatable,dimension(:) :: lon,lat
    integer                           :: i,j,xi,yi,tid
    
    do i=1,nlon
       do j=1,nlat
          tid=ocean(i,j)
          msk2d(i,j)=mask1d(tid)
       enddo
    enddo
    
  end subroutine mask_MAPL_2d

  !------------------------------------------------------------------------
  ! further mask the ocn mask from MAPL_Tripolar.nc based on the land-ocean rst eg: M6TP0072x0036-Pfafstetter.rst
  subroutine mask_MAPL_bcs(rst_ocn_lnd,mask_mapl,mask,nlon,nlat,nl,nt)
    
    integer,intent(in)  :: nlon,nlat,nl,nt
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600  
    !nl: number of lnd tile
    !nt: number of total tile
    integer,intent(in)  :: rst_ocn_lnd(nlon,nlat) !tile idx from land-ocean rst eg: M6TP0072x0036-Pfafstetter.rst
    integer,intent(in)  :: mask_mapl(nlon,nlat) !ocn mask map from subroutine mask_MAPL_2d
    integer,intent(out) :: mask(nlon,nlat) !ocn mask map masked further by land-ocean rst 
    
    mask=0
    !tile idx<=nl are land tiles
    !tile idx==nt are glacier tiles
    !tile idx==nt-1 are lake tiles
    !rest are all ocean tiles
    where(rst_ocn_lnd>nl.and.rst_ocn_lnd/=nt.and.rst_ocn_lnd/=nt-1.and.mask_mapl==1)mask=1
    
  end subroutine mask_MAPL_bcs

  !------------------------------------------------------------------------
  !find the ocean boundary cells (that are next to non-ocean cell) on the 30s map based on the ocn mask from mask_MAPL_bcs
  subroutine ocean_boundary(mask,boundary,nlon,nlat)

    integer, intent(in)  :: nlon,nlat
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600      
    integer, intent(in)  :: mask(nlon,nlat) !ocn mask from subroutine mask_MAPL_bcs
    integer, intent(out) :: boundary(nlon,nlat) !map of ocean boundary cells

    real*8,  allocatable :: lon(:),lat(:)
    integer              :: xi,yi,id
    integer              :: xp1,xm1,yp1,ym1
    
    boundary=mask
    boundary=-9999
    
    do xi=2,nlon-1
       do yi=2,nlat-1
          id=mask(xi,yi)
          if(id==1)then
             boundary(xi,yi)=0 
             if(mask(xi+1,yi)==1.and.&
                  mask(xi+1,yi-1)==1.and.&
                  mask(xi  ,yi-1)==1.and.&
                  mask(xi-1,yi-1)==1.and.&
                  mask(xi-1,yi)==1.and.&
                  mask(xi-1,yi+1)==1.and.&
                  mask(xi  ,yi+1)==1.and.&
                  mask(xi+1,yi+1)==1)then
                boundary(xi,yi)=-9999
             endif
          endif
       enddo
    enddo
    
  end subroutine ocean_boundary

  !------------------------------------------------------------------------
  ! counting the number of ocean boundary cells
  subroutine ocean_boundary_num(boundary,nlon,nlat,nsh)
    
    integer, intent(in)  :: nlon,nlat
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600          
    integer, intent(in)  :: boundary(nlon,nlat) !map of ocean boundary cells   
    integer, intent(out) :: nsh   !number of ocean boundary cells 
    
    integer              :: i,xi,yi,k
    
    k=0
    do xi=1,nlon
       do yi=1,nlat
          if(boundary(xi,yi)==0)then
             k=k+1
          endif
       enddo
    enddo
    nsh=k

  end subroutine ocean_boundary_num

  !------------------------------------------------------------------------
  ! list the lat and lon of ocean boundary cells
  subroutine ocean_boundary_points(boundary,lon30s,lat30s,lonsh,latsh,nlon,nlat,nsh)

    integer,intent(in) :: nlon,nlat,nsh
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600   
    !nsh: number of ocean boundary cells    
    integer,intent(in) :: boundary(nlon,nlat) !map of ocean boundary cells
    real*8,intent(in)  :: lon30s(nlon),lat30s(nlat) !lon and lat value arrays of the 30s map
    real*8,intent(out) :: lonsh(nsh),latsh(nsh) !lists of the lat and lon values of the ocean boundary cells
    integer i,xi,yi,k
    
    k=0
    do xi=1,nlon
       do yi=1,nlat
          if(boundary(xi,yi)==0)then
             k=k+1
             lonsh(k)=lon30s(xi)
             latsh(k)=lat30s(yi)
          endif
       enddo
    enddo
  end subroutine ocean_boundary_points

  !------------------------------------------------------------------------
  ! move the outlet locations to the nearest ocean boundary cell 
  subroutine move_to_ocean(lonsi,latsi,lons,lats,mask,lonsh,latsh,lons_adj,lats_adj,ns,nlon,nlat,nsh)

    integer, intent(in)  :: ns,nlon,nlat,nsh
    !ns: number of the outlets
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600   
    !nsh: number of ocean boundary cells        
    integer, intent(in)  :: lonsi(ns),latsi(ns) !lon and lat idx of the outlets before moving
    real*8,  intent(in)  :: lons(ns),lats(ns) !lon and lat values of the outlets before moving
    integer, intent(in)  :: mask(nlon,nlat) !!ocn mask from subroutine mask_MAPL_bcs
    real*8,  intent(in)  :: lonsh(nsh),latsh(nsh) !lists of the lat and lon values of the ocean boundary cells
    real*8,  intent(out) :: lons_adj(ns),lats_adj(ns) !lon and lat values of the outlets after moving
    
    real,allocatable :: dist(:)
    
    integer :: i,j
    real :: dy,dy2,dx,dx2,dxA,dxB,dist_temp 
    
    allocate(dist(ns))
    do i=1,ns
       IF(mask(lonsi(i),latsi(i))==0)THEN
          dist(i)=1.e12
          do j=1,nsh
             dy=abs(lats(i)-latsh(j))
             dy2=dy*dy
             dxA=abs(lons(i)-lonsh(j)) 
             dxB=360.-dxA
             dx=min(dxA,dxB)
             dx2=dx*dx
             dist_temp=sqrt(dx2+dy2)
             if(dist_temp<dist(i))then
                dist(i)=dist_temp
                lons_adj(i)=lonsh(j)
                lats_adj(i)=latsh(j)
             endif
          enddo
       ELSE
          lons_adj(i)=lons(i)
          lats_adj(i)=lats(i)
          dist(i)=0.
       ENDIF
    enddo
    deallocate(dist)
    
  end subroutine move_to_ocean

  !------------------------------------------------------------------------
  !convert the lon and lat values of the outlets to lon and lat idx on the 30s map
  subroutine sinkxy_ocean(lons,lats,lon30s,lat30s,loni,lati,ns,nlon,nlat)
    
    integer, intent(in)                :: ns,nlon,nlat
    !ns: number of the outlets
    !nlon: number of lon in 30s map, eg 43200
    !nlat: number of lat in 30s map, eg 21600       
    real*8,  intent(in)                :: lons(ns),lats(ns) !lon and lat values of the outlets
    real*8,  intent(in)                :: lon30s(nlon),lat30s(nlat) !lon and lat value arrays of the 30s map
    integer, intent(out)               :: loni(ns),lati(ns) !lon and lat idx of the outlets

    real*8,  allocatable, dimension(:) :: lat_dis,lon_dis
    integer                            :: i,temp(1)
    
    allocate(lat_dis(nlat),lon_dis(nlon))
    
    do i=1,ns
       lat_dis=abs(lat30s-lats(i))
       temp=minloc(lat_dis)
       lati(i)=temp(1)
    enddo
    do i=1,ns
       lon_dis=abs(lon30s-lons(i))
       temp=minloc(lon_dis)
       loni(i)=temp(1)
    enddo
    
    deallocate(lat_dis,lon_dis)
    
  end subroutine sinkxy_ocean

  !------------------------------------------------------------------------
  ! put the list of lon and lat idx of the outlets back to the 30s map lons(nx,ny),lats(nx,ny)
  subroutine update_outlets(loni_ocn,lati_ocn,ns_map,lons,lats,nx,ny,ns)

    integer,intent(in)    :: nx,ny,ns
    !nx: number of lon in 30s map, eg 43200
    !ny: number of lat in 30s map, eg 21600  
    !ns: number of the outlets
    integer,intent(in)    :: loni_ocn(ns),lati_ocn(ns) !lon and lat idx of the outlets
    integer,intent(in)    :: ns_map(nx,ny) !it stores the outlets idx on the 30s map
    integer,intent(inout) :: lons(nx,ny),lats(nx,ny) !map of lon/lat idx of outlets for each land/Greenland cells 
    
    integer               :: i,j,lonc,latc,ind
    
    do i=1,nx
       do j=1,ny
          if(lons(i,j)/=-999.and.lats(i,j)/=-999.)then
             
             lonc=lons(i,j)
             latc=lats(i,j)      
             ind=ns_map(lonc,latc)
             if(ind<1.or.ind>ns)then
                print *,"ns_map is Incorrect, ind=",ind
                stop
             endif
             lons(i,j)=loni_ocn(ind)
             lats(i,j)=lati_ocn(ind)
             
          endif
       enddo
    enddo
    
  end subroutine update_outlets

  !------------------------------------------------------------------------

end program Runoff

! ============================ EOF =====================================================
