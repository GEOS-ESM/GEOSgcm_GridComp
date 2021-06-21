#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

PROGRAM mkIrrigTiles

  use MAPL
  use MAPL_ConstantsMod,               ONLY:   &
       MAPL_RADIUS, MAPL_PI               
  use gFTL_StringVector  
  use pFIO
  use rmTinyCatchParaMod,              ONLY: NC_VarID, RegridRaster
  use CubedSphere_GridMod
  use easeV2_conv
  use MAPL_SortMod
  
  implicit none
  include 'netcdf.inc'
  character*40                    :: BCSNAME, TILFILE, IMxJM
  character*1                     :: opt
  integer                         :: n, iargc, NXT, ITILES, NTILES, RC, NC_RASTER, NR_RASTER
  character*256                   :: arg
  real                            :: thresh
  integer, pointer, dimension (:) :: tile_id, irr_crp_tile, irr_pad
  character*200                   :: filename
  type(Netcdf4_FileFormatter)     :: InFmt
  
  ! Read arguments
  
  nxt = 1
  call getarg(nxt,arg)

  do while(arg(1:1)=='-')

     opt=arg(2:2)
     if(len(trim(arg))==2) then
        nxt = nxt + 1
        call getarg(nxt,arg)
     else
        arg = arg(3:)
     end if

     select case (opt)
     case ('t')
        TILFILE = trim(arg)
     case ('b')
        BCSNAME = trim(arg)
     case ('r')
        IMxJM   = trim(arg)
     case ('x')
        read(arg,'(i6)') nc_raster
     case ('y')
        read(arg,'(i6)') nr_raster
     case ('p')
        read (arg,* ) RC
        thresh = RC / 100.
     case default
        print *, "USAGE   : bin/mkIrrigTiles.x -b BCSNAME -t TILFILE -r IMxJM -p thresh%"
        call exit(1)
     end select

     nxt = nxt + 1
     call getarg(nxt,arg)

  end do

  call create_irrig_mask
  call write_raster_file('irrigation'//trim(IMxJM)//'.dat' )
  call write_tilfile    ('irrigation'//trim(IMxJM)//'.dat' )
  call write_clim_files ('green_clim'//trim(IMxJM)//'.data')
  call write_clim_files ('lai_clim'//trim(IMxJM)//'.data'  )
  call write_clim_files ('lnfm_clim'//trim(IMxJM)//'.data' )
  call write_clim_files ('ndvi_clim'//trim(IMxJM)//'.data' )
  call write_clim_files ('visdf'//trim(IMxJM)//'.dat'      )
  call write_clim_files ('nirdf'//trim(IMxJM)//'.dat'      )

  call write_nc4_files   ('clsm/catch_params.nc4'  )
  call write_nc4_files   ('clsm/catchcn_params.nc4')
  call write_nc4_files   ('irrigation'//trim(IMxJM)//'.dat')
  call write_nc4_files   ('vegdyn'//trim(IMxJM)//'.dat')
  call reset_irrig_fracs ('IRRIG/'//trim(BCSNAME)//'/irrigation'//trim(IMxJM)//'.dat')

  call write_tables
  
  
contains

  SUBROUTINE create_irrig_mask
    
    implicit none
    integer, allocatable, dimension (:) :: tile_id_tmp, irr_not_tmp, mos_typ, irr_pad_tmp
    real, allocatable, dimension (:)    :: IRRIGFRAC, PADDYFRAC
    real,    allocatable, dimension (:) :: tile_lon, tile_lat
    real                                :: minlon,maxlon,minlat,maxlat
    integer                             :: tindex1,pfaf1, v
    
    filename = trim(BCSNAME)//'/clsm/catchment.def'
    open (10, file = trim (filename), form ='formatted',action = 'read', status = 'old')
    open (11, file = trim(BCSNAME)//'/clsm/mosaic_veg_typs_fracs', &
         form ='formatted',action = 'read', status = 'old')
    read (10, *) NTILES    
    allocate (tile_lon (1:NTILES))
    allocate (tile_lat (1:NTILES))
    allocate (mos_typ  (1:NTILES))
    do n = 1, NTILES
       read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
       tile_lon(n) = (minlon + maxlon)/2.
       tile_lat(n) = (minlat + maxlat)/2.
       read (11,*) tindex1,pfaf1,mos_typ(n)
    end do
    close(10, status = 'keep')
    close(11, status = 'keep')

    allocate (tile_id_tmp (1:2*NTILES))
    allocate (irr_not_tmp (1:2*NTILES))
    allocate (irr_pad_tmp (1:2*NTILES))    
    allocate (IRRIGFRAC   (1:NTILES))
    allocate (PADDYFRAC   (1:NTILES))
    
    filename =  trim(BCSNAME)//'/irrigation'//trim(IMxJM)//'.dat'
    call InFmt%open(trim(filename), pFIO_READ,rc=rc) ; VERIFY_(RC)
    call MAPL_VarRead (InFmt,'IRRIGFRAC',IRRIGFRAC, rc=rc) ; VERIFY_(RC)
    call MAPL_VarRead (InFmt,'PADDYFRAC',PADDYFRAC, rc=rc) ; VERIFY_(RC)
    call inFmt%close(rc=rc)  ; VERIFY_(RC)

    filename = 'IRRIG/'//trim(BCSNAME)//'/clsm/irrig2bcs_mapping'
    open (10, file = trim (filename),form ='formatted',action = 'write', status = 'unknown') 
    tile_id_tmp = 0
    irr_not_tmp = 0
    irr_pad_tmp = 0  
    ITILES      = 1
    do n = 1, NTILES
       tile_id_tmp (ITILES) = n
       write (10,*) ITILES, n, 0, 0
       if((IRRIGFRAC(n) + PADDYFRAC(n)) >= thresh) then
          v = grassNeighbor (n, mos_typ, tile_lat, tile_lon)
          if(IRRIGFRAC(n) > 0) then
             ITILES = ITILES + 1
             tile_id_tmp (ITILES) = n          
             irr_not_tmp (ITILES) = v
             irr_pad_tmp (ITILES) = 1
             write (10,*) ITILES, n, v, irr_pad_tmp (ITILES)
          endif
          if(PADDYFRAC(n) > 0) then
             ITILES = ITILES + 1
             tile_id_tmp (ITILES) = n          
             irr_not_tmp (ITILES) = v
             irr_pad_tmp (ITILES) = irr_pad_tmp (ITILES-1) + 2
             write (10,*) ITILES, n, v, irr_pad_tmp (ITILES)             
          endif
       endif
       ITILES = ITILES + 1
    end do
    close(10, status = 'keep')
    
    ITILES = ITILES - 1
    allocate (tile_id      (1:ITILES))
    allocate (irr_crp_tile (1:ITILES))
    allocate (irr_pad      (1:ITILES))
    tile_id     (:) = tile_id_tmp (1:ITILES)
    irr_crp_tile(:) = irr_not_tmp (1:ITILES)
    irr_pad     (:) = irr_pad_tmp (1:ITILES)
    deallocate (tile_id_tmp, irr_not_tmp, irr_pad_tmp)
    deallocate (tile_lon, tile_lat, mos_typ)

  END SUBROUTINE create_irrig_mask

  ! ------------------------------------------------------------------------------

  integer function grassNeighbor (tid_in, mos_typ, tile_lat, tile_lon)
     
    implicit none
    real, dimension (:), intent(in)     :: tile_lat, tile_lon
    integer, intent (in)                :: tid_in
    integer, dimension (NTILES)         :: mos_typ
    integer                             :: i, nplus    
    logical                             :: tile_found
    logical, allocatable, dimension (:) :: mask
    integer, allocatable, dimension (:) :: sub_tid
    real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist
    real                                :: dw, min_lon, max_lon, min_lat, max_lat
    integer, allocatable, dimension (:) :: TILEID
      
    allocate (mask   (1:  NTILES))
    allocate (TILEID (1:  NTILES))
    forall (i=1:NTILES) TILEID (i) = i
    
    dw = 0.5
    grassNeighbor = -9999
    
    ZOOMOUT : do  
       
       tile_found = .false. 
       
       ! Min/Max lon/lat of the working window
       ! -------------------------------------
       
       min_lon = MAX(tile_lon (tid_in) - dw, -180.)
       max_lon = MIN(tile_lon (tid_in) + dw,  180.)
       min_lat = MAX(tile_lat (tid_in) - dw,  -90.)
       max_lat = MIN(tile_lat (tid_in) + dw,   90.) 
       
       mask = .false.
       mask =  ((tile_lat >= min_lat .and. tile_lat <= max_lat).and.(tile_lon >= min_lon .and. tile_lon <= max_lon).and.(mos_typ == 4))
       nplus =  count(mask = mask)
       
       if(nplus < 0) then
          dw = dw + 0.5
          CYCLE
       endif
       
       allocate (sub_tid (1:nplus))
       allocate (sub_lon (1:nplus))
       allocate (sub_lat (1:nplus))
       allocate (rev_dist  (1:nplus))
       
       sub_tid = PACK (TILEID  , mask= mask) 
       sub_lon = PACK (tile_lon, mask= mask)
       sub_lat = PACK (tile_lat, mask= mask)
         
       ! compute distance from the tile
       
       sub_lat = sub_lat * MAPL_PI/180.
       sub_lon = sub_lon * MAPL_PI/180.
       
       SEEK : if(grassNeighbor < 0) then
          
          rev_dist  = 1.e20
          
          do i = 1,nplus
             
             rev_dist(i) = haversine(to_radian(tile_lat(tid_in)), to_radian(tile_lon(tid_in)), &
                  sub_lat(i), sub_lon(i))
             
          end do
          
          FOUND : if(minval (rev_dist) < 1.e19) then
             if(mos_typ(sub_tid(minloc(rev_dist,1))) == 4) then
                grassNeighbor = sub_tid(minloc(rev_dist,1)) 
                tile_found = .true.
             endif
          endif FOUND
          
       endif SEEK
       
       deallocate (sub_tid, sub_lon, sub_lat, rev_dist)
       
       if(tile_found) GO TO 100
       
       ! if not increase the window size
       dw = dw + 0.5
       
    end do ZOOMOUT
    
100 continue
    
    deallocate (mask)
    
  end function grassNeighbor
  
  ! *****************************************************************************
  
  function to_radian(degree) result(rad)
    
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
  end function haversine
  
  ! ------------------------------------------------------------------------------

  SUBROUTINE reset_irrig_fracs (infile)

    implicit none
    character(*), intent (in)        :: infile
    real, allocatable, dimension (:) :: PADDYFRAC, IRRIGFRAC
    real, allocatable, dimension (:) :: CROPIRRIGFRAC
    type(FileMetadata)               :: Cfg
    integer                          :: j,dim1
    type(Variable), pointer          :: myVariable
    character(len=:), pointer        :: dname

    print *, 'RESET ', trim (infile)
    
    allocate (PADDYFRAC (1:ITILES))
    allocate (IRRIGFRAC (1:ITILES))
    allocate (CROPIRRIGFRAC (1:ITILES))

    call InFmt%open(trim(infile), pFIO_WRITE,rc=rc) ; VERIFY_(RC)
    Cfg  = InFmt%read(rc=rc)
    call MAPL_VarRead (InFmt,'PADDYFRAC', PADDYFRAC,rc=rc); VERIFY_(RC)
    call MAPL_VarRead (InFmt,'IRRIGFRAC', IRRIGFRAC,rc=rc); VERIFY_(RC)   

    myVariable => cfg%get_variable("CROPIRRIGFRAC")
    dname => myVariable%get_ith_dimension(2)
    dim1 = cfg%get_dimension(dname)

    do j = 1, dim1
       call MAPL_VarRead ( InFmt,"CROPIRRIGFRAC",CROPIRRIGFRAC,offset1=j,rc=rc) ; VERIFY_(RC)
       if (j ==3) then
          ! Paddy
          where (irr_pad > 1)
             CROPIRRIGFRAC = 1.
          elsewhere
             CROPIRRIGFRAC = 0.
          endwhere
       else
          ! Irrig crops
          where (irr_pad == 1)
             CROPIRRIGFRAC = CROPIRRIGFRAC/IRRIGFRAC
          elsewhere
             CROPIRRIGFRAC = 0.
          endwhere          
       endif
       call MAPL_VarWrite (InFmt,"CROPIRRIGFRAC",CROPIRRIGFRAC,offset1=j,rc=rc) ; VERIFY_(RC)
    end do

    where (irr_pad == 1)
       ! irrigated crop tiles
       IRRIGFRAC = 1.
       PADDYFRAC = 0.
    endwhere
    where (irr_pad > 1)
       ! paddy tiles
       PADDYFRAC = 1.
       IRRIGFRAC = 0.
    endwhere

    call MAPL_VarWrite (InFmt,'PADDYFRAC', PADDYFRAC)
    call MAPL_VarWrite (InFmt,'IRRIGFRAC', IRRIGFRAC)

    call inFmt%close (rc=rc)
    
  END SUBROUTINE reset_irrig_fracs
  
  ! ------------------------------------------------------------------------------

  SUBROUTINE write_nc4_files (infile)

    implicit none
    character(*), intent (in)        :: infile
    type(Netcdf4_FileFormatter)      :: OutFmt
    type(FileMetadata)               :: InCfg,OutCfg
    integer                          :: dim1,dim2, ndims, i, j
    type(StringVariableMap), pointer :: variables
    type(Variable), pointer          :: var
    type(StringVariableMapIterator)  :: var_iter
    type(StringVector), pointer      :: var_dimensions
    character(len=:), pointer        :: vname,dname
    logical                          :: irrg, catch, vegdyn
    real,allocatable                 :: tmp_var(:),  tmp_var2(:)

    print *, trim (infile)
    
    call InFmt%open(trim(BCSNAME)//'/'//trim(infile), pFIO_READ,rc=rc) ; VERIFY_(RC)
    InCfg  = InFmt%read(rc=rc)                           ; VERIFY_(RC)
    OutCfg = InCfg
    call OutCfg%modify_dimension('tile', ITILES, rc=rc)  ; VERIFY_(RC)
    call OutFmt%create('IRRIG/'//trim(BCSNAME)//'/'//trim(infile),rc=rc)
    call OutFmt%write(OutCfg,rc=rc)

    allocate (tmp_var (1: NTILES))
    allocate (tmp_var2(1: ITILES))
    variables => InCfg%get_variables()
    var_iter = variables%begin()
    
    do while (var_iter /= variables%end())
          
       vname => var_iter%key()
       var => var_iter%value()
       var_dimensions => var%get_dimensions()
          
       ndims = var_dimensions%size()
       if (trim(vname) =='CROPCLASSNAME') then
          call var_iter%next()
          cycle
       endif
          
       if (ndims == 1) then
          
          call MAPL_VarRead (InFmt,vname,tmp_var, rc=rc) ; VERIFY_(RC)
          tmp_var2 = tmp_var(tile_id)
          call MAPL_VarWrite(OutFmt,vname,tmp_var2, rc = rc) ; VERIFY_(RC)

          if(index(infile,"vegdyn")/=0) then
             if(trim(vname) == 'ITY') then
                tmp_var2 = tmp_var(tile_id)
                print *, 'Switching ITY to grass (4) on irrig tiles in : ', trim (infile)
                where (irr_crp_tile > 0)
                   tmp_var2 = 4.
                end where
                call MAPL_VarWrite(OutFmt,vname,tmp_var2)
             endif
             
             if(trim(vname) == 'Z2CH') then
                tmp_var2 = tmp_var(tile_id)
                print *, 'Setting Z2CH to 0.6 on irrig tiles in : ', trim (infile)
                where (irr_crp_tile > 0)
                   tmp_var2 = 0.6
                end where
                call MAPL_VarWrite(OutFmt,vname,tmp_var2)
             endif
          endif
          
          if(index(infile,"catch_")/=0) then
             if(trim(vname) == 'OLD_ITY') then
                tmp_var2 = tmp_var(tile_id)
                print *, 'Switching OLD_ITY to grass (4) on irrig tiles in : ', trim (infile)
                where (irr_crp_tile > 0)
                   tmp_var2 = 4.
                end where
                call MAPL_VarWrite(OutFmt,vname,tmp_var2)
             endif
             
             if(trim(vname) == 'BF3') then
                tmp_var2 = tmp_var(tile_id)
                print *, 'Switching BF3 to 25 on paddy tiles in : ', trim (infile)
                where (irr_pad >= 1)
                   tmp_var2 = 25.
                end where
                call MAPL_VarWrite(OutFmt,vname,tmp_var2)
             endif
          endif
          
          if(index(infile,"irrigation")/=0) then
             if((trim(vname) == 'PADDYFRAC').OR.(trim(vname) == 'IRRIGFRAC')) then
                tmp_var2 = tmp_var(tile_id)
                print *, 'Ensure PADDYFRAC/IRRIGFRAC on non-irrigated tiles are zeros in : ', trim (infile)
                where (irr_crp_tile == 0)
                   tmp_var2 = 0.
                end where
                call MAPL_VarWrite(OutFmt,vname,tmp_var2)
             endif
          endif
          
       else if (ndims == 2) then
             
          dname => var%get_ith_dimension(2)
          dim1=InCfg%get_dimension(dname)
          do j=1,dim1
             
             call MAPL_VarRead ( InFmt,vname,tmp_var,offset1=j, rc=rc) ; VERIFY_(RC)
             call MAPL_VarWrite(OutFmt,vname,tmp_var(tile_id),offset1=j, rc=rc) ; VERIFY_(RC)

             if(index(infile,"catchcn_")/=0) then
                if(trim(vname) == 'ITY') then
                   if (j == 1) print *, 'Switching ITY to 18 19 16 17 on irrig tiles in : ', trim (infile)
                   tmp_var2 = tmp_var(tile_id)
                   if(j == 1) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 18.
                      end where
                   endif
                   if(j == 2) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 19.
                      end where
                   endif
                   if(j == 3) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 16.
                      end where
                   endif
                   if(j == 4) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 17.
                      end where
                   endif
                   call MAPL_VarWrite(OutFmt,vname,tmp_var2,offset1=j)
                endif
                
                if(trim(vname) == 'FVG') then
                   if (j == 1)print *, 'Switching FVG to 0. 0.95 0. 0.05 on irrig tiles in : ', trim (infile)
                   tmp_var2 = tmp_var(tile_id)
                   if(j == 1) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 0.
                      end where
                   endif
                   if(j == 2) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 0.95
                      end where
                   endif
                   if(j == 3) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 0.
                      end where
                   endif
                   if(j == 4) then
                      where (irr_crp_tile > 0)
                         tmp_var2 = 0.05
                      end where
                   endif
                   call MAPL_VarWrite(OutFmt,vname,tmp_var2,offset1=j)
                endif
             endif

             if(index(infile,"irrigation_")/=0) then
                if(trim(vname) == 'CROPIRRIGFRAC') then
                   if (j == 1) print *, 'Set CROPIRRIGFRAC  on non-irrigated tiles to zeros in : ', trim (infile)
                   tmp_var2 = tmp_var(tile_id)
                   where (irr_crp_tile == 0)
                      tmp_var2 = 0.
                   end where
                   call MAPL_VarWrite(OutFmt,vname,tmp_var2,offset1=j)
                endif
             endif
             
          enddo
             
       else if (ndims == 3) then
             
          dname => var%get_ith_dimension(2)
          dim1=InCfg%get_dimension(dname)
          dname => var%get_ith_dimension(3)
          dim2=InCfg%get_dimension(dname)
          do i=1,dim2
             do j=1,dim1
                call MAPL_VarRead ( InFmt,vname,tmp_var   ,offset1=j,offset2=i, rc=rc) ; VERIFY_(RC)
                call MAPL_VarWrite(OutFmt,vname,tmp_var(tile_id),offset1=j,offset2=i, rc=rc) ; VERIFY_(RC)
             enddo
          enddo          
       end if
       call var_iter%next()
    enddo

    call inFmt%close (rc=rc)
    call OutFmt%close(rc=rc)
    deallocate (tmp_var, tmp_var2)
    
  END SUBROUTINE write_nc4_files

  ! ------------------------------------------------------------------------------
  
  SUBROUTINE write_clim_files (infile)

    implicit none
    character(*), intent (in)        :: infile
    real, dimension (:), allocatable :: vec_in
    real    :: yr,mn,dd,hh,mm,ss,yr1,mn1,dd1,hh1,mm1,ss1,maxcat,nx
    integer :: istat, n

    print *, trim (infile)

    allocate (vec_in (1: NTILES))
    
    open (10, FILE = trim(BCSNAME)//'/'//trim(infile),form = 'unformatted', status = 'old', action = 'read')
    open (20, file = 'IRRIG/'//trim(BCSNAME)//'/'//trim(infile), form = 'unformatted', status = 'unknown', action = 'write')
    do n = 1, 100
       read(10,iostat=istat) yr,mn,dd,hh,mm,ss,yr1,mn1,dd1,hh1,mm1,ss1,maxcat,nx 
       if(.not. IS_IOSTAT_END(istat)) then
          write(20) yr,mn,dd,hh,mm,ss,yr1,mn1,dd1,hh1,mm1,ss1,real(itiles),nx
          read (10) vec_in
          write(20) vec_in (tile_id)
       endif
    end  do
    close (10, status = 'keep')
    close (20, status = 'keep')
    deallocate (vec_in)
    
  END SUBROUTINE write_clim_files

  ! ------------------------------------------------------------------------------

  SUBROUTINE write_tables

    implicit none

    integer, dimension (:,:), allocatable :: iRtable, iWtable
    real,    dimension (:,:), allocatable :: rRtable, rWtable
    integer :: i, n, icol, rcol
    real,allocatable                      :: ITY(:), BF3(:)

    allocate (ITY (1: ITILES))
    call InFmt%open('IRRIG/'//trim(BCSNAME)//'/clsm/catch_params.nc4', pFIO_READ,rc=rc) ; VERIFY_(RC)
    call MAPL_VarRead (InFmt,'OLD_ITY',ITY, rc=rc) ; VERIFY_(RC)
    call inFmt%close (rc=rc)
    
    ! (1) catchment.def
    ! -----------------
    icol = 2
    rcol = 5
    open (10, FILE = trim(BCSNAME)//'/clsm/catchment.def', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/catchment.def', form = 'formatted', status = 'unknown', action = 'write')
    read  (10,*) n
    write (11,*) ITILES
    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))

    do n = 1, NTILES
        read (10,'(i10,i8,5(2x,f9.4))') iRtable(n,:),rRtable(n,:)
    end do
    
    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do

    do n = 1, ITILES
       write (11,'(i10,i8,5(2x,f9.4))') n, iWtable(n,2),rWtable(n,:)
    end do
    
    deallocate (iRtable, iWtable,rRtable, rWtable)
    close (10, status = 'keep')
    close (11, status = 'keep')
    
    ! (2) soil_param.dat
    ! ------------------
    icol = 4
    rcol = 16
    open (10, FILE = trim(BCSNAME)//'/clsm/soil_param.dat', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/soil_param.dat', form = 'formatted', status = 'unknown', action = 'write')

    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))

    do n = 1, NTILES
       read (10,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4, f8.4)')iRtable(n,:),rRtable(n,:)
    end do

    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do

    do n = 1, ITILES
       write (11,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4, f8.4)')n,iWtable (n,2:) ,rWtable(n,:)
    end do
    
    deallocate (iRtable, iWtable,rRtable, rWtable)    
    close (10, status = 'keep')
    close (11, status = 'keep')
    
    ! (3) ar.new
    ! ----------
    icol = 2
    rcol = 12
    open (10, FILE = trim(BCSNAME)//'/clsm/ar.new', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/ar.new', form = 'formatted', status = 'unknown', action = 'write')

    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))

    do n = 1, NTILES
       read (10,'(i10,i8,f5.2,11(2x,e14.7))') iRtable(n,:),rRtable(n,:)
    end do

    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do
    
    do n = 1, ITILES
       write (11,'(i10,i8,f5.2,11(2x,e14.7))')n,iWtable (n,2) ,rWtable(n,:)
    end do
    
    deallocate (iRtable, iWtable,rRtable, rWtable)
    close (10, status = 'keep')
    close (11, status = 'keep')
    
    ! (4) bf.dat
    ! ----------
    icol = 2
    rcol = 4
    open (10, FILE = trim(BCSNAME)//'/clsm/bf.dat', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/bf.dat', form = 'formatted', status = 'unknown', action = 'write')
    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))
    allocate (BF3     (1:ITILES))
    
    do n = 1, NTILES
       read (10,'(i10,i8,f5.2,3(2x,e13.7))') iRtable(n,:),rRtable(n,:)
    end do

    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do
    
    BF3 = rWtable (:,rcol)
    where (irr_pad >= 1)
       BF3 = 25.
    end where
    rWtable (:,rcol) = BF3
    do n = 1, ITILES
       write (11,'(i10,i8,f5.2,3(2x,e13.7))') n, iWtable(n,2),rWtable(n,:)
    end do
    
    deallocate (iRtable, iWtable,rRtable, rWtable, BF3)
    close (10, status = 'keep')
    close (11, status = 'keep')
    
    ! (5) ts.dat
    ! ----------
    icol = 2
    rcol = 5
    open (10, FILE = trim(BCSNAME)//'/clsm/ts.dat', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/ts.dat', form = 'formatted', status = 'unknown', action = 'write')
    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))

    do n = 1, NTILES
       read (10,'(i10,i8,f5.2,4(2x,e13.7))') iRtable(n,:),rRtable(n,:)
    end do

    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do
    
    do n = 1, ITILES
       write (11,'(i10,i8,f5.2,4(2x,e13.7))')n, iWtable(n,2),rWtable(n,:)
    end do
    deallocate (iRtable, iWtable,rRtable, rWtable)
    close (10, status = 'keep')
    close (11, status = 'keep')
    
    ! (6) tau_param.dat
    ! -----------------
    icol = 2
    rcol = 4
    open (10, FILE = trim(BCSNAME)//'/clsm/tau_param.dat', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/tau_param.dat', form = 'formatted', status = 'unknown', action = 'write')
    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))

    do n = 1, NTILES
       read (10,'(i10,i8,4f10.7)')iRtable(n,:),rRtable(n,:)
    end do

    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do

    do n = 1, ITILES
       write (11,'(i10,i8,4f10.7)')n, iWtable(n,2),rWtable(n,:)
    end do
    
    deallocate (iRtable, iWtable,rRtable, rWtable)
    close (10, status = 'keep')
    close (11, status = 'keep')
    
    ! (7) mosaic_veg_typs_fracs
    ! -------------------------
    icol = 4
    rcol = 4
    open (10, FILE = trim(BCSNAME)//'/clsm/mosaic_veg_typs_fracs', form = 'formatted', status = 'old', action = 'read')
    open (11, FILE = 'IRRIG/'//trim(BCSNAME)//'/clsm/mosaic_veg_typs_fracs', form = 'formatted', status = 'unknown', action = 'write')
    allocate (iRtable (1:NTILES, 1:icol))
    allocate (iWtable (1:ITILES, 1:icol))
    allocate (rRtable (1:NTILES, 1:rcol))
    allocate (rWtable (1:ITILES, 1:rcol))

    do n = 1, NTILES
        read (10,'(i10,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')iRtable(n,:),rRtable(n,:)
    end do

    do i = 1, icol
       iWtable (:,i) = iRtable (tile_id,i) 
    end do
    do i = 1, rcol
       rWtable (:,i) = rRtable (tile_id,i) 
    end do

    do n = 1, ITILES
       write (11,'(i10,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')n, iWtable(n,2),INT(ITY(n)),iWtable(n,4),rWtable(n,:)
    end do
    
    deallocate (iRtable, iWtable,rRtable, rWtable)
    close (10, status = 'keep')
    close (11, status = 'keep')
    
  END SUBROUTINE write_tables

  ! ------------------------------------------------------------------------------

  SUBROUTINE write_tilfile (infile)

    implicit none

    character(*), intent (in) :: infile
    integer      :: NT, NPF, NC, NR, NT_NEW,NG, IDUM, i, N, icol, rcol, non_land
    character*20 :: cdum
    integer, dimension (:,:), allocatable :: iRtable, iWtable
    real,    dimension (:,:), allocatable :: rRtable, rWtable
    integer*8, allocatable, dimension (:) :: rSRTM, wSRTM
    real, allocatable, dimension (:) :: IRRIGFRAC, PADDYFRAC
    
    call InFmt%open(trim(BCSNAME)//'/'//trim(infile), pFIO_READ,rc=rc) ; VERIFY_(RC)
    allocate (IRRIGFRAC (1:NTILES))
    allocate (PADDYFRAC (1:NTILES))
    call MAPL_VarRead (InFmt,'PADDYFRAC', PADDYFRAC, rc=rc) ; VERIFY_(RC)
    call MAPL_VarRead (InFmt,'IRRIGFRAC', IRRIGFRAC, rc=rc) ; VERIFY_(RC)
    call inFmt%close (rc=rc)
    
    open (10,file =  trim(BCSNAME)//'/'//trim(tilfile), form = 'formatted', action = 'read', status = 'old')
    open (11,file = 'IRRIG/'//trim(BCSNAME)//'/'//trim(tilfile), form = 'formatted', action = 'write', status = 'unknown')
    
    read (10, *) NT, NPF, NC, NR
    NT_NEW   = NT + ITILES - NTILES
    non_land = NT - NTILES
    write (11,'(4I10)')NT_NEW, NPF, NC, NR
    read (10, *) NG
    write(11, *) NG
    
    do n = 1, NG
       read (10, '(a)') cdum
       write(11, '(a)') trim (cdum)
       read (10, *) IDUM
       write(11, '(I10)') IDUM
       read (10, *) IDUM
       write(11, '(I10)') IDUM
    end do
        
    if(index(tilfile,"EASE")/=0) then
       
       icol = 5
       rcol = 3
       allocate (iRtable (1:NTILES, 1:icol))
       allocate (iWtable (1:ITILES, 1:icol))
       allocate (rRtable (1:NTILES, 1:rcol))
       allocate (rWtable (1:ITILES, 1:rcol))
       allocate (rSRTM   (1:NTILES))
       allocate (wSRTM   (1:ITILES))
       
       do n = 1,  ntiles
          read (10,'(i10,i9,2f10.4,2i5,f19.12,i10,i15,e13.4)') iRtable (n,1),iRtable (n,2),rRtable(n,1),rRtable(n,2),iRtable (n,3),iRtable (n,4),rRtable(n,3),iRtable (n,5),rSRTM(n)
       end do
       
       do i = 1, icol
          iWtable (:,i) = iRtable (tile_id,i) 
       end do
       do i = 1, rcol
          rWtable (:,i) = rRtable (tile_id,i) 
       end do
       wSRTM = rSRTM(tile_id)
       ! rWtable(n,3) - split cell fractions
       do n = 1, ITILES
          if( irr_crp_tile(n) > 0) then
             ! either irrig crop or paddy thus non-irrig tile fraction 
             if(irr_pad(n) <= 2) rWtable(n-1,3) = (1. - IRRIGFRAC(tile_id(n)) - PADDYFRAC(tile_id(n))) *  rWtable(n-1,3)
            
             ! irrigated crop tile
             if(irr_pad(n) == 1) rWtable(n,3)   = IRRIGFRAC(tile_id(n)) * rWtable(n,3)

             ! irrigated paddy tile
             if(irr_pad(n) >= 2) rWtable(n,3)   = PADDYFRAC(tile_id(n)) * rWtable(n,3)
          endif
       end do
       do n = 1, ITILES
          write (11,'(i10,i9,2f10.4,2i5,f19.12,i10,i15,e13.4)') iWtable (n,1),iWtable (n,2),rWtable(n,1),rWtable(n,2),iWtable (n,3),iWtable (n,4),rWtable(n,3),iWtable (n,5),wSRTM(n)
       end do
       do n = 1, non_land
          read (10,'(i10,i9,2f10.4,2i5,f19.12,i10,i15,e13.4)') iRtable (1,1),iRtable (1,2),rRtable(1,1),rRtable(1,2),iRtable (1,3),iRtable (1,4),rRtable(1,3),iRtable (1,5),rSRTM(1)
          write (11,'(i10,i9,2f10.4,2i5,f19.12,i10,i15,e13.4)') iRtable (1,1),iRtable (1,2),rRtable(1,1),rRtable(1,2),iRtable (1,3),iRtable (1,4),rRtable(1,3),iRtable (1,5),rSRTM(1)          
       end do
       
    else
       
       icol = 7
       rcol = 5
       allocate (iRtable (1:NTILES, 1:icol))
       allocate (iWtable (1:ITILES, 1:icol))
       allocate (rRtable (1:NTILES, 1:rcol))
       allocate (rWtable (1:ITILES, 1:rcol))
       
       do n = 1,  ntiles  
          read(10,'(I10,3E20.12,9(2I10,E20.12,I10))') iRtable (n,1),rRtable(n,1),rRtable(n,2),rRtable(n,3),iRtable (n,2),iRtable (n,3),rRtable(n,4),iRtable (n,4),&
               iRtable (n,5),iRtable (n,6),rRtable(n,5),iRtable (n,7)
       end do

       do i = 1, icol
          iWtable (:,i) = iRtable (tile_id,i) 
       end do
       do i = 1, rcol
          rWtable (:,i) = rRtable (tile_id,i) 
       end do
       
       ! rWtable(n,4), rWtable(n,5)  split cell fractions
       do n = 1, ITILES
          if( irr_crp_tile(n) > 0) then
             ! either irrig crop or paddy thus non-irrig tile fraction and pfaf fraction
             if(irr_pad(n) <= 2) then
                rWtable(n-1,1) = (1. - IRRIGFRAC(tile_id(n))- PADDYFRAC(tile_id(n))) *  rWtable(n-1,1) ! area
                rWtable(n-1,4) = (1. - IRRIGFRAC(tile_id(n))- PADDYFRAC(tile_id(n))) *  rWtable(n-1,4) ! cell frac
                rWtable(n-1,5) = (1. - IRRIGFRAC(tile_id(n))- PADDYFRAC(tile_id(n))) *  rWtable(n-1,5) ! pfaf_frac
             endif
             ! irrigated crop tile
             if(irr_pad(n) == 1) then
                rWtable(n,1)   = IRRIGFRAC(tile_id(n)) * rWtable(n,1)
                rWtable(n,4)   = IRRIGFRAC(tile_id(n)) * rWtable(n,4)
                rWtable(n,5)   = IRRIGFRAC(tile_id(n)) * rWtable(n,5)               
             endif

             ! irrigated paddy tile
             if(irr_pad(n) >= 2) then
                rWtable(n,1)   = PADDYFRAC(tile_id(n)) * rWtable(n,1)
                rWtable(n,4)   = PADDYFRAC(tile_id(n)) * rWtable(n,4)
                rWtable(n,5)   = PADDYFRAC(tile_id(n)) * rWtable(n,5)                  
             endif
          endif
       end do

       do n = 1, ITILES
          write (11,'(I10,3E20.12,9(2I10,E20.12,I10))') iWtable (n,1),rWtable(n,1),rWtable(n,2),rWtable(n,3),iWtable (n,2),iWtable (n,3),rWtable(n,4),iWtable (n,4),&
               iWtable (n,5),iWtable (n,6),rWtable(n,5),iWtable (n,7)
       end do
       
       ! non-land

       do n = 1, non_land
          read (10,'(I10,3E20.12,9(2I10,E20.12,I10))') iRtable (1,1),rRtable(1,1),rRtable(1,2),rRtable(1,3),iRtable (1,2),iRtable (1,3),rRtable(1,4),iRtable (1,4),&
               iRtable (1,5),iRtable (1,6),rRtable(1,5),iRtable (1,7)
          write(11,'(I10,3E20.12,9(2I10,E20.12,I10))') iRtable (1,1),rRtable(1,1),rRtable(1,2),rRtable(1,3),iRtable (1,2),iRtable (1,3),rRtable(1,4),iRtable (1,4),&
               iRtable (1,5),iRtable (1,6),rRtable(1,5),iRtable (1,7)          
       end do
       
    endif

    close (10, status = 'keep')
    close (11, status = 'keep')
    
  END SUBROUTINE write_tilfile

  ! ----------------------------------------------------------------------------------

  SUBROUTINE write_raster_file (infile)

    implicit none
    character(*), intent (in) :: infile
    REAL,          PARAMETER     :: UNDEF = -9999.
    
    ! Global Irrigated Area data (GIA)
    ! --------------------------------
    
    integer,       parameter :: NX_GIA = 43200
    integer,       parameter :: NY_GIA = 21600,  NY_GIAData = 18000  
    character*300, parameter :: GIA_file = 'data/CATCH/IRRIGATION/global_irrigated_areas.nc4'

    ! GRIPC data
    ! ----------
    
    integer,       parameter :: NX_gripc = 86400
    integer,       parameter :: NY_gripc = 43200,  NY_GripcData = 36000
    character*300, parameter :: GRIPC_file = 'data/CATCH/IRRIGATION/irrigtype_salmon2013.flt'

    integer, target, allocatable, dimension (:,:) :: gripc, gia, tid_old, tid_new
    integer, pointer,dimension (:,:)              :: subset_tidold, subset_gripc, subset_tidnew
    integer, allocatable, dimension (:,:)         :: data_in
    real,    allocatable, dimension (:,:)         :: rdata_in
    integer :: i,j,n,k, r, status, NCID
    integer :: NT, NPF, NC, NR, NG, IM, JM, i_ease, j_ease
    character*20 :: GRIDNAME
    real*8        :: dxy, d2r, r2d, lats, lons
    real*8,   allocatable :: xs(:,:), ys(:,:)
    real          :: x,y, xout, yout
    character*5   :: MGRID
    integer       :: i1s, i2s, i1n, i2n, j1w, j2w, j1e, j2e, i1, i2, j1, j2, icol, rcol
    real, allocatable, dimension(:)    :: irrtile_lat, irrtile_lon, tile_area
    real, allocatable, dimension(:)    :: pix_cnt
    integer, dimension (:,:), allocatable :: iRtable, iWtable
    real,    dimension (:,:), allocatable :: rRtable, rWtable
    integer*8, allocatable, dimension (:) :: rSRTM, wSRTM
    real, allocatable, dimension (:) :: IRRIGFRAC, PADDYFRAC
    
    call InFmt%open(trim(BCSNAME)//'/'//trim(infile), pFIO_READ,rc=rc) ; VERIFY_(RC)
    allocate (IRRIGFRAC (1:NTILES))
    allocate (PADDYFRAC (1:NTILES))
    call MAPL_VarRead (InFmt,'PADDYFRAC', PADDYFRAC, rc=rc) ; VERIFY_(RC)
    call MAPL_VarRead (InFmt,'IRRIGFRAC', IRRIGFRAC, rc=rc) ; VERIFY_(RC)
    call inFmt%close (rc=rc)
    
    allocate (gia     (1:nc_raster, nr_raster))
    allocate (gripc   (1:nc_raster, nr_raster))

    ! GIA irrigated pixel mask at 30 arcsec
    ! -------------------------------------
    
    allocate(data_in(NX_GIA,NY_GIA)) 
    data_in = UNDEF   
    status = NF_OPEN (trim(GIA_file),NF_NOWRITE, ncid) ; VERIFY_(STATUS)    
    
    do j = NY_GIAData, 1, -1
       status = NF_GET_VARA_INT(NCID,NC_VarID(NCID,'IrrigClass') ,(/1,j/),(/NX_GIA, 1/), data_in (:,j + 3600 )) ; VERIFY_(STATUS)           
    end do
      
    status = NF_CLOSE(NCID) ; VERIFY_(STATUS)
    ! we select only 1,2, & 3
    where (data_in < 1)
       data_in = undef
    end where
    where (data_in > 3)
       data_in = undef
    end where

    gia = undef
    
    if (NC_RASTER /= NX_GIA) then
       call RegridRaster (data_in, gia)
    else
       gia = data_in
    endif
    DEALLOCATE (data_in)
    
    ! GRIPC paddy and crop mask at 30 arcsec
    ! --------------------------------------
    
    allocate( rdata_in(NX_gripc,NY_gripc)) 
    rdata_in = UNDEF
      
    open ( 10, file = trim(GRIPC_file), form = 'unformatted', access='direct', recl=(NX_gripc))
    
    !- Read input file::
    
    do j = 1, NY_gripcdata
       r = NY_gripc -j + 1
       read(10,rec=j) rdata_in(:, r)
       do i = 1, NX_gripc
          ! 2 crop, 3 paddy
          if( rdata_in(i, r) <= 1. ) rdata_in(i, r) = undef
          if( rdata_in(i, r) == 4. ) rdata_in(i, r) = undef
       end do
    end do
    close(10)

    if (NC_RASTER /= NX_gripc) then
       call RegridRaster (nint(rdata_in), gripc)
    else
       gripc = nint(rdata_in)
    endif
    DEALLOCATE (rdata_in)

    where (gia == undef)
       gripc = undef
    endwhere

    where ((gia > 0) .and. (gripc == undef))
       gripc = 2
    end where

    ! tile_id rasters
    ! ---------------

    allocate (tid_old (1:nc_raster, nr_raster))
    allocate (tid_new (1:nc_raster, nr_raster))
    tid_new = 0
        
    open (10,file=trim(BCSNAME)//'/rst/'//tilfile(1:index(tilfile,'.')-1)//'.rst',status='old',action='read',  &
         form='unformatted',convert='little_endian')
    open (11,file='IRRIG/'//trim(BCSNAME)//'/rst/'//tilfile(1:index(tilfile,'.')-1)//'.rst',status='unknown',action='write',  &
         form='unformatted',convert='little_endian')
    
    do j=1,nr_raster       
       ! read a row       
       read(10) tid_old (:,j)
    end do

    close (10, status = 'keep')    

    ! create raster
    ! -------------

      open (10,file =  trim(BCSNAME)//'/'//trim(tilfile), form = 'formatted', action = 'read', status = 'old')
    read (10, *) NT, NPF, NC, NR
    read (10, *) NG
    read (10, '(a)') GRIDNAME
    read (10, *) IM
    read (10, *) JM
    
    if(index(tilfile,"EASE")/=0) then
       
       open (12,file = 'IRRIG/'//trim(BCSNAME)//'/'//trim(adjustl(GRIDNAME))//'_'//trim(adjustl(GRIDNAME))//'-Irrigation.TIL', form = 'formatted', action = 'write', status = 'unknown')
       write (12, *) ITILES, NPF, NC, NR
       write (12, *) 2
       write (12, '(a)') trim(GRIDNAME)
       write (12, *) IM
       write (12, *) JM
       write (12, '(a)') trim(GRIDNAME)//'-Irr'
       write (12, *) IM
       write (12, *) JM

       if (NG == 2) then
          read (10, '(a)') GRIDNAME
          read (10, *) K
          read (10, *) K          
       end if

       icol = 5
       rcol = 3
       allocate (iRtable (1:NTILES, 1:icol))
       allocate (iWtable (1:ITILES, 1:icol))
       allocate (rRtable (1:NTILES, 1:rcol))
       allocate (rWtable (1:ITILES, 1:rcol))
       allocate (rSRTM   (1:NTILES))
       allocate (wSRTM   (1:ITILES))
       
       do n = 1,  ntiles
          read (10,'(i10,i9,2f10.4,2i5,f19.12,i10,i15,e13.4)') iRtable (n,1),iRtable (n,2),rRtable(n,1),rRtable(n,2),iRtable (n,3),iRtable (n,4),rRtable(n,3),iRtable (n,5),rSRTM(n)
       end do
       
       do i = 1, icol
          iWtable (:,i) = iRtable (tile_id,i) 
       end do
       do i = 1, rcol
          rWtable (:,i) = rRtable (tile_id,i) 
       end do
       wSRTM = rSRTM(tile_id)
       ! rWtable(n,3) - split cell fractions
       do n = 1, ITILES
          if( irr_crp_tile(n) > 0) then
             ! either irrig crop or paddy thus non-irrig tile fraction 
             if(irr_pad(n) <= 2) rWtable(n-1,3) = (1. - IRRIGFRAC(tile_id(n)) - PADDYFRAC(tile_id(n))) *  rWtable(n-1,3)
            
             ! irrigated crop tile
             if(irr_pad(n) == 1) rWtable(n,3)   = IRRIGFRAC(tile_id(n)) * rWtable(n,3)

             ! irrigated paddy tile
             if(irr_pad(n) >= 2) rWtable(n,3)   = PADDYFRAC(tile_id(n)) * rWtable(n,3)
          endif
       end do
       
    endif
    close (10, status = 'keep')

    ! Allocate and define the Cell vertices
    !--------------------------------------

    if(index(tilfile,"EASE")/=0) then
       
       allocate(xs(IM+1,JM+1), ys(IM+1,JM+1),stat=STATUS)  
       if(index(tilfile,"M09")/=0) MGRID = 'M09'
       if(index(tilfile,"M36")/=0) MGRID = 'M36'
       
       do  j = 1, JM
          do i = 1, IM
             x = real(i-1)        -0.5
              y = real(JM - j)+0.5
              call easeV2_inverse(MGRID, x, y, yout, xout)
              ys (i,j) = dble(yout)
              xs (i,j) = dble(xout)
           end do
        end do
        
        do  j = JM + 1, JM + 1
           do i = IM + 1, IM + 1
              x = real(i-1)         -0.5
              y =  -0.5
              call easeV2_inverse(MGRID, x, y, yout, xout)
              ys (i,j) = dble(yout)
              xs (i,j) = dble(xout)        
           end do
        end do

        where (ys > 90.)
           ys = 90.D0
        endwhere
        where (ys < -90.)
           ys = -90.D0
        endwhere
        where (xs > 180.)
           xs = 180.D0
        endwhere
        where (xs < -180.)
           xs = -180.D0
        endwhere
     else
        allocate(xs(IM+1,(IM+1)*6), ys(IM+1,(IM+1)*6),stat=STATUS)  
        call Get_CubedSphere_Grid(IM+1, (IM+1)*6, xs, ys, 0, .true.)
     endif

     dxy = 360.D0/dble(nc_raster)

     ! populate new tid raster
 
     if(index(tilfile,"EASE")/=0) then     
        do j=2,JM-1
           do i = 1,IM

              i1s = floor   ((xs (i,j) + 180.)/dxy) + 1
              i2s = ceiling ((xs (i+1,j) + 180.)/dxy) + 1
              i1n = floor   ((xs (i,j+1) + 180.)/dxy) + 1
              i2n = ceiling ((xs (i+1,j+1) + 180.)/dxy) + 1
              j1w = floor   ((ys (i,j) + 90.)/dxy) + 1  
              j2w = ceiling ((ys (i,j+1) + 90.)/dxy) + 1 
              j1e = floor   ((ys (i+1,j) + 90.)/dxy) + 1 
              j2e = ceiling ((ys (i+1,j+1) + 90.)/dxy) + 1
              
              i1 = max(1, minval ((/i1s,i1n/)))
              i2 = min(nc_raster, maxval ((/i2s,i2n/)))
              j1 = max(1, minval ((/j1w,j1e/)))
              j2 = min(nr_raster,maxval ((/j2w,j2e/)))
              if (i2 > i1) then
                 if (associated (subset_tidold)) NULLIFY (subset_tidold)
                 if (associated (subset_tidnew)) NULLIFY (subset_tidnew)
                 if (associated (subset_gripc))  NULLIFY (subset_gripc )             
                 subset_tidold => tid_old (i1:i2,j1:j2)
                 subset_tidnew => tid_new (i1:i2,j1:j2)
                 subset_gripc  => gripc   (i1:i2,j1:j2)
                 
                 ! update tile id raster
                 call update_raster (subset_tidold, subset_gripc, subset_tidnew)
              endif
           end do
        end do

     else
        ! cubed-sphere
        do k = 1,6
           do j=1,IM
              do i = 1,IM
                 
                 i1s = floor   ((xs (i,j + (k-1)*IM) + 180.)/dxy) + 1
                 i2s = ceiling ((xs (i+1,j + (k-1)*IM) + 180.)/dxy) + 1
                 i1n = floor   ((xs (i,j+1 + (k-1)*IM) + 180.)/dxy) + 1
                 i2n = ceiling ((xs (i+1,j + 1 + (k-1)*IM) + 180.)/dxy) + 1
                 j1w = floor   ((ys (i,j+ (k-1)*IM) + 90.)/dxy) + 1  
                 j2w = ceiling ((ys (i,j+1 + (k-1)*IM) + 90.)/dxy) + 1 
                 j1e = floor   ((ys (i+1,j + (k-1)*IM) + 90.)/dxy) + 1 
                 j2e = ceiling ((ys (i+1,j + 1 + (k-1)*IM) + 90.)/dxy) + 1
                 
                 i1 = minval ((/i1s,i1n/))
                 i2 = maxval ((/i2s,i2n/))
                 j1 = minval ((/j1w,j1e/))
                 j2 = maxval ((/j2w,j2e/))
                 if (associated (subset_tidold)) NULLIFY (subset_tidold)
                 if (associated (subset_tidnew)) NULLIFY (subset_tidnew)
                 if (associated (subset_gripc))  NULLIFY (subset_gripc )
                 subset_tidold => tid_old (i1:i2,j1:j2)
                 subset_tidnew => tid_new (i1:i2,j1:j2)
                 subset_gripc  => gripc   (i1:i2,j1:j2)

                 ! update tile id raster
                 call update_raster (subset_tidold, subset_gripc, subset_tidnew)

              end do
           end do
        end do
     endif

     ! write raster file
     print *, 'Writing irrigation raster ..'
     do j=1,nr_raster       
        ! read a row       
        write (11) tid_new (:,j)
     end do
    
     close (11, status = 'keep')    
     
     if(index(tilfile,"EASE")/=0) then
        
        d2r    = (dble(MAPL_PI)/180.D0)
        ! write new til file only for EASE grids
        allocate (tile_area   (1: ITILES))
        allocate (irrtile_lat (1: ITILES))
        allocate (irrtile_lon (1: ITILES))
        allocate (pix_cnt     (1: ITILES))

        pix_cnt     = 0.
        tile_area   = 0.
        irrtile_lat = 0.
        irrtile_lon = 0.
        
        do j=1,nr_raster
           lats = -90.D0 + dxy/2.D0 + (j-1)*dxy
           do i = 1, nc_raster
              if((tid_new(i,j) >= 1).and.(tid_new(i,j) <= ITILES)) then
                 lons = -180.D0 + dxy/2.D0 + (i-1)*dxy
                 irrtile_lat(tid_new(i,j)) = irrtile_lat(tid_new(i,j)) + real(lats)
                 irrtile_lon(tid_new(i,j)) = irrtile_lon(tid_new(i,j)) + real(lons)
                 tile_area  (tid_new(i,j)) = tile_area  (tid_new(i,j)) + &
                      real((sin(d2r*(lats+0.5*dxy)) - sin(d2r*(lats-0.5*dxy)))*(dxy*d2r))
                 pix_cnt (tid_new(i,j)) = pix_cnt (tid_new(i,j)) + 1.
              endif
           end do
        end do
        where (pix_cnt > 0.)
           irrtile_lat = irrtile_lat / pix_cnt
           irrtile_lon = irrtile_lon / pix_cnt
        end where

        do n = 1, ITILES
           if(tile_area(n) > 0) then
              write (12,'(I10,3E20.12,9(2I10,E20.12,I10))') 100,tile_area(n),irrtile_lon(n),irrtile_lat(n),iWtable (n,3),iWtable (n,4),rWtable(n,3),1,&
                   iWtable (n,2),1,1.,iWtable (n,2)
           else
              write (12,'(I10,3E20.12,9(2I10,E20.12,I10))') 100,tile_area(n),irrtile_lon(n-1),irrtile_lat(n-1),iWtable (n,3),iWtable (n,4),rWtable(n,3),1,&
                   iWtable (n,2),1,1.,iWtable (n,2)
           endif
        end do
        close (12, status = 'keep')
     endif
         
  END SUBROUTINE write_raster_file

  ! --------------------------------------------------------------------

  SUBROUTINE update_raster (subset_tidold, subset_gripc, subset_tidnew)

    implicit none
    integer, dimension (:,:), intent (in)   :: subset_tidold, subset_gripc
    integer, dimension (:,:), intent (inout):: subset_tidnew
    integer  :: tid_land(1), tid_crop(1), tid_pad(1),n, NPLUS, NPADDY, NCROP
    integer, allocatable, dimension (:) :: loc_int
    logical, allocatable, dimension (:) :: unq_mask

    NPLUS = count(subset_tidold >= 1 .and. subset_tidold <= NTILES)
    
    if (NPLUS > 0) then
       NPADDY = count(subset_gripc == 3)
       NCROP  = count(subset_gripc == 2)

       allocate (loc_int (1:NPLUS))
       allocate (unq_mask(1:NPLUS))
       loc_int = pack(subset_tidold,mask = (subset_tidold >= 1 .and. subset_tidold <= NTILES))
       call MAPL_Sort (loc_int)
       unq_mask = .true.
       do n = 2,NPLUS 
          unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))
       end do
       
       do n = 1,NPLUS
          if( unq_mask(n)) then
             tid_land = findloc (tile_id,loc_int(n),mask = (irr_pad == 0))
             if(tid_land(1) > 0) then
                where (subset_tidold == loc_int(n) .and. subset_gripc < 2)
                   subset_tidnew = tid_land(1)
                endwhere
             endif
             if(NPADDY > 0) then
                tid_pad = findloc (tile_id,loc_int(n),mask = (irr_pad >= 2))
                if(tid_pad(1)  > 0) then
                    where (subset_tidold == loc_int(n) .and. subset_gripc == 3)
                      subset_tidnew = tid_pad(1)                      
                   endwhere
                endif
             endif

             if(NCROP > 0) then
                tid_crop= findloc (tile_id,loc_int(n),mask = (irr_pad == 1))
                if(tid_crop(1)  > 0) then
                    where (subset_tidold == loc_int(n) .and. subset_gripc == 2)
                      subset_tidnew = tid_crop(1)
                   endwhere
                endif
             endif            
          endif
       end do
       DEALLOCATE (loc_int, unq_mask)
    endif
  END SUBROUTINE update_raster
 
END PROGRAM mkIrrigTiles
