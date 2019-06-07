PROGRAM mk_RouteRestarts

  use routing_model,                    ONLY:     &
       SEARCH_DNST

  USE catch_constants, ONLY:          &
       N_All => N_Pfaf_Catchs !, N_LND => N_Pfaf_LandCatchs

  use MAPL_IOMod

  implicit none

  character*256 :: Usage="mk_RouteRestarts OutTileFile RestartTime"
  character*256 :: arg(2), inFile, OutTileFile, OutRstFile, RouteDataFile
  character*10  :: RestartTime, CDUM
  CHARACTER*2   :: MM
  integer       :: I,N, STATUS, NCFID,NCOUT, Pfaf_land,K, ThisTile,N_LND
  integer       :: IM, JM, NC, NR, NT, i1, i2, pfaf_index, dnst_pfaf_index, N_Land, N_Ocean, N_Lakes, N_Ice
  real          :: min_lon, max_lon, min_lat, max_lat, mean_elevation, cat_area, length, &
       ElevDiff, UP_lon, UP_lat, mouth_lon, mouth_lat, VDUM, dx, dy
  real, dimension (:), allocatable     :: AREACAT,LENGSC, WSTREAM,WRIVER, tmp_var, DN_lon, DN_lat,cum_area
  integer, allocatable, dimension (:)  :: DNST,  DNST2, Pfaf_all,type
  integer, target, allocatable, dimension (:) :: i_glb,j_glb, ORiverMouth,LRiverMouth
  integer, pointer, dimension (:)             :: i_ocean,j_ocean, i_lakes, j_lakes, i_land, j_land
  logical, allocatable, dimension (:)  :: mask
  integer, allocatable, dimension (:,:):: tile_id, ocean_id, lake_id
  integer(kind =8)                     :: pfaf_code

  INCLUDE 'netcdf.inc'

  RouteDataFile = '/gpfsm/dnb42/projects/p16/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/SRTM-TopoData/Pfafcatch-routing.dat'

  I = command_argument_count()
  
  if( I /=2 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     stop
  end if

  do n = 1,I
     call getarg(n,arg(n))
  end do
  read(arg(1),'(a)') OutTileFile
  read(arg(2),'(a)') RestartTime

  ! STEP 1: Copy restarts from offline simulation

  MM = RestartTime (5: 6)

  inFile = '/discover/nobackup/rreichle/l_data/LandRestarts_for_Regridding/route/' &
               //'route_internal_rst.YYYY'//MM//'01'

  print '(a150)','ROUTE_INTERNAL_RST(restarts) : '//trim(inFile)

  allocate (AREACAT (1:N_All))
  allocate (LENGSC  (1:N_All))
  allocate (DNST    (1:N_All))
  allocate (WSTREAM (1:N_All))
  allocate (WRIVER  (1:N_All))
  allocate (tmp_var (1:N_All))

  STATUS = NF_OPEN   (trim(inFile),NF_NOWRITE,NCFID)
  STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'AREACAT'   ), (/1/), (/N_All/),AREACAT )
  STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'LENGSC'    ), (/1/), (/N_All/),LENGSC  )
  STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'DNSTR'     ), (/1/), (/N_All/),tmp_var )
  STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WSTREAM'   ), (/1/), (/N_All/),WSTREAM )
  STATUS = NF_GET_VARA_REAL(NCFID,VarID(NCFID,'WRIVER'    ), (/1/), (/N_All/),WRIVER  )
  
  STATUS   = NF_CLOSE (NCFID) 

  DNST    = NINT (tmp_var(:))
  
  ! Reset DNST to by pass submerged catchments: 
  ! Out of 291284 (N_Pfaf_Catchs in catch_constants.f90) only 290188 (N_Pfaf_LandCatchs in catch_constants.f90)
  ! are overland Pfafstetetter catchments. The rest is submerged under lakes or large river segments.
  ! Here, we check if the assigned downstream catchment is submerged and if so we bypass the submerged catchments(s)
  ! and link to the first available over land catchment downstream of the catchment in question.

  allocate   (DNST2 (1:N_All), Pfaf_all (1:N_All))

  Pfaf_all(:) = -9999
  DNST2       = -1 

  ! Read land only catchment information
  i1 = index(OutTileFile,'/')
  i2 = index(OutTileFile,'.',back=.true.)

  OPEN (10, FILE = TRIM(OutTileFile(1:i1))//'til'//'/Pfafstetter.til', &
       FORM='FORMATTED', STATUS = 'OLD', ACTION = 'READ')
  READ (10,*)  N_LND
  DO N = 1,4            
     READ (10,'(A)')            
  END DO
  
  DO N = 1, N_LND
     READ (10,*)K,VDUM,VDUM, VDUM, Pfaf_land, K
     if (Pfaf_land <= N_ALL) Pfaf_all (Pfaf_land) = Pfaf_land
  END DO
  
  close (10, status = 'keep')

  DO N = 1, N_All

     ! STEP 2: Bypass submerged catchments and relink to the next downstream land catchment 
    
     if(Pfaf_all(N) >= 1) call SEARCH_DNST (N, N_All, DNST, Pfaf_all, DNST2(N))

  END DO

  ! STEP 3: Find ocean/lake tile ids at river mouths

  ! (3a) Reading til file
  print '(a100)','ROUTE_INTERNAL_RST (OutTileFile) : '// trim(OutTileFile)

  open (10, file=trim(OutTileFile), form='formatted', status ='old',action = 'read') 
  read (10,*) NT,NC,NR
  read (10,*) I
  read (10,'(a)') CDUM
  read (10,*) IM
  read (10,*) JM
  read (10,'(a)') CDUM
  read (10,*) I
  read (10,*) I

  allocate (tile_id (1:nc, 1:nr))
  allocate (ocean_id(1:im, 1:jm))
  allocate (lake_id (1:im, 1:jm))
  ocean_id = -9999
  lake_id  = -9999
  
  allocate (mask (1:NT))
  allocate (type (1:NT))
  allocate (i_glb(1:NT))
  allocate (j_glb(1:NT))  

  DO I = 1, NT

     read (10, *) type (i), cat_area,  min_lon, max_lon, i_glb(i),j_glb(i)

  END DO

  close (10, status = 'keep')

  mask = .false.
  mask = (type == 100)
  N_Land = count (mask = mask)
  allocate (i_land (1:N_Land))
  allocate (j_land (1:N_Land)) 
  i_land = pack(i_glb, mask = mask)
  j_land = pack(j_glb, mask = mask) 

  mask = .false.
  mask = (type == 0)
  N_Ocean = count (mask = mask)
  allocate (i_ocean (1:N_Ocean))
  allocate (j_ocean (1:N_Ocean)) 
  i_ocean = pack(i_glb, mask = mask)
  j_ocean = pack(j_glb, mask = mask) 

  do i = 1,N_ocean
     ocean_id (i_ocean(i), j_ocean(i)) = i  ! + N_Land
  end do

  mask = .false.
  mask = (type == 19)
  N_Lakes = count (mask = mask)  
  allocate (i_lakes (1:N_Lakes))
  allocate (j_lakes (1:N_Lakes)) 
  i_lakes = pack (i_glb, mask = mask)
  j_lakes = pack (j_glb, mask = mask) 

  do i = 1,N_lakes
     lake_id (i_lakes(i), j_lakes(i)) = i  ! + N_ocean + N_Land
  end do

  mask = .false.
  mask = (type == 20)
  N_ice = count (mask = mask)    

  deallocate (mask, i_lakes, j_lakes, i_ocean, j_ocean)

 ! print *,N_land, N_ocean, N_lakes, N_ice
 ! print *, NT, N_land + N_ocean + N_lakes + N_ice

  ! (3b) Reading rst file
  i1 = index(OutTileFile,'/')
  i2 = index(OutTileFile,'.',back=.true.)
 
  OutRstFile = OutTileFile(1:i1)//'rst'//OutTileFile(i1:i2)//'rst'
  print '(a100)','ROUTE_INTERNAL_RST (OutRstFile) : '// trim(OutRstFile)
  open (10, file=trim(OutRstFile), form='unformatted', status ='old',action = 'read', convert='little_endian') 
  do i =1,nr
     read (10) tile_id(:,i)
  end do
  close (10, status = 'keep') 

  ! (3c) Reading Pfafcatch-routing.dat for lat/lon of downstream confluence
  print '(a150)','ROUTE_INTERNAL_RST (RouteDataFile) : '// trim(RouteDataFile)
  open (10, file=trim(RouteDataFile), form='formatted', status ='old',action = 'read') 
  read (10,*) I1
  allocate (dn_lat(1:i1))
  allocate (dn_lon(1:i1))
  allocate (cum_area(1:i1))
  allocate (ORiverMouth (1:i1))
  allocate (LRiverMouth (1:i1))
  ORiverMouth = -9999
  LRiverMouth = -9999

  do i = 1,i1
  read(10,'(i8,i15,4(1x,f9.4),1x,e10.3,4(1x,e9.3),I8,6(1x,f9.4))') & 
     pfaf_index, pfaf_code, min_lon, max_lon, min_lat, max_lat, mean_elevation,           &
     cat_area,cum_area(i), length, ElevDiff, dnst_pfaf_index, DN_lon(i), DN_lat(i), UP_lon,   &
     UP_lat, mouth_lon, mouth_lat
  end do
  close (10, status = 'keep')

  ! (3d) mapping river mouths
  dx = 360./real(NC)
  dy = 180./real(NR)

  do i = 1, size (dn_lat)
     if(DNST2(i) == -1) then

        IM =  floor((DN_lon(i) + 180.)/dx) + 1
        JM =  floor((DN_lat(i) +  90.)/dy) + 1 
        ThisTile = tile_id(IM, JM)

        im = i_glb(ThisTile)
        jm = j_glb(ThisTile)
        
        if(ocean_id(im,jm) > 0) then
           ! Ocean
           ORiverMouth (i) =  ocean_id(im,jm)
           ! print *, 'Ocean tile at  :',i,ORiverMouth (i),DN_lat(i),DN_lon(i)
        elseif((ocean_id(im,jm) <= 0).and.(lake_id(im,jm) > 0)) then
           ! Lake
           LRiverMouth (i) =  lake_id(im,jm)
           ! print *, 'Lake tile at  :',i,  LRiverMouth (i) ,DN_lat(i),DN_lon(i)       
        else
           !  print *, 'No ocean or lake for :',i, cum_area(i),DN_lat(i),DN_lon(i)
           ! stop
        endif
        
     endif
  end do

  ! STEP 4: Now write route_internal_rst

  call create_route_internal_nc4 (N_All,NCOUT)   
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'AREACAT'   ), (/1/), (/N_All/),AREACAT )
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'LENGSC'    ), (/1/), (/N_All/),LENGSC  )
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'DNSTR'     ), (/1/), (/N_All/),REAL(DNST2))
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'WSTREAM'   ), (/1/), (/N_All/),WSTREAM )
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'WRIVER'    ), (/1/), (/N_All/),WRIVER  )
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'ORIVERMOUTH'), (/1/),(/N_All/),REAL (ORIVERMOUTH))
  STATUS = NF_PUT_VARA_REAL(NCOUT,VarID(NCOUT,'LRIVERMOUTH'), (/1/),(/N_All/),REAL (LRIVERMOUTH))
  status = NF_CLOSE (NCOUT)  

  deallocate (AREACAT,LENGSC, WSTREAM,WRIVER, tmp_var)
  deallocate (DNST, DNST2, Pfaf_all)

contains
  
  ! ----------------------------------------------------------------------
  
  integer function VarID (NCFID, VNAME) 
    
    integer, intent (in)      :: NCFID
    character(*), intent (in) :: VNAME
    integer                   :: status
    
    STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
    IF (STATUS .NE. NF_NOERR) &
         CALL HANDLE_ERR(STATUS, trim(VNAME))  
    
  end function VarID
  
  ! -----------------------------------------------------------------------
  
  SUBROUTINE HANDLE_ERR(STATUS, Line)
    
    INTEGER,      INTENT (IN) :: STATUS
    CHARACTER(*), INTENT (IN) :: Line
    
    IF (STATUS .NE. NF_NOERR) THEN
       PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
       STOP 'Stopped'
    ENDIF
    
  END SUBROUTINE HANDLE_ERR

  ! ---------------------------------------------------------------------------------------------

  SUBROUTINE create_route_internal_nc4 (catch, NCFOutID)

    implicit none

    integer, intent (in)       :: catch
    integer, intent (inout)    :: NCFOutID    
    integer :: CatchID, TimID, VID, status 

    status = NF_CREATE ('OutData/route_internal_rst', NF_NETCDF4, NCFOutID)
    status = NF_DEF_DIM(NCFOutID, 'tile', catch, CatchID)
    status = NF_DEF_DIM(NCFOutID, 'time'   , 1 , TimID)
    
    status = NF_DEF_VAR(NCFOutID, 'time'  , NF_DOUBLE, 1 , TimID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
         LEN_TRIM('minutes since  2014-01-01 00:00:00'), trim('minutes since 2014-01-01 00:00:00'))   
    status = NF_DEF_VAR(NCFOutID, 'AREACAT'  , NF_FLOAT, 1 , CatchID  , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',&
         LEN_TRIM('AREA OF CATCHMENT'),trim('AREA OF CATCHMENT')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units', &
         LEN_TRIM('km+2'), trim('km+2')) 
    
    status = NF_DEF_VAR(NCFOutID, 'LENGSC' , NF_FLOAT, 1 , CatchID , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name', &
         LEN_TRIM('LENGTH OF CHANNEL SEGMENT'),          &
         trim('LENGTH OF CHANNEL SEGMENT')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',     &
         LEN_TRIM('km'), trim('km')) 
    
    status = NF_DEF_VAR(NCFOutID, 'DNSTR' , NF_FLOAT, 1 , CatchID , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name', &
         LEN_TRIM('INDEX OF DOWNTREAM CATCHMENT'),       &
         trim('INDEX OF DOWNTREAM CATCHMENT')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',     &
         LEN_TRIM('-'), trim('-')) 
    
    status = NF_DEF_VAR(NCFOutID, 'WSTREAM' , NF_FLOAT, 1 , CatchID , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name', &
         LEN_TRIM('VOLUME OF WATER IN LOCAL STREAM'),    &
         trim('VOLUME OF WATER IN LOCAL STREAM')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',     &
         LEN_TRIM('m+3'), trim('m+3'))    
    
    status = NF_DEF_VAR(NCFOutID, 'WRIVER' , NF_FLOAT, 1 , CatchID , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name', &
         LEN_TRIM('VOLUME OF WATER IN RIVER'),           &
         trim('VOLUME OF WATER IN RIVER')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',     &
         LEN_TRIM('m+3'), trim('m+3'))        
    
    status = NF_DEF_VAR(NCFOutID, 'ORIVERMOUTH' , NF_FLOAT, 1 , CatchID , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name', &
         LEN_TRIM('TileID of the ocean tile at the river mouth'),           &
         trim('TileID of the ocean tile at the river mouth')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',     &
         LEN_TRIM('-'), trim('-'))        
    
    status = NF_DEF_VAR(NCFOutID, 'LRIVERMOUTH' , NF_FLOAT, 1 , CatchID , vid)
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name', &
         LEN_TRIM('TileID of the lake tile at the river mouth'),           &
         trim('TileID of the lake tile at the river mouth')) 
    status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',     &
         LEN_TRIM('-'), trim('-'))        
    
    status = NF_ENDDEF(NCFOutID)  

  END SUBROUTINE create_route_internal_nc4
  
END PROGRAM mk_RouteRestarts
