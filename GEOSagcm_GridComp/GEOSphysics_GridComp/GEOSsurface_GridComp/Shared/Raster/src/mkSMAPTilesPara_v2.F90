PROGRAM mkSMAPTilesPara_v2
!     This program constructs land and lake tiles for the SMAP-EASE-M09 and M36 grids (just set MGRID) 
!         for CLSM implementation.
!     f90 -c create_smap_tiles.f90
!     f90 -c smapconv.f
!     f90 -o create_smap_tiles create_smap_tiles.o smapconv.o
!
      use easeV2_conv
      use rmTinyCatchParaMod
      use process_hres_data
      use MAPL_SortMod
      use MAPL_ConstantsMod

      implicit none

      integer i,j,ig,jg,i0,iop,n,d1,d2,j1,j2,i1,i2,ix, jx,icount,pcount
      integer :: NC = i_raster, NR = j_raster, NT = 16330000, ND = 10000, ND_raster = 10000
      
      integer, parameter :: SRTM_maxcat = 291284, nc_esa = 129600, nr_esa = 64800

      ! For regridding

      integer, allocatable, target, dimension (:,:) &
                                          :: geos_msk 
      REAL,    allocatable, DIMENSION (:) :: loc_val
      INTEGER, ALLOCATABLE, DIMENSION (:) :: density, loc_int
      logical, dimension (:), allocatable :: unq_mask   
      integer, pointer  , dimension (:,:) :: subset
      integer, pointer    , dimension (:) :: subset1, subset_smap
      real,    pointer    , dimension (:) :: subset2
      integer :: dx_esa, dy_esa, NBINS, NPLUS

      integer*8, allocatable, dimension (:) ::  SRTM_catid

      integer,allocatable, dimension (:,:), target :: tileid_index,catid_index
      integer,allocatable, dimension (:,:)         :: catid, iaster
      integer,allocatable, dimension (:)           :: land_id,water_id,ice_id
      integer,allocatable, dimension (:)           :: my_land, all_id
      real, allocatable, dimension (:)             :: smap_grid_area,tile_area,SRTM_CatchArea 
      integer*1,allocatable, dimension (:,:)       :: veg, i2aster
      real*4, dimension (:,:), allocatable         :: q0,raster
      REAL, dimension (:), allocatable             :: tile_ele, tile_area_land  

      INTEGER*8 :: PFAF_CODE  
      integer l,imn,imx,jmn,jmx,mval,l_index,i_index,w_index,typ,pfaf,cindex
      integer :: LakeType, IceType, OceanType
      character(3) :: easegrid      
      real :: clat, clon, r_smap, s_smap, smap_convert, da
      real :: fr_gcm
      integer :: ind_col, ind_row, status, ncid, nciv,nland_cells, DOM_INDX
      REAL (kind=8), PARAMETER :: RADIUS=6378137.0,pi=3.14159265358979323846
      character*100 :: veg_class (12)
      character*5 :: MGRID
      character*100 :: gfile,gtopo30
      integer :: nc_smap,nr_smap, N_args, iargc 
      real :: EASE_grid_area, CELL_km
      REAL :: dx,dy,d2r,lats,mnx,mxx,mny,mxy,sum1,sum2,jgv, VDUM,pix_area
      character(40) :: arg
      character*200 :: tmpstring, tmpstring1, tmpstring2	      
      logical :: regrid = .false.
      character*128          :: MaskFile
      include 'netcdf.inc'

      N_args = iargc()

      if(N_args < 1) then
        print *,'USAGE : bin/mkSMAPTiles -smap_grid MXX'
	print *,'Allowed SMAP grids are: M01 M03 M09 M36'
        stop
      end if

      i=0      

      do while ( i < N_args )

         i = i+1
         
         call getarg(i,arg)
         
         if     ( trim(arg) == '-smap_grid' ) then
            i = i+1
            call getarg(i,MGRID)
            
         else ! stop for any other arguments
            
            print *,'USAGE : bin/mkSMAPTiles -smap_grid MXX'
            print *,'Allowed SMAP grids are: M09 M36 Ml'
            stop
            
         endif
         
      end do
      
      call system('cd data/ ; ln -s /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/ CATCH')  
      call system('cd ..')
      
      
      ! Setting SMAP Grid specifications
      ! --------------------------------
      
      if (trim(MGRID) == 'M09') then
         
         CELL_km = 9.008055210146     ! nominal cell size in kilometers
         nc_smap = 3856
         nr_smap = 1624
         gfile = 'SMAP_EASEv2_'//trim(MGRID)//'_3856x1624'
         EASE_grid_area = CELL_km*CELL_km
         
      elseif(trim(MGRID) == 'M36') then
         
         CELL_km = 36.032220840584    ! nominal cell size in kilometers
         nc_smap = 964
         nr_smap = 406
         gfile = 'SMAP_EASEv2_'//trim(MGRID)//'_964x406'
         EASE_grid_area = CELL_km*CELL_km
         
         !    elseif(trim(MGRID) == 'M25') then
         ! 
         !     	 CELL_km = 3.00003356589     ! nominal cell size in kilometers		 
         !         nc_smap = 1383
         !         nr_smap = 586
         !         gfile = 'SMAP_EASE_M25_1383x586'
         !         EASE_grid_area = CELL_km*CELL_km
         
      else if (trim(MGRID) .eq. 'M03') then ! SMAP  3 km grid
         CELL_km = 3.0026850700487     ! nominal cell size in kilometers
         nc_smap = 11568
         nr_smap = 4872
         gfile = 'SMAP_EASEv2_M03_11568x4872'
         EASE_grid_area = CELL_km*CELL_km
         regrid = .true.
         NC = 21600
         NR = 10800
         NT = 500000000
         
      else if (trim(MGRID) .eq. 'M01') then ! SMAP  1 km grid
         CELL_km = 1.00089502334956     ! nominal cell size in kilometers
         nc_smap = 34704
         nr_smap = 14616
         gfile = 'SMAP_EASEv2_M01_34704x14616'
         EASE_grid_area = CELL_km*CELL_km
         regrid = .true.
         NC = 43200
         NR = 21600   
         NT = 1500000000
         
      else  !
         
         print *,'Unknown SMAP Grid stopping..'
         stop
         
      endif

      allocate(land_id    (1:NT))
      allocate(water_id   (1:NT))
      allocate(ice_id     (1:NT))
      land_id     = 0
      water_id    = 0
      ice_id      = 0             
      OceanType = 0
      IceType   =11
      LakeType  =10

      ND        = 10*10**(nint(log10(1.*nr_smap)))
      
      !   Check for the 10 arc-sec MaskFile
      ! -----------------------------------
      
      call getenv ("MASKFILE"        ,MaskFile        )

      print *, 'Using MaskFile ', trim(MaskFile)

      if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then         
         ! New ESA (Veg) + SRTM (catchments) based mask file
         ! is overlaid on SMAP 
         ! -------------------------------------------------
         
         nc = 43200  ! Number of rows in raster file
         nr = 21600  ! Number of columns in raster file

         regrid = .true.
         dx_esa = nc_esa / nc ! x-dimension (or # of ESA columns within the raster grid cell)
         dy_esa = nr_esa / nr ! y-dimension (or # of ESA rows within the raster grid cell)

         allocate(tileid_index(1:nc,1:nr))
         allocate(SRTM_catid  (1:SRTM_maxcat+2))
         allocate(catid_index (1:nc,1:nr))          
         allocate(veg         (1:nc,1:nr))
         allocate(geos_msk    (1:nc_esa,1:dy_esa))
         allocate(SRTM_CatchArea (1:SRTM_maxcat))

         OPEN (10, FILE = 'data/CATCH/SRTM-TopoData/Pfafcatch-routing.dat', &
              FORM = 'FORMATTED',STATUS='OLD',ACTION='READ') 

         READ (10,*) I
         DO N = 1, I
            READ (10, '(i8,i15,4(1x,f9.4),1x,e10.3,4(1x,e9.3),I8,6(1x,f9.4))')    &
                 DOM_INDX,PFAF_CODE,VDUM,VDUM,VDUM,VDUM,VDUM,                     &
                 SRTM_CatchArea (N)
         END DO
         CLOSE (10, STATUS='KEEP')

         dx  = 360._8/nc
         dy  = 180._8/nr
         d2r = PI/180._8
         da  = MAPL_radius*MAPL_radius*pi*pi*dx*dy/180./180./1000000.    
         
         tileid_index = 0        
         catid_index  = 0
         veg          = 0
         
         status    = NF_OPEN ('data/CATCH/GEOS5_10arcsec_mask.nc', NF_NOWRITE, ncid)
         status    = NF_GET_VARA_INT64 (ncid,3,(/1/),(/SRTM_maxcat/),SRTM_catid(1:SRTM_maxcat))  ! Read pfafstetter IDs
         if(status /=0) then
            PRINT *, NF_STRERROR(STATUS)
            print *, 'Problem with NF_OPEN',trim(MaskFile)
         endif

         SRTM_catid (SRTM_maxcat + 1) = 190000000
         SRTM_catid (SRTM_maxcat + 2) = 200000000 
         i1 = 0  ! count # of 30-arcsec pixels

         do j=1,nr

            clat = -90. + float(j-1)*dy + dy/2.

            status  = NF_GET_VARA_INT (ncid,4,(/1,(j-1)*dy_esa +1/),(/nc_esa,dy_esa/),geos_msk) ! Read 10-arcsec rows that lie within the raster row 'j'  
            
            if(status /=0) then
               PRINT *, NF_STRERROR(STATUS)
               print *, 'Problem with NF_GET_VARA_INT',trim(MaskFile),status
            endif
            
            do i = 1,nc    

               clon = -180. + float(i-1)*dx + dx/2.

               if (associated (subset)) NULLIFY (subset)
               subset => geos_msk ((i-1)*dx_esa + 1 : i*dx_esa, 1:dy_esa) ! rectangular array contains ESA pixels that lie within the raster grid cell at i,j
               if(maxval (subset) > SRTM_maxcat) then
                  where (subset == 190000000) subset = SRTM_maxcat + 1
                  where (subset == 200000000) subset = SRTM_maxcat + 2
               endif
               
               if (maxval (subset) > 0) then ! check whether there are Non-ocean ESA pixels 
                  ! catID of the middle pixel

                  veg (i,j) = 1 ! veg is set to land

                  NPLUS = count(subset >= 1 .and. subset <= SRTM_maxcat + 2) ! Count non-ocean ESA pixels within  
                  allocate (loc_int (1:NPLUS))
                  allocate (unq_mask(1:NPLUS))
                  loc_int = pack(subset,mask = (subset >= 1 .and. subset <= SRTM_maxcat + 2)) ! loc_int contains catch_indices of non-ocean ESA pixels 
                  call MAPL_Sort (loc_int)
                  unq_mask = .true.
                  do n = 2,NPLUS 
                     unq_mask(n) = .not.(loc_int(n) == loc_int(n-1)) ! count number of unique numbers in loc_int for binning
                  end do
                  NBINS = count(unq_mask)

                  if (NBINS > 1) then
                     allocate(loc_val (1:NBINS))
                     allocate(density (1:NBINS))
                     loc_val = 1.*pack(loc_int,mask =unq_mask)                               ! loc_val contains available non-ocean catch_indices within the i,j grid cell,
                                                                                             ! Those numbers will be used as bin values
                     call histogram (dx_esa*dy_esa, NBINS, density, loc_val, real(subset))   ! density is the pixel count for each bin value
                     catid_index (i,j) = loc_val (maxloc(density,1))                         ! picks maximum density as the dominant catchment_index at i,j
                     deallocate (loc_val, density)
                  else
                     catid_index (i,j) = loc_int (1)
                  endif
                  deallocate (loc_int, unq_mask)

                  if(catid_index (i,j) == SRTM_maxcat + 1) veg (i,j) = LakeType
                  if(catid_index (i,j) == SRTM_maxcat + 2) veg (i,j) = IceType
                  if((catid_index(i,j) >= 1).and.(catid_index  (i,j) <=  SRTM_maxcat)) i1 = i1 + 1

                  ! count in if this is i,j pixel is a land, lake or ice within ind_col,ind_row SMAP grid cell
                  
                  call easeV2_convert(trim(MGRID), clat, clon, r_smap, s_smap)
                  
                  ind_col = nint(r_smap) + 1 
                  ind_row = nint(s_smap) + 1
                  
                  if((ind_row.ge.1).and.(veg(i,j).ne.OceanType).and.(ind_row.le.nr_smap)) then
                     l=  ind_row*ND +  ind_col
                     
                     if(veg(i,j)==LakeType) then 
                        water_id(l) = 1
                     else if(veg(i,j)==IceType) then
                        ice_id  (l) = 1
                     else
                        land_id (l) = 1
                     endif
                  endif
               endif
            end do
         enddo

         status    = NF_CLOSE (ncid)  
         deallocate (geos_msk)

         print *,'Read ', trim (MaskFile) 
         print *,'Min and Max of tile indices:',minval(catid_index),maxval(catid_index)

      else
         
         ! Old IGBP (Veg) + HYDRO1k (catchments) based mask will
         ! Overlaid on SMAP mask
         ! -----------------------------------------------------
         
         allocate(iaster      (i_raster,j_raster)) 
         allocate(i2aster     (i_raster,j_raster))         
         allocate(veg         (1:nc,1:nr))
         allocate(catid       (1:nc,1:nr))
         allocate(catid_index (1:nc,1:nr))          
         allocate(tileid_index(1:nc,1:nr))

         dx  = 360._8/nc
         dy  = 180._8/nr
         d2r = PI/180._8
         da  = MAPL_radius*MAPL_radius*pi*pi*dx*dy/180./180./1000000.    
         
         tileid_index = 0        

         !  Simple Biosphere 2 Model Legend 
         !  Value Class Name 
         !  (ftp://edcftp.cr.usgs.gov/pub/data/glcc/globe/latlon/sib22_0.leg)
         !  the types vary 0-11 (array index minus 1) 
         
         veg_class(1)  = 'Ocean'
         veg_class(2)  = 'Broadleaf Evergreen Trees' 
         veg_class(3)  = 'Broadleaf Deciduous Trees' 
         veg_class(4)  = 'Broadleaf and Needleleaf Trees' 
         veg_class(5)  = 'Needleleaf Evergreen Trees' 
         veg_class(6)  = 'Needleleaf Deciduous Trees' 
         veg_class(7)  = 'Short Vegetation/C4 Grassland'
         veg_class(8)  = 'Shrubs with Bare Soil' 
         veg_class(9)  = 'Dwarf Trees and Shrubs' 
         veg_class(10) = 'Agriculture or C3 Grassland' 
         veg_class(11) = 'Water, Wetlands'
         veg_class(12) = 'Ice/Snow'
         
         ! reading SiB2 land cover classification data - the origin of the 
         ! 2.5'x2.5' vegetation raster file is global 1min IGBP data 
         ! (ftp://edcftp.cr.usgs.gov/pub/data/glcc/globe/latlon/sib22_0.leg)
         
         open (10,file='data/CATCH/sib22.5_v2.0.dat', &
              form='unformatted', &
              action='read', convert='big_endian',status='old')
         
         READ(10)i2aster
         
         close (10,status='keep')
         
         if(regrid) then
            call RegridRaster1 (i2aster,veg)
         else
            veg = i2aster
         endif
         
         deallocate (i2aster)
         
         !   reading 2.5'x2.5' global raster file of Pfafstetter Catchment IDs
         !   In this version, the dateline has been overlaid over the catchments those straddle 
         !   across. The numbers contain for
         !    1 global ocean catchment                : Pfafstetter ID 0
         !    36716 global land catchments            : Pfafstetter IDs 1000-5999900
         !    1 global inland water (lakes) catchment : Pfafstetter ID 6190000
         !    1 global ice catchment                  : Pfafstetter ID 6200000
         
         open (10,file='data/CATCH/global.cat_id.catch.DL', form='formatted', &
              action='read', status='old')!
         
         do j=1,j_raster
            read(10,*)(iaster(i,j),i=1,i_raster)
         end do
         
         close (10,status='keep')
         
         if(regrid) then
            call RegridRaster(iaster,catid)
         else
            catid =  iaster
         endif
         
         print *,'Read global.cat_id.catch.DL' 
         print *,'Min and Max of Pfafstetter IDs:', minval(catid),maxval(catid)
         
         ! reading the 2.5'x2.5' global raster file of tile indices for the 
         !  above Pfafstetter Catchments
         !  1 global ocean catchment                : tile_index 36719
         !  36716 global land catchments            : tile_index 1-36716
         !  1 global inland water (lakes) catchment : tile_index 36717
         !  1 global ice catchment                  : tile_index 36718
         ! ------------------------------------------------------------
         
         open (10,file='data/CATCH/'  &
              //'PfafstatterDL.rst', form='unformatted',        &
              action='read',convert='little_endian', status='old')
         
         do j=1,j_raster
            read(10)(iaster(i,j),i=1,i_raster)
         end do
         
         close (10,status='keep')
         
         if(regrid) then
            call RegridRaster(iaster,catid_index)
         else
            catid_index =  iaster
         endif
         
         deallocate (iaster)
         
         print *,'Read PfafstatterDL.rst' 
         print *,'Min and Max of tile indices:',minval(catid_index),maxval(catid_index)
 
         ! While looping through the nc x nr grid (tile raster), this section counts # of  
         ! SMAP grid cells that contain land, ice or water, seperately.
         ! Each SMAP grid cell is assigned with an ID = ind_row*ND +  ind_col. 
         ! This is just the prelimiminery assessment in the process of assigning separate  
         !     tiles for land, water and ice fractions within the SMAP Grid cell
         ! The program checks each nc x nr pixels whether there is a SMAP grid cell underneath, and counts
         ! number of water, land and ice pixels as seen on veg raster.
         ! -----------------------------------------------------------------------------------------------
         
         
         do i = 1 ,nc
            
            clon = -180. + float(i-1)*dx + dx/2.
            
            do j =nr ,1 ,-1
               
               clat = -90. + float(j-1)*dy + dy/2.
               call easeV2_convert(trim(MGRID), clat, clon, r_smap, s_smap)
               
               ind_col = nint(r_smap) + 1 
               ind_row = nint(s_smap) + 1
               
               if((ind_row.ge.1).and.(veg(i,j).ne.OceanType).and.(ind_row.le.nr_smap)) then
                  l=  ind_row*ND +  ind_col
                  
                  if(veg(i,j)==LakeType) then 
                     water_id(l) = 1
                  else if(veg(i,j)==IceType) then
                     ice_id  (l) = 1
                  else
                     land_id (l) = 1
                  endif
               endif
            end do
         end do

      endif
      
      ! Reading SRTM elevation data - to be consistent with AGCM
      ! --------------------------------------------------------     
      
      allocate(raster      (i_raster,j_raster))
      allocate(q0(nc,nr)) 
      
      gtopo30 = 'data/CATCH/srtm30_withKMS_2.5x2.5min.data'
     
      open (10,file=trim(gtopo30),form='unformatted',status='old',convert='little_endian')
      read (10) raster
      close (10,status='keep') 
      
      if(regrid) then
         call RegridRasterReal(raster,q0)
      else
         q0 =  raster
      endif
      
      deallocate (raster)
      
      print *,'# of Land  pixels in SMAP: ',sum (land_id)
      print *,'# of water pixels in SMAP: ',sum (water_id)
      print *,'# of ice   pixels in SMAP: ',sum (ice_id)

      l_index=0
      w_index=sum (land_id)
      i_index=sum (land_id) + sum (water_id)
      nland_cells = w_index


      allocate(tile_area (1:i_index + sum (ice_id)))
      allocate(smap_grid_area     (1:NT))
      allocate(tile_ele      (1:w_index))
      allocate(tile_area_land(1:w_index)) 
      allocate(my_land       (1:i_index + sum (ice_id)))
      allocate(all_id        (1:i_index + sum (ice_id)))

      land_id = 0
      water_id= 0
      ice_id  = 0

      my_land    = 0
      all_id     = 0
      smap_grid_area = 0. 
      tile_area_land = 0.       
      tile_ele       = 0.
      tile_area      = 0.

      ! While looping through the nc x nr grid, this section derives land, ice and water tiles.
      ! Each SMAP grid cell is assigned with an ID = ind_row*ND +  ind_col 
      !         ind_col, ind_row are overlying SMAP grid cell indices 
      ! Based on the above sums: 
      !         l_index Grid cells have land fractions (sum(land_id)) 
      !         w_index SMAP Grid cells have inland water fractions (sum(water_id))
      !         i_index SMAP Grid cells have ice fractions (sum(ice_id))
      ! hence, tile_index        1                     to l_index                     represent land tiles  
      !         tile_index       l_index +1            to l_index + w_index           represent water (lakes) tiles  
      !         tile_index       l_index + w_index +1  to l_index + w_index + i_index represent ice tiles
      ! global nc x nr array of tileid_index(nc,nr) contains corresponding tile_index values which 
      !        is derived in the below loop
 
      ND_raster = 10*10**(nint(log10(1.*NR)))
      i2 = 1

      do i = 1 ,nc
         
         clon = -180. + float(i-1)*dx + dx/2.
         
         do j =nr ,1 ,-1
            lats = -90._8 + (j - 0.5_8)*dy
            clat = -90. + float(j-1)*dy + dy/2.
            call easeV2_convert(trim(MGRID), clat, clon, r_smap, s_smap)
            
            ind_col = nint(r_smap) + 1 
            ind_row = nint(s_smap) + 1

            l=  ind_row*ND +  ind_col
            pix_area =(sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)

            if((ind_row.ge.1).and.(veg(i,j).ge.1).and.(ind_row.le.nr_smap)) then
                                           
               if(veg(i,j)==LakeType) then
                  if(water_id(l)==0) then
                     w_index = w_index + 1
                     water_id(l) = w_index  
                     tileid_index(i,j)= water_id(l)
                  else
                     tileid_index(i,j)= water_id(l)
                  endif
               endif
               
               if(veg(i,j)==IceType) then
                  if(ice_id(l)==0) then
                     i_index = i_index + 1
                     ice_id  (l) = i_index
                     tileid_index(i,j)= ice_id  (l)  !i_index
                  else
                     tileid_index(i,j)= ice_id  (l)  !i_index
                  endif
               endif
               
               if(veg(i,j).lt.LakeType) then
                  if(land_id(l)==0) then
                     l_index = l_index + 1     
                     land_id (l) = l_index        
                     tileid_index(i,j)= land_id (l) !1-l_index
                  else
                     tileid_index(i,j)= land_id (l) !1-l_index
                  endif
 
               endif
               tile_area(tileid_index(i,j))= tile_area(tileid_index(i,j)) + &
                    pix_area     
               my_land(tileid_index(i,j)) = l
               all_id (tileid_index(i,j)) = j*ND_raster + i 
            endif

            if((ind_row.ge.1).and.(ind_row.le.nr_smap)) then  
               smap_grid_area(l) = smap_grid_area(l) +             &
                    pix_area
            endif

            ! computing tile area/elevation
            ! -----------------------------

           if((tileid_index(i,j) > 0).and.(tileid_index(i,j) <= nland_cells))then
               tile_ele(tileid_index(i,j))       = tile_ele(tileid_index(i,j)) + q0(i,j) *  &	
                    pix_area
               tile_area_land(tileid_index(i,j)) = tile_area_land(tileid_index(i,j))     +  &	
                    pix_area
            endif  
         end do
      end do
      
      deallocate(land_id, q0)
      deallocate(water_id)
      deallocate(ice_id  )

      tile_ele = tile_ele/tile_area_land

      ! adjustment Global Mean Topography to 614.649 (615.662 GTOPO 30) m
      ! -----------------------------------------------------------------
      sum1=0.
      sum2=0. 

      do j=1,l_index	 
         sum1 = sum1 + tile_ele(j)*tile_area(j)
      enddo

      if(sum1/sum(tile_area(1:l_index)).ne. 614.649D0 ) then
         print *,'Global Mean Elevation (over land): ', sum1/sum(tile_area(1:l_index))
         tile_ele =tile_ele*(614.649D0 / (sum1/sum(tile_area(1:l_index))))				  	
         sum1=0.
         sum2=0. 
         do j=1,l_index	 
            sum1 = sum1 + tile_ele(j)*tile_area(j)
         enddo
         print *,'Global Mean Elevation after scaling to SRTM : ',sum1/sum(tile_area(1:l_index))
      endif
              print *,'Total Land Area :', sum(tile_area(1:l_index))* MAPL_RADIUS * MAPL_RADIUS/1000./1000., &
                   sum(tile_area_land(1:l_index))* MAPL_RADIUS * MAPL_RADIUS/1000./1000.

      print *,'Creating ... ', trim(gfile)//'rst'

     !-------------------------------------------

      open (10, file ='rst/'//trim(gfile)//'.rst',form='unformatted',status='unknown',  &
           action='write')
      
      do j=1,nr
         write(10)(tileid_index(i,j),i=1,nc)
      end do
      
      close (10,status='keep')

      print *,'Creating ... ', trim(gfile)//'til ,catchment.def'
 
    !-----------------------------------------------------------

      open (11,file='clsm/catchment.def',     &
           form='formatted',status='unknown')
      write(11,*)l_index

      open  (10, file ='til/'//trim(gfile)//'.til',form='formatted',status='unknown',action='write')
      write (10,*)i_index, nc, nr
      write (10,*)1
      write (10,*)'SMAP-EASEv2-'//trim(MGRID)
      write (10,*)nc_smap
      write (10,*)nr_smap
      write (10,*)'NO-OCEAN'
      write (10,*) -9999
      write (10,*) -9999      

      do l=1,i_index

         ig    = my_land(l)-ND*(my_land(l)/ND)
         jg    = my_land(l)/ND
  
         cindex= catid_index(all_id(l)-ND_raster*(all_id(l)/ND_raster),all_id(l)/ND_raster)

         if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
            pfaf = cindex
         else
            pfaf = catid(all_id(l)-ND_raster*(all_id(l)/ND_raster),all_id(l)/ND_raster)
         endif

         if ((l > l_index).and.(l <= w_index)) typ =19
         if (l > w_index) typ = 20

         if (l <= l_index) then 
            typ = 100
            call easeV2_inverse (trim(MGRID), real(ig-1),real(jg-1), clat, clon) 
            
            mnx = clon - 180./real(nc_smap)
            mxx = clon + 180./real(nc_smap)
            
            jgv = real(jg-1) + 0.5
            
            call easeV2_inverse (trim(MGRID), real(ig-1),jgv, clat, clon) 

            mny = clat
         
            jgv = real(jg-1) - 0.5
         
            call easeV2_inverse (trim(MGRID), real(ig-1),jgv, clat, clon) 

            mxy = clat 

            write (11,'(i8,i8,5(2x,f9.4), i4)')l,pfaf,mnx,mxx,mny,mxy,tile_ele(l)

         endif

         call easeV2_inverse (trim(MGRID), real(ig-1),  real(jg-1), clat, clon)
         
         fr_gcm= tile_area(l)/smap_grid_area(jg*ND +  ig)

         if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
            write(10,'(i10,i9,2f10.4,2i5,f19.12,i10,i15,e13.4)') &
                 typ,pfaf,clon,clat,ig-1,jg-1,fr_gcm ,pfaf,SRTM_catid(cindex) 
         else
            write(10,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i8)') &
                 typ,pfaf,clon,clat,ig-1,jg-1,fr_gcm ,cindex 
         endif
         
      end do
      
      close (10,status='keep')      
      close (11,status='keep')          

      deallocate (tileid_index,catid_index,veg)
      deallocate (tile_area, smap_grid_area, tile_ele, tile_area_land, my_land, all_id)
 
      if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then         

         print *,'Creating SMAP-Catch_TransferData.nc files.'

         !---------------------------------------------------

         deallocate (SRTM_CatchArea, SRTM_catid)
         
      endif
      
      ! create Grid2Catch transfer file
      ! -------------------------------

      CALL CREATE_ROUT_PARA_FILE (NC, NR, trim(gfile), MGRID=MGRID)  
      
      ! now run mkCatchParam
      ! --------------------

      tmpstring1 = '-e EASE -g '//trim(gfile)
      write(tmpstring2,'(2(a2,x,i5,x))')'-x',nc,'-y',nr
      tmpstring = 'bin/mkCatchParam_openmp '//trim(tmpstring2)//' '//trim(tmpstring1)
      print *,trim(tmpstring)
      
      call system (tmpstring)   

   END PROGRAM mkSMAPTilesPara_v2

