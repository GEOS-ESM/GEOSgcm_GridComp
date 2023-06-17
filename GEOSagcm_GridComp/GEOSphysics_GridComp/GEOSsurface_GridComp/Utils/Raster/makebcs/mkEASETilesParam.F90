#define I_AM_MAIN
#include "MAPL_ErrLog.h"

PROGRAM mkEASETilesParam 
  
  ! This program constructs land, landice, and lake tiles for EASE grid tile spaces such as
  !  those used for the SMAP Level-4 products and other offline projects

  ! This program resulted from the merger and cleanup of mkSMAPTilesPara.F90
  !  and mkSMAPTilesPara_v2.F90 in September 2022.
  ! Before the merger and cleanup, the EASE grid parameters were hard-coded here.
  !  For EASEv2 M25, the outdated scale value was used here.
  ! The program was renamed to mkEASETileParam from mkSMAPTilesPara_v2
  !
  ! - wjiang + reichle, 21 Sep 2022
  
  ! additional major cleanup (reichle, 15 Jun 2023)
  ! - rewrote some code blocks for clarity and efficiency
  ! - removed obsolete variable declarations
  ! - removed obsolete use statements
  ! - added "only:" qualifier to some use statement
  ! - commented out obsolete code blocks
  ! - renamed some variables for clarity
  ! - removed repetition of identical operations
  ! - added comments
  ! - white-space changes for improved readability
  
      use EASE_conv,          only : EASE_extent, EASE_convert, EASE_inverse
      use rmTinyCatchParaMod, only : i_raster, j_raster, SRTM_maxcat 
      use rmTinyCatchParaMod, only : RegridRasterReal                  
      use rmTinyCatchParaMod, only : MAKE_BCS_INPUT_DIR 
      use process_hres_data,  only : histogram
      use MAPL_SortMod
      use MAPL_ConstantsMod
      use MAPL_ExceptionHandling
      use netcdf
      
      implicit none
      
      integer, parameter :: nc_esa = 129600       ! number of cols in 10-arcsec ESA mask file
      integer, parameter :: nr_esa =  64800       ! number of rows in 10-arcsec ESA mask file
      
      ! define tile types used for processing here (values may be from ESA mask?) 
      
      integer, parameter :: OceanType  =  0  
      integer, parameter :: LandType   =  1 ! land    type used for processing here; in GEOS, land    tiles are type=100
      integer, parameter :: LakeType   = 10 ! lake    type used for processing here; in GEOS, lake    tiles are type= 19
      integer, parameter :: IceType    = 11 ! landice type used for processing here; in GEOS, landice tiles are type= 20

      integer :: i, j, ig, jg, n
      integer :: NC, NR, N_ease_grid_cells, NDND, ND_raster
      
      ! For regridding

      integer,     allocatable, dimension(:,:), target   :: geos_msk 

      REAL,        allocatable, DIMENSION(:)             :: loc_val
      INTEGER,     ALLOCATABLE, DIMENSION(:)             :: density, loc_int
      logical,     allocatable, dimension(:)             :: unq_mask   
      integer,                  dimension(:,:), pointer  :: subset

      integer                                            :: dx_esa, dy_esa, NBINS, NPLUS

      integer*8,   allocatable, dimension(:)             :: SRTM_catid
      real(kind=8),allocatable, dimension(:)             :: SRTM_catid_r8

      integer,     allocatable, dimension(:,:), target   :: tileid_index, catid_index

      ! integer,     allocatable, dimension(:,:)           :: catid, iaster
      
      integer,     allocatable, dimension(:)             :: land_id, water_id, ice_id
      integer,     allocatable, dimension(:)             :: my_land, all_id
      real,        allocatable, dimension(:)             :: ease_grid_area, tile_area    !, SRTM_CatchArea 
      integer*1,   allocatable, dimension(:,:)           :: veg                          !, i2aster
      real*4,      allocatable, dimension(:,:)           :: q0, raster
      REAL,        allocatable, dimension(:)             :: tile_ele
      
      !INTEGER*8    :: PFAF_CODE  

      integer      :: l, l_index, i_index, w_index, typ, pfaf, cindex
      real         :: clat, clon, r_ease, s_ease
      real         :: fr_gcm
      integer      :: ind_col, ind_row, status, ncid, varid, nciv                        !, DOM_INDX
      integer      :: n_land, n_lake, n_landice, n_landlakelandice
      
      ! REAL (kind=8), PARAMETER :: RADIUS=6378137.0,pi=3.14159265358979323846

      !character*100          :: veg_class(12)

      character*200          :: gfile, gtopo30
      integer                :: nc_ease, nr_ease, N_args, command_argument_count 
      REAL                   :: dx, dy, d2r, lats, mnx, mxx, mny, mxy, jgv, pix_area     ! VDUM
      real                   :: dx_ease, mean_land_ele
      character(40)          :: arg, EASElabel_ 

      character(len=:), allocatable :: EASElabel 

      !character*200          :: tmpstring, tmpstring1, tmpstring2	      

      logical                :: regrid   = .false.
      character*128          :: MaskFile
      
      !logical                :: pfaf_til = .false.
      
      character*1            :: PF
      
      !character(len=6)       :: EASE_Version
      
      character(len=10)      :: nc_string, nr_string
      character(len=128)     :: usage1, usage2
      character(len=128)     :: Iam = "mkEASETilesParam"

      character(len=512)     :: fname_mask

      ! --------------------------------------------------------------------------------------

      call get_environment_variable( "MAKE_BCS_INPUT_DIR", MAKE_BCS_INPUT_DIR )

      usage1 = 'USAGE : bin/mkEASETilesParam.x -ease_label EASELabel                    '
      usage2 = '        where EASELabel = *EASEv[x]_M[yy]*, x={1,2}, yy={01,03,09,25,36}'

      N_args = command_argument_count()

      if(N_args < 1) then
        print *,trim(usage1)
        print *,trim(usage2)
        stop
      end if

      i=0      
      do while ( i < N_args )

         i = i+1
         
         call get_command_argument(i,arg)
         
         if     ( trim(arg) == '-ease_label' ) then
            i = i+1
            call get_command_argument(i,EASELabel_)

         ! WY noted: this may be used in the future for irrigation tiles
         !elseif ( trim(arg) == '-pfaf_til' ) then
         !   i = i+1
         !   call get_command_argument(i,PF)
         !   if (PF == 'T') pfaf_til = .true.

         else ! stop for any other arguments
            print *,trim(usage1)
            print *,trim(usage2)
            stop
         endif
         
      end do
      
      ! Get EASE Grid specifications
      ! --------------------------------
      
      EASElabel = trim(EASELabel_)
      
      call ease_extent( EASELabel, nc_ease, nr_ease )
      
      write(nc_string, '(i0)') nc_ease
      write(nr_string, '(i0)') nr_ease
      
      gfile = trim(EASElabel)//'_'//trim(nc_string)//'x'//trim(nr_string)
      
!      ! get appropriate number for length of land_id, water_id, ice_id vectors
!
!      if      (index(EASELabel,'M01') /=0) then ! EASE   1 km grid
!         N_ease_grid_cells = 1500000000
!      else if (index(EASELabel,'M03') /=0) then ! EASE   3 km grid
!         N_ease_grid_cells = 1500000000
!      else if (index(EASELabel,'M09') /=0) then ! EASE   9 km grid
!         N_ease_grid_cells = 16330000
!      else if (index(EASELabel,'M25') /=0) then ! EASE  25 km grid
!         N_ease_grid_cells = 16330000
!      else if (index(EASELabel,'M36') /=0) then ! EASE  36 km grid
!         N_ease_grid_cells = 16330000
!      else
!         print *,"ERROR: Unknown EASELabel, stopping."
!         stop
!      endif

      N_ease_grid_cells = nc_ease * nr_ease
      
      allocate(land_id (1:N_ease_grid_cells))
      allocate(water_id(1:N_ease_grid_cells))
      allocate(ice_id  (1:N_ease_grid_cells))

      land_id     =  0
      water_id    =  0
      ice_id      =  0             

      ! Step size for conversion of 2-dim indexing into 1-dim indexing:  l = ind_row*NDND + ind_col
      ! If conversion was simply to get from 2-dim to 1-dim indexing for EASE grid cells, 
      !  NDND=nc_ease would suffice
      ! Here, NDND is computed as a power of 10, leaving room for "empty" indices that may or may not be needed later.
      ! - reichle, 15 Jun 2023

      !NDND        = 10*10**(nint(log10(1.*nr_ease)))
      
      NDND = nc_ease

      !   Check for the 10 arc-sec MaskFile
      ! -----------------------------------
      
      call get_environment_variable ("MASKFILE"        ,MaskFile        )

      print *, 'Using MaskFile ', trim(MaskFile)
      
      !   This section was used to make Irrigated Tiles 
      !if(pfaf_til)  then

      !   nc = 43200  ! Number of rows in raster file
      !   nr = 21600
      !   call mkEASEv2Raster
         
      !else
      !   if((trim(MGRID) == 'M09').or.(trim(MGRID) == 'M36'))call write_tilfile 
      !endif
      
      if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then         
         
         ! New ESA (Veg) + SRTM (catchments) based mask file
         ! is overlaid on the EASE grid
         ! -------------------------------------------------
         
         ! ESA/SRTM mask file is 10-arcsec resolution;  coarsen here to 30-arcsec raster grid:

         ! NOTE: Only coastlines (and lake shorelines?) are at 10-arcsec resolution.
         !       The original Pfaf catchments are at 1-arcmin resolution (~2 km).

         nc = 43200  ! Number of rows    in raster file
         nr = 21600  ! Number of columns in raster file

         regrid = .true.

         dx_esa = nc_esa / nc ! =3  (# of ESA columns within a raster grid cell in x-dim)
         dy_esa = nr_esa / nr ! =3  (# of ESA rows    within a raster grid cell in y-dim)

         allocate(tileid_index(1:nc,1:nr))
         allocate(SRTM_catid    (1:SRTM_maxcat+2)               )
         allocate(SRTM_catid_r8 (1:SRTM_maxcat+2), source = 0.d0)
         allocate(catid_index (1:nc,1:nr))          
         allocate(veg         (1:nc,1:nr))
         allocate(geos_msk    (1:nc_esa,1:dy_esa))

         ! the following block is not needed (perhaps for routing)
         
!         allocate(SRTM_CatchArea(1:SRTM_maxcat))
!
!         OPEN (10, FILE = trim(MAKE_BCS_INPUT_DIR)//'/land/topo/v1/SRTM-TopoData/Pfafcatch-routing.dat', &
!              FORM = 'FORMATTED',STATUS='OLD',ACTION='READ') 
!
!         READ (10,*) I
!         DO N = 1, I
!            READ (10, '(i8,i15,4(1x,f9.4),1x,e10.3,4(1x,e9.3),I8,6(1x,f9.4))')    &
!                 DOM_INDX,PFAF_CODE,VDUM,VDUM,VDUM,VDUM,VDUM,                     &
!                 SRTM_CatchArea (N)
!         END DO
!         CLOSE (10, STATUS='KEEP')

         dx  = 360._8/nc               ! raster grid spacing (dlon)
         dy  = 180._8/nr               ! raster grid spacing (dlat)

         d2r = MAPL_PI_R8/180._8       ! degree-to-radians conversion factor -- d2r declared as REAL ?!?!?!

         !da  = MAPL_radius*MAPL_radius*pi*pi*dx*dy/180./180./1000000.    
         
         tileid_index = 0        
         catid_index  = 0
         veg          = OceanType      ! initialize to ocean (type=0)
         
         ! read list of Pfafstetter catchment IDs ('PfafID') from 10-arcsec mask file into variable 'SRTM_catid[_r8]'
         ! [vector of length SRTM_maxcat]

         fname_mask = trim(MAKE_BCS_INPUT_DIR) // '/shared/mask/GEOS5_10arcsec_mask.nc'

         print *, 'Opening ', trim(fname_mask)

         status    = NF90_OPEN( fname_mask, NF90_NOWRITE, ncid )
         if(status /=0) then
            PRINT *, trim(NF90_STRERROR(STATUS))
            print *, 'Problem with NF90_OPEN():  ', trim(fname_mask), '  ', ncid, status
         else
            print *, 'ncid=', ncid
         endif
         
         status    = nf90_inq_varid( ncid, name='PfafID', varid=varid )
         if(status /=0) then
            PRINT *, trim(NF90_STRERROR(STATUS))
            print *, 'Problem with NF90_INQ_VARID():  PfafID  ', ncid, varid, status
         endif

         status    = nf90_get_var( ncid, varid, SRTM_catid_r8, (/1/),(/SRTM_maxcat/) )   
         if(status /=0) then
            PRINT *, trim(NF90_STRERROR(STATUS))
            print *, 'Problem with NF90_GET_VAR():  PfafID  ', ncid, varid, SRTM_maxcat, status
         endif

         SRTM_catid = int8(SRTM_catid_r8)                         ! convert data to integer*8    -- contains 12-digit Pfaf code
         SRTM_catid (SRTM_maxcat + 1) = 190000000                 ! append ID for Lake type
         SRTM_catid (SRTM_maxcat + 2) = 200000000                 ! append ID for Landice type
         
         ! -----------------------------------------------
         !
         ! first loop through 30-arcsec raster grid cells 
         ! - aggregate mask from 10-arcsec resolution to 30-arcsec raster grid by determining dominant tile type
         ! - for each EASE grid cell, determine tile types based on aggregated (30-arcsec) mask

         !i1 = 0  ! count # of 30-arcsec pixels - NOT USED

         status = nf90_inq_varid(ncid, name='CatchIndex', varid=varid)
         if(status /=0) then
            PRINT *, trim(NF90_STRERROR(STATUS))
            print *, 'Problem with NF90_INQ_VARID():  CatchIndex  ', ncid, varid, status
         endif

         do j=1,nr
            
            clat = -90. + float(j-1)*dy + dy/2.              ! center lat of raster grid cell (*,j)
                        
            ! read slice of variable 'CatchIndex' from 10-arcsec mask file into variable 'geos_msk' 
            ! [2d-array: 129600-by-3]

            status = NF90_GET_VAR( ncid, varid, geos_msk, (/1,(j-1)*dy_esa+1/), (/nc_esa,dy_esa/) ) ! Read 10-arcsec rows that lie within the raster row 'j'  
            !status  = NF_GET_VARA_INT (ncid,4,(/1,(j-1)*dy_esa +1/),(/nc_esa,dy_esa/),geos_msk) ! Read 10-arcsec rows that lie within the raster row 'j'  
            
            if(status /=0) then
               PRINT *, trim(NF90_STRERROR(STATUS))
               print *, 'Problem with NF_GET_VAR():  CatchIndex  ', ncid, varid, j, dy_esa, nc_esa, status
            endif
            
            do i = 1,nc    

               clon = -180. + float(i-1)*dx + dx/2.          ! center lon of raster grid cell (i,*)

               ! extract [3-by-3] subset of ESA/SRTM 10-arcsec mask file that corresponds to 30-arcsec raster grid cell (i,j)

               if (associated (subset)) NULLIFY (subset)

               subset => geos_msk ((i-1)*dx_esa + 1 : i*dx_esa, 1:dy_esa) ! rectangular array contains ESA pixels that lie within the raster grid cell at i,j
               
               if(maxval (subset) > SRTM_maxcat) then
                  where (subset == 190000000) subset = SRTM_maxcat + 1    ! Lake type      (convert ID from 190000000 to SRTM_maxcat+1)
                  where (subset == 200000000) subset = SRTM_maxcat + 2    ! Landice type   (convert ID from 200000000 to SRTM_maxcat+2)
               endif

               if (maxval(subset) > 0) then  ! check if there are Non-ocean ESA pixels 

                  ! raster grid cell has at least one 10-arcsec pixel that is land or landice or lake
                  !
                  ! now find out the *dominant* Pfafstetter ID (could be Lake or Landice) within the raster grid cell
                  
                  NPLUS = count( subset>=1 .and. subset<=SRTM_maxcat+2 )                     ! Count non-ocean ESA pixels within raster grid cell  
                  allocate(loc_int (1:NPLUS))
                  allocate(unq_mask(1:NPLUS))
                  loc_int = pack( subset, mask = ( subset>= 1 .and. subset<=SRTM_maxcat+2) ) ! loc_int contains catch_indices of non-ocean ESA pixels 
                  call MAPL_Sort(loc_int)
                  unq_mask = .true.
                  do n=2,NPLUS 
                     unq_mask(n) = .not. (loc_int(n)==loc_int(n-1))                          ! count number of unique numbers in loc_int for binning
                  end do
                  NBINS = count(unq_mask)
                  
                  if (NBINS > 1) then
                     allocate(loc_val(1:NBINS))
                     allocate(density(1:NBINS))
                     loc_val = 1.*pack(loc_int,mask =unq_mask)                               ! loc_val contains available non-ocean catch_indices within the i,j grid cell,
                                                                                             ! Those numbers will be used as bin values
                     call histogram( dx_esa*dy_esa, NBINS, density, loc_val, real(subset) )  ! density is the pixel count for each bin value
                     catid_index(i,j) = loc_val(maxloc(density,1))                           ! picks maximum density as the dominant catchment_index at i,j
                     deallocate(loc_val, density)
                  else
                     catid_index(i,j) = loc_int (1)                                          
                  endif
                  deallocate(loc_int, unq_mask)

                  ! now catid_index(i,j) = dominant PfafID (or LakeID or LandiceID) in raster grid cell (i,j)
                                    
                  if (catid_index (i,j) <= SRTM_maxcat    ) veg(i,j) = LandType
                  if (catid_index (i,j) == SRTM_maxcat + 1) veg(i,j) = LakeType
                  if (catid_index (i,j) == SRTM_maxcat + 2) veg(i,j) = IceType

                  !if( (catid_index(i,j)>= 1) .and. (catid_index(i,j)<=SRTM_maxcat) ) i1 = i1 + 1        ! i1 = counter for raster grid cells with land as dominant type -- NOT USED?

                  ! set land_id or water_id or ice_id of EASE grid cell to 1 accordingly
                  ! --> EASE grid cell may contain more than one tile type
                  
                  ! count in if this is i,j pixel is a land, lake or ice within ind_col,ind_row EASE grid cell
                  
                  ! get 1-based ind_col and ind_row indices of EASE grid cell that contains raster grid cell (i,j)
                  
                  call EASE_convert(EASELabel, clat, clon, r_ease, s_ease)
                  
                  ind_col = nint(r_ease) + 1 
                  ind_row = nint(s_ease) + 1     ! can be negative or greater than nr_ease (lat near N/S pole)
                  
                  if( (veg(i,j).ne.OceanType) .and. (ind_row.ge.1) .and. (ind_row.le.nr_ease) ) then
                     
                     ! raster grid cell has tile type other than ocean and is located within the lat band covered
                     !  by the EASE grid (approx 85S-85N): 1 <= ind_row <= nr_ease
                     
                     ! 1-dim indexing used here appears wasteful; why not ll=(ind_row-1)*nc_ease+ind_col, with NT=nc_ease*nr_ease ??
                     
                     !l = ind_row*NDND + ind_col
                     
                     l = (ind_row-1)*NDND + ind_col  ! 1-dim index for all EASE grid cells
                     
                     if(     veg(i,j)==LakeType) then 
                        water_id(l) = 1
                     else if(veg(i,j)==IceType ) then
                        ice_id  (l) = 1
                     else
                        land_id (l) = 1
                     endif
                  endif
               endif
               
            end do   ! i=1,nc
         enddo       ! j=1,nr
         
         status    = NF90_CLOSE (ncid)  
         deallocate (geos_msk)

         print *,'Done reading ', trim(MaskFile) 
         print *,'Min and Max of tile indices:', minval(catid_index), maxval(catid_index)

      else

         print *,'MaskFile = ', trim(MaskFile)
         print *,'ERROR: Selected mask file not supported for creating tiles on the EASE grid, stopping.'
         stop

!         ! Old IGBP (Veg) + HYDRO1k (catchments) based mask will
!         ! Overlaid on EASE mask
!         ! -----------------------------------------------------
!         
!         allocate(iaster      (i_raster,j_raster)) 
!         allocate(i2aster     (i_raster,j_raster))         
!         allocate(veg         (1:nc,1:nr))
!         allocate(catid       (1:nc,1:nr))
!         allocate(catid_index (1:nc,1:nr))          
!         allocate(tileid_index(1:nc,1:nr))
!
!         dx  = 360._8/nc
!         dy  = 180._8/nr
!         d2r = MAPL_PI_R8/180._8
!         !da  = MAPL_radius*MAPL_radius*pi*pi*dx*dy/180./180./1000000.    
!         
!         tileid_index = 0        
!
!         !  Simple Biosphere 2 Model Legend 
!         !  Value Class Name 
!         !  (ftp://edcftp.cr.usgs.gov/pub/data/glcc/globe/latlon/sib22_0.leg)
!         !  the types vary 0-11 (array index minus 1) 
!         
!         veg_class(1)  = 'Ocean'
!         veg_class(2)  = 'Broadleaf Evergreen Trees' 
!         veg_class(3)  = 'Broadleaf Deciduous Trees' 
!         veg_class(4)  = 'Broadleaf and Needleleaf Trees' 
!         veg_class(5)  = 'Needleleaf Evergreen Trees' 
!         veg_class(6)  = 'Needleleaf Deciduous Trees' 
!         veg_class(7)  = 'Short Vegetation/C4 Grassland'
!         veg_class(8)  = 'Shrubs with Bare Soil' 
!         veg_class(9)  = 'Dwarf Trees and Shrubs' 
!         veg_class(10) = 'Agriculture or C3 Grassland' 
!         veg_class(11) = 'Water, Wetlands'
!         veg_class(12) = 'Ice/Snow'
!         
!         ! reading SiB2 land cover classification data - the origin of the 
!         ! 2.5'x2.5' vegetation raster file is global 1min IGBP data 
!         ! (ftp://edcftp.cr.usgs.gov/pub/data/glcc/globe/latlon/sib22_0.leg)
!         
!         open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/veg/pft/v1/sib22.5_v2.0.dat', &
!              form='unformatted', &
!              action='read', convert='big_endian',status='old')
!         
!         READ(10)i2aster
!         
!         close (10,status='keep')
!         
!         if(regrid) then
!            call RegridRaster1 (i2aster,veg)
!         else
!            veg = i2aster
!         endif
!         
!         deallocate (i2aster)
!         
!         !   reading 2.5'x2.5' global raster file of Pfafstetter Catchment IDs
!         !   In this version, the dateline has been overlaid over the catchments those straddle 
!         !   across. The numbers contain for
!         !    1 global ocean catchment                : Pfafstetter ID 0
!         !    36716 global land catchments            : Pfafstetter IDs 1000-5999900
!         !    1 global inland water (lakes) catchment : Pfafstetter ID 6190000
!         !    1 global ice catchment                  : Pfafstetter ID 6200000
!         
!         open (10,file= trim(MAKE_BCS_INPUT_DIR)//'/shared/mask/global.cat_id.catch.DL', form='formatted', &
!              action='read', status='old')!
!         
!         do j=1,j_raster
!            read(10,*)(iaster(i,j),i=1,i_raster)
!         end do
!         
!         close (10,status='keep')
!         
!         if(regrid) then
!            call RegridRaster(iaster,catid)
!         else
!            catid =  iaster
!         endif
!         
!         print *,'Read global.cat_id.catch.DL' 
!         print *,'Min and Max of Pfafstetter IDs:', minval(catid),maxval(catid)
!         
!         ! reading the 2.5'x2.5' global raster file of tile indices for the 
!         !  above Pfafstetter Catchments
!         !  1 global ocean catchment                : tile_index 36719
!         !  36716 global land catchments            : tile_index 1-36716
!         !  1 global inland water (lakes) catchment : tile_index 36717
!         !  1 global ice catchment                  : tile_index 36718
!         ! ------------------------------------------------------------
!         
!         open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/topo/'  &
!              //'PfafstatterDL.rst', form='unformatted',        &
!              action='read',convert='little_endian', status='old')
!         
!         do j=1,j_raster
!            read(10)(iaster(i,j),i=1,i_raster)
!         end do
!         
!         close (10,status='keep')
!         
!         if(regrid) then
!            call RegridRaster(iaster,catid_index)
!         else
!            catid_index =  iaster
!         endif
!         
!         deallocate (iaster)
!         
!         print *,'Read PfafstatterDL.rst' 
!         print *,'Min and Max of tile indices:',minval(catid_index),maxval(catid_index)
! 
!         ! While looping through the nc x nr grid (tile raster), this section counts # of  
!         ! EASE grid cells that contain land, ice or water, seperately.
!         ! Each EASE grid cell is assigned with an ID = ind_row*NDND +  ind_col. 
!         ! This is just the prelimiminery assessment in the process of assigning separate  
!         !     tiles for land, water and ice fractions within the EASE Grid cell
!         ! The program checks each nc x nr pixels whether there is a EASE grid cell underneath, and counts
!         ! number of water, land and ice pixels as seen on veg raster.
!         ! -----------------------------------------------------------------------------------------------
!         
!         
!         do i = 1 ,nc
!            
!            clon = -180. + float(i-1)*dx + dx/2.
!            
!            do j =nr ,1 ,-1
!               
!               clat = -90. + float(j-1)*dy + dy/2.
!               call EASE_convert(EASELabel, clat, clon, r_ease, s_ease)
!               
!               ind_col = nint(r_ease) + 1 
!               ind_row = nint(s_ease) + 1
!               
!               if((ind_row.ge.1).and.(veg(i,j).ne.OceanType).and.(ind_row.le.nr_ease)) then
!                  l=  ind_row*NDND +  ind_col
!                  
!                  if(veg(i,j)==LakeType) then 
!                     water_id(l) = 1
!                  else if(veg(i,j)==IceType) then
!                     ice_id  (l) = 1
!                  else
!                     land_id (l) = 1
!                  endif
!               endif
!            end do
!         end do

         
      endif    ! (GEOS5_10arcsec_mask)
      
      ! --------------------------------------------------------     
      !
      ! Read SRTM elevation data - to be consistent with AGCM
      
      allocate(raster(i_raster,j_raster))     ! 2.5-min raster grid    ( 8640-by-4320 )
      allocate(q0    (nc,      nr)      )     ! 30-arcsec raster grid  (43200-by-21600)
      
      gtopo30 = trim(MAKE_BCS_INPUT_DIR)//'/land/topo/v1/srtm30_withKMS_2.5x2.5min.data'
     
      open (10,file=trim(gtopo30),form='unformatted',status='old',convert='little_endian')
      read (10) raster
      close (10,status='keep') 
      
      ! remap SRTM elevation data from 2.5-min to 30-arcsec raster grid

      if(regrid) then
         call RegridRasterReal(raster,q0)
      else
         q0 = raster
      endif
      
      deallocate (raster)
      
      ! ---------------------------------------------------------
      !
      ! determine number of land, lake, and landice tiles 

      n_land    = sum(land_id)     ! number of land    tiles
      n_lake    = sum(water_id)    ! number of lake    tiles
      n_landice = sum(ice_id)      ! number of landice tiles

      n_landlakelandice = n_land + n_lake + n_landice

      print *,'# of Land              tiles: ', n_land
      print *,'# of Lake              tiles: ', n_lake
      print *,'# of Landice           tiles: ', n_landice
      print *,'# of Land+Lake+Landice tiles: ', n_landlakelandice

      l_index     = 0                   ! start index for land    EASE cells
      w_index     = n_land              ! start index for lake    EASE cells
      i_index     = w_index + n_lake    ! start index for landice EASE cells

      allocate(ease_grid_area(1:N_ease_grid_cells))

      allocate(tile_area     (1:n_landlakelandice))
      allocate(my_land       (1:n_landlakelandice))
      allocate(all_id        (1:n_landlakelandice))

      allocate(tile_ele      (1:n_land))      

      ! ===========================================================================
      !
      ! prepare for second loop through raster grid cells

      land_id        = 0
      water_id       = 0
      ice_id         = 0

      my_land        = 0
      all_id         = 0

      ease_grid_area = 0. 
      tile_ele       = 0.
      tile_area      = 0.

      ! While looping through the nc x nr grid, this section derives land, ice and water tiles.
      ! Each EASE grid cell is assigned with an ID = ind_row*ND +  ind_col 
      !         ind_col, ind_row are overlying EASE grid cell indices 
      ! Based on the above sums: 
      !         l_index Grid cells have land fractions (sum(land_id)) 
      !         w_index EASE Grid cells have inland water fractions (sum(water_id))
      !         i_index EASE Grid cells have ice fractions (sum(ice_id))
      ! hence, tile_index        1                     to l_index                     represent land tiles  
      !         tile_index       l_index +1            to l_index + w_index           represent water (lakes) tiles  
      !         tile_index       l_index + w_index +1  to l_index + w_index + i_index represent ice tiles
      ! global nc x nr array of tileid_index(nc,nr) contains corresponding tile_index values which 
      !        is derived in the below loop
 
      ND_raster = 10*10**(nint(log10(1.*NR)))

      !i2 = 1   -- NOT USED

      do i = 1 ,nc
         
         clon = -180. + float(i-1)*dx + dx/2.               ! center lon of raster grid cell (*,j)
         
         do j =nr ,1 ,-1

            lats = -90._8 + (j - 0.5_8)*dy                  ! center lat of raster grid cell (*,j) -- lats declared REAL ?!?!?!  same as clat ?!?!?!
            clat = -90. + float(j-1)*dy + dy/2.             ! center lat of raster grid cell (*,j)

            ! get 1-based ind_col and ind_row indices of EASE grid cell that contains raster grid cell (i,j)
            
            call EASE_convert(EASELabel, clat, clon, r_ease, s_ease)  
            
            ind_col = nint(r_ease) + 1 
            ind_row = nint(s_ease) + 1     ! can be negative or greater than nr_ease (lat near N/S pole)
            
            if( (ind_row.ge.1) .and. (ind_row.le.nr_ease) ) then
               
               ! raster grid cell (i,j) is located within the lat band covered by the (cylindrical) EASE grid (approx 85S-85N)
               
               !l=  ind_row*NDND +  ind_col    ! 1-dim index for EASE grid cells
               
               l = (ind_row-1)*NDND + ind_col   ! 1-dim index for all EASE grid cells
               
               pix_area = ( sin(d2r*(lats+0.5*dy)) - sin(d2r*(lats-0.5*dy)) )*(dx*d2r)   ! area of 30-arcsec raster grid cell
            
               if (veg(i,j).ne.OceanType) then
                  
                  ! set tile ID and compute tile area
                  
                  select case (veg(i,j))
                     
                  case (LakeType)   ! raster grid cell (i,j) is Lake 
                     
                     if(water_id(l)==0) then
                        ! raster grid cell is first Lake type seen for this EASE grid cell, tile ID has not yet been set
                        ! recall: tile IDs for Lake: [(n_land+1):(n_land+n_lake)], and w_index was initalized to n_land above 
                        w_index           = w_index + 1                 
                        water_id(l)       = w_index     ! needed so if condition will be false for next raster grid cell of same type
                        tileid_index(i,j) = w_index  
                     else
                        ! no action needed (tile ID does not depend on number of contributing Lake raster grid cells)
                     endif
                     
                  case (IceType)    ! raster grid cell (i,j) is Landice
                     
                     if(ice_id(l)==0) then
                        ! raster grid cell is first Landice type seen for this EASE grid cell, tile ID has not yet been set
                        ! recall: tile IDs for Landice: [(n_land+n_lake+1):n_landlakelandice], and i_index was initalized to n_land+n_lake above 
                        i_index           = i_index + 1
                        ice_id(l)         = i_index     ! needed so if condition will be false for next raster grid cell of same type
                        tileid_index(i,j) = i_index
                     else
                        ! no action needed (tile ID does not depend on number of contributing Landice raster grid cells)                     
                     endif
                     
                  case (LandType)   ! raster grid cell (i,j) is Land
                     
                     if(land_id(l)==0) then
                        ! raster grid cell is first Land type seen for this EASE grid cell, tile ID has not yet been set
                        ! recall: tile IDs for Land: [1:n_land], and l_index was initalized to 0 above 
                        l_index           = l_index + 1     
                        land_id(l)        = l_index     ! needed so if condition will be false for next raster grid cell of same type
                        tileid_index(i,j) = l_index
                     else
                        ! no action needed (tile ID does not depend on number of contributing Land raster grid cells)                     
                     endif
                                    
                     ! sum up area and (area-weighted) elevation (only over raster grid cells of type land!)
                     
                     tile_ele(      tileid_index(i,j)) = tile_ele(      tileid_index(i,j)) + q0(i,j) * pix_area  ! q0 = elevation
                     
                     ! tile_area_land should be obsolete because identical to tile_area(1:n_land)
                     !tile_area_land(tileid_index(i,j)) = tile_area_land(tileid_index(i,j)) +           pix_area
                     
                  case default
                     
                     print *,'ERROR: unknown tile type value in veg(i,j): ', veg(i,j), '   STOPPING.'
                     stop
                     
                  end select
                  
                  ! sum up area of raster grid cells contributing to each tile (land or lake or landice)
                  
                  tile_area(tileid_index(i,j)) = tile_area(tileid_index(i,j)) + pix_area     
                  
                  my_land(tileid_index(i,j)) = l                    ! store 1-dim index of all EASE grid cells
                  all_id( tileid_index(i,j)) = j*ND_raster + i      ! store 1-dim index of all 30-arcsec raster grid cells; BUG????

                  ! BUG??? all_id stores only the 1-dim index of the 30-arcsec raster grid cells that is the 
                  !        last one to contribute to the EASE tile specified by tileid_index(i,j)
                  !        This does not seem to be the desired dominant ID across the EASE grid cell
                  
               endif  ! raster grid cell has tile type other than ocean 
               
               ! compute total area of raster grid cells of any tile type that contribute to EASE grid cell
               
               ease_grid_area(l) = ease_grid_area(l) + pix_area

            endif    ! raster grid cell is located within lat band of EASE grid
               
         end do
      end do
      
      deallocate(land_id, q0)
      deallocate(water_id)
      deallocate(ice_id  )

      tile_ele = tile_ele/tile_area(1:n_land)   ! finalize tile elevation

      ! adjustment Global Mean Topography to 614.649 (615.662 GTOPO 30) m
      ! -----------------------------------------------------------------
      mean_land_ele=0.

      do j=1,l_index	 
         mean_land_ele = mean_land_ele + tile_ele(j)*tile_area(j)
      enddo

      mean_land_ele = mean_land_ele/sum(tile_area(1:l_index))

      if(mean_land_ele .ne. 614.649D0 ) then
         print *,'Global Mean Elevation (over land): ', mean_land_ele
         tile_ele = tile_ele*(614.649D0/mean_land_ele)				  	
         ! verify that adjust elevation is correct
         mean_land_ele=0.
         do j=1,l_index	 
            mean_land_ele = mean_land_ele + tile_ele(j)*tile_area(j)
         enddo
         print *,'Global Mean Elevation after scaling to SRTM : ',mean_land_ele/sum(tile_area(1:l_index))
      endif
      print *,'Total Land Area :', sum(tile_area(1:l_index))* MAPL_RADIUS * MAPL_RADIUS/1000./1000.  ! , &
           !sum(tile_area_land(1:l_index))* MAPL_RADIUS * MAPL_RADIUS/1000./1000.

      !-------------------------------------------
      !     
      ! write *.rst
      
      print *,'Writing ... ', trim(gfile)//'.rst'

      open (10, file ='rst/'//trim(gfile)//'.rst',form='unformatted',status='unknown',  &
           action='write')
      
      do j=1,nr
         write(10)(tileid_index(i,j),i=1,nc)
      end do
      
      close (10,status='keep')

      !-----------------------------------------------------------
      
      ! write catchment.def and *.til files
      
      print *,'Writing ... ', trim(gfile)//'.til  and  catchment.def'
      
      open (11,file='clsm/catchment.def',     &
           form='formatted',status='unknown')
      write(11,*)l_index

      open  (10, file ='til/'//trim(gfile)//'.til',form='formatted',status='unknown',action='write')
      write (10,*) i_index,SRTM_maxcat, nc, nr 
      write (10,*) 1
      write (10,*) EASELabel
      write (10,*) nc_ease
      write (10,*) nr_ease

      !     write (10,*)'NO-OCEAN'
      !     write (10,*) -9999
      !     write (10,*) -9999      
      
      dx_ease = 180./real(nc_ease)

      do l=1,i_index

         !ig    = my_land(l)-NDND*(my_land(l)/NDND)
         !jg    = my_land(l)/NDND
  
         ! get row and column indices from 1-dim index  (that is, invert l=(ind_row-1)*NDND+ind_col)
         
         jg = (my_land(l)-1)/NDND + 1          ! = ind_row (1-based)   [note integer division]
         
         ig = my_land(l) - NDND*(jg-1)         ! = ind_col (1-based)
       
         ! extract original PfafID from catid_index
         !
         ! BUG??? catid_index is in 30-arcsec raster space and contains the dominant PfafID
         !        (or LakeID, or LandiceID) within the 30-arcsec raster grid cell,
         !        where dominant is w.r.t. the 10-arcsec mask with PfafIDs
         !        It seems that here we just pick the PfafID associated with a single 
         !        30-arcsec raster grid cell within the EASE tile, which happens to be the 
         !        last of the 30-arcsec raster grid cells that contribute to the EASE tile
         ! 
         


 
         cindex= catid_index(all_id(l)-ND_raster*(all_id(l)/ND_raster),all_id(l)/ND_raster)

         if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
            pfaf = cindex
         else
            print *,'MaskFile = ', trim(MaskFile)
            print *,'ERROR: Selected mask file not supported for creating tiles on the EASE grid, stopping.'
            stop
            ! pfaf = catid(all_id(l)-ND_raster*(all_id(l)/ND_raster),all_id(l)/ND_raster)
         endif
         
         if ((l>l_index) .and. (l<=w_index)) typ = 19     ! Lake tile
         
         if ( l>w_index                    ) typ = 20     ! Landice tile
         
         if (l<=l_index) then 

            typ = 100                                     ! Land tile

            ! get min/max lat/lon of EASE grid cell
            ! BUG: This is *not* the desired min/max lat/lon of the land tile!!!
            
            call EASE_inverse( EASELabel, real(ig-1), real(jg-1), clat, clon ) 
            
            mnx = clon - dx_ease
            mxx = clon + dx_ease
            
            jgv = real(jg-1) + 0.5
            
            call EASE_inverse( EASELabel, real(ig-1), jgv, clat, clon ) 

            mny = clat
         
            jgv = real(jg-1) - 0.5
         
            call EASE_inverse( EASELabel, real(ig-1), jgv, clat, clon ) 

            mxy = clat 
            
            ! write tile properties into catchment.def file

            write (11,'(i10,i8,5(2x,f9.4), i4)') l, pfaf, mnx, mxx, mny, mxy, tile_ele(l)

         endif

         ! get area fraction of tile within EASE grid cell
         ! NOTE: the area of the EASE grid cell here has to be the sum of the areas of the
         !       contributing raster grid cells, which is *not* the same for all EASE grid cells
         
         call EASE_inverse( EASELabel, real(ig-1), real(jg-1), clat, clon )
         
         fr_gcm = tile_area(l) / ease_grid_area((jg-1)*NDND+ig)
         
         ! write tile properties into *.til file
         
         if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
            
            ! Note: If-condition uses MaskFile as a proxy for format & specs of *.til file (columns written)
            
            ! pfaf       = running *index* of Pfafstetter catchment (1:SRTM_maxcat)
            ! SRTM_catid = 12-digit Pfafstetter code (encoding network/routing info)
            
            write(10,'(i10,i9,2f10.4,2i6,f19.12,i10,i15,e13.4)')                         &
                 typ, pfaf, clon, clat, ig-1, jg-1, fr_gcm, pfaf, SRTM_catid(cindex) 
         else

            ! write statement below was used with MaskFile prior to availability of SRTM-based Pfafstetter catchments
            ! obsolete with EASE grid
            ! - reichle, 15 Jun 2023

            print *,'MaskFile = ', trim(MaskFile)
            print *,'ERROR: Selected mask file not supported for creating tiles on the EASE grid, stopping.'
            stop
!
!            write(10,'(i10,i9,2f10.4,2i5,f19.12,i10,e13.4,i8)') &
!                 typ,pfaf,clon,clat,ig-1,jg-1,fr_gcm ,cindex 
         endif
         
      end do
      
      close(10,status='keep')      
      close(11,status='keep')          
      
      deallocate( tileid_index, catid_index,veg )
      deallocate( tile_area, ease_grid_area, tile_ele, my_land, all_id )
      
      ! Commented out "empty" if-block. -rreichle, 15 Jun 2023
!
!      if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then         
!
!         print *,'Creating SMAP-Catch_TransferData.nc files.'
!
!         !---------------------------------------------------
!
!         deallocate (SRTM_CatchArea, SRTM_catid, SRTM_catid_r8)
!         
!      endif
      
      ! create Grid2Catch transfer file
      ! -------------------------------

      ! CALL CREATE_ROUT_PARA_FILE (NC, NR, trim(gfile), MGRID=MGRID)  
      
      ! now run mkCatchParam
      ! --------------------

      ! WY Note: now mkCatchParam is run in the make_bcs script, not here
      !          and nthread will be reset to run mkCatchParam

      ! tmpstring1 = '-e EASE -g '//trim(gfile)//' -v '//trim(LBCSV)
      ! write(tmpstring2,'(2(a2,x,i5,x))')'-x',nc,'-y',nr
      ! tmpstring = 'bin/mkCatchParam.x '//trim(tmpstring2)//' '//trim(tmpstring1)
      ! print *,trim(tmpstring)
      
      ! call execute_command_line (tmpstring)


!!! commented out. It may be used in the future for irrigation tiles
!!!    contains
!!!
!!!      ! -------------------------------------------------------------------------------
!!!
!!!      SUBROUTINE mkEASEv2Raster
!!!
!!!        implicit none
!!!
!!!        integer       :: i, j, i_ease, j_ease
!!!        real*8,   allocatable :: xs(:,:), ys(:,:)
!!!        real          :: x,y, xout, yout
!!!        
!!!        allocate (xs ( nc_ease+1, nr_ease+1))
!!!        allocate (ys ( nc_ease+1, nr_ease+1))
!!!        
!!!        do  j = 1, nr_ease+1
!!!           do i = 1, nc_ease+1
!!!              x = real(i-1)        -0.5
!!!              y = real(nr_ease - j)+0.5
!!!              call EASE_inverse(MGRID, x, y, yout, xout)
!!!              ys (i,j) = dble(yout)
!!!              xs (i,j) = dble(xout)
!!!           end do
!!!        end do
!!!
!!!        call  LRRasterize(EASElabel,xs,ys,nc=nc,nr=nr,xmn = xs(1,1), xmx= xs(nc_ease+1, nr_ease+1), &
!!!                       ymn=ys(1,1), ymx = ys(nc_ease+1, nr_ease+1), Here=.false., Verb=.false.)       
!!!
!!!        stop
!!!      end SUBROUTINE mkEASEv2Raster
!!!
!!!      ! ------------------------------------------------------------
!!!      
!!!      SUBROUTINE write_tilfile 
!!!
!!!        implicit none
!!!
!!!        character*200 :: infile
!!!        integer      :: NT, NF, NC, NR, NPF, NG, IDUM, i, N, icol, rcol
!!!        character*20 :: cdum
!!!        integer, dimension (:,:), allocatable :: iRtable
!!!        real,    dimension (:,:), allocatable :: rRtable
!!!
!!!        infile = 'til/'//trim(EASElabel)//'_'//trim(EASElabel)//'-Pfafstetter.'
!!!        
!!!        open (10,file =  trim(infile)//'ind', form = 'formatted', action = 'read', status = 'old')
!!!        open (11,file =  trim(infile)//'TIL', form = 'formatted', action = 'write', status = 'unknown')
!!!
!!!        read (10, *) NT, NF, NC, NR
!!!        write (11,'(4I10)')NT, NF, NC, NR
!!!        read (10, *) NG
!!!        write(11, *) NG
!!!        
!!!        do n = 1, NG
!!!           read (10, '(a)') cdum
!!!           write(11, '(a)') trim (cdum)
!!!           read (10, *) IDUM
!!!           write(11, '(I10)') IDUM
!!!           read (10, *) IDUM
!!!           write(11, '(I10)') IDUM
!!!        end do
!!!        
!!!        icol = 7
!!!        rcol = 5
!!!        allocate (iRtable (1, 1:icol))
!!!        allocate (rRtable (1, 1:rcol))
!!!        
!!!        do n = 1,  nt
!!!           read(10,'(I10,3E20.12,9(2I10,E20.12,I10))') iRtable (1,1),rRtable(1,1), &
!!!                rRtable(1,2),rRtable(1,3),iRtable (1,2),iRtable (1,3),rRtable(1,4),iRtable (1,4),&
!!!                iRtable (1,5),iRtable (1,6),rRtable(1,5),iRtable (1,7)
!!!           write(11,'(I10,3E20.12,9(2I10,E20.12,I10))') iRtable (1,1),rRtable(1,1), &
!!!                rRtable(1,2),rRtable(1,3),iRtable (1,2)-1,nr_ease - iRtable (1,3),rRtable(1,4),iRtable (1,4),&
!!!                iRtable (1,5),iRtable (1,6),rRtable(1,5),iRtable (1,7)
!!!        end do
!!!     
!!!        close (10, status = 'keep')
!!!        close (11, status = 'keep')
!!!    
!!!   END SUBROUTINE write_tilfile
 END PROGRAM

