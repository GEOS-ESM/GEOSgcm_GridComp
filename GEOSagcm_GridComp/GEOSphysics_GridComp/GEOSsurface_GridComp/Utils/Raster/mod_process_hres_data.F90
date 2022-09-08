#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

!
! A Collection subroutine that helps process MODIS Albedo, GEOLAND2 LAI, and
! NGDC-HWSD-STATSGO merged soil data on their native grids 3-23-2012
!   Contact: Sarith Mahanama  sarith.p.mahanama@nasa.gov
!   Email  : sarith.p.mahanama@nasa.gov
!
! CHANGE LOG:
!
! jkolassa, reichle, May 2022: 
! The bcs file "CLM4.5_veg_typs_fracs" was not used in CatchmentCNCLM45 and is no longer
! produced by make_bcs.
! Separate mappings from ESA GlobCover to CatchmentCNCLM40 and CatchmentCNCLM45 PFTs
! were initially implemented because the underlying CLM4.0 and CLM4.5 models have different
! plant functional types and distributions. Ultimately, the decision was made to use the same 
! (CLM4.0-based) PFT distribution for both CatchmentCNCLM40 and CatchmentCNCLM45, and the 
! obsolete mapping of ESA GlobCover data to CLM4.5 PFTs (subroutine ESA2CLM_45) was removed.


MODULE process_hres_data
use rmTinyCatchParaMod
use MAPL_SortMod
use date_time_util  
use leap_year
use MAPL_ConstantsMod
use lsm_routines, ONLY: sibalb

#if defined USE_EXTERNAL_FINDLOC
use findloc_mod, only: findloc
#endif

implicit none

include 'netcdf.inc'	

private

public :: soil_para_hwsd,hres_lai,hres_gswp2, merge_lai_data, grid2tile_modis6
public :: modis_alb_on_tiles_high,modis_scale_para_high,hres_lai_no_gswp
public :: histogram, create_mapping, esa2mosaic , esa2clm
public :: grid2tile_ndep_t2m_alb, CREATE_ROUT_PARA_FILE, map_country_codes, get_country_codes
public :: CLM45_fixed_parameters, CLM45_clim_parameters, gimms_clim_ndvi, grid2tile_glass,  open_landparam_nc4_files

! Below structure is used to regrid high resolution data to high resolution tile raster

integer, parameter   :: N_tiles_per_cell = 9
integer  , parameter :: nc_esa = 129600, nr_esa = 64800
real, parameter      :: pi= MAPL_PI,RADIUS=MAPL_RADIUS
integer, parameter   :: N_GADM = 256 + 1, N_STATES = 50

real, parameter      :: SOILDEPTH_MIN_HWSD = 1334.   ! minimum soil depth for HWSD soil parameters

type :: do_regrid
   integer                               :: NT
   integer, dimension (N_tiles_per_cell) :: TID
   integer, dimension (N_tiles_per_cell) :: count
end type do_regrid
type, public :: regrid_map
   integer :: nc_data = 1
   integer :: nr_data = 1
   integer, allocatable, dimension (:,:)   :: ij_index
   type(do_regrid), pointer, dimension (:) :: map
end type regrid_map

contains
!
! ---------------------------------------------------------------------
!

  SUBROUTINE ESA2CLM (nc, nr, gfile)

    implicit none

    integer  , intent (in) :: nc, nr
    character (*)          :: gfile
    
    integer  , parameter   :: N_lon_clm = 1152, N_lat_clm = 768, lsmpft = 17
    integer*2, allocatable, target, dimension (:,:) :: esa_veg
    integer*2, pointer    , dimension (:,:) :: subset
    integer  , allocatable, dimension (:)   :: tile_id, i_esa2clm, j_esa2clm
    integer :: i,j, k,n, status, ncid, varid, maxcat, dx,dy, esa_type, tid, cid, ii, jj   
    real    :: dx_clm, dy_clm, x_min_clm (N_lon_clm), y_min_clm (N_lat_clm), clm_fracs(lsmpft)
    real    :: minlon,maxlon,minlat,maxlat,tile_lat, scale, ftot
    integer :: cpt1, cpt2, cst1, cst2  ! CLM-carbon types
    real    :: cpf1, cpf2, csf1, csf2  ! CLM-carbon fractions
    DOUBLE PRECISION,  allocatable, dimension (:) :: lon_esa, lat_esa
    DOUBLE PRECISION       :: EDGEN, EDGEE, EDGES, EDGEW 
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: PCT_PFT_DBL
    REAL, ALLOCATABLE, DIMENSION (:,:,:) :: PCTPFT 
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    integer :: NBINS, NPLUS
    integer, allocatable, dimension (:,:) :: clm_veg
    integer :: esa_clm_veg (2)
    real    :: esa_clm_frac(2)
    logical :: file_exists
    REAL, ALLOCATABLE, DIMENSION (:,:) :: NITYP,NFVEG

    ! Reading CLM pft data file
    !--------------------------

    ALLOCATE (PCTPFT      (1:N_lon_clm, 1:N_lat_clm, 1:lsmpft))
    ALLOCATE (PCT_PFT_DBL (1:N_lon_clm, 1:N_lat_clm, 1:lsmpft))
    status  = NF_OPEN ('data/CATCH/surfdata_0.23x0.31_simyr2000_c100406.nc', NF_NOWRITE, ncid)   
    status  = NF_GET_VARA_DOUBLE (ncid,1,(/1/),(/1/),EDGEN) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_DOUBLE (ncid,2,(/1/),(/1/),EDGEE) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_DOUBLE (ncid,3,(/1/),(/1/),EDGES) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_DOUBLE (ncid,4,(/1/),(/1/),EDGEW) ; VERIFY_(STATUS)
    status  = NF_INQ_VARID (ncid,'PCT_PFT',VarID) ; VERIFY_(STATUS)

    do k = 1, lsmpft
       status  = NF_GET_VARA_DOUBLE (ncid,VarID,(/1,1,k/),(/N_lon_clm, N_lat_clm, 1/),PCT_PFT_DBL(:,:,k)) ; VERIFY_(STATUS)
    end do

    status = NF_CLOSE(ncid)

    ! change type 6 to 10 for Australia only gkw: to remove CLM artificial tree line, and stay true to ESA
    ! ----------------------------------------------------------------------------------------------------

    PCT_PFT_DBL(360:494,215:341,11) = PCT_PFT_DBL(360:494,215:341,11) + PCT_PFT_DBL(360:494,215:341, 7)
    PCT_PFT_DBL(360:494,215:341, 7) = 0.

    ! CLM description (17)                                    CatchmentCNCLM description (19)        
    ! --------------------                                    ------------------------------ 

    ! 'BARE'   1  	bare                                     (does not have bare soil)
    ! 'NLEt'   2 	needleleaf evergreen temperate tree    1
    ! 'NLEB'   3 	needleleaf evergreen boreal tree       2
    ! 'NLDB'   4  	needleleaf deciduous boreal tree       3
    ! 'BLET'   5 	broadleaf evergreen tropical tree      4
    ! 'BLEt'   6 	broadleaf evergreen temperate tree     5
    ! 'BLDT'   7 	broadleaf deciduous tropical tree      6
    ! 'BLDt'   8 	broadleaf deciduous temperate tree     7
    ! 'BLDB'   9 	broadleaf deciduous boreal tree        8
    ! 'BLEtS' 10 	broadleaf evergreen temperate shrub    9
    ! 'BLDtS' 11 	broadleaf deciduous temperate shrub   10  broadleaf deciduous temperate shrub [moisture +  deciduous]
    ! 'BLDtSm'  	broadleaf deciduous temperate shrub   11  broadleaf deciduous temperate shrub [moisture stress only]
    ! 'BLDBS' 12 	broadleaf deciduous boreal shrub      12
    ! 'AC3G'  13 	arctic c3 grass                       13
    ! 'CC3G'  14 	cool c3 grass                         14  cool c3 grass [moisture +  deciduous]
    ! 'CC3Gm'           cool c3 grass                         15  cool c3 grass [moisture stress only]
    ! 'WC4G'  15 	warm c4 grass                         16
    ! 'WC4Gm'   	warm c4 grass                         17
    ! 'CROP'  16 	crop                                  18  crop [moisture +  deciduous]
    ! 'CROPm'   	crop                                  19  crop [moisture stress only]
    !         17        water

    dx_clm = 360./N_lon_clm
    dy_clm = 180./N_lat_clm

    do i = 1, N_lon_clm 
       x_min_clm (i) = (i-1)*dx_clm + EDGEW - 180.
    end do

    do i = 1,  N_lat_clm
       y_min_clm (i) = (i-1)*dy_clm  + EDGES
    end do

    PCTPFT (1:N_lon_clm/2             ,:,:) =  REAL (PCT_PFT_DBL(N_lon_clm/2 + 1: N_lon_clm,:,:))
    PCTPFT (N_lon_clm/2 + 1: N_lon_clm,:,:) =  REAL (PCT_PFT_DBL(1:N_lon_clm/2             ,:,:))

    DEALLOCATE (PCT_PFT_DBL)

    ! Find primary and secondary types in the CLM data file
    ! -----------------------------------------------------

    ! allocate (clm_veg (1:N_lon_clm,1:N_lat_clm,1:2))
    !
    ! do j = 1, N_lat_clm 
    !    do i = 1, N_lon_clm  
    !       if(maxval(PCT_PFT(i,j,:)) > 0.) then
    !          clm_fracs = PCT_PFT(i,j,:)
    !          if (maxval (clm_fracs) == 100.) then 
    !             clm_veg(i,j,:) = maxloc (clm_fracs)
    !          else 
    !             clm_veg(i,j,0)             = maxloc (clm_fracs)		
    !             clm_fracs (clm_veg(i,j,0)) = 0.
    !             clm_veg(i,j,1)             = maxloc (clm_fracs)	 
    !          endif
    !       else 
    !          clm_veg(i,j,:) = 17
    !       endif
    !    end do
    ! end do
    
    ! Reading ESA vegetation types
    !-----------------------------

    allocate (esa_veg (1:nc_esa, 1: nr_esa))
    allocate (lon_esa (1:nc_esa))
    allocate (lat_esa (1:nr_esa))

    status    = NF_OPEN ('data/CATCH/ESA_GlobalCover.nc', NF_NOWRITE, ncid)   

    if(status /=0) then
       PRINT *, NF_STRERROR(STATUS)
       print *, 'Problem with NF_OPEN','ESA_GlobalCover.nc'
       stop
    endif

       status  = NF_GET_VARA_DOUBLE (ncid,1,(/1/),(/nr_esa/),lat_esa)
       status  = NF_GET_VARA_DOUBLE (ncid,2,(/1/),(/nc_esa/),lon_esa)

    do j = 1,nr_esa
       status  = NF_GET_VARA_INT2 (ncid,3,(/1,j/),(/nc_esa,1/),esa_veg(:,j))
       if(status /=0) then
          PRINT *, NF_STRERROR(STATUS)
          print *, 'Problem with NF_GET ESA_GlobalCover.nc : ', STATUS
          stop
       endif
    end do

    status = NF_CLOSE(ncid)

    ! Find I,J of overlying CLM grid cells for each ESA pixel 
    !--------------------------------------------------------
    allocate (i_esa2clm (1:nc_esa))
    allocate (j_esa2clm (1:nr_esa))

    do i = 1, N_lon_clm 
       where ((real(lon_esa) >=  x_min_clm(i)).and.(real(lon_esa) <  (x_min_clm(i) + dx_clm))) i_esa2clm= i
    end do

    i_esa2clm(129545:nc_esa) = 1

    do j = 1, N_lat_clm 
       where ((real(lat_esa) >=  y_min_clm(j)).and.(real(lat_esa) <  (y_min_clm(j) + dy_clm))) j_esa2clm= j
    end do

    !
    ! Reading number of tiles
    ! -----------------------

    open (10, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
         action =  'read')
    
    read (10, *) maxcat

    close (10, status = 'keep')
    
    !
    ! Loop through tile_id raster
    ! ___________________________


    allocate (tile_id (1:nc             ))   
    allocate (clm_veg (1:maxcat,1:lsmpft))
    clm_veg = 0.

    dx = nc_esa / nc
    dy = nr_esa / nr

    open (10,file=trim(gfile)//'.rst',status='old',action='read',  &
          form='unformatted',convert='little_endian')

    do j=1,nr
       
       ! read a row
       
       read(10)tile_id(:)
       
       do i = 1,nc

          ii = i_esa2clm ((i-1)*dx + dx/2)
          jj = j_esa2clm ((j-1)*dy + dy/2)

          if((tile_id (i) >= 1).and.(tile_id(i)  <= maxcat)) then

             if (associated (subset)) NULLIFY (subset)
             subset => esa_veg((i-1)*dx +1 :i*dx, (j-1)*dy +1:j*dy)
             NPLUS = count(subset >= 1 .and. subset <= 230)

             if(NPLUS > 0)  then
                allocate (loc_int (1:NPLUS))
                allocate (unq_mask(1:NPLUS))
                loc_int = pack(subset,mask = (subset >= 1 .and. subset <= 230))
                call MAPL_Sort (loc_int)
                unq_mask = .true.
                do n = 2,NPLUS 
                   unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))
                end do
                NBINS = count(unq_mask)
                
                allocate(loc_val (1:NBINS))
                allocate(density (1:NBINS))
                loc_val = 1.*pack(loc_int,mask =unq_mask)
                call histogram (size(subset,1)*size(subset,2), NBINS, density, loc_val, real(subset))   
                
                do k = 1, nbins
                   
                   if (density (k) > 0) then
                      
                      esa_type = int (loc_val(k))
                      
                      ! if (esa_type ==  10)  clm_veg (tile_id(i), 17) = 1.* density(k)   ! lakes inland water
                      
                      if (esa_type ==  11)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type  11: Post-flooding or irrigated croplands 
                      if (esa_type ==  14)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type  14: Rainfed croplands 
                      if (esa_type ==  20)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type  20: Mosaic Cropland (50-70%) / Vegetation (grassland, shrubland, forest) (20-50%) 
                      if (esa_type == 190)  clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   ! ESA type 190:	Artificial surfaces and associated areas (urban areas >50%) 
                      
                      ! if (esa_type == 200)  clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 1.* density(k)   ! ESA type 200:	Bare areas
                      ! if (esa_type == 210)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ocean
                      ! if (esa_type == 220)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ice  
                      ! gkw: bare soil excluded! only considering vegetated land                   
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  30) then 
                         ! ESA type  30: Mosaic Vegetation (grassland, shrubland, forest) (50-70%) / Cropland (20-50%) 
                         clm_veg (tile_id(i),  16) = clm_veg (tile_id(i),  16) + 0.5* density(k)
                         if(sum(PCTPFT(ii,jj,2:15)) > 0.) then
                            do n = 2,  15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.5* density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:15))
                            enddo
                         else 
                            clm_veg (tile_id(i),  16) = clm_veg (tile_id(i),  16) + 1.0* density(k)
                         endif
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  40) then
                         ! ESA type  40:	Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,5:6)) > 0.) then
                            do n = 5, 6 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,5:6))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  5) = clm_veg (tile_id(i),  5) + 1.0* density(k)
                            else
                               clm_veg (tile_id(i),  6) = clm_veg (tile_id(i),  6) + 1.0* density(k)
                            endif
                         endif
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if ((esa_type ==  50) .or. (esa_type ==  60)) then    
                         ! ESA type  50:	Closed (>40%) broadleaved deciduous forest (>5m) 
                         ! ESA type  60:	Open (15-40%) broadleaved deciduous forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,7:9)) > 0.) then
                            do n = 7, 9 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:9))
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 1.0* density(k)
                            else
                               if(abs(y_min_clm(jj) + 0.5*dy_clm) <  60.) clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 1.0* density(k)
                               if(abs(y_min_clm(jj) + 0.5*dy_clm) >= 60.) clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 1.0* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  70) then    
                         ! ESA type  70:	Closed (>40%) needleleaved evergreen forest (>5m)
                         
                         if(sum(PCTPFT(ii,jj,2:3)) > 0.) then
                            do n = 2, 3 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:3))	
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  90) then 
                         !ESA type  90:	Open (15-40%) needleleaved deciduous or evergreen forest (>5m) 
                         
                         if(sum(PCTPFT(ii,jj,2:4)) > 0.) then
                            do n = 2, 4 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n)   + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,2:4))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 100) then 
                         !  ESA type 100:	Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m) 
                         
                         if((sum(PCTPFT(ii,jj,2:4)) + sum(PCTPFT(ii,jj,7:9))) > 0.) then
                            do n = 2, 9 
                               if((n /= 5) .and. (n /= 6)) clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/(sum(PCTPFT(ii,jj,2:4)) + sum(PCTPFT(ii,jj,7:9)))	
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.5* density(k)
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 0.5* density(k)			
                            elseif (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.5* density(k)
                               clm_veg (tile_id(i),  2) = clm_veg (tile_id(i),  2) + 0.5* density(k)
                            else 
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.5* density(k)
                               clm_veg (tile_id(i),  3) = clm_veg (tile_id(i),  3) + 0.5* density(k)			
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 110) then   
                         ! ESA type 110:	Mosaic Forest/Shrubland (50-70%) / Grassland (20-50%) 
                         
                         if(sum(PCTPFT(ii,jj,7:12))  > 0.) then
                            do n = 7, 12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) +  0.6*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:12)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.3* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.3* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.3* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.3* density(k)
                            else
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.3* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.3* density(k)			
                            end if
                         end if
                         
                         if(sum(PCTPFT(ii,jj,13:15))  > 0.) then
                            do n =13, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.4*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 0.4* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.4* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.4* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------

                      if (esa_type == 120) then   
                         ! ESA type 120:	Mosaic Grassland (50-70%) / Forest/Shrubland (20-50%) 
                         
                         if(sum(PCTPFT(ii,jj,7:12)) > 0.) then
                            do n = 7, 12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) +  0.4*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,7:12)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  7) = clm_veg (tile_id(i),  7) + 0.2* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.2* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 0.2* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.2* density(k)
                            else 
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 0.2* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.2* density(k)			
                            end if
                         end if
                         
                         if(sum(PCTPFT(ii,jj,13:15)) > 0.) then
                            do n =13, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.6*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15)) 
                            enddo
                         else
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 0.6* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.6* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.6* density(k)
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 130) then
                         ! 	Closed to open (>15%) shrubland (<5m) 
                         
                         if(sum(PCTPFT(ii,jj,10:12)) > 0.) then
                            do n = 10,12 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:12))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  11) = clm_veg (tile_id(i),  11) + 1.0* density(k)
                            else 	
                               clm_veg (tile_id(i),  12) = clm_veg (tile_id(i),  12) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 140) then
                         ! ESA type 140:	Closed to open (>15%) grassland 
                         
                         if(sum(PCTPFT(ii,jj,13:15)) > 0.) then
                            do n = 13,15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,13:15))
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i),  15) = clm_veg (tile_id(i),  15) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then	
                               clm_veg (tile_id(i),  14) = clm_veg (tile_id(i),  14) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i),  13) = clm_veg (tile_id(i),  13) + 1.0* density(k)	
                            end if
                         end if
                      end if
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 150) then
                         ! ESA type 150:	Sparse (<15%) vegetation (woody vegetation, shrubs, grassland) 
                         
                         if(sum(PCTPFT(ii,jj,10:15)) > 0.) then
                            do n = 10, 15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.0*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:15)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.5* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.5* density(k)			
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 0.5* density(k)
                               clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 0.5* density(k)
                            else 
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 0.5* density(k)
                               clm_veg (tile_id(i), 12) = clm_veg (tile_id(i), 12) + 0.5* density(k)			
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if((esa_type == 160) .or. (esa_type == 170)) then  
                         ! ESA type 160:	Closed (>40%) broadleaved forest regularly flooded - Fresh water      ! ESA type 170:	Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded
                         
                         if(sum(PCTPFT(ii,jj,5:9)) > 0.) then
                            do n = 5,9 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,5:9)) 
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 23.5) then
                               clm_veg (tile_id(i),  5) = clm_veg (tile_id(i),  5) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 60.) then
                               clm_veg (tile_id(i),  8) = clm_veg (tile_id(i),  8) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i),  9) = clm_veg (tile_id(i),  9) + 1.0* density(k)	
                            end if
                         end if
                      endif
                      
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type == 180) then
                         ! ESA type 180:	Closed to open (>15%) vegetation (grassland, shrubland, woody vegetation) on regularly flooded or waterlogged soil - Fresh, brackish or saline water 
                         
                         if(sum(PCTPFT(ii,jj,10:15)) > 0.) then
                            do n = 10,15 
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 1.*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,10:15))	
                            enddo
                         else 
                            if(abs(y_min_clm(jj) + 0.5*dy_clm) < 30.) then
                               clm_veg (tile_id(i), 15) = clm_veg (tile_id(i), 15) + 1.0* density(k)
                            else if (abs(y_min_clm(jj) + 0.5*dy_clm) < 55.) then
                               clm_veg (tile_id(i), 14) = clm_veg (tile_id(i), 14) + 1.0* density(k)	
                            else 	
                               clm_veg (tile_id(i), 13) = clm_veg (tile_id(i), 13) + 1.0* density(k)	
                            end if
                         end if
                      endif
                   endif
                enddo
                deallocate (loc_int,unq_mask,loc_val,density)
             endif
          end if
       enddo
    end do
    
    
    deallocate (tile_id, PCTPFT,esa_veg,lon_esa,lat_esa,i_esa2clm,j_esa2clm)  
    close (10,status='keep')    

    !
    ! Now create CLM-carbon_veg_fracs file
    ! ------------------------------------

    open (10,file='clsm/CLM_veg_typs_fracs',  &
         form='formatted',status='unknown')
    open (11, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
         action =  'read')
    
    read (11, *) maxcat   

    inquire(file='clsm/catchcn_params.nc4', exist=file_exists)
    if(file_exists) then
       status = NF_OPEN ('clsm/catchcn_params.nc4', NF_WRITE, ncid) ; VERIFY_(STATUS)
       allocate (NITYP (1:MAXCAT, 1:4))
       allocate (NFVEG (1:MAXCAT, 1:4))    
    endif

    do k = 1, maxcat

       read (11,'(i10,i8,5(2x,f9.4))') tid,cid,minlon,maxlon,minlat,maxlat
       tile_lat = (minlat + maxlat)/2.
       scale = (ABS (tile_lat) - 32.)/10.
       scale = min (max(scale,0.),1.)

       esa_clm_veg = 0
       esa_clm_frac= 0.

       clm_fracs = clm_veg (k,:)
             
       if (sum (clm_fracs) == 0.) then ! gkw: no vegetated land found; set to BLDtS
          esa_clm_veg (1) = 11              ! broadleaf deciduous shrub 
          esa_clm_frac(1) = 100.
       else
          esa_clm_veg (1) = maxloc(clm_fracs,1)
          esa_clm_frac(1) = maxval(clm_fracs) 
       endif

       clm_fracs (esa_clm_veg (1)) = 0.

       if (sum (clm_fracs) == 0.) then ! gkw: no vegetated secondary type found, set to primary with zero fraction
          esa_clm_veg (2) = esa_clm_veg (1)
          esa_clm_frac(1) = 100.
          esa_clm_frac(2) = 0.
       else
          esa_clm_veg (2) = maxloc(clm_fracs,1)
          esa_clm_frac(1) = 100.*clm_veg (k,esa_clm_veg (1))/(clm_veg (k,esa_clm_veg (1)) + clm_veg (k,esa_clm_veg (2)))
          esa_clm_frac(2) = 100. - esa_clm_frac(1)
       end if

! Now splitting CLM types for CNCLM  model
! --------------------------------------------
 
! CLM types 2- 10,12,13 are not being split.
! .............................................
     
       if ((esa_clm_veg (1) >= 2).and.(esa_clm_veg (1) <= 10)) then
          CPT1 = esa_clm_veg (1) - 1
          CPT2 = esa_clm_veg (1) - 1
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 2).and.(esa_clm_veg (2) <= 10)) then
          CST1 = esa_clm_veg (2) - 1
          CST2 = esa_clm_veg (2) - 1
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! .............................................

       if ((esa_clm_veg (1) >= 12).and.(esa_clm_veg (1) <= 13)) then
          CPT1 = esa_clm_veg (1)
          CPT2 = esa_clm_veg (1)
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 12).and.(esa_clm_veg (2) <= 13)) then
          CST1 = esa_clm_veg (2)
          CST2 = esa_clm_veg (2)
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! Now splitting
! .............

       if (esa_clm_veg (1) == 11) then
          CPT1 = 10
          CPT2 = 11
          CPF1 = esa_clm_frac(1) * scale
          CPF2 = esa_clm_frac(1) * (1. - scale)
       endif

       if (esa_clm_veg (2) == 11) then
          CST1 = 10
          CST2 = 11
          CSF1 = esa_clm_frac(2) * scale       
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

! .............

       if (esa_clm_veg (1) == 14) then
          CPT1 = 14
          CPT2 = 15
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 14) then
          CST1 = 14
          CST2 = 15
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

! .............

       if (esa_clm_veg (1) == 15) then
          CPT1 = 16
          CPT2 = 17
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 15) then
          CST1 = 16
          CST2 = 17
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale)
       endif
! .............

       if (esa_clm_veg (1) == 16) then
          CPT1 = 18
          CPT2 = 19
          CPF1 = esa_clm_frac(1) * scale        
          CPF2 = esa_clm_frac(1) * (1. - scale) 
       endif

       if (esa_clm_veg (2) == 16) then
          CST1 = 18
          CST2 = 19
          CSF1 = esa_clm_frac(2) * scale        
          CSF2 = esa_clm_frac(2) * (1. - scale) 
       endif

       ! fractions must sum to 1
       ! -----------------------
       ftot = cpf1 + cpf2 + csf1 + csf2

       if(ftot /= 100.) then
          cpf1 = 100. * cpf1 / ftot  
          cpf2 = 100. * cpf2 / ftot 
          csf1 = 100. * csf1 / ftot 
          csf2 = 100. * csf2 / ftot 
       endif
    
       write (10,'(2I10,4I3,4f7.2,2I3,2f7.2)')     &
            tid,cid,cpt1, cpt2, cst1, cst2, cpf1, cpf2, csf1, csf2, &
            esa_clm_veg (1), esa_clm_veg (2), esa_clm_frac(1), esa_clm_frac(2)

       if (allocated (NITYP)) NITYP (k, :) = (/REAL(cpt1), REAL(cpt2), REAL(cst1), REAL(cst2)/)
       if (allocated (NFVEG)) NFVEG (k, :) = (/cpf1, cpf2, csf1, csf2/)

    end do

    if(file_exists) then

       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ITY'    ) ,(/1,1/),(/maxcat,1/), NITYP (:, 1)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ITY'    ) ,(/1,2/),(/maxcat,1/), NITYP (:, 2)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ITY'    ) ,(/1,3/),(/maxcat,1/), NITYP (:, 3)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ITY'    ) ,(/1,4/),(/maxcat,1/), NITYP (:, 4)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'FVG'    ) ,(/1,1/),(/maxcat,1/), NFVEG (:, 1)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'FVG'    ) ,(/1,2/),(/maxcat,1/), NFVEG (:, 2)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'FVG'    ) ,(/1,3/),(/maxcat,1/), NFVEG (:, 3)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'FVG'    ) ,(/1,4/),(/maxcat,1/), NFVEG (:, 4)) ; VERIFY_(STATUS)
       DEALLOCATE (NITYP, NFVEG)
       STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
   
    endif
    
    close (10, status = 'keep')
    close (11, status = 'keep')    

  END SUBROUTINE ESA2CLM
!
! ---------------------------------------------------------------------
!
  SUBROUTINE ESA2MOSAIC (nc, nr, gfile)
    
    implicit none

    integer  , intent (in) :: nc, nr
    character (*)          :: gfile
    integer  , parameter   :: nc_esa = 129600, nr_esa = 64800
    integer*2, allocatable, target, dimension (:,:) :: esa_veg
    integer*2, pointer    , dimension (:,:) :: subset
    integer  , allocatable, dimension (:)   :: tile_id, ityp
    integer :: i,j, k, status, ncid, maxcat, dx,dy, esa_type, tid, cid
    integer :: mos1, mos2
    real    :: mfrac, sfrac, tfrac, tem (6)
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    real   , allocatable       :: veg (:,:)
    integer :: NBINS, NPLUS
    real, pointer, dimension (:)  :: z2, z0
    real, dimension (6) :: VGZ2 = (/35.0, 20.0, 17.0, 0.6, 0.5, 0.6/) ! Dorman and Sellers (1989)
    logical :: file_exists

    ! Reading ESA vegetation types
    !-----------------------------

    allocate (esa_veg (1:nc_esa, 1: nr_esa))

    status    = NF_OPEN ('data/CATCH/ESA_GlobalCover.nc', NF_NOWRITE, ncid)   

    if(status /=0) then
       PRINT *, NF_STRERROR(STATUS)
       print *, 'Problem with NF_OPEN','ESA_GlobalCover.nc'
       stop
    endif

    do j = 1,nr_esa
       status  = NF_GET_VARA_INT2 (ncid,3,(/1,j/),(/nc_esa,1/),esa_veg(:,j))
       if(status /=0) then
          PRINT *, NF_STRERROR(STATUS)
          print *, 'Problem with NF_GET ESA_GlobalCover.nc : ', STATUS
          stop
       endif
    end do
    status = NF_CLOSE(ncid)

!
! Reading number of tiles
! -----------------------

    open (10, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
         action =  'read')
    
    read (10, *) maxcat
    
    close (10, status = 'keep')

!
! Loop through tile_id raster
! ___________________________


    allocate (tile_id (1:nc))   
    allocate(veg(1:maxcat,1:6))
    veg = 0.

    dx = nc_esa / nc
    dy = nr_esa / nr

    open (10,file=trim(gfile)//'.rst',status='old',action='read',  &
          form='unformatted',convert='little_endian')
    
    do j=1,nr

       ! read a row

       read(10)tile_id(:)

       do i = 1,nc
          if((tile_id (i) >= 1).and.(tile_id(i)  <= maxcat)) then
             if (associated (subset)) NULLIFY (subset)
             subset => esa_veg((i-1)*dx +1 :i*dx, (j-1)*dy +1:j*dy)
             
             NPLUS = count(subset >= 1 .and. subset <= 230)
             
             if(NPLUS > 0)  then
                allocate (loc_int (1:NPLUS))
                allocate (unq_mask(1:NPLUS))
                loc_int = pack(subset,mask = (subset >= 1 .and. subset <= 230))
                call MAPL_Sort (loc_int)
                unq_mask = .true.
                
                do k = 2,NPLUS 
                   unq_mask(k) = .not.(loc_int(k) == loc_int(k-1))
                end do
                NBINS = count(unq_mask)
                
                allocate(loc_val (1:NBINS))
                allocate(density (1:NBINS))
                loc_val = 1.*pack(loc_int,mask =unq_mask)
                call histogram (size(subset,1)*size(subset,2), NBINS, density, loc_val, real(subset))   
                
                do k = 1, nbins

                   if (density (k) > 0) then
                      esa_type = int (loc_val(k))
                      !                 if (esa_type ==  10)  veg (tile_id(i),10) = 1.* density (k)   ; lakes inland water
                      if (esa_type ==  10)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)   ! inconsistent mask/veg  
                      if (esa_type ==  11)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  14)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  20)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  30)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type ==  40)  veg (tile_id(i), 1) = veg (tile_id(i), 1) + 1.* density (k)
                      if (esa_type ==  50)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 1.* density (k)
                      if (esa_type ==  60)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 1.* density (k)
                      if (esa_type ==  70)  veg (tile_id(i), 3) = veg (tile_id(i), 3) + 1.* density (k)
                      if (esa_type ==  90)  veg (tile_id(i), 3) = veg (tile_id(i), 3) + 1.* density (k)
                      if (esa_type == 100)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 0.5* density (k)
                      if (esa_type == 100)  veg (tile_id(i), 3) = veg (tile_id(i), 3) + 0.5* density (k)
                      if (esa_type == 110)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 0.3* density (k)
                      if (esa_type == 110)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 0.3* density (k)
                      if (esa_type == 110)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 0.4* density (k)
                      if (esa_type == 120)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 0.2* density (k)
                      if (esa_type == 120)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 0.2* density (k)
                      if (esa_type == 120)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 0.6* density (k)
                      if (esa_type == 130)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 1.* density (k)
                      if (esa_type == 140)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      
                      if((j > NINT(real(nr)*(40./180.))).and.(j < NINT(real(nr)*(140./180.)))) then
                         if (esa_type == 150)  veg (tile_id(i),5) = veg (tile_id(i),5) + 0.5* density (k)
                         if (esa_type == 150)  veg (tile_id(i),4) = veg (tile_id(i),4) + 0.5* density (k)
                      else
                         if (esa_type == 150)  veg (tile_id(i),6) = veg (tile_id(i),6) + 0.5* density (k)
                         if (esa_type == 150)  veg (tile_id(i),4) = veg (tile_id(i),4) + 0.5* density (k)
                      end if
                      
                      if((j > NINT(real(nr)*(70./180.))).and.(j < NINT(real(nr)*(110./180.)))) then 
                         if (esa_type == 160)  veg (tile_id(i), 1) = veg (tile_id(i), 1) + 1.* density (k) 
                      else
                         if (esa_type == 160)  veg (tile_id(i), 2) = veg (tile_id(i), 2) + 1.* density (k)
                      end if
                      
                      if (esa_type == 170)  veg (tile_id(i), 1) = veg (tile_id(i), 1) + 1.* density (k)
                      if (esa_type == 180)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type == 190)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)
                      if (esa_type == 200)  veg (tile_id(i), 5) = veg (tile_id(i), 5) + 1.* density (k)
                      if (esa_type == 210)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)  ! inconsistent mask/veg  
                      if (esa_type == 220)  veg (tile_id(i), 4) = veg (tile_id(i), 4) + 1.* density (k)  ! inconsistent mask/veg  
                      !                         if (esa_type == 210)  veg (tile_id(i),11) = 1.* density (k)      ; ocean
                      !                         if (esa_type == 220)  veg (tile_id(i), 9) = 1.* density (k)     ; ice     
                   endif
                enddo
                deallocate (loc_int,unq_mask,loc_val,density)
             endif
          endif
       end do
    end do

    deallocate (tile_id)   
    close (10,status='keep')    

! Canopy height and ASCAT roughness length

    call ascat_r0 (nc,nr,gfile, z0)

    if(jpl_height) then
       call jpl_canoph (nc,nr,gfile, z2)
    else
       allocate (z2(1:maxcat))       
    endif
!
! Now create mosaic_veg_fracs file
! --------------------------------

    allocate (ityp (1:maxcat))       

    open (10,file='clsm/mosaic_veg_typs_fracs',  &
         form='formatted',status='unknown')
    open (11, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
         action =  'read')
    
    read (11, *) maxcat   
 
    do k = 1, maxcat

       read (11,'(i10,i8,5(2x,f9.4))') tid,cid
       tem = 0.
       tem(1:6)=veg (k,1:6)

       if(sum(tem).gt.0)then

          mfrac = -10.
          sfrac = -10.
          mos1  = 100
          mos2  = 100         

          do i = 1,6
             if(mfrac.le.tem(i))then
                sfrac = mfrac
                mos2  = mos1
                mfrac = tem(i)
                mos1  = i
             elseif(sfrac.le.tem(i)) then
                if(tem(i).lt.mfrac)then
                   sfrac = tem(i)
                   mos2  = i
                endif
             endif             
          end do

          mfrac = max (mfrac,0.)
          sfrac = max (sfrac,0.)
          tfrac = mfrac + sfrac
          mfrac = mfrac / tfrac
          sfrac = sfrac / tfrac

          if (mos1 == 100) then
             mos1 = 4
             mos2 = 4
             mfrac= 1.
             sfrac= 0. 
          endif

          if (sfrac == 0.) mos2 = mos1 ! No secondary type
          if(.not.jpl_height) z2(k) = VGZ2(mos1)
          ityp (k) = mos1
          write (10,'(i10,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')     &
            tid,cid,mos1,mos2,100.*mfrac,100.*sfrac, z2(k), z0 (k)
          
       endif
    end do

    close (10,status='keep')    
    close (11,status='keep') 

    inquire(file='clsm/catch_params.nc4', exist=file_exists)

    if(file_exists) then
       status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid                                ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'OLD_ITY'    ) ,(/1/),(/maxcat/), real(ityp)) ; VERIFY_(STATUS)
       STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
    endif
    
    inquire(file='clsm/vegdyn.data', exist=file_exists)

    if(file_exists) then
       status = NF_OPEN ('clsm/vegdyn.data', NF_WRITE, ncid                                 ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ITY'    ) ,(/1/),(/maxcat/), real(ityp)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'Z2CH'   ) ,(/1/),(/maxcat/), z2        ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ASCATZ0') ,(/1/),(/maxcat/), Z0        ) ; VERIFY_(STATUS)
       STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
    else
       open (20,file='clsm/vegdyn.data',status='unknown',action='write',form='unformatted', &
            convert='little_endian')   
       write (20) real(ityp)
       write (20) z2 (:)
       write (20) z0 (:)
       close (20)     
   endif
       
   deallocate (veg, z0, z2, ityp)
   
  END SUBROUTINE ESA2MOSAIC

!
!----------------------------------------------------------------------
!
  SUBROUTINE HISTOGRAM (NLENS, NBINS, density, loc_val, x, BIN)

    implicit none

    integer, intent (in) :: NBINS, NLENS
    real,    intent (in) :: x (NLENS)
    integer, intent (out):: density (NBINS)
    real,    intent (inout) :: loc_val (NBINS)
    real,    intent (in), optional :: bin
    real :: xdum(NLENS), xl, xu, min_value
    integer :: n

    if(present (bin)) min_value =  real(floor(minval(x)))

    DO N = 1, NBINS
       if(present (bin)) then
          xl = (N - 1)*BIN + min_value
          loc_val (n) = xl
          xu = xl + bin
          XDUM = 0.
          where((x >= xl).and.(x < xu))XDUM = 1
       else
          XDUM = 0.
          where(x == loc_val (n)) XDUM = 1
       endif
          density(n) = int(sum(XDUM))       
    END DO

END SUBROUTINE HISTOGRAM

!
!----------------------------------------------------------------------
!

 SUBROUTINE create_mapping (nc,nr,nc_data,nr_data,rmap, gfile)

   implicit none

   integer, intent (in) :: nc,nr,nc_data,nr_data
   type (regrid_map), intent (inout) :: rmap
   integer :: i,j,n, i1,i2,j1,j2,ncatch, nbins, status, NPLUS,pix_count
   REAL,    allocatable, DIMENSION (:) :: loc_val
   INTEGER, ALLOCATABLE, DIMENSION (:) :: density, loc_int
   logical, dimension (:), allocatable :: unq_mask    
   integer, allocatable, target, dimension (:,:) :: tile_id
   integer, pointer , dimension (:,:) :: subset
   real   :: dx_data, dy_data,dx_geos, dy_geos, lon1,lon2,lat1,lat2
   character (*), intent(in) :: gfile
   integer, pointer  :: iraster  (:,:)

! Reading rst file

   open (10,file=trim(gfile)//'.rst',status='old',action='read',  &
        form='unformatted',convert='little_endian')
   allocate (tile_id    (1:nc,1:nr))         
   
   do j=1,nr
      read(10)tile_id(:,j)
   end do
   close (10,status='keep')

! Reading number of catchments

   open (10,file='clsm/catchment.def',status='old',action='read',  &
        form='formatted')
   read (10, *) ncatch
   close (10, status = 'keep')
   
   dx_data = 360./real(nc_data)
   dy_data = 180./real(nr_data)
   dx_geos = 360./real(nc)
   dy_geos = 180./real(nr)

   if((nc_data  >= nc).and.(nr_data >= nr)) then

      allocate(iraster(nc_data,nr_data),stat=STATUS); VERIFY_(STATUS)
      call RegridRaster(tile_id,iraster)
      NPLUS = count(iraster >= 1 .and. iraster <= ncatch)
      allocate (rmap%ij_index(1:nc_data, 1:nr_data), source = 0)
      allocate (rmap%map (1:NPLUS))     
      rmap%map%NT = 0
      pix_count = 1
      do j = 1,nr_data
         do i =  1,nc_data
            if((iraster (i,j) >=1).and.(iraster (i,j) <=ncatch)) then
               rmap%map(pix_count)%NT = 1
               rmap%map(pix_count)%TID  (rmap%map(pix_count)%NT) = iraster (i,j)
               rmap%map(pix_count)%count(rmap%map(pix_count)%NT) = 1
               rmap%ij_index(i,j) = pix_count
               pix_count = pix_count + 1
            endif
         end do
      end do
      deallocate (iraster) ; VERIFY_(STATUS)

   else
      NPLUS = count(tile_id >= 1 .and. tile_id <= ncatch)
      allocate (rmap%ij_index(1:nc_data, 1:nr_data), source = 0)
      allocate (rmap%map (1:NPLUS))     
      rmap%map%NT = 0
      pix_count   = 1
      do j = 1,nr_data
         
         lat1 = -90. + (j-1)*dy_data
         lat2 = lat1 + dy_data
         j1   = floor  ((-90. + (j-1)*dy_data + 90.)/dy_geos) + 1
         j2   = ceiling((-90. + (j)*dy_data + 90.  )/dy_geos) 
         
         do i =  1,nc_data
            
            lon1 = -180. + (i-1)*dx_data
            lon2 = lon1 + dx_data 
            i1   = floor  ((-180. + (i-1)*dx_data + 180.)/dx_geos) + 1
            i2   = ceiling((-180. + (i)*dx_data + 180.  )/dx_geos) 
            
            if(j2 > j1 .or. i2 > i1)  then
               subset => tile_id (i1:i2,j1:j2)
               NPLUS = count(subset >= 1 .and. subset <= ncatch)
               if(NPLUS > 0)  then
                  allocate (loc_int (1:NPLUS))
                  allocate (unq_mask(1:NPLUS))
                  loc_int = pack(subset,mask = (subset >= 1 .and. subset <= ncatch))
                  call MAPL_Sort (loc_int)
                  unq_mask = .true.
                  do n = 2,NPLUS 
                     unq_mask(n) = .not.(loc_int(n) == loc_int(n-1))
                  end do
                  NBINS = count(unq_mask)
                  
                  allocate(loc_val (1:NBINS))
                  allocate(density (1:NBINS))
                  loc_val = 1.*pack(loc_int,mask =unq_mask)
                  call histogram (size(subset,1)*size(subset,2), NBINS, density, loc_val, real(subset))   
                  
                  DO N =1,NBINS
                     if(density(n) > 0) then 
                        rmap%map(pix_count)%NT = rmap%map(pix_count)%NT + 1
                        if(rmap%map(pix_count)%NT > N_tiles_per_cell) then
                           print *,'N_tiles_per_cell exceeded :', rmap%map(pix_count)%NT
                           print *,i,j,i1,i2,j1,j2
                           print *,'NT',rmap%map(pix_count)%NT 
                           print *,rmap%map(pix_count)%TID
                           print *,rmap%map(pix_count)%count
                           stop
                        endif
                        rmap%map(pix_count)%TID  (rmap%map(pix_count)%NT) = NINT(loc_val(n))
                        rmap%map(pix_count)%count(rmap%map(pix_count)%NT) = density(n)
                        
                     endif
                  END DO
                  rmap%ij_index(i,j) = pix_count
                  pix_count = pix_count + 1
                  deallocate (loc_val, density)
                  deallocate (loc_int, unq_mask)
               endif
               NULLIFY (subset)
            else
               if((tile_id (i1,j1) > 0).and.(tile_id(i1,j1).le.ncatch)) then
                  rmap%map(pix_count)%NT       = 1
                  rmap%map(pix_count)%TID(1)   = tile_id (i1,j1)
                  rmap%map(pix_count)%COUNT(1) = 1
                  rmap%ij_index(i,j) = pix_count
                  pix_count = pix_count + 1
               endif
            endif
         end do
      end do
   end if
      
 END SUBROUTINE create_mapping

!
!----------------------------------------------------------------------
!
  SUBROUTINE merge_lai_data (MaskFile)
    implicit none
    type (date_time_type) :: bf_geol2_time,af_geol2_time,date_time_new,bf_lai_time,   &
       af_lai_time
    character (*) :: MaskFile
    integer :: n,k,ntiles,t,ierr
    integer, allocatable, dimension (:) :: pfaf
    ! South AMerica/ Africa/ Australia are from GEOLAND2
    integer :: i1,i2,i3,i4,i5,i6
    integer, parameter :: i1_hydr = 1011000, i2_hydr =  1999900  ! South America
    integer, parameter :: i3_hydr = 3021000, i4_hydr =  3990000  ! Africa 
    integer, parameter :: i5_hydr = 5000142, i6_hydr =  5999900  ! Australia
    integer, parameter :: i1_srtm = 229075 , i2_srtm =  267083   ! South America
    integer, parameter :: i3_srtm = 75369  , i4_srtm =  140751   ! Africa 
    integer, parameter :: i5_srtm = 267084 , i6_srtm =  SRTM_maxcat  ! Australia

    REAL, ALLOCATABLE, dimension (:) :: geol2_lai_bf,geol2_lai_af,geol2_lai, lai    
    real :: dum, gyr,gmn,gdy,gyr1,gmn1,gdy1, slice1,slice2
    real :: yr,mn,dy,yr1,mn1,dy1

    if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
       i1 = i1_srtm
       i2 = i2_srtm
       i3 = i3_srtm
       i4 = i4_srtm
       i5 = i5_srtm
       i6 = i6_srtm
    else
       i1 = i1_hydr
       i2 = i2_hydr
       i3 = i3_hydr
       i4 = i4_hydr
       i5 = i5_hydr
       i6 = i6_hydr
    endif

    open (10, file ='clsm/catchment.def',form='formatted',status='old',action='read')
    read (10,*) ntiles

    allocate (pfaf(1:ntiles))
    allocate (geol2_lai_bf(1:ntiles))
    allocate (geol2_lai_af(1:ntiles))
    allocate (geol2_lai   (1:ntiles))
    allocate (lai         (1:ntiles))  
           
    do n =1,ntiles
       read (10,*) k,pfaf(n)
    end do

    close (10,status='keep')

!
        open (41,file='clsm/lai.GEOLAND2_10-DayClim',  &
         form='unformatted',status='old',convert='little_endian',action='read')
        open (42,file='clsm/lai.MODIS_8-DayClim',      &
         form='unformatted',status='old',convert='little_endian',action='read')
        open (43,file='clsm/lai.dat',      &
         form='unformatted',status='unknown',convert='little_endian',action='write')

       read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1       
       read(41) geol2_lai_bf
       call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,bf_geol2_time)

       read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1
       read(41) geol2_lai_af
       call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,af_geol2_time)

       do t = 1, 48

          read(42) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1 
          read(42) lai 
          write(43) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(ntiles),1.
	  call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,date_time_new)

!          date_time_new%year   = nint(yr) + 2001 
!          date_time_new%month  = nint(mn)
!          date_time_new%day    = nint(dy)
!          date_time_new%hour   = 0            
!          date_time_new%min    = 0            
!          date_time_new%sec    = 0 
!          call get_dofyr_pentad(date_time_new)   
          if (datetime_le_refdatetime(date_time_new,af_geol2_time)) then
           
          else
             read(41,IOSTAT=ierr) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1
             if(ierr == 0) then
                geol2_lai_bf = geol2_lai_af
                read(41) geol2_lai_af
                bf_geol2_time = af_geol2_time
                call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,af_geol2_time)
             else
                print *,'END OF GEOL2 LAI FILE'
                stop
             endif
          endif

          if(t==1) then 
             date_time_new%year = date_time_new%year + 1
             geol2_lai_af       = geol2_lai_bf
             af_geol2_time      = bf_geol2_time
             af_geol2_time%year = af_geol2_time%year + 1

             do k = 1,34 
                read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1       
                read(41) geol2_lai_bf
                call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,bf_geol2_time)                
             end do
          endif

!          print *,t
!          print *,'DATE_TIME_NEW :',date_time_new
!          print *,'bf_geol2_time :',bf_geol2_time
!          print *,'af_geol2_time :',af_geol2_time

          call Time_Interp_Fac (date_time_new, bf_geol2_time, af_geol2_time, slice1, slice2)
          geol2_lai    = (slice1*geol2_lai_bf + slice2*geol2_lai_af)
          
          if(t == 1) then
             rewind(41)
             read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1       
             read(41) geol2_lai_bf
             call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,bf_geol2_time)
             
             read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1
             read(41) geol2_lai_af
             call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,af_geol2_time)
          endif

! replace South America with GEOLAND2

          DO n =1,ntiles             
             if((pfaf(n) >= i1).and.(pfaf(n) <= i2)) lai(n) = geol2_lai(n)
             if((pfaf(n) >= i3).and.(pfaf(n) <= i4)) lai(n) = geol2_lai(n)
             if((pfaf(n) >= i5).and.(pfaf(n) <= i6)) lai(n) = geol2_lai(n)
          end do
          write (43) lai(:)
       end do

       close (41,status = 'keep')
       close (42,status = 'keep')
       close (43,status = 'keep')

       deallocate (pfaf,geol2_lai_bf, geol2_lai_af,geol2_lai,lai)

  END SUBROUTINE merge_lai_data
!
!----------------------------------------------------------------------
!
  SUBROUTINE modis_scale_para_high (ease_grid,MA,gfile)
    
    implicit none
    type (date_time_type) :: gf_green_time,af_green_time,end_time, &
              bf_lai_time,af_lai_time,date_time_new,bf_modis_time,   &
              af_modis_time
    logical :: ease_grid
    character*6 :: MA
    CHARACTER*20 :: version,resoln,continent
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: tsteps,zth, slr,tarea
    INTEGER :: typ,j_dum,ierr,indr1,ip2
    character*100 :: path,fname,fout,metpath
    character (*) :: gfile
    integer :: n,maxcat,ip
    integer :: yy,j,month
    integer, allocatable, dimension (:) :: vegcls 
    real, allocatable, dimension (:) :: &
         modisvf, modisnf,albvf,albnf, &
         green,lai,lai_before,lai_after,grn_before,grn_after
    real, allocatable, dimension (:) :: &
         calbvf,calbnf, zero_array, one_array, albvr,albnr
    character*300 :: ifile1,ifile2,ofile
    integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    integer :: day, hour, min ,secs_in_day,k
    real :: yr,mn,dy,yr1,mn1,dy1,dum, slice1,slice2
    logical :: save_sib = .false.

    fname='clsm/catchment.def'
    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)maxcat
    allocate (albvf    (1:maxcat))
    allocate (albnf    (1:maxcat))
    allocate (calbvf   (1:maxcat))
    allocate (calbnf   (1:maxcat))
    allocate (modisvf  (1:maxcat))
    allocate (modisnf  (1:maxcat))
    allocate (lai      (1:maxcat))
    allocate (green    (1:maxcat))
    allocate (lai_before (1:maxcat))
    allocate (grn_before (1:maxcat))
    allocate (lai_after  (1:maxcat))
    allocate (grn_after  (1:maxcat))
    allocate (vegcls     (1:maxcat))
    allocate (zero_array (1:maxcat))
    allocate (one_array  (1:maxcat))
    allocate (albvr      (1:maxcat))
    allocate (albnr      (1:maxcat))
    close (10,status='keep')

    fname=trim(gfile)//'.til'
    open (10,file=fname,status='old',action='read',form='formatted')

    fname='clsm/mosaic_veg_typs_fracs'
    open (20,file=fname,status='old',action='read',form='formatted')

    read (10,*)ip
    read (10,*)j_dum

    do n = 1, j_dum
       read (10,'(a)')version
       read (10,*)nc_gcm
       read (10,*)nr_gcm
    end do    
    
    do n = 1,ip
      if (ease_grid) then     
         read(10,*,IOSTAT=ierr) typ !,pfs,lont,latt,ig,jg,fr_gcm
      else
         !read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
         !   typ,tarea,lont,latt,ig,jg,fr_gcm,indx_dum,pfs,j_dum,fr_cat,j_dum
         read(10,*,IOSTAT=ierr) typ
      endif
       if (typ == 100) then
          ip2 = n 
          !read (20,'(i10,i8,2(2x,i3),2(2x,f6.4))')     &
          !  indr1,indr1,vegcls(ip2),indr1,fr_gcm,fr_gcm
          read (20,*,IOSTAT=ierr) indr1,indr1,vegcls(ip2)
       endif
       if(ierr /= 0)write (*,*)'Problem reading', n, ease_grid
    end do
    close (10,status='keep')
    close (20,status='keep')

    albvf    =0.
    albnf    =0.
    calbvf   =0.
    calbnf   =0.
    modisvf  =0.
    modisnf  =0.
    zero_array = 0.
    one_array  = 1.
    albvr      = 0.
    albnr      = 0.

! MODIS Albedo files
    if(MA == 'MODIS1') then 
       open (10,file='clsm/AlbMap.WS.16-day.tile.0.3_0.7.dat',&
            form='unformatted',convert='little_endian', &
            action='read',status='old')   
       open (11,file='clsm/AlbMap.WS.16-day.tile.0.7_5.0.dat',&
            form='unformatted',convert='little_endian', &
            action='read',status='old') 
    endif

    if(MA == 'MODIS2') then
       open (10,file='clsm/AlbMap.WS.8-day.tile.0.3_0.7.dat',&
            form='unformatted',convert='little_endian', &
            action='read',status='old')   
       open (11,file='clsm/AlbMap.WS.8-day.tile.0.7_5.0.dat',&
            form='unformatted',convert='little_endian', &
            action='read',status='old') 
    endif

    read(10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read (10) modisvf (:)
    read (11) modisnf (:)

! SiB Albedo Parameterization files
    if (save_sib) then
       open (20,file='clsm/sib_visdf.dat',convert='little_endian', &
             action='write',status='unknown',form='unformatted')
       open (21,file='clsm/sib_nirdf.dat',convert='little_endian', &
             action='write',status='unknown',form='unformatted')  
       write(20) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
       write(21) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1. 
    endif

! MODIS scale parameter files
    open (30,file='clsm/visdf.dat',convert='little_endian', &
             action='write',status='unknown',form='unformatted')
    open (31,file='clsm/nirdf.dat',convert='little_endian', &
             action='write',status='unknown',form='unformatted')
    write(30) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
    write(31) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,date_time_new)

    bf_modis_time   = date_time_new

    read(10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1

    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_modis_time)
    end_time        = af_modis_time
    end_time%year   = NINT(yr1) + 2001 + 1 

    rewind (10)
    rewind (11)

    read(10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read (10) modisvf (:)
    read (11) modisnf (:)

    fname='clsm/lai.dat'
    open (40,file=fname,status='old',action='read',form='unformatted', &
         convert='little_endian')
    read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(40) lai_before
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,bf_lai_time)
 
    read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(40) lai_after
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)
 
   if(date_time_new%dofyr < bf_lai_time%dofyr) then
       do while ((date_time_new%dofyr >  af_lai_time%dofyr)) 
          lai_before = lai_after
          bf_lai_time = af_lai_time
          read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
          read(40) lai_after
          call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)
        end do
    endif

    fname='clsm/green.dat'
    open (41,file=fname,status='old',action='read',form='unformatted', &
         convert='little_endian')
    read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(41) grn_before
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,gf_green_time)
    
    read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(41) grn_after
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)

    calbvf   =0.
    calbnf   =0.
    albvf    =0.
    albnf    =0.
    tsteps   =0.    

    do while (datetime_le_refdatetime(date_time_new,end_time))
       
!       write (*,'(a48,i4.4,i2.2,i2.2)') '     Computing MODIS scale parameters for month: ', &
            
              
          if (datetime_le_refdatetime(date_time_new,af_lai_time)) then

          else

            read(40,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1 
            if(ierr == 0) then 
               lai_before = lai_after
               read(40) lai_after
               bf_lai_time = af_lai_time
               call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)
            else
               rewind(40) 
               read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
               read(40) lai_before
               call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,bf_lai_time)
               read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
               read(40) lai_after
               call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)

               if(date_time_new%dofyr < bf_lai_time%dofyr) then
                  do while ((date_time_new%dofyr >  af_lai_time%dofyr))
                     lai_before = lai_after
                     bf_lai_time = af_lai_time
                     read(40) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                     read(40) lai_after
                     call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)
                  end do
               endif
            endif
         endif
         call Time_Interp_Fac (date_time_new, bf_lai_time, af_lai_time, slice1, slice2)
         lai    = (slice1*lai_before + slice2*lai_after)

         if (datetime_le_refdatetime(date_time_new,af_green_time)) then
            
         else
            
            read(41,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
            if(ierr == 0) then 
               grn_before = grn_after
               gf_green_time = af_green_time
               read(41) grn_after
               call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)
            endif
         endif
!            else
!               rewind(41) 
!               read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
!               read(41) grn_before
!               gf_green_time%month = NINT(mn)
!               gf_green_time%day   = NINT(dy)
!               call get_dofyr_pentad(gf_green_time)
!               af_green_time%month = NINT(mn1)
!               af_green_time%day   = NINT(dy1)
!               call get_dofyr_pentad(af_green_time)    
!               if(date_time_new%dofyr < gf_green_time%dofyr) then
!                  do while ((date_time_new%dofyr >  af_green_time%dofyr)) 
!                     read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
!                     read(41) grn_before
!                     gf_green_time%year  = date_time_new%year
!                     gf_green_time%month = NINT(mn)
!                     gf_green_time%day   = NINT(dy)
!                     call get_dofyr_pentad(gf_green_time)
!                     af_green_time%year  = date_time_new%year
!                     if ((yr1-yr) == 1.)af_green_time%year = af_green_time%year+1
!                     af_green_time%month = NINT(mn1)
!                     af_green_time%day   = NINT(dy1)
!                     call get_dofyr_pentad(af_green_time)           
!                  end do
!                  read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
!                  read(41) grn_after
!               endif
!            endif
!         endif

         call Time_Interp_Fac (date_time_new, gf_green_time, af_green_time, slice1, slice2)
         green  = (slice1*grn_before + slice2*grn_after)
          
         call sibalb (                                 &
              MAXCAT,vegcls,lai,green, zero_array,     &
              one_array,one_array,one_array,one_array, &
              ALBVR, ALBNR, albvf, albnf)
         
          calbvf = calbvf + albvf
          calbnf = calbnf + albnf         
          tsteps = tsteps + 1.  
          call augment_date_time( 86400, date_time_new ) 

          if (datetime_le_refdatetime(date_time_new,af_modis_time)) then
             
          else
             bf_modis_time = af_modis_time
             calbvf = calbvf/tsteps
             calbnf = calbnf/tsteps

             modisvf = modisvf/(calbvf + 1.e-20)
             modisnf = modisnf/(calbnf + 1.e-20)

             do n =1, maxcat
!                if(modisvf(n).le.0)print *,'Negative MODISVF scale param at cell',n, modisvf(n)
!                if(modisnf(n).le.0)print *,'Negative MODISNF scale param at cell',n, modisnf(n)
!                if(modisvf(n).gt.100)print *,'Too large MODISVF scale param at cell',n, modisvf(n)
!                if(modisnf(n).gt.100)print *,'Too large MODISNF scale param at cell',n, modisnf(n)
                if(modisvf(n).le.0.) modisvf(n) = 1.
                if(modisnf(n).le.0.) modisnf(n) = 1.
                if(modisvf(n).gt.100)modisvf(n)= 1.
                if(modisnf(n).gt.100)modisnf(n)= 1.
            enddo

            if (save_sib) then
               write (20) calbvf (:)
               write (21) calbnf (:)               
            endif

            write (30) modisvf (:)
            write (31) modisnf (:)

	    read(10,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1       
  
            if(ierr == 0) then
	       
               read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
               read (10) modisvf (:)
               read (11) modisnf (:)  
               write(30) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
               write(31) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.

               if (save_sib) then
                  write(20) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
                  write(21) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1. 
               endif
               
 	       bf_modis_time = af_modis_time
               call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_modis_time)
               calbvf   =0.
               calbnf   =0.
               albvf    =0.
               albnf    =0.
               tsteps   =0.   
            endif
         endif
      end do
  
      deallocate (modisvf,modisnf,albvf,albnf)
      deallocate (green,lai)
      deallocate (vegcls)
      deallocate (calbvf,calbnf)
      deallocate (lai_before,grn_before, lai_after,grn_after)
      deallocate (zero_array, one_array, albvr, albnr)

      close (10, status='keep')
      close (11, status='keep')        
      close (30, status='keep')
      close (31, status='keep')
      if (save_sib) then
         close (20, status='keep')
         close (21, status='keep')
      endif
  
END SUBROUTINE modis_scale_para_high

!
! ---------------------------------------------------------------------------------------
! 
  SUBROUTINE modis_alb_on_tiles_high (nc_data,nr_data,rmap,MA,gfiler)
!
! Processing MODIS Albedo and creating 8-day climatological data 
!
  implicit none 
  integer, intent (in) :: nc_data,nr_data
  type (regrid_map), intent (in) :: rmap
  character*6 :: MA
  character(*)  :: gfiler
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj, pix_count
  integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
      time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, dimension (:,:) :: &
         net_data1,net_data2
  REAL, ALLOCATABLE, dimension (:) :: vec_AlbVis, count_AlbVis,vec_AlbNir, count_AlbNir
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  character *10 :: vname
  REAL :: sf

!
! Reading number of cathment-tiles from catchment.def file
! 
      fname='clsm/catchment.def' 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      close (10,status='keep')

      if(MA=='MODIS1') fname =trim(c_data)//'MODIS-Albedo/MODISalb.c004.v2.WS_H11V13.nc'
      if(MA=='MODIS2') fname =trim(c_data)//'MODIS-Albedo2/MCD43GF_wsa_H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS)
      allocate (MMDD      (0: n_tslices + 1))
      allocate (MMDD_next (0: n_tslices + 1))

      status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)

      if(nc_data/=i_highd .or. nr_data/=j_highd) then
         print *,'Inconsistent mapping and dimensions in modis_alb_on_tiles_high   -so stopping ...'
         stop
      end if

      mmdd(0) = mmdd(n_tslices)
      mmdd(n_tslices + 1)= mmdd(1)

      mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
      mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)


      allocate(net_data1 (1:nc_10,1:nr_10))
      allocate(net_data2 (1:nc_10,1:nr_10))       

       !
       ! reading Albedo data
       !

       if(MA == 'MODIS1') then 
          open (31,file='clsm/AlbMap.WS.16-day.tile.0.3_0.7.dat',        &
               form='unformatted',status='unknown',convert='little_endian')
          open (32,file='clsm/AlbMap.WS.16-day.tile.0.7_5.0.dat',        &
               form='unformatted',status='unknown',convert='little_endian')
       endif

       if(MA == 'MODIS2') then 
          open (31,file='clsm/AlbMap.WS.8-day.tile.0.3_0.7.dat',        &
               form='unformatted',status='unknown',convert='little_endian')
          open (32,file='clsm/AlbMap.WS.8-day.tile.0.7_5.0.dat',        &
               form='unformatted',status='unknown',convert='little_endian')
       endif

       allocate(vec_AlbVis(maxcat))
       allocate(count_AlbVis(1:maxcat))    
       allocate(vec_AlbNir(maxcat))
       allocate(count_AlbNir(1:maxcat))    

       do t =0,n_tslices+1
       
          time_slice = t
          yr = 1
          yr1= 1
          if(t == 0) then
             time_slice =  n_tslices
             yr         =  1 - 1
          endif

          if(t >= n_tslices) then 
             yr1 = 1 + 1
             if(t ==n_tslices + 1) then
                time_slice =  1
                yr = 1 + 1
             endif
          endif

          read(mmdd(t),'(i2.2,i2.2)') mn,dd
          read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1

          ! Reading, interpolating or aggregating on to catchment-tiles
          
          vec_AlbVis   =0.
          count_AlbVis = 0.
          vec_AlbNir   =0.
          count_AlbNir = 0. 

          do jx = 1,18	
             do ix = 1,36
                write (vv,'(i2.2)')jx
                write (hh,'(i2.2)')ix 
                if(MA=='MODIS1') fname =trim(c_data)//'MODIS-Albedo/MODISalb.c004.v2.WS_H'//hh//'V'//vv//'.nc'
                if(MA=='MODIS2') fname =trim(c_data)//'MODIS-Albedo2/MCD43GF_wsa_H'//hh//'V'//vv//'.nc'
                status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
                if(status == 0) then
                   status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                   status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
                   status = NF_GET_att_INT  (ncid,4,'UNDEF',d_undef); VERIFY_(STATUS)
                   status = NF_GET_att_REAL (ncid,4,'ScaleFactor',sf); VERIFY_(STATUS)
                   status = NF_GET_VARA_INT (ncid,4,(/1,1,time_slice/),(/nc_10,nr_10,1/),net_data1); VERIFY_(STATUS)
                   status = NF_GET_VARA_INT (ncid,5,(/1,1,time_slice/),(/nc_10,nr_10,1/),net_data2); VERIFY_(STATUS)

                   do j = jLL,jLL + nr_10 -1 
                      do i = iLL, iLL + nc_10 -1
                         pix_count = rmap%ij_index(i,j)
                         if (pix_count ==0) cycle
                         if(net_data1(i-iLL +1 ,j - jLL +1) > 0) then
                            
                            if(rmap%map(pix_count)%nt > 0) then
                               do n = 1, rmap%map(pix_count)%nt
                                  vec_AlbVis(rmap%map(pix_count)%tid(n))  = vec_AlbVis(rmap%map(pix_count)%tid(n)) +  &
                                       sf*net_data1(i-iLL +1 ,j - jLL +1)*rmap%map(pix_count)%count(n) 
                                  count_AlbVis(rmap%map(pix_count)%tid(n))= count_AlbVis(rmap%map(pix_count)%tid(n)) + &
                                        1.*rmap%map(pix_count)%count(n)
                               end do
                            endif
                         endif
                         if(net_data2(i-iLL +1 ,j - jLL +1) > 0) then
                            if(rmap%map(pix_count)%nt > 0) then
                               do n = 1, rmap%map(pix_count)%nt
                                  vec_AlbNir(rmap%map(pix_count)%tid(n))  = vec_AlbNir(rmap%map(pix_count)%tid(n)) +  &
                                       sf*net_data2(i-iLL +1 ,j - jLL +1)*rmap%map(pix_count)%count(n) 
                                  count_AlbNir(rmap%map(pix_count)%tid(n))= count_AlbNir(rmap%map(pix_count)%tid(n)) + &
                                        1.*rmap%map(pix_count)%count(n)
                               end do
                            endif
                          endif
                      enddo
                   enddo
                   status = NF_CLOSE(ncid)
                endif
             end do
          end do

          DO n =1,maxcat
             if(count_AlbVis(n)/=0.) vec_AlbVis(n)=vec_AlbVis(n)/count_AlbVis(n)
             if(count_AlbNir(n)/=0.) vec_AlbNir(n)=vec_AlbNir(n)/count_AlbNir(n)
          END DO
          write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
          write(31)  vec_AlbVis(:)
          write(32) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
          write(32)  vec_AlbNir(:)
       end do
       close(31,status='keep')
       close(32,status='keep')

       deallocate (net_data1,net_data2)
       deallocate (count_AlbVis, count_AlbNir)
       deallocate (vec_AlbVis, vec_AlbNir)
   
     END SUBROUTINE modis_alb_on_tiles_high
!
! ---------------------------------------------------------------------------------------
! 
  SUBROUTINE hres_lai (nx,ny,gfiler,c_data,lai_name,merge)
!
! Processing GEOLAND2/MODIS LAI and creating 10-day climatological data 
!
  implicit none 
  integer, intent (in) :: nx, ny 
  character(*)  :: gfiler,c_data,lai_name
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj,ierr
  integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
      time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2
  real :: dum, gyr,gmn,gdy,gyr1,gmn1,gdy1, slice1,slice2
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, dimension (:,:) :: &
         net_data1
  integer (kind=2) , allocatable, target, dimension (:,:) :: LAI_HIGH
  integer (kind=2), pointer, dimension (:,:) :: Raster
  REAL, ALLOCATABLE, dimension (:) :: vec_lai, count_lai
  REAL, ALLOCATABLE, dimension (:) :: gswp2_lai_bf,gswp2_lai_af,gswp2_lai 
  integer, allocatable, target, dimension (:,:) :: tile_id
  integer, pointer :: iRaster(:,:)
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  REAL :: sf
  logical :: first_entry = .true.
  type (date_time_type) :: bf_gswp2_time,af_gswp2_time,date_time_new,bf_lai_time,   &
       af_lai_time
  integer, intent(in), optional :: merge 
!
! Reading number of cathment-tiles from catchment.def file
!--------------------------------------------------------- 
      if (first_entry) then
          nullify(iraster) ; first_entry = .false.
      end if

      fname='clsm/catchment.def' 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      close (10,status='keep')

      fname =trim(c_data)//trim(lai_name)//'lai_clim.H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS)
      allocate (MMDD      (0: n_tslices + 1))
      allocate (MMDD_next (0: n_tslices + 1))

      status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)

       mmdd(0) = mmdd(n_tslices)
       mmdd(n_tslices + 1)= mmdd(1)

       mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
       mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)


       allocate(tile_id(1:nx,1:ny))
       allocate(net_data1 (1:nc_10,1:nr_10))
       
       fname=trim(gfiler)//'.rst'
       !          
       ! Reading tile-id raster file
       !
       open (10,file=fname,status='old',action='read',  &
            form='unformatted',convert='little_endian')
       
       do j=1,ny
          read(10)tile_id(:,j)
       end do
       
       close (10,status='keep')

       !
       ! writing GEOLAND2 LAI data
       !

       if(present(merge)) then
       open (31,file='clsm/lai.'//lai_name(1:index(lai_name,'/')-1),  &
            form='unformatted',status='unknown',convert='little_endian')
       else
       open (31,file='clsm/lai.dat',  &
            form='unformatted',status='unknown',convert='little_endian')
       endif

       allocate(vec_lai(maxcat))
       allocate(lai_high(1:i_highd,1:j_highd))  
       allocate(count_lai(1:maxcat))    
       allocate(gswp2_lai_bf (1:maxcat))    
       allocate(gswp2_lai_af (1:maxcat))    
       allocate(gswp2_lai    (1:maxcat))   

       !
       ! reading GSWP2 LAI data
       !
       
       open (41,file='clsm/lai.gswp2',  &
         form='unformatted',status='old',convert='little_endian',action='read')
       read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1
       read(41) gswp2_lai_bf
       call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,bf_gswp2_time)

       read(41) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1
       read(41) gswp2_lai_af
       call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,af_gswp2_time)

       do t =0,n_tslices+1
       
          time_slice = t
          yr = 1
          yr1= 1
          if(t == 0) then
             time_slice =  n_tslices
             yr         =  1 - 1
          endif

          if(t >= n_tslices) then 
             yr1 = 1 + 1
             if(t ==n_tslices + 1) then
                time_slice =  1
                yr = 1 + 1
             endif
          endif

          read(mmdd(t),'(i2.2,i2.2)') mn,dd
          read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1
 
          lai_high = -9999

          do jx = 1,18	
             do ix = 1,36
                write (vv,'(i2.2)')jx
                write (hh,'(i2.2)')ix 
                fname = trim(c_data)//trim(lai_name)//'lai_clim.H'//hh//'V'//vv//'.nc'
                status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
                if(status == 0) then
                   status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                   status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
                   status = NF_GET_att_INT  (ncid,4,'UNDEF',d_undef); VERIFY_(STATUS)
                   status = NF_GET_att_REAL (ncid,4,'ScaleFactor',sf); VERIFY_(STATUS)
                   status = NF_GET_VARA_INT (ncid, 4,(/1,1,time_slice/),(/nc_10,nr_10,1/),net_data1); VERIFY_(STATUS)
                   
                   do j = jLL,jLL + nr_10 -1 
                      do i = iLL, iLL + nc_10 -1 
                         if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) &
                              lai_high(i,j) = net_data1(i-iLL +1 ,j - jLL +1)
                      enddo
                   enddo
                   status = NF_CLOSE(ncid)
                endif
             end do
          end do
          
          ! Regridding 
          
          nx_adj = nx
          ny_adj = ny
          
          regrid = nx/=i_highd .or. ny/=j_highd
          
          if(regrid) then
             if(nx > i_highd) then 
                allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
                call RegridRaster2(lai_high,raster)	
                iRaster => tile_id
                if(ny < j_highd) then
                   print *,'nx > i_highd and ny < j_highd'
                   stop 
                endif
             else
                if(.not. associated(iraster)) then
                   allocate(iraster(i_highd,j_highd),stat=STATUS); VERIFY_(STATUS)  
                endif

!		if( associated(iraster)) deallocate(iraster)
!	       allocate(iraster(i_highd,j_highd),stat=STATUS); VERIFY_(STATUS)    
               call RegridRaster(tile_id,iraster)	
                raster  => lai_high
                nx_adj = i_highd
                ny_adj = j_highd
                
                if(ny > j_highd) then
                   print *,'nx < i_highd and ny > j_highd'
                   stop 
                endif
             endif
          else
             raster  => lai_high
             iRaster => tile_id
          end if
          
          ! Interpolation or aggregation on to catchment-tiles
          
          vec_lai =0.
          count_lai = 0.
          
          do j=1,ny_adj
             do i=1,nx_adj
                if((iRaster(i,j).gt.0).and.(iRaster(i,j).le.maxcat)) then
                   if ((raster(i,j).ge.0)) then
                      vec_lai(iRaster(i,j)) = &
                           vec_lai(iRaster(i,j)) + sf*raster(i,j)
                      count_lai(iRaster(i,j)) = &
                           count_lai(iRaster(i,j)) + 1. 
                   endif
                endif
             end do
          end do

          write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
	  call Get_MidTime(real(yr),real(mn),real(dd),real(yr1),real(mn1),real(dd1),date_time_new)
!          date_time_new%year   = yr + 2001
!          date_time_new%month  = mn
!          date_time_new%day    = dd
!          date_time_new%hour   = 0            
!          date_time_new%min    = 0            
!          date_time_new%sec    = 0 
!          call get_dofyr_pentad(date_time_new)             

          if (datetime_le_refdatetime(date_time_new,af_gswp2_time)) then
             
          else
             read(41,IOSTAT=ierr) gyr,gmn,gdy,dum,dum,dum,gyr1,gmn1,gdy1
             if(ierr == 0) then
                gswp2_lai_bf = gswp2_lai_af
                read(41) gswp2_lai_af
                bf_gswp2_time = af_gswp2_time
                call Get_MidTime(gyr,gmn,gdy,gyr1,gmn1,gdy1,af_gswp2_time)
             else
                print *,'END OF GSWP2 LAI FILE'
                stop
             endif
          endif

          call Time_Interp_Fac (date_time_new, bf_gswp2_time, af_gswp2_time, slice1, slice2)
          gswp2_lai    = (slice1*gswp2_lai_bf + slice2*gswp2_lai_af)
          
!          print *, 'Merging GEOLAND2-AVHRR'
!          print *,  bf_gswp2_time
!          print *,  date_time_new
!          print *,  af_gswp2_time
!          print *,  slice1, slice2
!          print *, maxval(gswp2_lai), minval(gswp2_lai)

          DO n =1,maxcat
             if(count_lai(n)/=0.) vec_lai(n)= vec_lai(n)/count_lai(n)
             if(vec_lai(n)==0.) vec_lai(n)  = gswp2_lai(n)
          END DO
          
          write(31)  vec_lai(:)
       end do
       close(31,status='keep')
       close(41,status='keep')

       deallocate (net_data1)
       deallocate (LAI_HIGH)
       deallocate (count_lai)
       deallocate (vec_lai, iRaster)
       deallocate (gswp2_lai_bf,gswp2_lai_af,gswp2_lai, tile_id)

  END SUBROUTINE hres_lai
!
! ---------------------------------------------------------------------------------------
! 
  SUBROUTINE grid2tile_modis6 (nc_data,nr_data,ncol,nrow,gfiler,lai_name)
!
! Processing GEOLAND2/MODIS LAI and creating 10-day climatological data 
!
  implicit none 
  integer, intent (in) :: nc_data,nr_data, ncol,nrow
  real, parameter :: dxy = 1.
  integer :: QSize
  character(*)  :: gfiler,lai_name
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj,ierr,nx,ny
  integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
      time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2,tindex1,pfaf1
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, target,  dimension (:,:) :: net_data1
  integer, pointer, dimension (:,:) :: QSub
  real,    pointer, dimension (:,:)    :: subset
  REAL, ALLOCATABLE, dimension (:)      :: vec_lai, count_lai,tile_lon, tile_lat &
       , x, y !, distance
  real, allocatable, target, dimension (:,:) :: lai_grid
  INTEGER ::imn,imx,jmn,jmx,mval,d1,d2,l
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  REAL :: sf, dum,dist_save,tile_distance,minlat,maxlat,minlon,maxlon
  logical :: first_entry = .true.
  type (date_time_type) :: date_time_new,bf_lai_time,   &
       af_lai_time
  integer, dimension (:,:), allocatable, target :: tile_id
  integer ::  tileid_tile
  real    :: dxm, dym
! Reading rst file
!-----------------
   open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
        form='unformatted',convert='little_endian')
   allocate (tile_id    (1:ncol,1:nrow))         
   
   do j=1,nrow
      read(10)tile_id(:,j)
   end do
   close (10,status='keep')

   dxm = real(nc_data) /real(ncol) 
   dym = real(nr_data) /real(nrow)

   if ((mod( nc_data, ncol) /= 0).OR. (mod( nc_data, ncol) /= 0)) then
      print *, 'For now, 86400 should be evenly divisible by NC Talk to Sarith'
      stop
   endif
!
! Reading number of cathment-tiles from catchment.def file
!_________________________________________________________ 
!
      fname='clsm/catchment.def' 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      allocate (tile_lon(1:maxcat)) 
      allocate (tile_lat(1:maxcat)) 
      
      do n = 1, maxcat
         read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
         tile_lon(n) = (minlon + maxlon)/2.
         tile_lat(n) = (minlat + maxlat)/2.
      end do      
      close (10,status='keep')

      fname =trim(c_data)//trim(lai_name)//'lai_clim.H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS) 
      allocate (MMDD      (0: n_tslices + 1))
      allocate (MMDD_next (0: n_tslices + 1))

      status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)
 
      if(nc_data/=i_highd .or. nr_data/=j_highd) then
         print *,'Inconsistent mapping and dimensions in  hres_lai_no_gswp  -so stopping ...'
         stop
      end if
      
      mmdd(0) = mmdd(n_tslices)
      mmdd(n_tslices + 1)= mmdd(1)
      
      mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
      mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)
      
      
      allocate(net_data1   (1:nc_10,1:nr_10))
    
      ! writing MODIS6
      !
      open (31,file='clsm/lai.dat',  &
           form='unformatted',status='unknown',convert='little_endian')
     
      allocate (vec_lai     (maxcat))
      allocate (count_lai (1:maxcat))
 
!      allocate (vec_fill    (maxcat))
!      allocate (distance    (maxcat))
!      allocate (vec_lai_save(maxcat))
!       vec_fill = 0

      nx = nint (360./dxy)
      ny = nint (180./dxy)
      allocate (x(1:nx))
      allocate (y(1:ny))

      FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
      FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy

      allocate (lai_grid (1 : nx, 1 : ny)) 
      
      QSize = nint(dxy*nc_data/360.)
!      allocate (QSub (1:QSize,1:QSize))

      do t =0,n_tslices+1
         
         time_slice = t
         yr = 1
         yr1= 1
         if(t == 0) then
            time_slice =  n_tslices
            yr         =  1 - 1
         endif
         
         if(t >= n_tslices) then 
            yr1 = 1 + 1
            if(t ==n_tslices + 1) then
               time_slice =  1
               yr = 1 + 1
            endif
         endif
         
         read(mmdd(t),'(i2.2,i2.2)') mn,dd
         read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1
         
         ! Reading Interpolation or aggregation on to catchment-tiles
         
         vec_lai   = -9999.
         count_lai = 0.
         lai_grid  = -9999

         do jx = 1,18	
            do ix = 1,36
               write (vv,'(i2.2)')jx
               write (hh,'(i2.2)')ix 
               fname = trim(c_data)//trim(lai_name)//'lai_clim.H'//hh//'V'//vv//'.nc'
               status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
               if(status == 0) then
                  status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                  status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
                  status = NF_GET_att_INT  (ncid,4,'UNDEF',d_undef); VERIFY_(STATUS)
                  status = NF_GET_att_REAL (ncid,4,'ScaleFactor',sf); VERIFY_(STATUS)
                  status = NF_GET_VARA_INT (ncid, 4,(/1,1,time_slice/),(/nc_10,nr_10,1/),net_data1); VERIFY_(STATUS)
                  
                  do j = jLL,jLL + nr_10 -1 
                     do i = iLL, iLL + nc_10 -1 
                        if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) then
                           tileid_tile = tile_id (ceiling(i/dxm), ceiling (j/dym))
                           if((tileid_tile >= 1).and.(tileid_tile <= maxcat)) then
                                 if(vec_lai(tileid_tile) == -9999.) vec_lai(tileid_tile) = 0.                                 
                                 vec_lai(tileid_tile)   = vec_lai(tileid_tile) + &
                                      sf*net_data1(i-iLL +1 ,j - jLL +1)
                                 count_lai(tileid_tile) = &
                                      count_lai(tileid_tile) + 1.                                     
                           endif
                        endif
                     enddo
                  enddo

! After experimenting with few finer methods, in order to reduce the time taken by the gap filling procedure,
! creating a 0.25-degree gridded data set from finer LAI data and use it for filling the gaps seems the most practical/manageble method.
!---------------------------------------------------------------------------------------------------------------------------------------
                  do j = ceiling(1.*jLL/QSize),ceiling(1.*jLL/QSize) -1 + nr_10/QSize
                     do i = ceiling(1.*iLL/QSize),ceiling(1.*iLL/QSize) -1 + nc_10/QSize 
                         QSub  => net_data1((i-1)*QSize+2-iLL :i*QSize-iLL+1, (j-1)*QSize+2-jLL :j*QSize-jLL+1) 
                         if(maxval (QSub) > 0) lai_grid(i,j) = sf*sum(QSub, QSub>0)/(max(1,count(QSub>0)))
                     enddo
                  enddo                  
                  status = NF_CLOSE(ncid)
               endif
             end do
          end do
          
          NULLIFY (QSub)
          
          write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
          
          where (count_lai > 0.) vec_lai = vec_lai/count_lai

! Filling gaps
!-------------
          DO n =1,maxcat
             if(count_lai(n)==0.)  then 
                
                DO i = 1,nx - 1
                   if ((tile_lon(n) >= x(i)).and.(tile_lon(n) < x(i+1))) ix = i
                end do
                DO i = 1,ny -1
                   if ((tile_lat(n) >= y(i)).and.(tile_lat(n) < y(i+1))) jx = i
                end do
                
                l = 1
                do 
                  imx=ix + l
                  imn=ix - l
                  jmn=jx - l
                  jmx=jx + l
                  imn=MAX(imn,1)
                  jmn=MAX(jmn,1)
                  imx=MIN(imx,nx)
                  jmx=MIN(jmx,ny)
                  d1=imx-imn+1
                  d2=jmx-jmn+1
                  subset => lai_grid(imn: imx,jmn:jmx)

                  if(maxval(subset) > 0.) then 
                     vec_lai (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                     exit
                  endif
                  l = l + 1
                  NULLIFY (subset)
                end do
             endif
          END DO
          write(31)  vec_lai(:)
       end do
       close(31,status='keep')
       
       deallocate (net_data1, tile_id)
       deallocate (count_lai)
       deallocate (vec_lai) 
       deallocate (tile_lat,tile_lon)

     END SUBROUTINE grid2tile_modis6

!
! ---------------------------------------------------------------------------------------
! 
  SUBROUTINE hres_lai_no_gswp (nc_data,nr_data,rmap,gfiler,lai_name, merge)
!
! Processing GEOLAND2/MODIS LAI and creating 10-day climatological data 
!
  implicit none 
  integer, intent (in) :: nc_data,nr_data
  real, parameter :: dxy = 1.
  integer :: QSize
  type (regrid_map), intent (in) :: rmap
  character(*)  :: gfiler,lai_name
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj,ierr,nx,ny
  integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
      time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2,tindex1,pfaf1
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, target,  dimension (:,:) :: net_data1
  integer, pointer, dimension (:,:) :: QSub
  real,    pointer, dimension (:,:)    :: subset
  REAL, ALLOCATABLE, dimension (:)      :: vec_lai, count_lai,tile_lon, tile_lat &
       , x, y !, distance
  real, allocatable, target, dimension (:,:) :: lai_grid
  INTEGER ::imn,imx,jmn,jmx,mval,d1,d2,l,pix_count
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  REAL :: sf, dum,dist_save,tile_distance,minlat,maxlat,minlon,maxlon
  logical :: first_entry = .true.
  type (date_time_type) :: date_time_new,bf_lai_time,   &
       af_lai_time
  integer, intent(in), optional :: merge 

! Reading rst file
!-----------------
!   open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
!        form='unformatted',convert='little_endian')
!   allocate (tile_id    (1:nx,1:ny))         
!   
!   do j=1,ny
!      read(10)tile_id(:,j)
!   end do
!   close (10,status='keep')
!
!
! Reading number of cathment-tiles from catchment.def file
!_________________________________________________________ 
!
      fname='clsm/catchment.def' 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      allocate (tile_lon(1:maxcat)) 
      allocate (tile_lat(1:maxcat)) 
      
      do n = 1, maxcat
         read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
         tile_lon(n) = (minlon + maxlon)/2.
         tile_lat(n) = (minlat + maxlat)/2.
      end do      
      close (10,status='keep')

      fname =trim(c_data)//trim(lai_name)//'lai_clim.H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS) 
      allocate (MMDD      (0: n_tslices + 1))
      allocate (MMDD_next (0: n_tslices + 1))

      status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)
 
      if(nc_data/=i_highd .or. nr_data/=j_highd) then
         print *,'Inconsistent mapping and dimensions in  hres_lai_no_gswp  -so stopping ...'
         stop
      end if
      
      mmdd(0) = mmdd(n_tslices)
      mmdd(n_tslices + 1)= mmdd(1)
      
      mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
      mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)
      
      
      allocate(net_data1 (1:nc_10,1:nr_10))
      
      !
      ! writing MODIS/GEOLAND2 LAI data
      !
      
      if(present(merge)) then
         open (31,file='clsm/lai.'//lai_name(1:index(lai_name,'/')-1),  &
              form='unformatted',status='unknown',convert='little_endian')
      else
         open (31,file='clsm/lai.dat',  &
              form='unformatted',status='unknown',convert='little_endian')
      endif
      
      allocate (vec_lai     (maxcat))
      allocate (count_lai (1:maxcat))
 
!      allocate (vec_fill    (maxcat))
!      allocate (distance    (maxcat))
!      allocate (vec_lai_save(maxcat))
!       vec_fill = 0

      nx = nint (360./dxy)
      ny = nint (180./dxy)
      allocate (x(1:nx))
      allocate (y(1:ny))

      FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
      FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy

      allocate (lai_grid (1 : nx, 1 : ny)) 
      
      QSize = nint(dxy*nc_data/360.)
!      allocate (QSub (1:QSize,1:QSize))

      do t =0,n_tslices+1
         
         time_slice = t
         yr = 1
         yr1= 1
         if(t == 0) then
            time_slice =  n_tslices
            yr         =  1 - 1
         endif
         
         if(t >= n_tslices) then 
            yr1 = 1 + 1
            if(t ==n_tslices + 1) then
               time_slice =  1
               yr = 1 + 1
            endif
         endif
         
         read(mmdd(t),'(i2.2,i2.2)') mn,dd
         read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1
         
         ! Reading Interpolation or aggregation on to catchment-tiles
         
         vec_lai   = -9999.
         count_lai = 0.
         lai_grid  = -9999

         do jx = 1,18	
            do ix = 1,36
               write (vv,'(i2.2)')jx
               write (hh,'(i2.2)')ix 
               fname = trim(c_data)//trim(lai_name)//'lai_clim.H'//hh//'V'//vv//'.nc'
               status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
               if(status == 0) then
                  status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                  status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
                  status = NF_GET_att_INT  (ncid,4,'UNDEF',d_undef); VERIFY_(STATUS)
                  status = NF_GET_att_REAL (ncid,4,'ScaleFactor',sf); VERIFY_(STATUS)
                  status = NF_GET_VARA_INT (ncid, 4,(/1,1,time_slice/),(/nc_10,nr_10,1/),net_data1); VERIFY_(STATUS)
                  
                  do j = jLL,jLL + nr_10 -1 
                     do i = iLL, iLL + nc_10 -1 
                        if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) then
                           pix_count = rmap%ij_index(i,j)
                           if (pix_count ==0) cycle
                           if(rmap%map(pix_count)%nt > 0) then
                              do n = 1, rmap%map(pix_count)%nt
                                 if(vec_lai(rmap%map(pix_count)%tid(n)) == -9999.) vec_lai(rmap%map(pix_count)%tid(n)) = 0.                                 
                                 vec_lai(rmap%map(pix_count)%tid(n))   = vec_lai(rmap%map(pix_count)%tid(n)) + &
                                      sf*net_data1(i-iLL +1 ,j - jLL +1)*rmap%map(pix_count)%count(n)
                                 count_lai(rmap%map(pix_count)%tid(n)) = &
                                      count_lai(rmap%map(pix_count)%tid(n)) + 1.*rmap%map(pix_count)%count(n)                                     
                              end do
                           endif
                        endif
                     enddo
                  enddo

! After experimenting with few finer methods, in order to reduce the time taken by the gap filling procedure,
! creating a 0.25-degree gridded data set from finer LAI data and use it for filling the gaps seems the most practical/manageble method.
!---------------------------------------------------------------------------------------------------------------------------------------
                  do j = ceiling(1.*jLL/QSize),ceiling(1.*jLL/QSize) -1 + nr_10/QSize
                     do i = ceiling(1.*iLL/QSize),ceiling(1.*iLL/QSize) -1 + nc_10/QSize 
                         QSub  => net_data1((i-1)*QSize+2-iLL :i*QSize-iLL+1, (j-1)*QSize+2-jLL :j*QSize-jLL+1) 
                         if(maxval (QSub) > 0) lai_grid(i,j) = sf*sum(QSub, QSub>0)/(max(1,count(QSub>0)))
                     enddo
                  enddo                  
                  status = NF_CLOSE(ncid)
               endif
             end do
          end do
          
          NULLIFY (QSub)
          
          write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
          
          where (count_lai > 0.) vec_lai = vec_lai/count_lai

! Filling gaps
!-------------
          DO n =1,maxcat
             if(count_lai(n)==0.)  then 
                
                DO i = 1,nx - 1
                   if ((tile_lon(n) >= x(i)).and.(tile_lon(n) < x(i+1))) ix = i
                end do
                DO i = 1,ny -1
                   if ((tile_lat(n) >= y(i)).and.(tile_lat(n) < y(i+1))) jx = i
                end do
                
                l = 1
                do 
                  imx=ix + l
                  imn=ix - l
                  jmn=jx - l
                  jmx=jx + l
                  imn=MAX(imn,1)
                  jmn=MAX(jmn,1)
                  imx=MIN(imx,nx)
                  jmx=MIN(jmx,ny)
                  d1=imx-imn+1
                  d2=jmx-jmn+1
                  subset => lai_grid(imn: imx,jmn:jmx)

                  if(maxval(subset) > 0.) then 
                     vec_lai (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                     exit
                  endif
                  l = l + 1
                  NULLIFY (subset)
                end do

! Another Method in which search for a neighboring value while looping through nc_data*nr_data
!
!               
!               DO i = 1,nc_data - 1
!                  if ((tile_lon(n) >= x(i)).and.(tile_lon(n) < x(i+1))) ix = i
!               end do
!               DO i = 1,nr_data -1
!                  if ((tile_lat(n) >= y(i)).and.(tile_lat(n) < y(i+1))) jx = i
!               end do
!               
!               l = 1
!               do 
!                 imx=ix + l
!                 imn=ix - l
!                 jmn=jx - l
!                 jmx=jx + l
!                 imn=MAX(imn,1)
!                 jmn=MAX(jmn,1)
!                 imx=MIN(imx,nc_data)
!                 jmx=MIN(jmx,nr_data)
!                 d1=imx-imn+1
!                 d2=jmx-jmn+1
!                 ALLOCATE(subset(1:d1,1:d2))
!                 subset = -9999
!                 
!                 do j = 1,d2
!                    do i = 1,d1
!                       if (rmap(imn + i -1,jmn + j -1)%nt > 0) subset(i,j)=rmap(imn + i -1,jmn + j -1)%tid(1)
!                    end do
!                 end do
!                 
!                 mval = maxval(subset)
!                 deallocate (subset)
!
!                 if((mval > 0).and.(vec_lai_save(mval) > 0.)) then
!                    vec_lai (n) = vec_lai_save (mval)
!                    print *, count_lai(n),mval, vec_lai_save (mval)
!                    exit
!                 endif                  
!                 l = l + 1
!               end do
!                               
! The OLDEST METHOD - in which process tile space
!                if((vec_fill(n) > 0).and.(vec_lai_save(vec_fill(n)) > 0.)) then
!                   vec_lai (n) = vec_lai_save (vec_fill(n))
!                else
!
!                   distance = 1000000.
!                   where ((abs(tile_lat - tile_lat(n)) < 20.).and.                   &
!                           (abs(tile_lon - tile_lon(n)) < 10.))                      &
!                           distance =                                                &
!                           (tile_lon - tile_lon(n)) * (tile_lon - tile_lon(n)) +     &
!                           (tile_lat - tile_lat(n)) * (tile_lat - tile_lat(n))
!                   distance (n) = 1000000.
!                   k = minloc(distance,dim=1)
!
!!                   do i = 1,maxcat
!!                      if((i /= n).and.(abs(tile_lat(i) - tile_lat(n)) < 20.).and.  &
!!                           (abs(tile_lon(i) - tile_lon(n)) < 10.)) then
!!                         if(vec_lai_save(i).gt.0.) then                         
!!                            tile_distance = (tile_lon(i) - tile_lon(n)) * (tile_lon(i) - tile_lon(n)) + &
!!                                 (tile_lat(i) - tile_lat(n)) * (tile_lat(i) - tile_lat(n))
!!                            if(tile_distance < dist_save) then
!!                               k = i
!!                               dist_save = tile_distance
!!                            endif
!!                         endif
!!                      endif
!!                   enddo
!
!                   vec_lai (n) = vec_lai_save (k)
!                   vec_fill(n) = k
!                endif
             endif
          END DO
          write(31)  vec_lai(:)
       end do
       close(31,status='keep')
       
       deallocate (net_data1)
       deallocate (count_lai)
       deallocate (vec_lai) 
       deallocate (tile_lat,tile_lon)

     END SUBROUTINE hres_lai_no_gswp
!
! ---------------------------------------------------------------------------------------
! 
  SUBROUTINE hres_gswp2 (nc_data,nr_data,rmap, gfiler,lai_name,merge)
!
! Processing GSWP2 30sec LAI and grnFrac climatological data 
!
  implicit none 
  integer, intent (in) :: nc_data, nr_data 
  character(*)  :: gfiler,lai_name
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj,ierr
  integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
      time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2
  type (regrid_map), intent (in) :: rmap
  real :: dum, gyr,gmn,gdy,gyr1,gmn1,gdy1, slice1,slice2
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, target, dimension (:,:) :: &
         net_data1
  REAL, ALLOCATABLE, dimension (:) :: vec_lai, count_lai
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  REAL :: sf,minlat,maxlat,minlon,maxlon
  logical :: first_entry = .true.
  type (date_time_type) :: date_time_new,bf_lai_time,   &
       af_lai_time
  integer, intent(in), optional :: merge
  real, parameter :: dxy = 1.
  integer         :: nx, ny, QSize, pix_count
  REAL, ALLOCATABLE, dimension (:) :: x,y,tile_lon, tile_lat
  real, allocatable, target, dimension (:,:) :: data_grid
  integer, pointer, dimension (:,:) :: QSub
  INTEGER ::imn,imx,jmn,jmx,mval,d1,d2,l,tindex1,pfaf1 
  real,    pointer, dimension (:,:)    :: subset

  if(trim(lai_name) == 'lai'  ) vid = 4
  if(trim(lai_name) == 'green') vid = 5

    
  ! For Gap filling
  ! ---------------

  nx = nint (360./dxy)
  ny = nint (180./dxy)
  allocate (x(1:nx))
  allocate (y(1:ny))
  
  FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
  FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy
  
  allocate (data_grid (1 : nx, 1 : ny))   
  
  QSize = nint(dxy*nc_data/360.)
  
! Reading number of cathment-tiles from catchment.def file
! -------------------------------------------------------- 

  fname='clsm/catchment.def' 
  open (10,file=fname,status='old',action='read',form='formatted')
  read(10,*) maxcat
  allocate (tile_lon(1:maxcat)) 
  allocate (tile_lat(1:maxcat)) 
  
  do n = 1, maxcat
     read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
     tile_lon(n) = (minlon + maxlon)/2.
     tile_lat(n) = (minlat + maxlat)/2.
  end do
  
  close (10,status='keep')
  
  fname =trim(c_data)//'GSWP2_30sec_VegParam/GSWP2_VegParam_H11V13.nc'
  status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
  status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
  status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
  status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
  status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
  status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
  status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
  status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS)
  allocate (MMDD      (0: n_tslices + 1))
  allocate (MMDD_next (0: n_tslices + 1))
  
  status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
  status = NF_CLOSE(ncid); VERIFY_(STATUS)
  
  mmdd(0) = mmdd(n_tslices)
  mmdd(n_tslices + 1)= mmdd(1)
  
  mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
  mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)
  
  allocate(net_data1 (1:nc_10,1:nr_10))
  
  ! writing GSWP2 data
  ! ------------------

  if(present(merge)) then
     open (31,file='clsm/lai.gswp2',  &
          form='unformatted',status='unknown',convert='little_endian')
  else
     open (31,file='clsm/'//trim(lai_name)//'.dat',  &
          form='unformatted',status='unknown',convert='little_endian')
  endif
  
  allocate(vec_lai   (1:maxcat))
  allocate(count_lai (1:maxcat))    
  
  do t =0,n_tslices+1
     
     time_slice = t
     yr = 1
     yr1= 1
     if(t == 0) then
        time_slice =  n_tslices
        yr         =  1 - 1
     endif
     
     if(t >= n_tslices) then 
        yr1 = 1 + 1
        if(t ==n_tslices + 1) then
           time_slice =  1
           yr = 1 + 1
        endif
     endif
     
     read(mmdd(t),'(i2.2,i2.2)') mn,dd
     read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1
     
     vec_lai   = -9999.
     count_lai = 0.
     data_grid  = -9999

     do jx = 1,18	
        do ix = 1,36
           write (vv,'(i2.2)')jx
           write (hh,'(i2.2)')ix 
           fname = trim(c_data)//'GSWP2_30sec_VegParam/GSWP2_VegParam_H'//hh//'V'//vv//'.nc'
           status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
           if(status == 0) then
              status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
              status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
              status = NF_GET_att_INT  (ncid,vid,'UNDEF',d_undef); VERIFY_(STATUS)
              status = NF_GET_att_REAL (ncid,vid,'ScaleFactor',sf); VERIFY_(STATUS)
              status = NF_GET_VARA_INT (ncid, vid,(/1,1,time_slice/),(/nc_10,nr_10,1/),net_data1); VERIFY_(STATUS)
              
              do j = jLL,jLL + nr_10 -1 
                 do i = iLL, iLL + nc_10 -1 
                    if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) then
                       pix_count = rmap%ij_index(i,j)
                       if (pix_count == 0) cycle
                       if(rmap%map(pix_count)%nt > 0) then
                          do n = 1, rmap%map(pix_count)%nt
                             if(vec_lai(rmap%map(pix_count)%tid(n)) == -9999.) vec_lai(rmap%map(pix_count)%tid(n)) = 0.                                 
                             vec_lai(rmap%map(pix_count)%tid(n))   = vec_lai(rmap%map(pix_count)%tid(n)) + &
                                  sf*net_data1(i-iLL +1 ,j - jLL +1)*rmap%map(pix_count)%count(n)
                             count_lai(rmap%map(pix_count)%tid(n)) = &
                                  count_lai(rmap%map(pix_count)%tid(n)) + 1.*rmap%map(pix_count)%count(n)                                     
                          end do
                       endif
                    endif
                 enddo
              enddo

              ! After experimenting with few finer methods, in order to reduce the time taken by the gap filling procedure,
              ! creating a 1.-degree gridded data set from finer LAI data and use it for filling the gaps seems the most practical/manageble method.
              !---------------------------------------------------------------------------------------------------------------------------------------
              
              do j = ceiling(1.*jLL/QSize),ceiling(1.*jLL/QSize) -1 + nr_10/QSize
                 do i = ceiling(1.*iLL/QSize),ceiling(1.*iLL/QSize) -1 + nc_10/QSize 
                    QSub  => net_data1((i-1)*QSize+2-iLL :i*QSize-iLL+1, (j-1)*QSize+2-jLL :j*QSize-jLL+1) 
                    if(maxval (QSub) > 0) data_grid(i,j) = sf*sum(QSub, QSub>0)/(max(1,count(QSub>0)))
                 enddo
              enddo
              
              status = NF_CLOSE(ncid)
           endif
        end do
     end do
     
     write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
     where (count_lai > 0.) vec_lai = vec_lai/count_lai
 
     ! Filling gaps
     !-------------
     DO n =1,maxcat
        if(count_lai(n)==0.)  then 
           
           DO i = 1,nx - 1
              if ((tile_lon(n) >= x(i)).and.(tile_lon(n) < x(i+1))) ix = i
           end do
           DO i = 1,ny -1
              if ((tile_lat(n) >= y(i)).and.(tile_lat(n) < y(i+1))) jx = i
           end do
           
           l = 1
           do 
              imx=ix + l
              imn=ix - l
              jmn=jx - l
              jmx=jx + l
              imn=MAX(imn,1)
              jmn=MAX(jmn,1)
              imx=MIN(imx,nx)
              jmx=MIN(jmx,ny)
              d1=imx-imn+1
              d2=jmx-jmn+1
              subset => data_grid(imn: imx,jmn:jmx)
              
              if(maxval(subset) > 0.) then 
                 vec_lai (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                 exit
              endif
              l = l + 1
              NULLIFY (subset)
           end do
        endif
     end do     
     write(31)  vec_lai(:)
  end do

  close(31,status='keep')
  
  deallocate (net_data1)
  deallocate (count_lai)
  deallocate (vec_lai)

  END SUBROUTINE hres_gswp2

!----------------------------------------------------------------------  

  SUBROUTINE soil_para_hwsd (nx,ny,gfiler)

! Processing NGDC-HWSD-STATSGO merged soil properties with Woesten Soil
! Parameters and produces tau_param.dat and soil_param.dat files
 
      implicit none	    
      integer, intent (in) :: nx, ny 
      character(*)  :: gfiler
      real, dimension (:), allocatable ::           &
      	    a_sand,a_clay,a_silt,a_oc,a_bee,a_psis, &
            a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, &
            atau_2cm,btau_2cm 
      integer, dimension (100,3) :: table_map
      integer, dimension (3) :: nsoil_pcarbon 
      type (mineral_perc) :: min_percs
 
      integer :: n,maxcat,i,j,k,ktop,ncid,i_highd,j_highd,nx_adj,ny_adj
      integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,d_undef,   &
                 i1,i2,icount
      character*100 :: fname,fout
      character*10 :: string
      character*2 :: VV,HH

      logical, allocatable, dimension(:,:) :: land_pixels
      integer, allocatable, dimension (:,:) :: &
         net_data1,net_data2,net_data3,net_data4,net_data5,net_data6 ,net_data7 
      integer (kind=2) , allocatable, target, dimension (:,:) :: SOIL_HIGH,  &
          sand_top,clay_top,oc_top,sand_sub,clay_sub,oc_sub, grav_grid
      integer (kind=2), pointer, dimension (:,:) :: Raster, &  
           Raster1,Raster2,Raster3,Raster4,Raster5,Raster6
      integer (kind=4), allocatable, dimension (:) :: tileid_vec,arrayA,arrayB          
      integer (kind=2), allocatable, dimension (:) ::  &
         data_vec1, data_vec2,data_vec3, data_vec4,data_vec5, data_vec6
      REAL, ALLOCATABLE, dimension (:) :: soildepth, grav_vec,soc_vec,poc_vec  
!             ncells_top,ncells_top_pro,ncells_sub_pro                          ! ncells_* not used
      integer(kind=2) , allocatable, dimension (:) :: ss_clay,    &
             ss_sand,ss_clay_all,ss_sand_all,ss_oc_all
      REAL, ALLOCATABLE :: count_soil(:)
      integer, allocatable, target, dimension (:,:) :: tile_id
      integer, pointer :: iRaster(:,:)
      integer :: tindex, pfafindex,fac,o_cl,o_clp,fac_surf   !,vtype
      real,dimension(4) :: cFamily
      real   ,dimension(5) :: cF_lim
      logical :: first_entry = .true.
      logical :: regrid,write_debug
      INTEGER, allocatable, dimension (:) :: soil_class_top,soil_class_com
      REAL :: sf,factor,wp_wetness,fac_count,this_cond
      logical                            :: CatchParamsNC_file_exists
      REAL, ALLOCATABLE, DIMENSION (:,:) :: parms4file

      ! PEATCLSM:
      REAL, PARAMETER :: PEATMAP_THRESHOLD_1 = 0.5  ! for converting PEATMAP area fraction into peat/non-peat (on raster grid)
      REAL, PARAMETER :: PEATMAP_THRESHOLD_2 = 0.5  ! for aggregation from raster grid cells to tiles

      REAL, DIMENSION (:), POINTER       :: PMAP
      REAL, ALLOCATABLE, DIMENSION (:,:) :: PMAPR

! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1, li, ui, t_count
!
integer, dimension(:), allocatable :: low_ind, upp_ind
!
! ------------------------------------------------------------------
        
  ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
  !
  ! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION
  !
  !$ running_omp = .true.         ! conditional compilation
  !
  ! ECHO BASIC OMP VARIABLES
  !
  !$OMP PARALLEL DEFAULT(NONE) SHARED(running_omp,n_threads) 
  !
  !$OMP SINGLE
  !
  !$ n_threads = omp_get_num_threads()
  !
  !$ write (*,*) 'running_omp = ', running_omp
  !$ write (*,*)
  !$ write (*,*) 'parallel OpenMP with ', n_threads, 'threads'
  !$ write (*,*)
  !$OMP ENDSINGLE
  !
  !$OMP CRITICAL
  !$ write (*,*) 'thread ', omp_get_thread_num(), ' alive'
  !$OMP ENDCRITICAL
  !
  !$OMP BARRIER
  !
  !$OMP ENDPARALLEL

   if (first_entry) then
      nullify(iraster) ; first_entry = .false.
   endif

   ! define orgC content thresholds for orgC classes 1-4 (low, medium, high, peat)

      cF_lim(1) =  0.
      cF_lim(2) =  0.4                             ! 0.365    ! 0.3
      cF_lim(3) =  0.64                            ! 0.585    ! 4.0	   
      cF_lim(4) =  15./1.72   ! 15./1.72=8.72      ! 9.885    ! 8.5
      cF_lim(5) =  100.0

      ! define number of mineral classes in each orgC class

      nsoil_pcarbon(1) = 84 ! 84
      nsoil_pcarbon(2) = nsoil_pcarbon(1) + 84 ! 84
      nsoil_pcarbon(3) = nsoil_pcarbon(2) + 84 ! 57

      ! Read number of catchment-tiles (maxcat) from catchment.def file
      
      fname='clsm/catchment.def'
      
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      
      close (10,status='keep')
          
      ! Read tile-id raster file

      allocate(tile_id(1:nx,1:ny))
      
      fname=trim(gfiler)//'.rst'
      
      open (10,file=fname,status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,ny
         read(10)tile_id(:,j)
      end do
      
      close (10,status='keep')

      ! read soil depth data from GSWP2_soildepth_H[xx]V[yy].nc
      !
      ! get info common to all H[xx]V[yy] rectangles:
      
      fname =trim(c_data)//'SOIL-DATA/GSWP2_soildepth_H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      !status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)  ! cannot be needed here
      !status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)  ! cannot be needed here
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)
      
      ! GSWP2_soildepth_H[xx]V[yy].nc as of 29 Apr 2022:
      !
      ! Associated 43200-by-21600 global grid (1/120=0.083333 deg lat/lon):
      !
      ! N_lon_global = i_highd = 43200
      ! N_lat_global = j_highd = 21600
      !
      ! i_ind_offset_LL = iLL = 42001 
      ! j_ind_offset_LL = jLL = 19201 
      ! 
      ! Each file contains data for one rectangle of size 1200-by-1200, which is
      ! assumed to be the same for each H[xx]V[yy] rectangle
      !
      ! N_lon = nc_10 = 1200 
      ! N_lat = nr_10 = 1200

      allocate(soil_high(1:i_highd,1:j_highd))  
      allocate(net_data1 (1:nc_10,1:nr_10))

      soil_high = -9999
      do jx = 1,18
      	 do ix = 1,36
	    write (vv,'(i2.2)')jx
	    write (hh,'(i2.2)')ix 
	    fname = trim(c_data)//'SOIL-DATA/GSWP2_soildepth_H'//hh//'V'//vv//'.nc'
            status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
	    if(status == 0) then
		status = NF_GET_att_INT  (ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                status = NF_GET_att_INT  (ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
		status = NF_GET_att_INT  (ncid,4,'UNDEF',d_undef); VERIFY_(STATUS)
		status = NF_GET_att_REAL (ncid,4,'ScaleFactor',sf); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 4,(/1,1/),(/nc_10,nr_10/),net_data1); VERIFY_(STATUS)

		do j = jLL,jLL + nr_10 -1 
                   do i = iLL, iLL + nc_10 -1 
                    if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            soil_high(i,j) = net_data1(i-iLL +1 ,j - jLL +1)
                   enddo
                enddo
	        status = NF_CLOSE(ncid)
	    endif
	 end do 
      end do

      deallocate (net_data1)

      ! Regridding 
      
      ! *EASE* grid bcs use mask file GEOS5_10arcsec_mask.nc or GEOS5_10arcsec_mask_freshwater-lakes.nc,
      ! and the make_bcs script assigns NX=43200, NY=21600, which are passed into the present subroutine
      ! via command-line arguments of mkCatchParam.x.  That is, should have regrid=.false. for *EASE*
      ! grid tile space.

      nx_adj = nx
      ny_adj = ny

      regrid = nx/=i_highd .or. ny/=j_highd

      if(regrid) then
          if(nx > i_highd) then 
               allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(soil_high,raster)	
	       iRaster => tile_id
	       if(ny < j_highd) then
	          print *,'nx > i_highd and ny < j_highd'
		  stop 
	       endif
	  else
               if( .not.associated(iraster) ) then 
                   allocate(iraster(i_highd,j_highd),stat=STATUS); VERIFY_(STATUS)
               endif
               call RegridRaster(tile_id,iraster)	
               raster  => soil_high
	       nx_adj = i_highd
               ny_adj = j_highd

	       if(ny > j_highd) then
	          print *,'nx < i_highd and ny > j_highd'
		  stop 
	       endif
          endif	    	              
      else
          raster  => soil_high
	  iRaster => tile_id
      end if
      
      ! Interpolate/aggregate soil depth from raster grid to catchment-tiles
      
      allocate(soildepth(1:maxcat))
      allocate(count_soil(1:maxcat))  
      
      soildepth  = 0.      ! 1-d tile space
      count_soil = 0.      ! 1-d tile space
 
      do j=1,ny_adj
         do i=1,nx_adj
            if((iRaster(i,j).gt.0).and.(iRaster(i,j).le.maxcat)) then
               if ((raster(i,j).gt.0)) then
                  soildepth(iRaster(i,j)) = &
                       soildepth(iRaster(i,j)) + sf*raster(i,j)  ! integer "raster" --> real "soildepth"
                  count_soil(iRaster(i,j)) = &
                       count_soil(iRaster(i,j)) + 1. 
               endif
            endif
         end do
      end do

      DO n =1,maxcat
	 	if(count_soil(n)/=0.) soildepth(n)=soildepth(n)/count_soil(n)	
	          soildepth(n) = max(soildepth(n),SOILDEPTH_MIN_HWSD)
!                  soildepth(n) = soildepth(n) + 2000.
!                  soildepth(n) = min(soildepth(n),8000.) 
       END DO

      deallocate (SOIL_HIGH)
      !deallocate (count_soil)                 ! do not deallocate, needed again shortly
      NULLIFY(Raster)
      
      ! ---------------------------------------------------------------------------------
      !
      ! Read NGDC-HWSD-STATSGO merged soil texture from SoilProperties_H[xx]V[yy].nc'
      !
      ! get info common to all H[xx]V[yy] rectangles (could in theory differ from that
      !   of soildepth data read above but is the same as of 29 Apr 2022).

      fname =trim(c_data)//'SOIL-DATA/SoilProperties_H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      !status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)  ! cannot be needed here
      !status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)  ! cannot be needed here
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS) 
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS) 
      status = NF_CLOSE(ncid)      
      
      ! SoilProperties_H[xx]V[yy].nc as of 29 Apr 2022:
      !
      ! Associated 43200-by-21600 global grid (1/120=0.083333 deg lat/lon):
      !
      ! N_lon_global = i_highd = 43200
      ! N_lat_global = j_highd = 21600
      !
      ! i_ind_offset_LL = iLL = 42001 
      ! j_ind_offset_LL = jLL = 19201 
      ! 
      ! Each file contains soil texture data for one rectangle of size 1200-by-1200, which is
      ! assumed to be the same for each H[xx]V[yy] rectangle
      !
      ! N_lon = nc_10 = 1200 
      ! N_lat = nr_10 = 1200

      !regrid = nx/=i_highd .or. ny/=j_highd       ! not needed here, done below
      
      allocate(net_data1 (1:nc_10,1:nr_10))
      allocate(net_data2 (1:nc_10,1:nr_10))
      allocate(net_data3 (1:nc_10,1:nr_10))
      allocate(net_data4 (1:nc_10,1:nr_10))
      allocate(net_data5 (1:nc_10,1:nr_10))
      allocate(net_data6 (1:nc_10,1:nr_10))
      allocate(net_data7 (1:nc_10,1:nr_10))

      allocate(sand_top (1:i_highd,1:j_highd))  
      allocate(clay_top (1:i_highd,1:j_highd))  
      allocate(oc_top   (1:i_highd,1:j_highd))  
      allocate(sand_sub (1:i_highd,1:j_highd))  
      allocate(clay_sub (1:i_highd,1:j_highd))  
      allocate(oc_sub   (1:i_highd,1:j_highd))  
      allocate(grav_grid(1:i_highd,1:j_highd))   

      sand_top = -9999   ! integer*2
      clay_top = -9999   ! integer*2
      oc_top   = -9999   ! integer*2
      sand_sub = -9999   ! integer*2
      clay_sub = -9999   ! integer*2
      oc_sub   = -9999   ! integer*2
      grav_grid= -9999   ! integer*2

      do jx = 1,18
      	 do ix = 1,36
	    write (vv,'(i2.2)')jx
	    write (hh,'(i2.2)')ix 
	    fname = trim(c_data)//'SOIL-DATA/SoilProperties_H'//hh//'V'//vv//'.nc'
            status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
	    if(status == 0) then
		status = NF_GET_att_INT  (ncid, NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                status = NF_GET_att_INT  (ncid, NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
                ! assume UNDEF and ScaleFactor (sf) are the same for *all* variables read below
                ! (ok for SoilProperties_H[xx]V[yy].nc as of 29 Apr 2022).
		status = NF_GET_att_INT  (ncid, 4,'UNDEF',d_undef); VERIFY_(STATUS) 
		status = NF_GET_att_REAL (ncid, 4,'ScaleFactor',sf); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 4,(/1,1/),(/nc_10,nr_10/),net_data1); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 5,(/1,1/),(/nc_10,nr_10/),net_data2); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 6,(/1,1/),(/nc_10,nr_10/),net_data3); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 7,(/1,1/),(/nc_10,nr_10/),net_data4); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 8,(/1,1/),(/nc_10,nr_10/),net_data5); VERIFY_(STATUS)
		status = NF_GET_VARA_INT (ncid, 9,(/1,1/),(/nc_10,nr_10/),net_data6); VERIFY_(STATUS)
                status = NF_GET_VARA_INT (ncid,10,(/1,1/),(/nc_10,nr_10/),net_data7); VERIFY_(STATUS)
		do j = jLL,jLL + nr_10 -1 
                   do i = iLL, iLL + nc_10 -1 
                    if(net_data1(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            clay_top(i,j) = net_data1(i-iLL +1 ,j - jLL +1)
                    if(net_data2(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            sand_top(i,j) = net_data2(i-iLL +1 ,j - jLL +1)
                    if(net_data3(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            oc_top  (i,j) = net_data3(i-iLL +1 ,j - jLL +1)
                    if(net_data4(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            clay_sub(i,j) = net_data4(i-iLL +1 ,j - jLL +1)
                    if(net_data5(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            sand_sub(i,j) = net_data5(i-iLL +1 ,j - jLL +1)
                    if(net_data6(i-iLL +1 ,j - jLL +1) /= d_undef) &
		            oc_sub  (i,j) = net_data6(i-iLL +1 ,j - jLL +1)
                    if(net_data7(i-iLL +1 ,j - jLL +1) /= d_undef) &
		           grav_grid(i,j) = net_data7(i-iLL +1 ,j - jLL +1)
                    enddo
                enddo
	    status = NF_CLOSE(ncid)
	    endif
	 end do 
       end do

      deallocate (net_data1)
      deallocate (net_data2)
      deallocate (net_data3)
      deallocate (net_data4)
      deallocate (net_data5)
      deallocate (net_data6)
      deallocate (net_data7)
      
       ! ----------------------------------------------------------------------------

      if(use_PEATMAP) then 
         print *, 'PEATMAP_THRESHOLD_1 : ', PEATMAP_THRESHOLD_1
         allocate(pmapr (1:i_highd,1:j_highd))
         status  = NF_OPEN ('data/CATCH/PEATMAP_mask.nc4', NF_NOWRITE, ncid)
         status  = NF_GET_VARA_REAL (ncid,NC_VarID(NCID,'PEATMAP'), (/1,1/),(/i_highd, j_highd/), pmapr) ; VERIFY_(STATUS)      

         ! move HWSD sub-surface peat to peat-rich mineral Group 3 because merged surface peat defines sub-surface peat

         where (oc_sub*sf >= cF_lim(4))
            oc_sub = NINT(8./sf)
         endwhere

         ! Hybridize: add OC 1km PEATMAP pixels to HWSD oc_top

         where (pmapr >= PEATMAP_THRESHOLD_1)
           oc_top = NINT(33.0/sf)
         endwhere  

         deallocate (pmapr)
         status = NF_CLOSE(ncid)
      endif
      
      ! ----------------------------------------------------------------------------

      ! Regridding 
      
      ! *EASE* grid bcs use mask file GEOS5_10arcsec_mask.nc or GEOS5_10arcsec_mask_freshwater-lakes.nc,
      ! and the make_bcs script assigns NX=43200, NY=21600, which are passed into the present subroutine
      ! via command-line arguments of mkCatchParam.x.  That is, should have regrid=.false. for *EASE*
      ! grid tile space.

      nx_adj = nx
      ny_adj = ny

      regrid = nx/=i_highd .or. ny/=j_highd       

      if(regrid) then
          if(nx > i_highd) then 
               allocate(raster1(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(clay_top,raster1)	

               allocate(raster2(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(sand_top,raster2)	

               allocate(raster3(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(oc_top,  raster3)	

               allocate(raster4(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(clay_sub,raster4)	

               allocate(raster5(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(sand_sub,raster5)	

               allocate(raster6(nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(oc_sub,  raster6)	

               allocate(raster (nx,ny),stat=STATUS); VERIFY_(STATUS)
               call RegridRaster2(grav_grid,raster)

	       iRaster => tile_id

	       if(ny < j_highd) then
	          print *,'nx > i_highd and ny < j_highd'
		  stop 
	       endif
	  else
	       nx_adj = i_highd
               ny_adj = j_highd
               if( .not.associated(iraster) ) then
                  allocate(iRaster(i_highd,j_highd),stat=STATUS); VERIFY_(STATUS)
               endif
               call RegridRaster(tile_id,iRaster)	

               raster1 => clay_top
               raster2 => sand_top
               raster3 => oc_top
               raster4 => clay_sub
               raster5 => sand_sub
               raster6 => oc_sub
               raster  => grav_grid

	       if(ny > j_highd) then
	          print *,'nx < i_highd and ny > j_highd'
		  stop 
	       endif
          endif	    	              
      else
	  iRaster => tile_id
          raster1 => clay_top
          raster2 => sand_top
          raster3 => oc_top
          raster4 => clay_sub
          raster5 => sand_sub
          raster6 => oc_sub
          raster  => grav_grid
      end if
      
      ! ----------------------------------------------------------------------------

      ! compute peat fraction on tile for CLM45+ (for fires?)
      
      allocate(pmap  (1:maxcat))
      !allocate(count_soil(1:maxcat))            ! already allocated above

      pmap       = 0.    ! 1-d tile space; peat fraction in tile based on oc_top
      count_soil = 0.    ! 1-d tile space 
      
      do j=1,ny_adj
         do i=1,nx_adj
            if((iRaster(i,j).gt.0).and.(iRaster(i,j).le.maxcat)) then
               count_soil(iRaster(i,j)) = count_soil(iRaster(i,j)) + 1. 
               if (raster3(i,j)*sf >= cF_lim(4)) then
                  pmap (iRaster(i,j)) = pmap(iRaster(i,j)) + 1                  
               endif
            endif
         end do
      end do

      where (count_soil > 0) pmap = pmap /count_soil
      
      !deallocate (count_soil)                 ! do not deallocate, needed again shortly

      ! ----------------------------------------------------------------------------
      
      ! get number of "land" pixels (i1) on raster grid

      allocate(land_pixels(1:size(iRaster,1),1:size(iRaster,2)))
      land_pixels = (iRaster >=1).and.(iRaster<=maxcat)
      i1 = count(land_pixels)   
      deallocate (land_pixels) 

      ! allocate 1-d arrays for all "land" pixels on raster grid

      allocate (tileid_vec(1:i1))
      allocate (data_vec1 (1:i1))
      allocate (data_vec2 (1:i1))
      allocate (data_vec3 (1:i1))
      allocate (data_vec4 (1:i1))
      allocate (data_vec5 (1:i1))
      allocate (data_vec6 (1:i1))

      ! allocate 1-d arrays for all "land" tiles

      allocate (grav_vec  (1:maxcat))
      allocate (soc_vec   (1:maxcat))
      allocate (poc_vec   (1:maxcat))
      !allocate (ncells_top  (1:maxcat))            ! ncells_* not used
      !allocate (ncells_top_pro  (1:maxcat))        ! ncells_* not used
      !allocate (ncells_sub_pro  (1:maxcat))        ! ncells_* not used
      !allocate(count_soil(1:maxcat))               
  
      count_soil = 0.
      grav_vec   = 0.
      soc_vec    = 0.            ! soil orgC (top layer 0-30)
      poc_vec    = 0.            ! soil orgC (profile layer 0-100)

      !ncells_top = 0.           ! ncells_* not used
      !ncells_top_pro = 0.       ! ncells_* not used
      !ncells_sub_pro = 0.       ! ncells_* not used

      n =1
      do j=1,ny_adj
         do i=1,nx_adj
            if((iRaster(i,j).ge.1).and.(iRaster(i,j).le.maxcat)) then

               ! map from 2-d raster array to 1-d raster vec

	       tileid_vec (n) =  iRaster(i,j)         ! iRaster => tile_id       int*4
	       data_vec1  (n) =  Raster1(i,j)         ! raster1 => clay_top      int*2
	       data_vec2  (n) =  Raster2(i,j)         ! raster2 => sand_top      int*2
	       data_vec3  (n) =  Raster3(i,j)         ! raster3 => oc_top        int*2
	       data_vec4  (n) =  Raster4(i,j)         ! raster4 => clay_sub      int*2
	       data_vec5  (n) =  Raster5(i,j)         ! raster5 => sand_sub      int*2
	       data_vec6  (n) =  Raster6(i,j)         ! raster6 => oc_sub        int*2

               ! BUG???  It is unclear why here grav_vec is filled in the order of "tile_id"
               ! while data_vec[x] is filled in the order of the long/lat grid.
               ! Not sure if grav_vec is processed correctly below!
               ! -reichle, 29 Apr 2022

               if ((raster(i,j).gt.0)) then 
                  grav_vec(iRaster(i,j)) = &
                       grav_vec(iRaster(i,j)) + sf*raster(i,j)  ! raster  => grav_grid    int*2
                  count_soil(iRaster(i,j)) = &
                       count_soil(iRaster(i,j)) + 1. 
               endif
	       n = n + 1
	     endif
	  end do
      end do    

      DO n =1,maxcat
	 	if(count_soil(n)/=0.) grav_vec(n)=grav_vec(n)/count_soil(n)	
      END DO

      deallocate (count_soil)
      NULLIFY(Raster,Raster1,Raster2,Raster3,Raster4,Raster5,Raster6)
      deallocate (clay_top,sand_top,oc_top,clay_sub,sand_sub,oc_sub,grav_grid) 
      deallocate (tile_id)

      ! sort 1-d land pixels vectors according to tile_id

      allocate (arrayA    (1:i1))             ! 1-d land pixels on raster grid
      allocate (arrayB    (1:i1))             ! 1-d land pixels on raster grid

      arrayA = tileid_vec
      arrayB = data_vec1
      call MAPL_Sort (arrayA, arrayB)
      data_vec1 = arrayB

      arrayA = tileid_vec
      arrayB = data_vec2
      call MAPL_Sort (arrayA, arrayB)
      data_vec2 = arrayB

      arrayA = tileid_vec
      arrayB = data_vec3
      call MAPL_Sort (arrayA, arrayB)
      data_vec3 = arrayB

      arrayA = tileid_vec
      arrayB = data_vec4
      call MAPL_Sort (arrayA, arrayB)
      data_vec4 = arrayB

      arrayA = tileid_vec
      arrayB = data_vec5
      call MAPL_Sort (arrayA, arrayB)
      data_vec5 = arrayB

      arrayA = tileid_vec
      arrayB = data_vec6
      call MAPL_Sort (arrayA, arrayB)
      data_vec6 = arrayB

      tileid_vec= arrayA

      deallocate (arrayA, arrayB)

      ! -------------------------------------------------------------------- 
      !
      ! Read Woesten soil parameters and CLSM tau parameters for soil classes (1:253)

	allocate(a_sand  (1:n_SoilClasses))
	allocate(a_clay  (1:n_SoilClasses))
	allocate(a_silt  (1:n_SoilClasses))
	allocate(a_oc    (1:n_SoilClasses))
	allocate(a_bee   (1:n_SoilClasses))
	allocate(a_psis  (1:n_SoilClasses))
	allocate(a_poros (1:n_SoilClasses))
	allocate(a_wp    (1:n_SoilClasses))
	allocate(a_aksat (1:n_SoilClasses))
	allocate(atau    (1:n_SoilClasses))
	allocate(btau    (1:n_SoilClasses))
	allocate(atau_2cm(1:n_SoilClasses))
	allocate(btau_2cm(1:n_SoilClasses))
        allocate(a_wpsurf(1:n_SoilClasses))
        allocate(a_porosurf(1:n_SoilClasses))

      ! SoilClasses-SoilHyd-TauParam.dat and SoilClasses-SoilHyd-TauParam.peatmap differ
      ! only in the parameters for the peat class #253.  The file *.peatmap contains
      ! the PEATCLSM parameters from Table 2 of Bechtold et al. 2019 (doi:10.1029/2018MS001574).
      !
      ! Note: K_s = COND*exp(-zks*gnu)   ==> with zks=2 and gnu=1, K_s = 0.135335*COND
      !
      !         K_s      COND     [m/s]
      ! NLv4  7.86e-7   5.81e-6
      ! NLv5  3.79e-6   2.80e-5   <== note *typo* in Table 2 of Bechtold et al. 2019, which erroneously lists K_s=2.8e-5

      if(use_PEATMAP) then 
         fname = trim(c_data)//'SoilClasses-SoilHyd-TauParam.peatmap'
      else
         fname = trim(c_data)//'SoilClasses-SoilHyd-TauParam.dat'
      endif

      table_map = 0                      ! 100-by-3 look-up table

      open (11, file=trim(fname), form='formatted',status='old', &
           action = 'read')
      read (11,'(a)')fout        ! read header line

      do n =1,n_SoilClasses

      	 read (11,'(4f7.3,4f8.4,e13.5,2f12.7,2f8.4,4f12.7)')a_sand(n),a_clay(n),a_silt(n),a_oc(n),a_bee(n),a_psis(n), &
              a_poros(n),a_wp(n),a_aksat(n),atau(n),btau(n),a_wpsurf(n),a_porosurf(n),atau_2cm(n),btau_2cm(n)

         ! assemble scalar structure that holds mineral percentages of soil class n

	 min_percs%clay_perc = a_clay(n)
	 min_percs%silt_perc = a_silt(n)
	 min_percs%sand_perc = a_sand(n)
  
         ! "soil_class" is an integer function (see rmTinyCatchParam.F90) that assigns 
         !    an integer (mineral) soil class [1-100] for a given mineral percentage triplet

         ! "table_map"  is a 2-d array (100-by-3) that maps between overall soil class (1:252) and 
         !   (mineral_class 1:84, orgC_class).   "table_map" has no entry for the peat class #253.

	 if( n <= nsoil_pcarbon(1))                              table_map(soil_class (min_percs),1) = n  
	 if((n >  nsoil_pcarbon(1)).and.(n <= nsoil_pcarbon(2))) table_map(soil_class (min_percs),2) = n  
         if((n >  nsoil_pcarbon(2)).and.(n <= nsoil_pcarbon(3))) table_map(soil_class (min_percs),3) = n 

      end do   ! n=1,n_SoilClasses

      close (11,status='keep') 

      ! ------------------------------------------------------------
      !
      !  When Woesten soil parameters are not available for a particular soil class,
      !  as defined by "tiny" triangles in HWSD soil triangle, Woesten soil
      !  parameters from the nearest available "tiny" triangle will be substituted.
      !  For "tiny" triangles, see Fig 1b of De Lannoy et al. 2014 (doi:10.1002/2014MS000330).      

      do n =1,10
	  do k=1,n*2 -1

     	     min_percs%clay_perc = 100. -((n-1)*10 + 5)
	     min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
	     min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc

	     i = soil_class (min_percs)

	     if(table_map (i,1)== 0) then   
	       j = center_pix (a_clay(1:nsoil_pcarbon(1)),a_sand(1:nsoil_pcarbon(1)),                       &
	           min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.) 

	 	   min_percs%clay_perc = a_clay(j)
	 	   min_percs%silt_perc = a_silt(j)
	 	   min_percs%sand_perc = a_sand(j)		   
	       	   table_map (i,1)= table_map (soil_class (min_percs),1)
	     endif

     	     min_percs%clay_perc = 100. -((n-1)*10 + 5)
	     min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
	     min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc

	     if(table_map (i,2)== 0) then   
	       j = center_pix(a_clay(nsoil_pcarbon(1)+1 : nsoil_pcarbon(2)),         &
	                      a_sand(nsoil_pcarbon(1)+1 : nsoil_pcarbon(2)),         &
	           min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.) 
	 	   min_percs%clay_perc = a_clay(j + nsoil_pcarbon(1))
	 	   min_percs%silt_perc = a_silt(j + nsoil_pcarbon(1))
	 	   min_percs%sand_perc = a_sand(j + nsoil_pcarbon(1))		   
	       	   table_map (i,2)= table_map (soil_class (min_percs),2)	         
	     endif

     	     min_percs%clay_perc = 100. -((n-1)*10 + 5)
	     min_percs%sand_perc = 100. -  min_percs%clay_perc -2.-(k-1)*5.
	     min_percs%silt_perc = 100. -  min_percs%clay_perc - min_percs%sand_perc

	     if(table_map (i,3)== 0) then   
	       j = center_pix (a_clay(nsoil_pcarbon(2)+1 : nsoil_pcarbon(3)),         &
                               a_sand(nsoil_pcarbon(2)+1 : nsoil_pcarbon(3)),         &
	           min_percs%clay_perc,min_percs%sand_perc,min_percs%silt_perc,.true.) 
	 	   min_percs%clay_perc = a_clay(j + nsoil_pcarbon(2))
	 	   min_percs%silt_perc = a_silt(j + nsoil_pcarbon(2))
	 	   min_percs%sand_perc = a_sand(j + nsoil_pcarbon(2))		   
	       	   table_map (i,3)= table_map (soil_class (min_percs),3)	         	         
	     endif
	  end do         	  
      end do
!
! Now deriving soil types based on NGDC-HWSD-STATSGO merged soil property maps
!
      allocate (soil_class_top (1:maxcat))
      allocate (soil_class_com (1:maxcat))
      soil_class_top =-9999
      soil_class_com =-9999
      
      allocate(low_ind(n_threads))
      allocate(upp_ind(n_threads))
      low_ind(1)         = 1
      upp_ind(n_threads) = maxcat

      if (running_omp)  then
       do i=1,n_threads-1  
         upp_ind(i)   = low_ind(i) + (maxcat/n_threads) - 1 
         low_ind(i+1) = upp_ind(i) + 1
      end do 
     end if

!$OMP PARALLELDO DEFAULT(NONE)                          &
!$OMP SHARED( n_threads, low_ind, upp_ind, tileid_vec,  &
!$OMP         sf,data_vec1,data_vec2,data_vec3,         &
!$OMP         data_vec4,data_vec5,data_vec6,cF_lim,     &
!$OMP         table_map,soil_class_top,soil_class_com,  &
!$OMP         soc_vec,poc_vec,use_PEATMAP)              &
!ncells_* not used !$OMP         soc_vec,poc_vec,ncells_top,ncells_top_pro,&
!ncells_* not used !$OMP         ncells_sub_pro,use_PEATMAP)               &
!$OMP PRIVATE(n,i,j,k,icount,t_count,i1,i2,ss_clay,     &
!$OMP         ss_sand,ss_clay_all,ss_sand_all,          &
!$OMP         ss_oc_all,cFamily,factor,o_cl,o_clp,ktop, &
!$OMP         min_percs, fac_count, write_debug)

    ! loop through tiles (split into two loops for OpenMP)

    DO t_count = 1,n_threads
      DO n = low_ind(t_count),upp_ind(t_count)

	write_debug = .false.

!	if (n==171010)  write_debug = .true.

        ! initialize "icount" when starting loop through n at low_ind(t_count)
        ! recall: tileid_vec is a 1-d vector that covers all land pixels on the raster grid that
        !         contains the (sorted) tile IDs, with matching parameter vectors data_vec[x]

        if(n==low_ind(t_count)) then
             icount = 1
             ! Not sure what the following loops do.  Why not check backwards from low_ind(t_count)??
             do k=1,low_ind(t_count) - 1
                do while (tileid_vec(icount)== k)   
                   icount = icount + 1
                end do
             end do
        endif

        ! ------------------------------------------------------------------
        !
        ! determine the land raster grid cells i1:i2 that make up tile n
        
        ! NOTE change in meaning of "i1":
        !
        ! before: i1 = total no. of land pixels on the raster grid
        ! now:    i1 = starting index of land raster grid cells (within 1-d vector) that make up tile n (?)
             
        i1 = icount  
        
	loop: do while (tileid_vec(icount)== n)
	   if(icount <= size(tileid_vec,1)) icount = icount + 1
           if(icount > size(tileid_vec,1)) exit loop
	end do loop 

	i2 = icount -1
 	i = i2 - i1 + 1                ! number of land raster grid cells that make up tile n (?)


        ! -------------------------------------------------------------------
        ! 
        ! prep data
  
        allocate(ss_clay    (1:2*i))   ! for top layer (0-30)   -- why allocate 1:2*i and not 1:i??
        allocate(ss_sand    (1:2*i))   ! for top layer (0-30)   -- why allocate 1:2*i and not 1:i??

        allocate(ss_clay_all(1:2*i))   ! for top (0-30) and sub (30-100) layers
        allocate(ss_sand_all(1:2*i))   ! for top (0-30) and sub (30-100) layers
        allocate(ss_oc_all  (1:2*i))   ! for top (0-30) and sub (30-100) layers
	  
        ss_clay    = 0    ! int*2  -- why only clay and sand for top layer and not orgC ??  
        ss_sand    = 0	  ! int*2
        
        ss_clay_all= 0    ! int*2
        ss_sand_all= 0    ! int*2
        ss_oc_all  = 0    ! int*2

	ss_clay_all (1:i)     = data_vec1(i1:i2)  ! put top layer info into first i elements (1:i)
	ss_sand_all (1:i)     = data_vec2(i1:i2)
	ss_oc_all   (1:i)     = data_vec3(i1:i2)	

	ss_clay_all (1+i:2*i) = data_vec4(i1:i2)  ! put sub layer info into next i elements (i+1:2*i)
	ss_sand_all (1+i:2*i) = data_vec5(i1:i2)
	ss_oc_all   (1+i:2*i) = data_vec6(i1:i2)  ! <-- oc_sub	


        ! -----------------------------------------------------------------------
        !
        ! determine aggregate/dominant orgC *top* layer soil class ("o_cl") of tile n

	cFamily = 0.
!!        factor  = 1.

	do j=1,i
           if(j <= i) factor = 1.
	   if((ss_oc_all(j)*sf >=  cF_lim(1)).and. (ss_oc_all(j)*sf < cF_lim(2))) cFamily(1) = cFamily(1) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(2)).and. (ss_oc_all(j)*sf < cF_lim(3))) cFamily(2) = cFamily(2) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(3)).and. (ss_oc_all(j)*sf < cF_lim(4))) cFamily(3) = cFamily(3) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(4) ))                                   cFamily(4) = cFamily(4) + factor
	end do

	if (sum(cFamily) == 0.) o_cl  = 1    ! default is o_cl=1 (if somehow no grid cell has top-layer orgC >=0.)

!!        if (.not. use_PEATMAP) then
           
           ! assign dominant *top* layer org soil class (even if only a minority of the contributing 
           !   raster grid cells is peat)

           if (sum(cFamily)  > 0.) o_cl  = maxloc(cFamily, dim = 1)

!!        else

        if (use_PEATMAP) then
           
           ! PEATMAP: tile has *top* layer peat class only if more than 50% of the contributing 
           !   raster grid cells are peat (may loose some peat tiles w.r.t. non-PEATMAP bcs version)
           
           if (cFamily(4)/real(i) > PEATMAP_THRESHOLD_2) then 
              o_cl  = 4
           else
              if (sum(cFamily(1:3)) > 0.) o_cl  = maxloc(cFamily(1:3), dim = 1)  ! o_cl = 1, 2, or 3
           endif

        endif           
        

        ! determine aggregate/dominant orgC *profile* (0-100) soil class ("o_clp") of tile n,
        ! weight factor=1. for top (0-30) layer and weight factor=2.33 for sub (30-100) layer
        
	cFamily = 0.

	do j=1,2*i
	   if(j <= i) factor = 1.
	   if(j  > i) factor = 2.33
	   if((ss_oc_all(j)*sf >=  cF_lim(1)).and. (ss_oc_all(j)*sf < cF_lim(2))) cFamily(1) = cFamily(1) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(2)).and. (ss_oc_all(j)*sf < cF_lim(3))) cFamily(2) = cFamily(2) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(3)).and. (ss_oc_all(j)*sf < cF_lim(4))) cFamily(3) = cFamily(3) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(4) ))                                   cFamily(4) = cFamily(4) + factor
	end do
 
        ! NOTE: For PEATMAP, oc_sub was cut back to 8./sf above:
        !       "! move HWSD sub-surface peat to peat-rich mineral Group 3 because merged surface peat defines sub-surface peat"
        !       "where (oc_sub*sf >= cF_lim(4))                                                                                "
        !       "    oc_sub = NINT(8./sf)                                                                                      "
        !       "endwhere                                                                                                      "
        !       For PEATMAP, the sub-layer weight of 2.33 should only count towards cFamily(1:3), and in most cases the 
        !       maxloc statement below should therefore result in o_clp = 1, 2, or 3 only.  However, if the top-layer orgC
        !       is peat for most contributing raster grid cells and the sub-layer orgC values are relatively evenly spread 
        !       over orgC classes 1, 2, and 3, then maxloc(cFamily) can result in o_clp=4.

	if (sum(cFamily) == 0.) o_clp = 1
	if (sum(cFamily)  > 0.) o_clp = maxloc(cFamily, dim = 1)        

        ! ----------------------------------------------------------------------------------------
        ! 
        ! Determine *top* layer mineral/organic soil class of tile n 

        if(o_cl == 4) then 

           ! Top-layer soil class of tile n is peat. 
           ! Compute average top-layer orgC (only across raster grid cells whose top layer is peat).

           soil_class_top(n) = n_SoilClasses
	   ktop = 0
	   do j=1,i
             ! avg only across contributing raster grid cells that are peat
	     if(ss_oc_all(j)*sf >= cF_lim(4)) then           
	        soc_vec (n) = soc_vec(n) + ss_oc_all(j)*sf
		ktop = ktop + 1
             endif
	   end do
	   if(ktop.ne.0) soc_vec (n)   = soc_vec(n)/ktop
	   !ncells_top(n) = 100.*float(ktop)/float(i)            ! ncells_* not used

        else 
                      
           ! Top-layer soil class of tile n is mineral.
           ! Compute average top-layer orgC (only across raster grid cells within same orgC class)
           ! and collect all clay/sand pairs of raster grid cells within same orgC class. 

            !k = 1        !cleanup k counter
	    !ktop = 1     !cleanup k counter
	    ktop = 0      !cleanup k counter

	    do j=1,i      ! loop only through top-layer elements of ss_*_all

                ! avg only across contributing raster grid cells with orgC class as that assigned to tile n 
	        if((ss_oc_all(j)*sf >= cF_lim(o_cl)).and.(ss_oc_all(j)*sf < cF_lim(o_cl + 1))) then 

                   if((ss_clay_all(j)*sf >= 0.).and.(ss_sand_all(j)*sf >= 0.)) then    ! avoiding no-data-values

                           ktop = ktop + 1      !cleanup k counter
    			   ss_clay (ktop) = ss_clay_all(j)                    
			   ss_sand (ktop) = ss_sand_all(j)                    

                           ! adjust clay and sand content if outside joint physical bounds
                           if((ss_clay (ktop) + ss_sand (ktop)) > 9999) then  ! note: 9999 = 99.99%  (scale factor = 0.01)
                              if(ss_clay (ktop) >= ss_sand (ktop)) then
                                 ss_sand (ktop) = 10000 - ss_clay (ktop)
                               else
                                 ss_clay (ktop) = 10000 - ss_sand (ktop)
                               endif
                            endif
			   soc_vec (n) = soc_vec(n) + ss_oc_all(j)*sf    ! sum up top-layer orgC
                           !k = k + 1            !cleanup k counter
                           !ktop = ktop + 1      !cleanup k counter
		   endif
                endif	
	    end do
	    
	    !k = k - 1           !cleanup k counter
	    !ktop = ktop -1      !cleanup k counter

	    if(ktop.ne.0) soc_vec (n) = soc_vec(n)/ktop     ! normalize top-layer orgC

	    !ncells_top(n) = 100.*float(ktop)/float(i)            ! ncells_* not used

            ! debugging output
            if (write_debug) write(80+n,*)ktop,o_cl
            if(ktop > 0) then 
               if (write_debug) write (80+n,*)ss_clay(1:ktop)
               if (write_debug) write (80+n,*)ss_sand(1:ktop)
            endif

            ! Determine the raster grid cell j that has (top-layer) clay/sand content closest
            ! to the average (top-layer) clay/sand across all raster grid cells within the 
            ! dominant orgC class.

  	    j = center_pix_int0(sf, ktop,ktop, ss_clay(1:ktop),ss_sand(1:ktop))

            ! Assign soil class of raster grid cell j to tile n

            if(j >=1) then 
               min_percs%clay_perc = ss_clay(j)*sf
               min_percs%sand_perc = ss_sand(j)*sf
               min_percs%silt_perc = 100. - ss_clay(j)*sf - ss_sand(j)*sf
               soil_class_top (n) = table_map(soil_class (min_percs),o_cl)   
            endif
            
            ! debugging output
            if (write_debug) write(80+n,*)j

        endif

        ! debugging output
        if (write_debug) write(80+n,*)soil_class_top (n) 

        ! -------------------------------------------------------------------------------
        ! 
        ! determine aggregate sand/clay/orgC for *profile* layer of tile n 

        if(o_clp == 4) then 

           ! Profile-layer soil class of tile n is peat. 
           ! Compute average profile-layer orgC (only across raster grid cells and layers that are peat)
           
           soil_class_com(n) = n_SoilClasses
	   fac_count = 0.
	   k =0
	   ktop =0
	   do j=1,2*i
	     if(ss_oc_all(j)*sf >= cF_lim(4)) then
                 if(j <= i) factor = 1.                ! top layer contribution  1   <= j <=i
	         if(j  > i) factor = 2.33              ! sub layer contribution  i+1 <= j <=2*i
		 if(j  > i) k = k + 1                  ! sub layer counter
		 if(j <= i) ktop = ktop + 1            ! top layer counter
	         poc_vec (n) = poc_vec(n) + ss_oc_all(j)*sf*factor     ! weighted sum of orgC
		 fac_count = fac_count + factor                        ! sum of weights
             endif
	   end do
	   if(fac_count.ne.0) poc_vec (n) = poc_vec (n)/fac_count   ! normalize
           !ncells_sub_pro(n) = 100.*float(k)/float(i)              ! ncells_* not used
           !ncells_top_pro(n) = 100.*float(ktop)/float(i)           ! ncells_* not used
        else

           ! Profile-layer soil class of tile n is mineral.
           ! Compute average profile-layer orgC (only across raster grid cells within same orgC class)
           ! and collect all clay/sand pairs of raster grid cells within same orgC class. 

            !k = 1        !cleanup k counter
	    !ktop = 1     !cleanup k counter
            k = 0         !cleanup k counter
	    ktop = 0      !cleanup k counter

            ss_clay=0
            ss_sand=0
	    fac_count = 0.

	    do j=1,2*i    ! loop through both top (1<=j<=i) layer and sub (i+1<=j<=2*i) layer elements
        
                ! avg only across contributing raster grid cells and layers with orgC class as that assigned to tile n 
	        if((ss_oc_all(j)*sf >=  cF_lim(o_clp)).and.(ss_oc_all(j)*sf < cF_lim(o_clp + 1))) then 

                   if((ss_clay_all(j)*sf >= 0.).and.(ss_sand_all(j)*sf >= 0.)) then    ! avoiding no-data-values 

  	              if(j <= i) factor = 1.        ! top layer contribution
	              if(j  > i) factor = 2.33      ! sub layer contribution

	              poc_vec (n) = poc_vec(n) + ss_oc_all(j)*sf*factor   ! weighted sum of orgC
		      fac_count = fac_count + factor

                      k = k + 1                  ! counter for top and sub contributions        !cleanup k counter  
                      
                      if (j<=i) ktop = ktop + 1  ! counter for top contributions only           !cleanup k counter


!obsolete20220502  The code within the if-then and if-else statements below was nearly identical,
!obsolete20220502  except for the omission of the ktop counter from the else block.
!obsolete20220502
!obsolete20220502		      if(j <= i) then

                      ss_clay (k) = ss_clay_all(j)
                      ss_sand (k) = ss_sand_all(j)

                      ! adjust clay and sand content if outside joint physical bounds
                      if((ss_clay (k) + ss_sand (k)) > 9999) then  ! note: 9999 = 99.99%  (scale factor = 0.01)
                         if(ss_clay (k) >= ss_sand (k)) then
                            ss_sand (k) = 10000 - ss_clay (k)
                         else
                            ss_clay (k) = 10000 - ss_sand (k)
                         endif
                      endif
                      !k = k + 1           !cleanup k counter
                      !ktop = ktop + 1     !cleanup k counter

!obsolete20220502		       else
!obsolete20220502    			   ss_clay (k) = ss_clay_all(j)
!obsolete20220502			   ss_sand (k) = ss_sand_all(j)
!obsolete20220502                           if((ss_clay (k) + ss_sand (k)) > 9999) then
!obsolete20220502                              if(ss_clay (k) >= ss_sand (k)) then
!obsolete20220502                                ss_sand (k) = 10000 - ss_clay (k)
!obsolete20220502                              else
!obsolete20220502                                ss_clay (k) = 10000 - ss_sand (k)
!obsolete20220502                              endif
!obsolete20220502                           endif
!obsolete20220502                           !k = k + 1          !cleanup k counter                         
!obsolete20220502		       endif   
                   endif
                endif	
	    end do
	    
            !k = k - 1            !cleanup k counter
	    !ktop = ktop -1       !cleanup k counter

	    if(fac_count.ne.0) poc_vec (n) = poc_vec(n)/fac_count     ! normalize profile-layer orgC

	    !ncells_top_pro(n) = 100.*float(ktop)/float(i)              ! ncells_* not used
	    !ncells_sub_pro(n) = 100.*float(k-ktop)/float(i)            ! ncells_* not used

            ! debugging output
            if (write_debug) write (80+n,*)ktop,k,o_cl
            if (write_debug) write (80+n,*)ss_clay(1:k)
            if (write_debug) write (80+n,*)ss_sand(1:k)

            ! Determine the raster grid cell and layer j that has clay/sand content closest
            ! to the average (profile) clay/sand across all raster grid cells within the 
            ! dominant orgC class.

	    j = center_pix_int0 (sf, ktop,k, ss_clay(1:k),ss_sand(1:k))

            ! Assign soil class of raster grid cell and layer j to tile n
     
            if(j >=1) then 
               min_percs%clay_perc = ss_clay(j)*sf
               min_percs%sand_perc = ss_sand(j)*sf
               min_percs%silt_perc = 100. - ss_clay(j)*sf - ss_sand(j)*sf
               soil_class_com (n) = table_map(soil_class (min_percs),o_clp)  
            endif

            ! debugging output
            if (write_debug) write(80+n,*) j
            if (write_debug) write(80+n,*) soil_class_com (n) 
            if (write_debug) close(80+n)          

        endif

	deallocate (ss_clay,ss_sand,ss_clay_all,ss_sand_all,ss_oc_all)

      END DO
      END DO              ! loop through tiles
!$OMP ENDPARALLELDO

!      call process_peatmap (nx, ny, gfiler, pmap)

      ! -----------------------------------------------------------------------------
      !
      ! apply final touches and write output files: 
      ! - soil_param.first
      ! - tau_param.dat
      ! - catch_params.nc4 [soil hydraulic and srfexc-rzexc time scale parameters ONLY;
      !                     parameters from ar.new, bf.dat, and ts.dat parameters will be 
      !                     added to catch_params.nc4 by subroutine create_model_para_woesten()]
      
      inquire(file='clsm/catch_params.nc4', exist=CatchParamsNC_file_exists)

      if(CatchParamsNC_file_exists) then
         status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid) ; VERIFY_(STATUS)
         allocate (parms4file (1:maxcat, 1:10))
      endif

      fname ='clsm/soil_param.first'
      open (11,file=trim(fname),form='formatted',status='unknown',action = 'write')

      fname ='clsm/tau_param.dat'
      open (12,file=trim(fname),form='formatted',status='unknown',action = 'write')

      ! open catchment.def for reading tile index and Pfafstetter index
    
      fname='clsm/catchment.def'
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat     ! re-read header line

!obsolete20220502      fname ='clsm/mosaic_veg_typs_fracs'
!obsolete20220502      open (13,file=trim(fname),form='formatted',status='old',action = 'read')

      do n = 1, maxcat

!obsolete20220502         read (13,*) tindex,pfafindex,vtype

         ! fill gaps from neighbor for rare missing values caused by inconsistent masks

         if ((soil_class_top (n) == -9999).or.(soil_class_com (n) == -9999)) then

            ! if com-layer has data, the issue is only with top-layer

            if(soil_class_com (n) >= 1) soil_class_top (n) = soil_class_com (n)

            ! if there is nothing, look for the neighbor 
            ! 
            !  ^
            !  |
            !  | The comment above seems wrong; could have soil_class_top(n)>=1, unless
            !      earlier soil_class_com was set equal to soil_class_top whenever
            !      soil_class_top was available and soil_class_com was not. 

            if (soil_class_com (n) == -9999) then

               ! Look for neighbor j (regardless of soil_class_top) and set both
               ! soil_class_com(n) and soil_class_top(n) equal to the neighbor's 
               ! soil_class_com(j).

               do k = 1, maxcat
                  j  = 0
                  i1 = n - k
                  i2 = n + k
                  if(i1 >=     1) then
                     if (soil_class_com (i1) >=1) j = i1  ! tentatively use "lower" neighbor unless out of range
                  endif

                  if(1 <= i2 .and. i2 <=maxcat) then
                     if (soil_class_com (i2) >=1) j = i2  ! "upper" neighbor prevails unless out of range
                  endif

                  if (j > 0) then
                     soil_class_com (n) = soil_class_com (j)
                     !soil_class_top (n) = soil_class_com (n)    
                     soil_class_top (n) = soil_class_com (j)   ! should be faster/safer than usin gsoil_class_com(n)
                     grav_vec(n)        = grav_vec(j)
                     soc_vec(n)         = soc_vec (j)
                     poc_vec(n)         = poc_vec (j)
                  endif

                  if (soil_class_com (n) >=1) exit
               end do
            endif

         endif

         fac_surf = soil_class_top(n)
	 fac      = soil_class_com(n)

         if(use_PEATMAP) then
            ! the maximum peat soil depth is set to the value Michel used to derive parameters (5000.) 
            if (fac_surf == 253)  soildepth(n) = 5000. ! max(soildepth(n),5000.)
            ! reset subsurface to peat if surface soil type is peat
            if (fac_surf == 253)  fac      = 253
         endif

         wp_wetness = a_wp(fac) /a_poros(fac)
         
         this_cond  = a_aksat(fac)/exp(-1.0*zks*gnu)
         
         ! read tile index and Pfafstetter index from catchment.def
         
         read (10,*) tindex,pfafindex     

         ! write soil_param.first

         write (11,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4, f8.4)')tindex,pfafindex,      &
               fac_surf, fac, a_bee(fac),a_psis(fac),a_poros(fac),&
               this_cond,wp_wetness,soildepth(n),                 &
               grav_vec(n),soc_vec(n),poc_vec(n), &
               a_sand(fac_surf),a_clay(fac_surf),a_sand(fac),a_clay(fac), &
	       a_wpsurf(fac_surf)/a_porosurf(fac_surf),a_porosurf(fac_surf), pmap(n)

         ! write tau_param.dat
	       	    
         write (12,'(i10,i8,4f10.7)')tindex,pfafindex, &
	       atau_2cm(fac_surf),btau_2cm(fac_surf),atau(fac_surf),btau(fac_surf)  

         ! write catch_params.nc [soil hydraulic and srfexc-rzexc time scale parameters]

         if (allocated (parms4file)) then

            parms4file (n, 1) = a_bee(fac)
            parms4file (n, 2) = this_cond                  ! a_aksat(fac)/exp(-1.0*zks*gnu)
            parms4file (n, 3) = a_poros(fac)
            parms4file (n, 4) = a_psis(fac)
            parms4file (n, 5) = wp_wetness
            parms4file (n, 6) = soildepth(n)
            parms4file (n, 7) = atau_2cm(fac_surf)
            parms4file (n, 8) = btau_2cm(fac_surf)
            parms4file (n, 9) = atau(fac_surf)
            parms4file (n,10) = btau(fac_surf) 
  
  	 endif
      end do

      ! add "header" line to the bottom of soil_param.first
      
      write (11,'(a)')'                    '
      write (11,'(a)')'FMT=i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4,f8.4'
      write (11,'(a)')'TileIndex PfafID SoilClassTop SoilClassProfile BEE PSIS POROS Ks_at_SURF WPWET SoilDepth %Grav %OCTop %OCProf %Sand_top %Clay_top %Sand_prof %Clay_prof WPWET_SURF POROS_SURF PMAP'

      close (10, status = 'keep')	            
      close (11, status = 'keep')	            
      close (12, status = 'keep')	            

!obsolete20220502      close (13, status = 'keep')

      deallocate (data_vec1, data_vec2,data_vec3, data_vec4,data_vec5, data_vec6)
      deallocate (tileid_vec)
      deallocate (a_sand,a_clay,a_silt,a_oc,a_bee,a_psis,       &
            a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, &
            atau_2cm,btau_2cm)
      deallocate (soildepth, grav_vec,soc_vec,poc_vec,soil_class_top,soil_class_com)        
             !ncells_top,ncells_top_pro,ncells_sub_pro,soil_class_top,soil_class_com)            ! ncells_* not used

      ! write catch_params.nc4 [soil hydraulic and srfexc-rzexc time scale parameters]
      
      if(CatchParamsNC_file_exists) then
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ) ,(/1/),(/maxcat/), parms4file (:, 1)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'COND' ) ,(/1/),(/maxcat/), parms4file (:, 2)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'POROS') ,(/1/),(/maxcat/), parms4file (:, 3)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ) ,(/1/),(/maxcat/), parms4file (:, 4)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'WPWET') ,(/1/),(/maxcat/), parms4file (:, 5)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR') ,(/1/),(/maxcat/), parms4file (:, 6)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ATAU2') ,(/1/),(/maxcat/), parms4file (:, 7)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BTAU2') ,(/1/),(/maxcat/), parms4file (:, 8)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ATAU5') ,(/1/),(/maxcat/), parms4file (:, 9)) ; VERIFY_(STATUS) 
         status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BTAU5') ,(/1/),(/maxcat/), parms4file (:,10)) ; VERIFY_(STATUS) 
         STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
         DEALLOCATE (parms4file)
      endif

    END SUBROUTINE soil_para_hwsd

    ! --------------------------------------------------------------------------------------------------------

!obsolete20220502    INTEGER FUNCTION center_pix_int (sf,ktop, ktot, x,y,x0,y0,z0,ext_point)
!obsolete20220502      
!obsolete20220502      implicit none
!obsolete20220502      
!obsolete20220502      integer (kind =2), dimension (:), intent (in) :: x,y
!obsolete20220502      integer, intent (in) :: ktop,ktot
!obsolete20220502      real, intent (in) :: sf
!obsolete20220502      real :: xi,xj,yi,yj,xx0,yy0,zz0
!obsolete20220502      real, allocatable, dimension (:,:) :: length_m
!obsolete20220502      real, allocatable, dimension (:) :: length
!obsolete20220502      real, intent (inout) :: x0,y0,z0
!obsolete20220502      integer :: i,j,npix
!obsolete20220502      logical, intent(in) :: ext_point
!obsolete20220502      real :: zi, zj
!obsolete20220502      
!obsolete20220502      allocate (length_m (1:ktot,1:ktot))
!obsolete20220502      allocate (length   (1:ktot))
!obsolete20220502      length_m =0.
!obsolete20220502      length   =0.
!obsolete20220502      
!obsolete20220502      center_pix_int = -9999
!obsolete20220502      if(ktot /= 0) then
!obsolete20220502         do i = 1,ktot
!obsolete20220502            xi = sf*x(i)
!obsolete20220502            yi = sf*y(i)
!obsolete20220502            zi = 100. - xi - yi
!obsolete20220502            if (.not. ext_point) then
!obsolete20220502               x0 = xi
!obsolete20220502               y0 = yi
!obsolete20220502               z0 = zi
!obsolete20220502            endif
!obsolete20220502            
!obsolete20220502            do j = 1,ktot
!obsolete20220502               xj = sf*x(j)
!obsolete20220502               yj = sf*y(j)
!obsolete20220502               zj = 100. - xj - yj
!obsolete20220502               xx0= xj - x0
!obsolete20220502               yy0= yj - y0
!obsolete20220502               zz0= zj - z0
!obsolete20220502               
!obsolete20220502               if(ktot > ktop) then 
!obsolete20220502                  if(j <= ktop) then
!obsolete20220502                     length_m (i,j) = (xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5
!obsolete20220502                  else
!obsolete20220502                     length_m (i,j) = 2.33*((xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5)
!obsolete20220502                  endif
!obsolete20220502               else
!obsolete20220502                  length_m (i,j) = (xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5
!obsolete20220502               endif
!obsolete20220502            end do
!obsolete20220502            length (i) = sum(length_m (i,:))
!obsolete20220502         end do
!obsolete20220502         
!obsolete20220502         center_pix_int = minloc(length,dim=1)
!obsolete20220502      endif
!obsolete20220502      
!obsolete20220502    END FUNCTION center_pix_int
!obsolete20220502      
!obsolete20220502    !
!obsolete20220502

    ! ====================================================================
    !
    
    INTEGER FUNCTION center_pix_int0 (sf,ktop, ktot, x,y)
      
      implicit none

      ! In a nutshell, given a list of clay/sand pairs, this function determines 
      ! the element (pair) in this list that is closest to the average clay/sand 
      ! across all pairs.  
      !
      ! The input list of clay/sand can consist of only top (0-30) layer clay/sand
      ! pairs, or of pairs of clay/sand pairs for the top (0-30) and sub (30-70) 
      ! layers.  In the latter case, a weighted average is computed.
      !
      ! This is to ensure that ultimately the clay/sand values assigned to a tile
      ! represent an actual soil class.
      !
      ! sf = 0.01 (integer to real scale factor)
      ! ktop = # of pixels in top layer
      ! ktot = total # of pixels, top + subsurface combined
      ! x (clay), y (sand)
      integer (kind =2), dimension (:), intent (in) :: x,y
      integer,                          intent (in) :: ktop,ktot
      real,                             intent (in) :: sf

      real :: xi,xj,yi,yj
      real :: length
      
      integer :: i,j,npix
      real :: zi, zj, mindist,xc,yc,zc
      
      length          = 0.
      
      center_pix_int0 = -9999
      
      ! compute average clay/sand
      
      if(ktot /= 0) then
         ! There should be some data pixels
         if(ktot > ktop) then 
            ! Have both layers
            if(ktop > 0) then
               ! There are data in top layer
               xc = sf*0.3*sum(real(x(1:ktop)))/real(ktop) + sf*0.7*sum(real(x(ktop+1 : ktot)))/real(ktot - ktop)  
               yc = sf*0.3*sum(real(y(1:ktop)))/real(ktop) + sf*0.7*sum(real(y(ktop+1 : ktot)))/real(ktot - ktop)
            else
               ! There are no data in top layer
               xc = sf*sum(real(x(1:ktot)))/real(ktot)  
               yc = sf*sum(real(y(1:ktot)))/real(ktot)         
            endif
         else
            ! working on Top layer alone
            xc = sf*sum(real(x(1:ktot)))/real(ktot)  
            yc = sf*sum(real(y(1:ktot)))/real(ktot)
         endif
         zc = 100. - xc - yc              ! silt [percent]
      endif
      
      mindist=100000.*100000.
      
      do i = 1,ktot
         xi = sf*x(i)
         yi = sf*y(i)
         zi = 100. - xi - yi
         length = (xi-xc)**2+(yi-yc)**2+(zi-zc)**2
         if(mindist>length)then
            mindist=length
            center_pix_int0=i
         end if
      end do
      !print *,ktop,ktot,center_pix_int0
      
    END FUNCTION center_pix_int0
    
    ! --------------------------------------------------------------------------------------

! this subroutine seems obsolete, commented out for now - reichle, 9 Feb 2022

!   SUBROUTINE process_peatmap (nc, nr, gfiler, pmap)
!     
!     implicit none
!     integer  , parameter                         :: N_lon_pm = 43200, N_lat_pm = 21600
!     integer, intent (in)                         :: nc, nr
!     real, pointer, dimension (:), intent (inout) :: pmap
!     character(*), intent (in)                    :: gfiler
!     integer                                      :: i,j, status, varid, ncid
!     integer                                      :: NTILES        
!     REAL, ALLOCATABLE, dimension (:)             :: count_pix
!     REAL, ALLOCATABLE, dimension (:,:)           :: data_grid, pm_grid
!     INTEGER, ALLOCATABLE, dimension (:,:)        :: tile_id
!     character*100                                :: fout    
!     
!     ! Reading number of tiles
!     ! -----------------------
!     
!     open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
!     
!     read (20, *) NTILES
!     
!     close (20, status = 'keep')
!     
!     ! READ PEATMAP source data files and regrid
!     ! -----------------------------------------
!     
!     status  = NF_OPEN ('data/CATCH/PEATMAP_mask.nc4', NF_NOWRITE, ncid)
!     
!     allocate (pm_grid   (1 : NC      , 1 : NR))
!     allocate (data_grid (1 : N_lon_pm, 1 : N_lat_pm)) 
!     
!     status  = NF_INQ_VARID (ncid,'PEATMAP',VarID) ; VERIFY_(STATUS)
!     status  = NF_GET_VARA_REAL (ncid,VarID, (/1,1/),(/N_lon_pm, N_lat_pm/), data_grid) ; VERIFY_(STATUS)
!     
!     call RegridRasterReal(data_grid, pm_grid)
!     
!     status = NF_CLOSE(ncid)
!     
!     ! Grid to tile
!     ! ------------
!     
!     ! Reading tile-id raster file
!     
!     allocate(tile_id(1:nc,1:nr))
!     
!     open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
!          form='unformatted',convert='little_endian')
!     
!     do j=1,nr
!        read(10)tile_id(:,j)
!     end do
!     
!     close (10,status='keep')     
!     
!     allocate (pmap      (1:NTILES))
!     allocate (count_pix (1:NTILES))
!     
!     pmap      = 0.
!     count_pix = 0.
!     
!     do j = 1,nr
!        do i = 1, nc
!           if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then                
!              if(pm_grid(i,j) > 0.)  pmap (tile_id(i,j)) = pmap (tile_id(i,j)) + pm_grid(i,j)
!              count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
!           endif
!        end do
!     end do
!     
!     where (count_pix >   0.) pmap = pmap/count_pix
!     
!     deallocate (count_pix)
!     deallocate (pm_grid)
!     deallocate (tile_id)
!     
!   END SUBROUTINE process_peatmap
    
! ====================================================================

  SUBROUTINE grid2tile_ndep_t2m_alb (irst,jrst,gfiler)  

    implicit none

    !PFT Description
    !  1 needleleaf evergreen temperate tree
    !  2 needleleaf evergreen boreal tree
    !  3 needleleaf deciduous boreal tree
    !  4 broadleaf evergreen tropical tree
    !  5 broadleaf evergreen temperate tree
    !  6 broadleaf deciduous tropical tree
    !  7 broadleaf deciduous temperate tree
    !  8 broadleaf deciduous boreal tree
    !  9 broadleaf evergreen temperate shrub
    ! 10 broadleaf deciduous temperate shrub [moisture + deciduous]
    ! 11 broadleaf deciduous temperate shrub [moisture stress only]
    ! 12 broadleaf deciduous boreal shrub
    ! 13 arctic c3 grass
    ! 14 cool c3 grass [moisture + deciduous]
    ! 15 cool c3 grass [moisture stress only]
    ! 16 warm c4 grass [moisture + deciduous]
    ! 17 warm c4 grass [moisture stress only]
    ! 18 crop          [moisture + deciduous]
    ! 19 crop          [moisture stress only]

    integer  ,intent (in)       :: irst, jrst
    character (*), intent (in)  :: gfiler    

    integer, parameter :: nveg =  4   ! number of veg types
    integer, parameter :: npft = 19   ! number of PFT
    
    integer, parameter :: iclm = 1152 ! lon dimension CLM NDEP data
    integer, parameter :: jclm =  768 ! lat dimension CLM NDEP data

    integer, parameter :: iprn =  360 ! lon dimension Princeton 2m T data
    integer, parameter :: jprn =  180 ! lat dimension Princeton 2m T data

    integer, parameter :: imra =  576 ! lon dimension MERRA-2 2m T data
    integer, parameter :: jmra =  361 ! lat dimension MERRA-2 2m T data

    integer, parameter :: ialb = 7200 ! lon dimension MODIS soil background albedo data
    integer, parameter :: jalb = 3600 ! lat dimension MODIS soil background albedo data

!    integer, parameter :: irst = 43200 ! lon dimension of tile raster file
!    integer, parameter :: jrst = 21600 ! lat dimension of tile raster file

    logical, parameter :: dir_access_files = .false.

    integer, dimension (:,:), allocatable :: tile_id ! tile raster file 

    real, allocatable :: ndep_tile(:), t2mp_tile(:), t2mm_tile(:), alb_tile(:,:,:)
    real, allocatable :: data_grid (:,:), vector(:) 
    integer, allocatable :: icount(:)

    character :: ctype*1, cband*1

    real :: rdum, ftot, xg, yg, fill, alonw, alats, alone, alatn, rlonw, rlats, rlone, rlatn, xx, yy
    integer :: i, j, n, im, jm, lwi, idum, ntiles, nland, nv, ix, jx, itype, iband, isum, ntl, np, jalbx, ialbx, ncid, status
    logical :: file_exists

    ! read nland from catchment.def
    ! -----------------------------
 
    open (8, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
         action =  'read')
    
    read (8,*) nland

    close(8, status = 'keep')

    ! Read tile-id raster file; used for mapping gridded fields to tile space
    ! -----------------------------------------------------------------------

    allocate (tile_id(1:irst,1:jrst))
    allocate(vector(nland))
    allocate(icount(nland))
    
    open(8,file=trim(gfiler)//'.rst' ,status='old',action='read',form='unformatted')
    do j=1,jrst
      read(8) tile_id(:,j)
    end do
    close (8)

    !=====================================================================================================
    ! The below correction was moved to esa2clm - SM
    !
    ! create PFT data for offline
    ! ---------------------------
    !    open(8,file='clsm/CLM_veg_typs_fracs',form='formatted',status='old')
    !        
    !    allocate (fveg(nland,nveg))
    !    allocate (ityp(nland,nveg))
    !
    !    ! loop over tiles
    !    ! ---------------
    !
    !    if (dir_access_files) open(9,file='clsm/pft.dat',form='unformatted',convert='big_endian', &
    !            status='unknown',access='direct',recl=8)
    !
    !    do n = 1,nland
    !      read(8,*) idum,idum,ityp(n,:),fveg(n,:)
    !      fveg(n,:) = 0.01*fveg(n,:) ! convert from percent to fraction
    !
    !      ! fractions must sum to 1
    !      ! -----------------------
    !      ftot = sum(fveg(n,:))
    !      if(ftot /= 1.) fveg(n,:) = fveg(n,:) / ftot
    !
    !      ! prevent trivial fractions gkw: the carbon model ignores fractions that are less than or equal to 1e-4
    !      ! -------------------------
    !      if(fveg(n,1) <= 1.e-4) then
    !         fveg(n,2) = fveg(n,2) + fveg(n,1)
    !         fveg(n,1) = 0.
    !      endif
    !      
    !      if(fveg(n,2) <= 1.e-4) then
    !         fveg(n,1) = fveg(n,1) + fveg(n,2)
    !         fveg(n,2) = 0.
    !      endif
    !      
    !      if(fveg(n,3) <= 1.e-4) then
    !         if(fveg(n,4) > 1.e-4) then
    !            fveg(n,4) = fveg(n,4) + fveg(n,3)
    !         else if(fveg(n,2) > 1.e-4) then
    !            fveg(n,2) = fveg(n,2) + fveg(n,3)
    !         else if(fveg(n,1) > 1.e-4) then
    !            fveg(n,1) = fveg(n,1) + fveg(n,3)
    !         else
    !            print *, 'fveg3:',n,ityp(n,:),fveg(n,:)
    !            stop 'fveg3'
    !         endif
    !         fveg(n,3) = 0.
    !      endif
    !      
    !      if(fveg(n,4) <= 1.e-4) then
    !         if(fveg(n,3) > 1.e-4) then
    !            fveg(n,3) = fveg(n,3) + fveg(n,4)
    !         else if(fveg(n,2) > 1.e-4) then
    !            fveg(n,2) = fveg(n,2) + fveg(n,4)
    !         else if(fveg(n,1) > 1.e-4) then
    !            fveg(n,1) = fveg(n,1) + fveg(n,4)
    !         else
    !            print *, 'fveg4:',n,ityp(n,:),fveg(n,:)
    !            stop 'fveg4'
    !         endif
    !         fveg(n,4) = 0.
    !      endif
    !      
    !      if(abs(sum(fveg(n,:))-1.) > 1.e-6) stop 'fracs/=1'
    
    !      if (dir_access_files) write(9,rec=n) ityp(n,:),fveg(n,:)
    !
    !80    format('pft:',i8,2f10.4,4i3,4f7.4)
    !      
    !   end do ! end tile loop
    !
    !   close(8)   
    !   if (dir_access_files) close(9)
   
   !=====================================================================================================
      
   ! nitrogen deposition
   ! -------------------

   allocate(data_grid(iclm,jclm))
   allocate(ndep_tile(nland))
   
   open(8,file='data/CATCH/CLSM-CN/ndep_clm_simyr2000_0.23x0.31_c091106.gdat', &
        form='unformatted',status='old')
   read(8) data_grid
   close(8)
   
   ! regridding to raster grid irst x jrst
   ! -------------------------------------

   xx = iclm/real(irst)
   yy = (jclm-1)/real(jrst)  ! gkw: subtract 1, since 1 & jclm are centered at pole (dlat=180/(jclm-1))
   
   vector = 0.
   icount = 0 

   do j = 1,jrst
      jx = (j-1)*yy + 1 + 0.5 ! add half because CLM data is centered at south pole
      if(jx<1 .or. jx>jclm) stop 'jclm'
      do i = 1,irst

         if(tile_id(i,j)>0 .and. tile_id(i,j)<=nland) then 
            ix = (i-1)*xx + 1 + 0.5 ! add half because CLM data is centered on dateline
            ix = ix + iclm/2
            if(ix > iclm) ix = ix - iclm  ! shift 180 degrees; data starts at 0 lon
            if(ix<1 .or. ix>iclm) stop 'iclm'
            
            if(data_grid(ix,jx) >= 0.) then

               ! aggregation on to catchment-tiles
               ! ---------------------------------

               vector(tile_id(i,j)) = vector(tile_id(i,j)) + data_grid(ix,jx)
               icount(tile_id(i,j)) = icount(tile_id(i,j)) + 1 

            endif
         endif
      end do
   end do
   
   where (icount > 0) ndep_tile = (vector/icount)* (1.e9 / (86400. * 365.)) ! g/m2/yr --> ng/m2/s (for offline; GEOS5 will use g/m2/s)
   
   if (dir_access_files) then 
      ! write tile-space data
      ! ---------------------
      open(9,file='clsm/ndep.dat',form='unformatted',convert='big_endian', &
           status='unknown',access='direct',recl=1)
      do n = 1,nland
         write(9,rec=n) ndep_tile(n)
      end do
      close(9)
   endif
   deallocate(data_grid)
   
   !=====================================================================================================
   
   
   ! annual mean 2m air temperature climatology: Sheffield Princeton 1948-2012
   ! -------------------------------------------------------------------------
   allocate(data_grid(iprn,jprn))
   allocate(t2mp_tile(nland))
   
   open(8,file='data/CATCH/CLSM-CN/princeton_annual_mean_T2m_1948-2012.gdat', &
        form='unformatted',status='old')
   read(8) data_grid
   close(8)
   
   ! regridding to raster grid irst x jrst
   ! -------------------------------------
   xx = iprn/real(irst)
   yy = jprn/real(jrst)
   
   vector = 0.
   icount = 0

   do j = 1,jrst
      jx = (j-1)*yy + 1
      if(jx<1 .or. jx>jprn) stop 'jprn'
      do i = 1,irst

         if(tile_id(i,j)>0 .and. tile_id(i,j)<=nland) then 
            ix = (i-1)*xx + 1
            ix = ix + iprn/2              ! shift 180 degrees; data starts at 0 lon
            if(ix > iprn) ix = ix - iprn
            if(ix<1 .or. ix>iprn) stop 'iprn'
            if(data_grid(ix,jx) >= 0.) then   

               ! aggregation on to catchment-tiles
               ! ---------------------------------

               vector(tile_id(i,j)) = vector(tile_id(i,j)) + data_grid(ix,jx)
               icount(tile_id(i,j)) = icount(tile_id(i,j)) + 1 

            endif
         endif
      end do
   end do
   
   where (icount > 0) t2mp_tile = vector/icount
   
   if (dir_access_files) then 
      ! write tile-space data
      ! ---------------------
      open(9,file='clsm/cli_t2m_princeton.dat',form='unformatted',convert='big_endian', &
           status='unknown',access='direct',recl=1)
      do n = 1,nland
         write(9,rec=n) t2mp_tile(n)
      end do
      close(9)
   endif

   deallocate(data_grid)
   
   !=====================================================================================================


   ! annual mean 2m air temperature climatology: MERRA-2 1980-2014
   ! -------------------------------------------------------------
   allocate(data_grid(imra,jmra))
   allocate(t2mm_tile(nland))
   
   open(8,file='data/CATCH/CLSM-CN/MERRA2_annual_mean_T2m_1980-2014.gdat', &
        form='unformatted',status='old')
   read(8) data_grid
   close(8)
   
   ! regridding to raster grid irst x jrst
   ! -------------------------------------
   xx = imra/real(irst)
   yy = (jmra-1)/real(jrst)
   
   vector = 0.
   icount = 0

   do j = 1,jrst
      jx = (j-1)*yy + 1 + 0.5
      if(jx<1 .or. jx>jmra) stop 'jmra'
      do i = 1,irst

         if(tile_id(i,j)>0 .and. tile_id(i,j)<=nland) then 
            ix = (i-1)*xx + 1 + 0.5
            if(ix > imra) ix = ix - imra
            if(ix<1 .or. ix>imra) stop 'imra'
            if(data_grid (ix,jx) >= 0.) then
               
               ! aggregation on to catchment-tiles
               ! ---------------------------------

               vector(tile_id(i,j)) = vector(tile_id(i,j)) + data_grid(ix,jx)
               icount(tile_id(i,j)) = icount(tile_id(i,j)) + 1 

            endif
         endif
      end do
   end do

   where (icount > 0) t2mm_tile = vector/icount
   
   if (dir_access_files) then 
      ! write tile-space data
      ! ---------------------
      open(9,file='clsm/cli_t2m_merra2.dat',form='unformatted',convert='big_endian', &
           status='unknown',access='direct',recl=1)
      do n = 1,nland
         write(9,rec=n) t2mm_tile(n)
      end do
      close(9)
   endif

   deallocate(data_grid)
    
    !=====================================================================================================
    
    
    ! read soil background albedo if tile falls in MODIS grid cell, use that value  gkw: may want to interpolate or aggregate
    ! ----------------------------------------------------------------------------
    allocate(data_grid(ialb,jalb))
    allocate(alb_tile(nland,2,2))
    
    do itype = 1,2   
       do iband = 1,2
          
          if(itype == 1) then
             ctype = 'b' ! "b" (direct, black sky)
          else
             ctype = 'w' ! "w" (diffuse, white sky)
          endif
          
          if(iband == 1) then
             cband = '1' ! "1" (visible)
             fill = 0.10 ! fill value to use when albedo not defined over land
          else
             cband = '2' ! "2" (near IR)
             fill = 0.07
          endif
          
          open(8,file='data/CATCH/CLSM-CN/modis_'//ctype//'sa_soil_bb'//cband//'_cmg', &
               form='unformatted',status='old',access='direct',recl=ialb*jalb)
          read(8,rec=1) (data_grid(:,j), j = jalb,1,-1) ! data is from north to south
          where(data_grid <= 0.) data_grid = fill
          close(8)
          
          ! regridding to raster grid irst x jrst
          ! -------------------------------------
          xx = ialb/real(irst)
          yy = jalb/real(jrst)
          
          vector = 0.
          icount = 0
         
          do j = 1,jrst
             jx = (j-1)*yy + 1
             if(jx<1 .or. jx>jalb) stop 'jalb'
             do i = 1,irst

                if(tile_id(i,j)>0 .and. tile_id(i,j)<=nland) then 
                   ix = (i-1)*xx + 1
                   if(ix<1 .or. ix>ialb) stop 'ialb'
                   if(data_grid (ix,jx) >= 0.) then

                      ! aggregation on to catchment-tiles
                      ! ---------------------------------
                      
                      vector(tile_id(i,j)) = vector(tile_id(i,j)) + data_grid(ix,jx)
                      icount(tile_id(i,j)) = icount(tile_id(i,j)) + 1 

                   endif
                endif
             end do
          end do
          
          where (icount > 0) vector = vector/icount
          alb_tile(:,itype,iband) = vector (:)
          
          if (dir_access_files) then 
             ! write tile-space data
             ! ---------------------
             open(9,file='clsm/alb_'//ctype//'_'//cband//'_tile.dat', &
                  form='unformatted',convert='big_endian',access='direct',recl=nland)
             write(9,rec=1) alb_tile(:,itype,iband)
             close(9)
          endif
       end do ! end band loop
    end do   ! end type loop

    !=====================================================================================================
    ! Writing output file
    ! --------------------

    open (10, file = 'clsm/CLM_NDep_SoilAlb_T2m', form = 'formatted', status ='unknown', &
         action = 'write')
    
    do n = 1,nland
       write (10, '(f10.4,4f7.4,2f8.3)') ndep_tile(n), &
            alb_tile(n,1,1),alb_tile(n,2,1),     &
            alb_tile(n,1,2),alb_tile(n,2,2),     &
            t2mm_tile(n) ,t2mp_tile(n) 
       ! VISDR, VISDF, NIRDR, NIRDF
    end do

    inquire(file='clsm/catchcn_params.nc4', exist=file_exists)

    if(file_exists) then
       status = NF_OPEN ('clsm/catchcn_params.nc4', NF_WRITE, ncid ) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'NDEP'    ) ,(/1/),(/nland/), ndep_tile      ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BGALBVR' ) ,(/1/),(/nland/), alb_tile(:,1,1)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BGALBVF' ) ,(/1/),(/nland/), alb_tile(:,2,1)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BGALBNR' ) ,(/1/),(/nland/), alb_tile(:,1,2)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BGALBNF' ) ,(/1/),(/nland/), alb_tile(:,2,2)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'T2_M'    ) ,(/1/),(/nland/), t2mm_tile      ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'T2_S'    ) ,(/1/),(/nland/), t2mp_tile      ) ; VERIFY_(STATUS)
       STATUS = NF_CLOSE (NCID) ; VERIFY_(STATUS)
    endif

    close (10, status ='keep')

  end SUBROUTINE grid2tile_ndep_t2m_alb

!
! --------------------------------------------------------------------------------------
!

        SUBROUTINE CREATE_ROUT_PARA_FILE (NC, NR, gfile, MGRID, deltaXY)
          
          IMPLICIT NONE
          
          INTEGER,     INTENT (IN)                  :: NC, NR
          character*5, INTENT (IN), OPTIONAL        :: MGRID
          REAL,        INTENT (IN), OPTIONAL        :: deltaXY
          character(*),INTENT (IN)                  :: gfile
          real,   allocatable,  dimension (:)       :: pfaf_area
          integer,allocatable,  dimension (:)       :: pfaf_index
          INTEGER :: NBINS, NPLUS, PFAF,N,L, I1,I2,J1,J2,I,J,K,IL1,IL2,JL1,JL2,NC_RAT
          REAL    :: mnx,mxx,mny,mxy, lats, dxy30
          INTEGER, PARAMETER :: NC_SRTM = 21600, NR_SRTM = 10800 
          REAL    :: dx =360._8/NC_ESA,dy = 180._8/NR_ESA, d2r = PI/180._8
          integer :: max_pfaf_smap = 40
          INTEGER, TARGET,  ALLOCATABLE, DIMENSION (:,:) :: raster, tileid_index,&
               SUBSET_MSK
          INTEGER, POINTER,              DIMENSION (:,:) :: SUBSET_RST
          REAL,             ALLOCATABLE, DIMENSION (:,:) :: SUBSET_AREA
          REAL,             ALLOCATABLE, DIMENSION (:)   :: loc_val, loc_area
          INTEGER,          ALLOCATABLE, DIMENSION (:)   :: density, loc_int
          logical,          ALLOCATABLE, DIMENSION (:)   :: unq_mask   

          INCLUDE 'netcdf.inc'

          INTEGER :: CellID, MaxID, d2(2), STATUS, VID, NCID, NCID_MSK, NCAT
          integer, dimension(8)   :: date_time_values
          character (22)          :: time_stamp    

          ! Reading raster file

          allocate(raster               (1:nc,1:nr))

          open (10, file ='rst/'//trim(gfile)//'.rst',form='unformatted',status='old',  &
               action='read')
          
          do j=1,nr
             read(10)(raster (i,j),i=1,nc)
          end do

          close (10,status='keep')
 
          ! Creating SMAP-Catch_TransferData.nc that contains SMAP cells to Pfafstetter transfer infor

          open (10,file='clsm/catchment.def',form='formatted',status='old', action = 'read')
          
          read (10,*) NCAT

          if (PRESENT (MGRID)) then
             if (trim(MGRID) == 'M25') max_pfaf_smap = 30
             if (trim(MGRID) == 'M09') max_pfaf_smap = 12
             if (trim(MGRID) == 'M03') max_pfaf_smap = 5
          endif

          if (PRESENT (deltaXY)) then
             if (deltaXY <  0.125) max_pfaf_smap = 15
             if (deltaXY >= 0.125) max_pfaf_smap = 15
             if (deltaXY >= 0.25 ) max_pfaf_smap = 40
             if (deltaXY >= 0.5  ) max_pfaf_smap = 100
             if (deltaXY >= 1.0  ) max_pfaf_smap = 250
          endif

          status = NF_CREATE ('clsm/Grid2Catch_TransferData.nc', NF_NETCDF4, NCID)
          status = NF_DEF_DIM(NCID, 'N_GRID_CELLS'    , ncat,CellID)
          status = NF_DEF_DIM(NCID, 'MAX_CAT_PER_CELL', max_pfaf_smap ,MaxID )

          d2(1) = MaxID
          d2(2) = CellID

          status = NF_DEF_VAR(NCID, 'NCats_in_GRID', NF_INT  , 1 ,CellID, vid)
          status = NF_PUT_ATT_TEXT(NCID, vid, 'long_name',&
               LEN_TRIM('No. of watersheds contributed to the Grid cell'), &
               trim('No. of watersheds contributed to the Grid cell')) 
          status = NF_DEF_VAR(NCID, 'Pfaf_Index'   , NF_INT  , 2 ,d2    , vid)
          status = NF_PUT_ATT_TEXT(NCID, vid, 'long_name',&
               LEN_TRIM('Pfaf indices of those contributing watersheds'), &
               trim('Pfaf indices of those contributing watersheds')) 
          status = NF_DEF_VAR(NCID, 'Pfaf_Area '   , NF_FLOAT, 2 ,d2    , vid)
          status = NF_PUT_ATT_TEXT(NCID, vid, 'long_name', &
               LEN_TRIM('Area of watershed fraction'),&
               trim('Area of watershed fraction')) 
          status = NF_PUT_ATT_TEXT(NCID, vid, 'units',&
               LEN_TRIM('km2'), trim('km2')) 
!          status = NF_DEF_VAR(NCID, 'Pfaf_Frac '   , NF_FLOAT, 2 ,d2    , vid)
!          status = NF_PUT_ATT_TEXT(NCID, vid, 'long_name', &
!               LEN_TRIM('Fraction of Pfaf catchment contributed to the SMAP cell'),&
!               trim('Fraction of Pfaf catchment contributed to the SMAP cell')) 
!
!  Global attributes
!
          call date_and_time(VALUES=date_time_values)
          
          write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
               date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
               date_time_values(5),':',date_time_values(6),':',date_time_values(7)

          status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'CreatedBy', LEN_TRIM('Sarith Mahanama @ GMAO/GSFC/NASA'),  &
               trim('Sarith Mahanama @ GMAO/GSFC/NASA'))
          status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'Contact', LEN_TRIM('sarith.p.mahanama@nasa.gov'),          &
               trim('sarith.p.mahanama@nasa.gov'))
          status = NF_PUT_ATT_TEXT(NCID, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))
          status = NF_ENDDEF(NCID) 

          ! Now computing SMAP-cells to Pfafcatchment fractional areas 

          status    = NF_OPEN ('data/CATCH/GEOS5_10arcsec_mask.nc', NF_NOWRITE, ncid_msk)	
          nbins = 1

          allocate (pfaf_area (1:max_pfaf_smap))
          allocate (pfaf_index(1:max_pfaf_smap))
          
          dxy30  = 360._8/nc
          NC_RAT = nc_esa/nc

          DO N = 1, NCAT
             
             pfaf_index = 0
             pfaf_area  = 0.

             READ (10,'(i10,i8,5(2x,f9.4), i4)')l,pfaf,mnx,mxx,mny,mxy

             IL1 = FLOOR  ((180. + mnx)/DXY30 + 1.)
             IL2 = CEILING((180. + mxx)/DXY30 + 1.)
             JL1 = FLOOR  (( 90. + mny)/DXY30 + 1.)
             JL2 = CEILING(( 90. + mxy)/DXY30 + 1.)
             
             IF(IL2 > NC) IL2 = NC
             IF(JL2 > NR) JL2 = NR

             I1  = NC_RAT * IL1 - (NC_RAT -1)
             I2  = NC_RAT * IL2
             J1  = NC_RAT * JL1 - (NC_RAT -1)
             J2  = NC_RAT * JL2
            
             ALLOCATE (SUBSET_MSK (1 : I2 - I1 +1 , 1 : J2 - J1 + 1))
             ALLOCATE (SUBSET_AREA(1 : I2 - I1 +1 , 1 : J2 - J1 + 1))
             ALLOCATE (TILEID_INDEX(1 : I2 - I1 +1 , 1 : J2 - J1 + 1))

             DO J = J1, J2 
                lats = -90._8 + (j - 0.5_8)*dy
                SUBSET_AREA(:,J-J1 + 1) = (sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)
             END DO

             status   = NF_GET_VARA_INT (ncid_msk,4,(/I1,J1/),(/I2 - I1 +1,J2 - J1 + 1/),SUBSET_MSK)

             if (associated (subset_rst )) NULLIFY (subset_rst)
             SUBSET_RST => RASTER (IL1 : IL2, JL1 : JL2)

             call RegridRaster(SUBSET_RST, tileid_index)

             NPLUS = count((tileid_index ==N).AND.(SUBSET_MSK >= 1 .and. subset_MSK <= SRTM_maxcat))
             allocate (loc_int (1:NPLUS))
             allocate (loc_area(1:NPLUS))
             allocate (unq_mask(1:NPLUS))

             loc_int = pack(SUBSET_MSK  ,mask = ((tileid_index ==N).AND.(SUBSET_MSK >= 1 .and. subset_MSK <= SRTM_maxcat)))
             loc_area= pack(SUBSET_AREA ,mask = ((tileid_index ==N).AND.(SUBSET_MSK >= 1 .and. subset_MSK <= SRTM_maxcat)))

             call MAPL_Sort (loc_int, loc_area)

             unq_mask = .true.
             do K = 2,NPLUS 
                unq_mask(K) = .not.(loc_int(K) == loc_int(K-1)) ! count number of unique numbers in loc_int for binning
             end do
             NBINS = count(unq_mask)
             
             if (NBINS > max_pfaf_smap) then
                print *, 'NBINS exceeded max_pfaf_smap', NBINS, max_pfaf_smap
                STOP
             endif

             if (NBINS > 1) then
                L = 1
                pfaf_index(L) = loc_int (1)
                pfaf_area (L) = loc_area(1)
                DO K = 2,NPLUS
                   IF(.not.(loc_int(K) == loc_int(K-1))) L = L + 1
                   pfaf_index(L) = loc_int (K)
                   pfaf_area (L) = pfaf_area (L) + loc_area(K) * MAPL_RADIUS * MAPL_RADIUS/1000./1000.
                END DO
             else
                IF(NBINS == 1) THEN
                   pfaf_index(1) = loc_int (1)
                   pfaf_area (1) = sum (loc_area(1:NPLUS))
                   pfaf_area (1) = pfaf_area(1) * MAPL_RADIUS * MAPL_RADIUS/1000./1000.
                ELSE
                   PRINT *,'NO Catchments so skipping'
                   NBINS = 1
                   pfaf_index(1) = -1
                   pfaf_area (1)= -9999.
                ENDIF
             endif

             status = NF_PUT_VARA_INT (NCID, 1,(/N/),(/1/),nbins)
             status = NF_PUT_VARA_INT (NCID, 2,(/1,N/),(/nbins,1/),pfaf_index(1:nbins))
             status = NF_PUT_VARA_REAL(NCID, 3,(/1,N/),(/nbins,1/),pfaf_area (1:nbins))
             
             DEALLOCATE (SUBSET_MSK,SUBSET_AREA, loc_int,loc_area,unq_mask,tileid_index)
          END DO

        DEALLOCATE (RASTER) 
        status    = NF_CLOSE (ncid)  
        status    = NF_CLOSE (ncid_msk)
        close (10, status = 'keep')  

        END SUBROUTINE CREATE_ROUT_PARA_FILE

! -------------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE CLM45_fixed_parameters (nc,nr,gfiler)

      implicit none
      
      ! producing  CLM4.5 fixed parameters :


      ! 1) Population density /discover/nobackup/fzeng/clm4-to-clm4.5/data/firedata4.5/clmforc.Li_2012_hdm_0.5x0.5_AVHRR_simyr1850-2010_c130401.nc
      ! Use 2010       
      ! 2) /gpfsm/dnb31/fzeng/clm4-to-clm4.5/data/rawdata4.5
      ! mksrf_abm_0.5x0.5_AVHRR_simyr2000.c130201.nc
      ! mksrf_gdp_0.5x0.5_AVHRR_simyr2000.c130228.nc
      ! mksrf_peatf_0.5x0.5_AVHRR_simyr2000.c130228.nc
      ! one value per tile
      ! 3) field capacity one value per tile
     
      integer, intent (in)               :: nc, nr 
      character(*), intent (in)          :: gfiler
      integer  , parameter               :: N_lon_clm = 720, N_lat_clm = 360
      real, parameter                    :: dxy_clm = 0.5
      integer                            :: i,j, status, varid, ncid_hdm, ncid_abm, ncid_gdp, ncid_peatf
      integer                            :: NTILES, tid, cid, ABM_INT, sc_top, sc_com
      REAL, ALLOCATABLE, dimension (:)   :: hdm, abm, gdp, peatf
      REAL, ALLOCATABLE, dimension (:,:) :: hdm_grid, gdp_grid, peatf_grid, data_grid, count_pix
      INTEGER, ALLOCATABLE, dimension (:,:) :: tile_id, abm_grid, int_grid
      REAL                               :: hdm_r, gdp_r, peatf_r
      character*100                      :: fout
      real                               ::                     &
      	    a_sand,a_clay,a_silt,a_oc,a_bee,a_psis,             &
            a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, &
            atau_2cm,btau_2cm, field_cap (n_SoilClasses) 

      ! Reading number of tiles
      ! -----------------------
      
      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
      
      read (20, *) NTILES
      
      close (20, status = 'keep')

      ! READ CLM4.5 source data files and regrid
      ! ----------------------------------------

      status  = NF_OPEN ('data/CATCH/CLM45/clmforc.Li_2012_hdm_0.5x0.5_AVHRR_simyr1850-2010_c130401.nc', NF_NOWRITE, ncid_hdm  )
      status  = NF_OPEN ('data/CATCH/CLM45/mksrf_abm_0.5x0.5_AVHRR_simyr2000.c130201.nc'               , NF_NOWRITE, ncid_abm  )
      status  = NF_OPEN ('data/CATCH/CLM45/mksrf_gdp_0.5x0.5_AVHRR_simyr2000.c130228.nc'               , NF_NOWRITE, ncid_gdp  )
      status  = NF_OPEN ('data/CATCH/CLM45/mksrf_peatf_0.5x0.5_AVHRR_simyr2000.c130228.nc'             , NF_NOWRITE, ncid_peatf)
            
      allocate (hdm_grid   (1:NC,1:NR))
      allocate (abm_grid   (1:NC,1:NR))
      allocate (gdp_grid   (1:NC,1:NR))
      allocate (peatf_grid (1:NC,1:NR))
      allocate (data_grid (1 : N_lon_clm, 1 : N_lat_clm)) 
      allocate (int_grid  (1 : N_lon_clm, 1 : N_lat_clm)) 

      status  = NF_INQ_VARID (ncid_hdm,'hdm',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_REAL (ncid_hdm,VarID,(/1,1,161/),(/N_lon_clm, N_lat_clm, 1/),data_grid(:,:)) ; VERIFY_(STATUS)
      call RegridRasterReal(data_grid, hdm_grid)

      status  = NF_INQ_VARID (ncid_abm,'abm',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_INT (ncid_abm,VarID, (/1,1/),(/N_lon_clm, N_lat_clm/), int_grid) ; VERIFY_(STATUS)
      call RegridRaster (int_grid, abm_grid)

      status  = NF_INQ_VARID (ncid_gdp,'gdp',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_REAL (ncid_gdp,VarID, (/1,1/),(/N_lon_clm, N_lat_clm/), data_grid) ; VERIFY_(STATUS)
      call RegridRasterReal(data_grid, gdp_grid)

      status  = NF_INQ_VARID (ncid_peatf,'peatf',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_REAL (ncid_peatf,VarID, (/1,1/),(/N_lon_clm, N_lat_clm/), data_grid) ; VERIFY_(STATUS)
      call RegridRasterReal(data_grid, peatf_grid)

      status = NF_CLOSE(ncid_hdm  )
      status = NF_CLOSE(ncid_abm  )
      status = NF_CLOSE(ncid_gdp  )
      status = NF_CLOSE(ncid_peatf)

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do       

      close (10,status='keep')     
      
      allocate (hdm   (1:NTILES))
      allocate (abm   (1:NTILES))
      allocate (gdp   (1:NTILES))
      allocate (peatf (1:NTILES))
      allocate (count_pix (1:NTILES, 1:4))
      
      hdm       = 0.
      abm       = 0.
      gdp       = 0.
      peatf     = 0.
      count_pix = 0.

      do j = 1,nr
         do i = 1, nc
            if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

               ! peatf 0. < 1.
               if((peatf_grid(i,j) >= 0.).and.(peatf_grid(i,j) <= 1.)) then 
                 peatf (tile_id(i,j)) = peatf (tile_id(i,j)) + peatf_grid(i,j)
                 count_pix (tile_id(i,j), 1) = count_pix (tile_id(i,j), 1) + 1. 
               endif

               ! gdp 0. < 300.
               if((gdp_grid(i,j) >= 0.).and.(gdp_grid(i,j) <= 300.)) then 
                 gdp (tile_id(i,j)) = gdp (tile_id(i,j)) + gdp_grid(i,j)
                 count_pix (tile_id(i,j), 2) = count_pix (tile_id(i,j), 2) + 1. 
               endif

               ! abm 1 < 12
               if((abm_grid(i,j) >= 1).and.(abm_grid(i,j) <= 12)) then 
                 abm (tile_id(i,j)) = abm (tile_id(i,j)) + abm_grid(i,j)
                 count_pix (tile_id(i,j), 3) = count_pix (tile_id(i,j), 3) + 1. 
               endif

               ! hdm 0. < 20000.
               if((hdm_grid(i,j) >= 0.).and.(hdm_grid(i,j) <= 20000.)) then 
                 hdm (tile_id(i,j)) = hdm (tile_id(i,j)) + hdm_grid(i,j)
                 count_pix (tile_id(i,j), 4) = count_pix (tile_id(i,j), 4) + 1. 
               endif            
            endif
         end do
      end do
    
    ! Field Capacity
    ! --------------

    open (11, file='data/CATCH/SoilClasses-SoilHyd-TauParam.dat', form='formatted',status='old', &
           action = 'read')
    read (11,'(a)')fout
    do i =1,n_SoilClasses 
       read (11,'(4f7.3,4f8.4,e13.5,2f12.7,2f8.4,4f12.7)')a_sand,a_clay,a_silt,a_oc,a_bee,a_psis, &
            a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf,atau_2cm,btau_2cm, field_cap(i)       
    end do
    close (11,status='keep') 


    open (10,file='clsm/CLM4.5_abm_peatf_gdp_hdm_fc',  &
         form='formatted',status='unknown')     

    open (20, file = 'clsm/soil_param.dat', form = 'formatted', status = 'old', action =  'read')

    do i = 1, NTILES

       hdm_r   = 0.
       gdp_r   = 0.
       peatf_r = 0.
       abm_int = 7
       read (20,*) tid, cid, sc_top, sc_com

       if(count_pix(i,1) > 0.) peatf_r = peatf (i) / count_pix(i,1)
       if(count_pix(i,2) > 0.) gdp_r   = gdp   (i) / count_pix(i,2)
       if(count_pix(i,3) > 0.) abm_int = NINT(abm  (i) / count_pix(i,3))
       if(count_pix(i,4) > 0.) hdm_r   = hdm   (i) / count_pix(i,4)
             
       write (10,'(2I10, i3, f8.4, f8.2, f10.2, f8.4)' ) tid, cid, abm_int, peatf_r, gdp_r, hdm_r, field_cap(sc_com)  

    end do

    deallocate (hdm, abm, gdp, peatf)
    deallocate (hdm_grid, gdp_grid, peatf_grid, data_grid, count_pix)
    deallocate (tile_id, abm_grid)

    close (10, status = 'keep')
    close (20, status = 'keep')

  END SUBROUTINE CLM45_fixed_parameters

   ! ----------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE CLM45_clim_parameters (nc,nr,gfiler)

      implicit none
      ! Producing :  lightening frequency HRMC_COM_FR /gpfsm/dnb31/fzeng/clm4-to-clm4.5/data/firedata4.5/LISOTD_HRMC_V2.3.2014.hdf
      !  12 values per tile
      integer, intent (in)                   :: nc, nr 
      character(*), intent (in)              :: gfiler
      integer  , parameter                   :: N_lon_clm = 720, N_lat_clm = 360
      integer                                :: NTILES, status, varid, ncid
      real, dimension (:,:), allocatable     :: hrmc_grid, data_grid
      REAL, ALLOCATABLE, dimension (:)       :: hrmc, count_pix
      INTEGER, ALLOCATABLE, dimension (:,:)  :: tile_id
      integer                                :: yr,mn,yr1,mn1, k,t,i,j

     ! Reading number of tiles
      ! -----------------------
      
      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
      
      read (20, *) NTILES
      
      close (20, status = 'keep')

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file

      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do       

      close (10,status='keep')     
      
      ! READ CLM4.5 source data files and regrid
      ! ----------------------------------------

      status  = NF_OPEN ('data/CATCH/CLM45/LISOTD_HRMC_V2.3.2014.nc4', NF_NOWRITE, ncid)
      status  = NF_INQ_VARID (ncid,'HRMC_COM_FR',VarID) ; VERIFY_(STATUS)

      allocate (hrmc_grid   (1:NC,1:NR))
      allocate (data_grid (1 : N_lon_clm, 1 : N_lat_clm)) 
      allocate (hrmc (1:NTILES))
      allocate (count_pix (1:NTILES))

      ! writing tile-spaced output
      ! --------------------------

      open (31,file='clsm/lnfm.dat',status='unknown',action='write',form='unformatted', &
           convert='little_endian')
      
      do K=0,13
         yr = (k+11)/12
         mn = mod(k+11,12)+1
         yr1= (k+12)/12
         mn1= mod(k+12,12)+1
         write(31) float((/yr,mn,1,0,0,0,yr1,mn1,1,0,0,0,NTILES,1/))
         hrmc = 0.
         count_pix = 0.
         t = k
         if (t == 0 ) t = 12
         if (t == 13) t = 1
         status  = NF_GET_VARA_REAL (ncid,VarID, (/1,1,t/),(/N_lon_clm, N_lat_clm,1/), data_grid) ; VERIFY_(STATUS)
         call RegridRasterReal(data_grid, hrmc_grid)

         do j = 1,nr
            do i = 1, nc
               if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                  if((hrmc_grid(i,j) >= 0.).and.(hrmc_grid(i,j) <= 1.)) then 
                     hrmc (tile_id(i,j)) = hrmc (tile_id(i,j)) + hrmc_grid(i,j)
                     count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                  endif
               endif
            end do
         end do

         where (count_pix > 0.) hrmc = hrmc /count_pix
         write(31) hrmc
      end do

      close(31,status='keep')    

    END SUBROUTINE CLM45_clim_parameters

! ----------------------------------------------------------------------------------------------------------------------------

  SUBROUTINE grid2tile_glass (ncol,nrow,gfiler,lai_name)
!
! Processing GLASS LAI (AVHRR or MODIS) and creating 8-day climatological data 
!
  implicit none 
  integer  , parameter                   :: N_lon_glass = 7200, N_lat_glass = 3600
  integer, intent (in) :: ncol, nrow
  real, parameter :: dxy = 1.
  integer :: QSize
  character(*)  :: gfiler,lai_name
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj,ierr,nx,ny
  integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
      time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2,tindex1,pfaf1
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, target,  dimension (:,:) :: net_data1
  real,  pointer, dimension (:,:) :: QSub
  real,  pointer, dimension (:,:) :: subset
  REAL, ALLOCATABLE, dimension (:):: vec_lai, count_lai,tile_lon, tile_lat &
       , x, y !, distance
  real, allocatable, target, dimension (:,:) :: lai_grid, data_grid, data_grid2
  INTEGER ::imn,imx,jmn,jmx,mval,d1,d2,l, VarID
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  REAL :: sf, dum,dist_save,tile_distance,minlat,maxlat,minlon,maxlon
  logical :: first_entry = .true.
  type (date_time_type) :: date_time_new,bf_lai_time,   &
       af_lai_time, date_time_this
  integer, dimension (:,:), allocatable, target :: tile_id
  integer       ::  tileid_tile
  character*3   :: ddd

! Reading rst file
!-----------------
   open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
        form='unformatted',convert='little_endian')
   allocate (tile_id    (1:ncol,1:nrow))         
   
   do j=1,nrow
      read(10)tile_id(:,j)
   end do
   close (10,status='keep')

!
! Reading number of cathment-tiles from catchment.def file
!_________________________________________________________ 
!
      fname='clsm/catchment.def' 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      allocate (tile_lon(1:maxcat)) 
      allocate (tile_lat(1:maxcat)) 
      
      do n = 1, maxcat
         read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
         tile_lon(n) = (minlon + maxlon)/2.
         tile_lat(n) = (minlat + maxlat)/2.
      end do      
      close (10,status='keep')

      fname =trim(c_data)//'/MODIS_8-DayClim/MODIS_lai_clim.H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,3,string, n_tslices); VERIFY_(STATUS) 
      allocate (MMDD      (0: n_tslices + 1))
      allocate (MMDD_next (0: n_tslices + 1))

      status = NF_GET_VARA_text(ncid, 3,(/1,1/),(/4,n_tslices/),MMDD(1:n_tslices)); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)
       
      mmdd(0) = mmdd(n_tslices)
      mmdd(n_tslices + 1)= mmdd(1)
      
      mmdd_next(0:n_tslices - 1) =  mmdd(1:n_tslices)
      mmdd_next(n_tslices: n_tslices + 1) = mmdd (1:2)
     
      ! writing GLASS LAI
      !
      open (31,file='clsm/lai.dat',  &
           form='unformatted',status='unknown',convert='little_endian')
     
      allocate (vec_lai     (maxcat))
      allocate (count_lai (1:maxcat))

      nx = nint (360./dxy)
      ny = nint (180./dxy)
      allocate (x(1:nx))
      allocate (y(1:ny))

      FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
      FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy

      allocate (lai_grid (1 : nx, 1 : ny)) 
      
      QSize = nint(dxy*N_lon_glass/360.)
      allocate (QSub (1:QSize,1:QSize))
      allocate (net_data1 (1 : N_lon_glass, 1 : N_lat_glass)) 
      allocate (data_grid (1:NCOL,1:NROW))
      allocate (data_grid2 (1 : N_lon_glass, 1 : N_lat_glass)) 

      do t =0,n_tslices+1
         
         time_slice = t
         yr = 1
         yr1= 1
         if(t == 0) then
            time_slice =  n_tslices
            yr         =  1 - 1
         endif
         
         if(t >= n_tslices) then 
            yr1 = 1 + 1
            if(t ==n_tslices + 1) then
               time_slice =  1
               yr = 1 + 1
            endif
         endif
         
         read(mmdd(t),'(i2.2,i2.2)') mn,dd
         read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1
         
         date_time_this%year       = 2001
         date_time_this%month      = mn
         date_time_this%day        = dd
         date_time_this%hour   = 0            
         date_time_this%min    = 0            
         date_time_this%sec    = 0 
         call get_dofyr_pentad(date_time_this)                      

         write (ddd,'(i3.3)')  date_time_this%dofyr

         ! Reading Interpolation or aggregation on to catchment-tiles
         
         vec_lai   = -9999.
         count_lai = 0.
         lai_grid  = -9999

         status  = NF_OPEN (trim(c_data)//trim(lai_name)//ddd//'.nc4', NF_NOWRITE, ncid) ; VERIFY_(STATUS)
         status  = NF_INQ_VARID (ncid,'LAI',VarID) ; VERIFY_(STATUS)
         status  = NF_GET_VARA_INT(ncid,VarID, (/1,1/),(/N_lon_glass, N_lat_glass/), net_data1) ; VERIFY_(STATUS)

         call RegridRasterReal(0.01*real(net_data1), data_grid)
         data_grid2 = 0.01*real(net_data1)

         status = NF_CLOSE(ncid)

         do j = 1,nrow
            do i = 1, ncol
               if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.MAXCAT)) then        
                  if((data_grid(i,j) >= 0.).and.(data_grid(i,j) <= 10.)) then 
                     if(vec_lai(tile_id(i,j)) == -9999.) vec_lai(tile_id(i,j)) = 0. 
                     vec_lai (tile_id(i,j)) = vec_lai (tile_id(i,j)) + data_grid(i,j)
                     count_lai (tile_id(i,j)) = count_lai (tile_id(i,j)) + 1.                     
                  endif
               endif
            end do
         end do

         write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
         
         where (count_lai > 0.) vec_lai = vec_lai/count_lai

         ! After experimenting with few finer methods, in order to reduce the time taken by the gap filling procedure,
         ! creating a 0.25-degree gridded data set from finer LAI data and use it for filling the gaps seems the most practical/manageble method.
         !---------------------------------------------------------------------------------------------------------------------------------------

         iLL = 1
         jLL = 1
         do j = 1, N_lat_glass/QSize
            do i = 1,  N_lon_glass/QSize 
               QSub  => data_grid2((i-1)*QSize+2-iLL :i*QSize-iLL+1, (j-1)*QSize+2-jLL :j*QSize-jLL+1) 
               if(minval (QSub) <= 10.) lai_grid(i,j) = sum(QSub, QSub<=10.)/(max(1,count(QSub<=10.)))
            enddo
         enddo
                  
         NULLIFY (QSub)

! Filling gaps
!-------------
         DO n =1,maxcat
             if(count_lai(n)==0.)  then 
                
                DO i = 1,nx - 1
                   if ((tile_lon(n) >= x(i)).and.(tile_lon(n) < x(i+1))) ix = i
                end do
                DO i = 1,ny -1
                   if ((tile_lat(n) >= y(i)).and.(tile_lat(n) < y(i+1))) jx = i
                end do
                
                l = 1
                do 
                  imx=ix + l
                  imn=ix - l
                  jmn=jx - l
                  jmx=jx + l
                  imn=MAX(imn,1)
                  jmn=MAX(jmn,1)
                  imx=MIN(imx,nx)
                  jmx=MIN(jmx,ny)
                  d1=imx-imn+1
                  d2=jmx-jmn+1
                  subset => lai_grid(imn: imx,jmn:jmx)

                  if(maxval(subset) > 0.) then 
                     vec_lai (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                     exit
                  endif
                  l = l + 1
                  NULLIFY (subset)
                end do
             endif
          END DO
          write(31)  vec_lai(:)
       end do
       close(31,status='keep')
       
       deallocate (net_data1, tile_id)
       deallocate (count_lai)
       deallocate (vec_lai) 
       deallocate (tile_lat,tile_lon)

     END SUBROUTINE grid2tile_glass

   ! ----------------------------------------------------------------------------------------------------------------------------

    SUBROUTINE gimms_clim_ndvi (nc,nr,gfiler)
      
      implicit none
      ! Producing :  GIMMS NDVI 15-day climatology from 5 arcmin data
      !  24 values per tile
      integer, intent (in)                   :: nc, nr 
      character(*), intent (in)              :: gfiler
      integer  , parameter                   :: N_lon_gimms = 4320, N_lat_gimms = 2160
      integer                                :: NTILES, status, varid, ncid1, ncid2,ncid
      real, dimension (:,:), allocatable     :: ndvi_grid, data_grid 
      integer, dimension (:,:), allocatable     ::data_grid2
      REAL, ALLOCATABLE, dimension (:)       :: ndvi, count_pix
      INTEGER, ALLOCATABLE, dimension (:,:)  :: tile_id
      integer                                :: yr,mn,yr1,mn1, k,t,i,j,l
      integer, parameter :: scale_fac = 10000
      real,    parameter :: val_min = -0.3, val_max = 1.

     ! Reading number of tiles
      ! -----------------------
      
      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', action =  'read')
      
      read (20, *) NTILES
      
      close (20, status = 'keep')

      ! Grid to tile
      ! ------------
                
      ! Reading tile-id raster file
      
      allocate(tile_id(1:nc,1:nr))
      
      open (10,file=trim(gfiler)//'.rst',status='old',action='read',  &
           form='unformatted',convert='little_endian')
      
      do j=1,nr
         read(10)tile_id(:,j)
      end do
      
      close (10,status='keep')     
      
      ! READ GIMMS NDVI source data files and regrid
      ! ----------------------------------------
      
      status  = NF_OPEN ('data/CATCH/ndvi3g_geo_v1_YYYY_0106.nc4', NF_NOWRITE, ncid1) ; VERIFY_(STATUS)
      status  = NF_OPEN ('data/CATCH/ndvi3g_geo_v1_YYYY_0712.nc4', NF_NOWRITE, ncid2) ; VERIFY_(STATUS)
      status  = NF_INQ_VARID (ncid2,'ndvi',VarID) ; VERIFY_(STATUS)
      
      allocate (ndvi_grid   (1:NC,1:NR))
      allocate (data_grid (1 : N_lon_gimms, 1 : N_lat_gimms)) 
      allocate (data_grid2(1 : N_lon_gimms, 1 : N_lat_gimms)) 
      allocate (ndvi (1:NTILES))
      allocate (count_pix (1:NTILES))
      
      ! writing tile-spaced output
      ! --------------------------
      
      open (31,file='clsm/ndvi.dat',status='unknown',action='write',form='unformatted', &
           convert='little_endian')
      
      do K=0,13
         yr = (k+11)/12
         mn = mod(k+11,12)+1
         yr1= (k+12)/12
         mn1= mod(k+12,12)+1

         ndvi = 0.
         count_pix = 0.
         t = k
         if (k == 0 ) then 
            t = 12
            ncid = ncid2
            write(31) float((/yr,mn,16,0,0,0,yr1,mn1,1,0,0,0,NTILES,1/))

            status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,t/),(/N_lon_gimms, N_lat_gimms,1/), data_grid2) ; VERIFY_(STATUS)

            do j = 1,  N_lat_gimms
               data_grid (:,j) =   data_grid2 (:,N_lat_gimms - (j-1)) / real(scale_fac)
            end do
            
            call RegridRasterReal(data_grid, ndvi_grid)
            
            do j = 1,nr
               do i = 1, nc
                  if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                     if((ndvi_grid(i,j) >= val_min).and.(ndvi_grid(i,j) <= val_max)) then 
                        ndvi (tile_id(i,j)) = ndvi (tile_id(i,j)) + ndvi_grid(i,j)
                        count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                     endif
                  endif
               end do
            end do
            
            where (count_pix > 0.) ndvi = ndvi /count_pix
            write(31) ndvi
            
         elseif (k == 13) then

            t = 1
            ncid = ncid1
            write(31) float((/yr,mn,1,0,0,0,yr,mn,16,0,0,0,NTILES,1/))

            status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,t/),(/N_lon_gimms, N_lat_gimms,1/), data_grid2) ; VERIFY_(STATUS)

            do j = 1, N_lat_gimms
               data_grid (:,j) = data_grid2 (:,N_lat_gimms - (j-1)) / real(scale_fac)
            end do
            
            call RegridRasterReal(data_grid, ndvi_grid)
            
            do j = 1,nr
               do i = 1, nc
                  if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                     if((ndvi_grid(i,j) >= val_min).and.(ndvi_grid(i,j) <= val_max)) then 
                        ndvi (tile_id(i,j)) = ndvi (tile_id(i,j)) + ndvi_grid(i,j)
                        count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                     endif
                  endif
               end do
            end do
            
            where (count_pix > 0.) ndvi = ndvi /count_pix
            write(31) ndvi

         else

            do l = 1, 0 , -1
               t = k*2 - l
               if (k <= 6) ncid = ncid1
               if (k >= 7) ncid = ncid2
               if (k >= 7) t = t - 12
               if(l == 1) write(31) float((/yr,mn,1,0,0,0,yr,mn,16,0,0,0,NTILES,1/))
               if(l == 0) write(31) float((/yr,mn,16,0,0,0,yr1,mn1,1,0,0,0,NTILES,1/))

               ndvi = 0.
               count_pix = 0.

               status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,t/),(/N_lon_gimms, N_lat_gimms,1/), data_grid2) ; VERIFY_(STATUS)

               do j = 1,  N_lat_gimms
                  data_grid (:,j) =   data_grid2 (:,N_lat_gimms - (j-1)) / real(scale_fac)
               end do

               call RegridRasterReal(data_grid, ndvi_grid)
               
               do j = 1,nr
                  do i = 1, nc
                     if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then        
                        if((ndvi_grid(i,j) >= val_min).and.(ndvi_grid(i,j) <= val_max)) then 
                           ndvi (tile_id(i,j)) = ndvi (tile_id(i,j)) + ndvi_grid(i,j)
                           count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1.                     
                        endif
                     endif
                  end do
               end do
               
               where (count_pix > 0.) ndvi = ndvi /count_pix
               write(31) ndvi
               
            end do            
         endif
      end do

      close(31,status='keep')    

    END SUBROUTINE gimms_clim_ndvi

    ! --------------------------------------------------------------------------

    SUBROUTINE open_landparam_nc4_files(N_tile) 

      implicit none
      integer                 :: NCCatOUTID,  NCCatCNOUTID,  NCVegOUTID  
      integer                 :: STATUS, CellID1, CellID2, CellID3, SubID
      integer, intent (in)    :: N_tile
      integer, dimension(8)   :: date_time_values
      character (22)          :: time_stamp
      character (100)         :: MYNAME

      status = NF_CREATE ('clsm/catch_params.nc4'  , NF_NETCDF4, NCCatOUTID  ) ; VERIFY_(STATUS)
      status = NF_CREATE ('clsm/catchcn_params.nc4', NF_NETCDF4, NCCatCNOUTID) ; VERIFY_(STATUS)
      status = NF_CREATE ('clsm/vegdyn.data'       , NF_NETCDF4, NCVegOUTID  ) ; VERIFY_(STATUS)

      status = NF_DEF_DIM(NCCatOUTID  , 'tile' , N_tile, CellID1)
      status = NF_DEF_DIM(NCCatCNOUTID, 'tile' , N_tile, CellID2)
      status = NF_DEF_DIM(NCVegOUTID  , 'tile' , N_tile, CellID3)
      status = NF_DEF_DIM(NCCatCNOUTID, 'unknown_dim2' , 4, SubID)

      call DEF_VAR ( NCCatOUTID, CellID1,'OLD_ITY'   ,'vegetation_type.'            , '1'       )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARA1'      ,'shape_param_1'               ,'m+2 kg-1' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARA2'      ,'shape_param_2'               ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARA3'      ,'shape_param_3'               ,'m+2 kg-1' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARA4'      ,'shape_param_4'               ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARS1'      ,'wetness_param_1'             ,'m+2 kg-1' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARS2'      ,'wetness_param_2'             ,'m+2 kg-1' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARS3'      ,'wetness_param_3'             ,'m+4 kg-2' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARW1'      ,'min_theta_param_1'           ,'m+2 kg-1' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARW2'      ,'min_theta_param_2'           ,'m+2 kg-1' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARW3'      ,'min_theta_param_3'           ,'m+4 kg-2' )
      call DEF_VAR ( NCCatOUTID, CellID1,'ARW4'      ,'min_theta_param_4'           ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'ATAU2'     ,'2cm_water_transfer_param_5'  ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'ATAU5'     ,'5cm_water_transfer_param_5'  ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'BEE'       ,'clapp_hornberger_b'          ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'BF1'       ,'topo_baseflow_param_1'       ,'kg m-4'   )
      call DEF_VAR ( NCCatOUTID, CellID1,'BF2'       ,'topo_baseflow_param_2'       ,'m'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'BF3'       ,'topo_baseflow_param_3'       ,'log(m)'   )
      call DEF_VAR ( NCCatOUTID, CellID1,'BTAU2'     ,'2cm_water_transfer_param_6'  ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'BTAU5'     ,'5cm_water_transfer_param_6'  ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'COND'      ,'sfc_sat_hydraulic_conduct'   ,'m s-1'    )
      call DEF_VAR ( NCCatOUTID, CellID1,'GNU'       ,'vertical_transmissivity'     ,'m-1'      )
      call DEF_VAR ( NCCatOUTID, CellID1,'POROS'     ,'soil_porosity'               ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'PSIS'      ,'saturated_matric_potential'  ,'m'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'TSA1'      ,'water_transfer_param_1'      ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'TSA2'      ,'water_transfer_param_2'      ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'TSB1'      ,'water_transfer_param_3'      ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'TSB2'      ,'water_transfer_param_4'      ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'WPWET'     ,'wetness_at_wilting_point'    ,'1'        )
      call DEF_VAR ( NCCatOUTID, CellID1,'DP2BR'     ,'depth_to_bedrock'            ,'mm'       )

      call DEF_VAR (  NCVegOUTID, CellID3,'ITY'      ,'vegetation_type'             ,'1'        )
      call DEF_VAR (  NCVegOUTID, CellID3,'Z2CH'     ,'vegetation_height'           ,'m'        )
      call DEF_VAR (  NCVegOUTID, CellID3,'ASCATZ0'  ,'ASCAT_roughness_length'      ,'m'        )

      call DEF_VAR ( NCCatCNOUTID, CellID2,'BGALBNF' ,'MODIS soil albedo nir dif'   ,'1'        )
      call DEF_VAR ( NCCatCNOUTID, CellID2,'BGALBNR' ,'MODIS soil albedo nir dir'   ,'1'        )
      call DEF_VAR ( NCCatCNOUTID, CellID2,'BGALBVF' ,'MODIS soil albedo vis dif'   ,'1'        )
      call DEF_VAR ( NCCatCNOUTID, CellID2,'BGALBVR' ,'MODIS soil albedo vis dir'   ,'1'        )
      call DEF_VAR ( NCCatCNOUTID, CellID2,'T2_M'    ,'Clim 2m temperature (MERRA2)'    ,'K'    )
      call DEF_VAR ( NCCatCNOUTID, CellID2,'T2_S'    ,'Clim 2m temperature (Sheffield)' ,'K'    )
      call DEF_VAR ( NCCatCNOUTID, CellID2,'NDEP'    ,'CLM_nitrogen_deposition'     ,'g m-2 s-1')
      call DEF_VAR ( NCCatCNOUTID, CellID2,'FVG'     ,'vegetation_fraction'         ,'1'        ,SubID = SubID)
      call DEF_VAR ( NCCatCNOUTID, CellID2,'ITY'     ,'vegetation_type'             ,'1'        ,SubID = SubID)

      call date_and_time(VALUES=date_time_values)
          
      write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
           date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
           date_time_values(5),':',date_time_values(6),':',date_time_values(7)
!      call execute_command_line('setenv    MYNAME `finger $USER | cut -d: -f3 | head -1`')
!      call sleep (5)
      call get_environment_variable ("USER"        ,MYNAME        )
      status = NF_PUT_ATT_TEXT(NCCatOUTID  , NF_GLOBAL, 'CreatedBy', LEN_TRIM(MYNAME),  trim(MYNAME)      )
      status = NF_PUT_ATT_TEXT(NCCatOUTID  , NF_GLOBAL, 'Date'     , LEN_TRIM(time_stamp),trim(time_stamp))
      status = NF_PUT_ATT_TEXT(NCVegOUTID  , NF_GLOBAL, 'CreatedBy', LEN_TRIM(MYNAME),  trim(MYNAME)      )
      status = NF_PUT_ATT_TEXT(NCVegOUTID  , NF_GLOBAL, 'Date'     , LEN_TRIM(time_stamp),trim(time_stamp))
      status = NF_PUT_ATT_TEXT(NCCatCNOUTID, NF_GLOBAL, 'CreatedBy', LEN_TRIM(MYNAME),  trim(MYNAME)      )
      status = NF_PUT_ATT_TEXT(NCCatCNOUTID, NF_GLOBAL, 'Date'     , LEN_TRIM(time_stamp),trim(time_stamp))
      
      status = NF_ENDDEF(NCCatOUTID  )  
      status = NF_ENDDEF(NCVegOUTID  )  
      status = NF_ENDDEF(NCCatCNOUTID)  

      status    = NF_CLOSE (NCCatOUTID  )  
      status    = NF_CLOSE (NCVegOUTID  )  
      status    = NF_CLOSE (NCCatCNOUTID)  

    contains
  
      SUBROUTINE DEF_VAR (NCFID, CellID, VarName, long_name, units, SubID)
        
        implicit none
        integer, intent (in)           :: NCFID, CellID    
        character (*), intent (in)     :: VarName, long_name, units
        integer, intent (in), optional :: SubID
        integer                        :: STATUS, VID

        if(present (SubID)) then 
           status = NF_DEF_VAR(NCFID, trim(VarName) , NF_FLOAT, 2 ,(/CellID, SubID/), vid) ; VERIFY_(STATUS)
        else
           status = NF_DEF_VAR(NCFID, trim(VarName) , NF_FLOAT, 1 ,(/CellID/), vid) ; VERIFY_(STATUS)
        endif

        status = NF_PUT_ATT_TEXT(NCFID, vid, 'long_name', LEN_TRIM(long_name), trim(long_name)) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCFID, vid, 'units'    , LEN_TRIM(units)    , trim(units))     ; VERIFY_(STATUS)


      END SUBROUTINE DEF_VAR

    END SUBROUTINE open_landparam_nc4_files

    ! ----------------------------------------------------------------------------------------------


    SUBROUTINE map_country_codes (NC, NR,  gfiler)

      implicit none
      integer  , intent (in) :: nc, nr
      character (*)          :: gfiler

      integer, parameter :: GC = 43200
      integer, parameter :: GR = 21600
      INTEGER,      dimension (:), pointer :: index_RANGE 
      character*20, dimension (:), pointer :: ST_NAME     
      character*48, dimension (:), pointer :: CNT_NAME  
      
      integer :: CNT_CODE, ST_CODE
      integer :: i(GC),j(GR), k,n, status, ncid, varid, maxcat, I0(1), j0(1)
      INTEGER, TARGET, ALLOCATABLE, dimension (:,:):: ST_grid, cnt_grid
      real    :: lat_mn, lat_mx, lon_mn, lon_mx
      real (kind =8) :: XG(GC),YG(GR), y0, x0, dxy

      call get_country_codes (index_RANGE = index_RANGE, ST_NAME = ST_NAME, &
           CNT_NAME = CNT_NAME)
      
      ! Reading number of tiles
      ! -----------------------

      open (20, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
           action =  'read')
      
      read (20, *) maxcat
      

      ! READ country code source data files and regrid
      ! -----------------------------------------
      
      status  = NF_OPEN ('data/CATCH/GADM_Country_and_USStates_codes_1km.nc4', NF_NOWRITE, ncid)
      
      allocate (cnt_grid  (1 : GC, 1 : GR))
      allocate (st_grid   (1 : GC, 1 : GR))
      
      status  = NF_INQ_VARID (ncid,'UNIT_CODE',VarID) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,1/),(/GC, GR,1/), cnt_grid) ; VERIFY_(STATUS)
      status  = NF_GET_VARA_INT (ncid,VarID, (/1,1,2/),(/GC, GR,1/), st_grid) ; VERIFY_(STATUS)
      where (st_grid == 0) st_grid = 999
      status = NF_CLOSE(ncid)

      open (10,file='clsm/country_and_state_code.data',  &
         form='formatted',status='unknown')

      dxy = 360./GC
      do k = 1, GC 
         xg(k) = (k-1)*dxy -180. + dxy/2.
      end do
      do k = 1, GR 
         yg(k) = (k-1)*dxy -90. + dxy/2.
      end do      

      DO n = 1, MAXCAT
         read (20,*) i0,j0, lon_mn, lon_mx, lat_mn, lat_mx
         x0 = (lon_mn + lon_mx)/2.
         y0 = (lat_mn + lat_mx)/2.
         I = 0
         J = 0
         WHERE ((xg >= x0).and.(xg < x0 + dxy)) I = 1
         WHERE ((yg >= y0).and.(yg < y0 + dxy)) J = 1
         
         I0 =FINDLOC(I,1)
         J0 =FINDLOC(J,1)

         cnt_code = cnt_grid(I0(1), J0(1))
         st_code  = st_grid (I0(1), J0(1))

         if(cnt_code > 300) then
            CNT_CODE = 257
         endif

         if(st_code <= 50) then
            write (10, '(i10, 2I4, 1x, a48, a20)') n, cnt_code, st_code, CNT_NAME(FINDLOC(INDEX_RANGE, CNT_CODE)), ST_NAME (ST_CODE)
         else
            write (10, '(i10, 2I4, 1x, a48, a20)') n, cnt_code, st_code, CNT_NAME(FINDLOC(INDEX_RANGE, CNT_CODE)), 'OUTSIDE USA'
         endif
      
     END DO

     close (10, status = 'keep')
     close (20, status = 'keep')
   END SUBROUTINE map_country_codes

   ! -------------------------------------------------------------------------------------------
   
   SUBROUTINE get_country_codes (index_RANGE, ST_NAME, CNT_NAME, ST_NAME_ABR, CNT_NAME_ABR)

     implicit none
     
     INTEGER,      dimension (N_GADM  ), TARGET :: index_RANGE_DATA
     character*20, dimension (N_STATES), TARGET :: ST_NAME_DATA
     character*48, dimension (N_GADM  ), TARGET :: CNT_NAME_DATA
     INTEGER,      dimension (:), pointer, intent (inout), optional :: index_RANGE 
     character*20, dimension (:), pointer, intent (inout), optional :: ST_NAME     
     character*48, dimension (:), pointer, intent (inout), optional :: CNT_NAME    
     character*2,  dimension (:), pointer, intent (inout), optional :: ST_NAME_ABR 
     character*3,  dimension (:), pointer, intent (inout), optional :: CNT_NAME_ABR

     DATA ST_NAME_DATA /            &
          'AK  1 Alaska          ' ,&
          'AL  2 Alabama         ' ,&
          'AZ  3 Arizona         ' ,&
          'AR  4 Arkansas        ' ,&
          'CA  5 California      ' ,&
          'CO  6 Colorado        ' ,&
          'CT  7 Connecticut     ' ,&
          'DE  8 Delaware        ' ,&
          'FL  9 Florida         ' ,&
          'GA 10 Georgia         ' ,&
          'HI 11 Hawaii          ' ,&
          'IA 12 Iowa            ' ,&
          'ID 13 Idaho           ' ,&
          'IL 14 Illinois        ' ,&
          'IN 15 Indiana         ' ,&
          'KS 16 Kansas          ' ,&
          'KY 17 Kentucky        ' ,&
          'LA 18 Louisiana       ' ,&
          'MA 19 Massachusetts   ' ,&
          'MD 20 Maryland        ' ,&
          'ME 21 Maine           ' ,&
          'MI 22 Michigan        ' ,&
          'MN 23 Minnesota       ' ,&
          'MO 24 Missouri        ' ,&
          'MS 25 Mississippi     ' ,&
          'MT 26 Montana         ' ,&
          'NC 27 NorthCarolina   ' ,&
          'ND 28 NorthDakota     ' ,&
          'NE 29 Nebraska        ' ,&
          'NH 30 NewHampshire    ' ,&
          'NJ 31 NewJersey       ' ,&
          'NM 32 NewMexico       ' ,&
          'NV 33 Nevada          ' ,&
          'NY 34 NewYork         ' ,&
          'OH 35 Ohio            ' ,&
          'OK 36 Oklahoma        ' ,&
          'OR 37 Oregon          ' ,&
          'PA 38 Pennsylvania    ' ,&
          'RI 39 RhodeIsland     ' ,&
          'SC 40 SouthCarolina   ' ,&
          'SD 41 SouthDakota     ' ,&
          'TN 42 Tennessee       ' ,&
          'TX 43 Texas           ' ,&
          'UT 44 Utah            ' ,&
          'VA 45 Virginia        ' ,&
          'VT 46 Vermont         ' ,&
          'WA 47 Washington      ' ,&
          'WI 48 Wisconsin       ' ,&
          'WV 49 WestVirginia    ' ,&
          'WY 50 Wyoming         ' /
     
     DATA CNT_NAME_DATA  /                                     & 
          'ABW  14 Aruba                                   '  ,&
          'AFG   1 Afghanistan                             '  ,&
          'AGO   8 Angola                                  '  ,&
          'AIA   9 Anguilla                                '  ,&
          'ALA   3 Aland                                   '  ,&
          'ALB   4 Albania                                 '  ,&
          'AND   7 Andorra                                 '  ,&
          'ARE 241 United Arab Emirates                    '  ,&
          'ARG  12 Argentina                               '  ,&
          'ARM  13 Armenia                                 '  ,&
          'ASM   6 American Samoa                          '  ,&
          'ATA  10 Antarctica                              '  ,&
          'ATF  82 French Southern Territories             '  ,&
          'ATG  11 Antigua and Barbuda                     '  ,&
          'AUS  15 Australia                               '  ,&
          'AUT  16 Austria                                 '  ,&
          'AZE  17 Azerbaijan                              '  ,&
          'BDI  39 Burundi                                 '  ,&
          'BEL  23 Belgium                                 '  ,&
          'BEN  25 Benin                                   '  ,&
          'BES  29 Bonaire, Sint Eustatius and Saba        '  ,&
          'BFA  38 Burkina Faso                            '  ,&
          'BGD  20 Bangladesh                              '  ,&
          'BGR  37 Bulgaria                                '  ,&
          'BHR  19 Bahrain                                 '  ,&
          'BHS  18 Bahamas                                 '  ,&
          'BIH  30 Bosnia and Herzegovina                  '  ,&
          'BLM 190 Saint-Barthelemy                        '  ,&
          'BLR  22 Belarus                                 '  ,&
          'BLZ  24 Belize                                  '  ,&
          'BMU  26 Bermuda                                 '  ,&
          'BOL  28 Bolivia                                 '  ,&
          'BRA  33 Brazil                                  '  ,&
          'BRB  21 Barbados                                '  ,&
          'BRN  36 Brunei                                  '  ,&
          'BTN  27 Bhutan                                  '  ,&
          'BVT  32 Bouvet Island                           '  ,&
          'BWA  31 Botswana                                '  ,&
          'CAF  46 Central African Republic                '  ,&
          'CAN  42 Canada                                  '  ,&
          'CCK  52 Cocos Islands                           '  ,&
          'CHE 223 Switzerland                             '  ,&
          'CHL  48 Chile                                   '  ,&
          'CHN  49 China                                   '  ,&
          'CIV  57 Cote dIvoire                            '  ,&
          'CMR  41 Cameroon                                '  ,&
          'COD   0 Democratic Republic of the Congo        '  ,&
          'COG 185 Republic of Congo                       '  ,&
          'COK  55 Cook Islands                            '  ,&
          'COL  53 Colombia                                '  ,&
          'COM  54 Comoros                                 '  ,&
          'CPV  43 Cape Verde                              '  ,&
          'CRI  56 Costa Rica                              '  ,&
          'CUB  59 Cuba                                    '  ,&
          'CUW  60 Curacao                                 '  ,&
          'CXR  50 Christmas Island                        '  ,&
          'CYM  45 Cayman Islands                          '  ,&
          'CYP  61 Cyprus                                  '  ,&
          'CZE  62 Czech Republic                          '  ,&
          'DEU  86 Germany                                 '  ,&
          'DJI  65 Djibouti                                '  ,&
          'DMA  66 Dominica                                '  ,&
          'DNK  64 Denmark                                 '  ,&
          'DOM  67 Dominican Republic                      '  ,&
          'DZA   5 Algeria                                 '  ,&
          'ECU  68 Ecuador                                 '  ,&
          'EGY  69 Egypt                                   '  ,&
          'ERI  72 Eritrea                                 '  ,&
          'ESH 253 Western Sahara                          '  ,&
          'ESP 215 Spain                                   '  ,&
          'EST  73 Estonia                                 '  ,&
          'ETH  74 Ethiopia                                '  ,&
          'FIN  78 Finland                                 '  ,&
          'FJI  77 Fiji                                    '  ,&
          'FLK  75 Falkland Islands                        '  ,&
          'FRA  79 France                                  '  ,&
          'FRO  76 Faroe Islands                           '  ,&
          'FSM 146 Micronesia                              '  ,&
          'GAB  83 Gabon                                   '  ,&
          'GBR 242 United Kingdom                          '  ,&
          'GEO  85 Georgia                                 '  ,&
          'GGY  95 Guernsey                                '  ,&
          'GHA  87 Ghana                                   '  ,&
          'GIB  88 Gibraltar                               '  ,&
          'GIN  96 Guinea                                  '  ,&
          'GLP  92 Guadeloupe                              '  ,&
          'GMB  84 Gambia                                  '  ,&
          'GNB  97 Guinea-Bissau                           '  ,&
          'GNQ  71 Equatorial Guinea                       '  ,&
          'GRC  89 Greece                                  '  ,&
          'GRD  91 Grenada                                 '  ,&
          'GRL  90 Greenland                               '  ,&
          'GTM  94 Guatemala                               '  ,&
          'GUF  80 French Guiana                           '  ,&
          'GUM  93 Guam                                    '  ,&
          'GUY  98 Guyana                                  '  ,&
          'HKG 102 Hong Kong                               '  ,&
          'HMD 100 Heard Island and McDonald Islands       '  ,&
          'HND 101 Honduras                                '  ,&
          'HRV  58 Croatia                                 '  ,&
          'HTI  99 Haiti                                   '  ,&
          'HUN 103 Hungary                                 '  ,&
          'IDN 106 Indonesia                               '  ,&
          'IMN 110 Isle of Man                             '  ,&
          'IND 105 India                                   '  ,&
          'IOT  34 British Indian Ocean Territory          '  ,&
          'IRL 109 Ireland                                 '  ,&
          'IRN 107 Iran                                    '  ,&
          'IRQ 108 Iraq                                    '  ,&
          'ISL 104 Iceland                                 '  ,&
          'ISR 111 Israel                                  '  ,&
          'ITA 112 Italy                                   '  ,&
          'JAM 113 Jamaica                                 '  ,&
          'JEY 115 Jersey                                  '  ,&
          'JOR 116 Jordan                                  '  ,&
          'JPN 114 Japan                                   '  ,&
          'KAZ 117 Kazakhstan                              '  ,&
          'KEN 118 Kenya                                   '  ,&
          'KGZ 122 Kyrgyzstan                              '  ,&
          'KHM  40 Cambodia                                '  ,&
          'KIR 119 Kiribati                                '  ,&
          'KNA 193 Saint Kitts and Nevis                   '  ,&
          'KOR 213 South Korea                             '  ,&
          'KWT 121 Kuwait                                  '  ,&
          'LAO 123 Laos                                    '  ,&
          'LBN 125 Lebanon                                 '  ,&
          'LBR 127 Liberia                                 '  ,&
          'LBY 128 Libya                                   '  ,&
          'LCA 194 Saint Lucia                             '  ,&
          'LIE 129 Liechtenstein                           '  ,&
          'LKA 217 Sri Lanka                               '  ,&
          'LSO 126 Lesotho                                 '  ,&
          'LTU 130 Lithuania                               '  ,&
          'LUX 131 Luxembourg                              '  ,&
          'LVA 124 Latvia                                  '  ,&
          'MAC 132 Macao                                   '  ,&
          'MAF 191 Saint-Martin                            '  ,&
          'MAR 152 Morocco                                 '  ,&
          'MCO 148 Monaco                                  '  ,&
          'MDA 147 Moldova                                 '  ,&
          'MDG 134 Madagascar                              '  ,&
          'MDV 137 Maldives                                '  ,&
          'MEX 145 Mexico                                  '  ,&
          'MHL 140 Marshall Islands                        '  ,&
          'MKD 133 Macedonia                               '  ,&
          'MLI 138 Mali                                    '  ,&
          'MLT 139 Malta                                   '  ,&
          'MMR 154 Myanmar                                 '  ,&
          'MNE 150 Montenegro                              '  ,&
          'MNG 149 Mongolia                                '  ,&
          'MNP 168 Northern Mariana Islands                '  ,&
          'MOZ 153 Mozambique                              '  ,&
          'MRT 142 Mauritania                              '  ,&
          'MSR 151 Montserrat                              '  ,&
          'MTQ 141 Martinique                              '  ,&
          'MUS 143 Mauritius                               '  ,&
          'MWI 135 Malawi                                  '  ,&
          'MYS 136 Malaysia                                '  ,&
          'MYT 144 Mayotte                                 '  ,&
          'NAM 155 Namibia                                 '  ,&
          'NCL 159 New Caledonia                           '  ,&
          'NER 162 Niger                                   '  ,&
          'NFK 165 Norfolk Island                          '  ,&
          'NGA 163 Nigeria                                 '  ,&
          'NIC 161 Nicaragua                               '  ,&
          'NIU 164 Niue                                    '  ,&
          'NLD 158 Netherlands                             '  ,&
          'NOR 169 Norway                                  '  ,&
          'NPL 157 Nepal                                   '  ,&
          'NRU 156 Nauru                                   '  ,&
          'NZL 160 New Zealand                             '  ,&
          'OMN 170 Oman                                    '  ,&
          'PAK 171 Pakistan                                '  ,&
          'PAN 174 Panama                                  '  ,&
          'PCN 180 Pitcairn Islands                        '  ,&
          'PER 178 Peru                                    '  ,&
          'PHL 179 Philippines                             '  ,&
          'PLW 172 Palau                                   '  ,&
          'PNG 175 Papua New Guinea                        '  ,&
          'POL 181 Poland                                  '  ,&
          'PRI 183 Puerto Rico                             '  ,&
          'PRK 166 North Korea                             '  ,&
          'PRT 182 Portugal                                '  ,&
          'PRY 177 Paraguay                                '  ,&
          'PSE 173 Palestina                               '  ,&
          'PYF  81 French Polynesia                        '  ,&
          'QAT 184 Qatar                                   '  ,&
          'REU 186 Reunion                                 '  ,&
          'ROU 187 Romania                                 '  ,&
          'RUS 188 Russia                                  '  ,&
          'RWA 189 Rwanda                                  '  ,&
          'SAU 200 Saudi Arabia                            '  ,&
          'SDN 218 Sudan                                   '  ,&
          'SEN 201 Senegal                                 '  ,&
          'SGP 205 Singapore                               '  ,&
          'SGS 212 South Georgia and the South Sandwich Is '  ,&
          'SHN 192 Saint Helena                            '  ,&
          'SJM 220 Svalbard and Jan Mayen                  '  ,&
          'SLB 209 Solomon Islands                         '  ,&
          'SLE 204 Sierra Leone                            '  ,&
          'SLV  70 El Salvador                             '  ,&
          'SMR 198 San Marino                              '  ,&
          'SOM 210 Somalia                                 '  ,&
          'SPM 195 Saint Pierre and Miquelon               '  ,&
          'SRB 202 Serbia                                  '  ,&
          'SSD 214 South Sudan                             '  ,&
          'STP 199 Sao Tome and Principe                   '  ,&
          'SUR 219 Suriname                                '  ,&
          'SVK 207 Slovakia                                '  ,&
          'SVN 208 Slovenia                                '  ,&
          'SWE 222 Sweden                                  '  ,&
          'SWZ 221 Swaziland                               '  ,&
          'SXM 206 Sint Maarten                            '  ,&
          'SYC 203 Seychelles                              '  ,&
          'SYR 224 Syria                                   '  ,&
          'TCA 237 Turks and Caicos Islands                '  ,&
          'TCD  47 Chad                                    '  ,&
          'TGO 230 Togo                                    '  ,&
          'THA 228 Thailand                                '  ,&
          'TJK 226 Tajikistan                              '  ,&
          'TKL 231 Tokelau                                 '  ,&
          'TKM 236 Turkmenistan                            '  ,&
          'TLS 229 Timor-Leste                             '  ,&
          'TON 232 Tonga                                   '  ,&
          'TTO 233 Trinidad and Tobago                     '  ,&
          'TUN 234 Tunisia                                 '  ,&
          'TUR 235 Turkey                                  '  ,&
          'TUV 238 Tuvalu                                  '  ,&
          'TWN 225 Taiwan                                  '  ,&
          'TZA 227 Tanzania                                '  ,&
          'UGA 239 Uganda                                  '  ,&
          'UKR 240 Ukraine                                 '  ,&
          'UMI 244 United States Minor Outlying Islands    '  ,&
          'URY 245 Uruguay                                 '  ,&
          'USA 243 United States                           '  ,&
          'UZB 246 Uzbekistan                              '  ,&
          'VAT 248 Vatican City                            '  ,&
          'VCT 196 Saint Vincent and the Grenadines        '  ,&
          'VEN 249 Venezuela                               '  ,&
          'VGB  35 British Virgin Islands                  '  ,&
          'VIR 251 Virgin Islands, U.S.                    '  ,&
          'VNM 250 Vietnam                                 '  ,&
          'VUT 247 Vanuatu                                 '  ,&
          'WLF 252 Wallis and Futuna                       '  ,&
          'WSM 197 Samoa                                   '  ,&
          'XAD   2 Akrotiri and Dhekelia                   '  ,&
          'XCA  44 Caspian Sea                             '  ,&
          'XCL  51 Clipperton Island                       '  ,&
          'XKO 120 Kosovo                                  '  ,&
          'XNC 167 Northern Cyprus                         '  ,&
          'XPI 176 Paracel Islands                         '  ,&
          'XSP 216 Spratly Islands                         '  ,&
          'YEM 254 Yemen                                   '  ,&
          'ZAF 211 South Africa                            '  ,&
          'ZMB 255 Zambia                                  '  ,&
          'ZWE 256 Zimbabwe                                '  ,&
          'UNK 257 Unknown                                 '/
      
     DATA INDEX_RANGE_DATA / &
         14 ,&
          1 ,&
          8 ,&
          9 ,&
          3 ,&
          4 ,&
          7 ,&
        241 ,&
         12 ,&
         13 ,&
          6 ,&
         10 ,&
         82 ,&
         11 ,&
         15 ,&
         16 ,&
         17 ,&
         39 ,&
         23 ,&
         25 ,&
         29 ,&
         38 ,&
         20 ,&
         37 ,&
         19 ,&
         18 ,&
         30 ,&
        190 ,&
         22 ,&
         24 ,&
         26 ,&
         28 ,&
         33 ,&
         21 ,&
         36 ,&
         27 ,&
         32 ,&
         31 ,&
         46 ,&
         42 ,&
         52 ,&
        223 ,&
         48 ,&
         49 ,&
         57 ,&
         41 ,&
          0 ,&
        185 ,&
         55 ,&
         53 ,&
         54 ,&
         43 ,&
         56 ,&
         59 ,&
         60 ,&
         50 ,&
         45 ,&
         61 ,&
         62 ,&
         86 ,&
         65 ,&
         66 ,&
         64 ,&
         67 ,&
          5 ,&
         68 ,&
         69 ,&
         72 ,&
        253 ,&
        215 ,&
         73 ,&
         74 ,&
         78 ,&
         77 ,&
         75 ,&
         79 ,&
         76 ,&
        146 ,&
         83 ,&
        242 ,&
         85 ,&
         95 ,&
         87 ,&
         88 ,&
         96 ,&
         92 ,&
         84 ,&
         97 ,&
         71 ,&
         89 ,&
         91 ,&
         90 ,&
         94 ,&
         80 ,&
         93 ,&
         98 ,&
        102 ,&
        100 ,&
        101 ,&
         58 ,&
         99 ,&
        103 ,&
        106 ,&
        110 ,&
        105 ,&
         34 ,&
        109 ,&
        107 ,&
        108 ,&
        104 ,&
        111 ,&
        112 ,&
        113 ,&
        115 ,&
        116 ,&
        114 ,&
        117 ,&
        118 ,&
        122 ,&
         40 ,&
        119 ,&
        193 ,&
        213 ,&
        121 ,&
        123 ,&
        125 ,&
        127 ,&
        128 ,&
        194 ,&
        129 ,&
        217 ,&
        126 ,&
        130 ,&
        131 ,&
        124 ,&
        132 ,&
        191 ,&
        152 ,&
        148 ,&
        147 ,&
        134 ,&
        137 ,&
        145 ,&
        140 ,&
        133 ,&
        138 ,&
        139 ,&
        154 ,&
        150 ,&
        149 ,&
        168 ,&
        153 ,&
        142 ,&
        151 ,&
        141 ,&
        143 ,&
        135 ,&
        136 ,&
        144 ,&
        155 ,&
        159 ,&
        162 ,&
        165 ,&
        163 ,&
        161 ,&
        164 ,&
        158 ,&
        169 ,&
        157 ,&
        156 ,&
        160 ,&
        170 ,&
        171 ,&
        174 ,&
        180 ,&
        178 ,&
        179 ,&
        172 ,&
        175 ,&
        181 ,&
        183 ,&
        166 ,&
        182 ,&
        177 ,&
        173 ,&
         81 ,&
        184 ,&
        186 ,&
        187 ,&
        188 ,&
        189 ,&
        200 ,&
        218 ,&
        201 ,&
        205 ,&
        212 ,&
        192 ,&
        220 ,&
        209 ,&
        204 ,&
         70 ,&
        198 ,&
        210 ,&
        195 ,&
        202 ,&
        214 ,&
        199 ,&
        219 ,&
        207 ,&
        208 ,&
        222 ,&
        221 ,&
        206 ,&
        203 ,&
        224 ,&
        237 ,&
         47 ,&
        230 ,&
        228 ,&
        226 ,&
        231 ,&
        236 ,&
        229 ,&
        232 ,&
        233 ,&
        234 ,&
        235 ,&
        238 ,&
        225 ,&
        227 ,&
        239 ,&
        240 ,&
        244 ,&
        245 ,&
        243 ,&
        246 ,&
        248 ,&
        196 ,&
        249 ,&
         35 ,&
        251 ,&
        250 ,&
        247 ,&
        252 ,&
        197 ,&
          2 ,&
         44 ,&
         51 ,&
        120 ,&
        167 ,&
        176 ,&
        216 ,&
        254 ,&
        211 ,&
        255 ,&
        256 ,&
        257 /

     if(present(index_RANGE )) index_RANGE => index_RANGE_DATA 
     if(present(ST_NAME     )) ST_NAME     => ST_NAME_DATA     
     if(present(CNT_NAME    )) CNT_NAME    => CNT_NAME_DATA    
     if(present(ST_NAME_ABR )) ST_NAME_ABR => ST_NAME_DATA (:)(1:2) 
     if(present(CNT_NAME_ABR)) CNT_NAME_ABR=> CNT_NAME_DATA(:)(1:3) 
     
   END SUBROUTINE get_country_codes

  END MODULE  process_hres_data

! ----------------------------------------------------------------------------------------------


