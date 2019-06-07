#include "Raster.h"

!
! A Collection subroutine that helps process MODIS Albedo, GEOLAND2 LAI, and
! NGDC-HWSD-STATSGO merged soil data on their native grids 3-23-2012
!   Contact: Sarith Mahanama  sarith.p.mahanama@nasa.gov
!   Email  : sarith.p.mahanama@nasa.gov

MODULE process_hres_data
use rmTinyCatchParaMod
use MAPL_SortMod
use date_time_util  
use leap_year
use MAPL_ConstantsMod

implicit none

include 'netcdf.inc'	

private

public :: soil_para_hwsd,hres_lai,hres_gswp2, merge_lai_data, grid2tile_modis6
public :: modis_alb_on_tiles_high,modis_scale_para_high,hres_lai_no_gswp
public :: histogram, regrid_map, create_mapping, esa2mosaic , esa2clm, ESA2CLM_45
public :: grid2tile_ndep_t2m_alb, CREATE_ROUT_PARA_FILE
public :: CLM45_fixed_parameters, CLM45_clim_parameters, gimms_clim_ndvi

! Below structure is used to regrid high resolution data to high resolution tile raster

integer, parameter   :: N_tiles_per_cell = 8
integer  , parameter :: nc_esa = 129600, nr_esa = 64800, SRTM_maxcat = 291284
real, parameter      :: pi= RASTER_PI,RADIUS=MAPL_RADIUS

type :: regrid_map

   integer                               :: NT
   integer, dimension (N_tiles_per_cell) :: TID
   integer, dimension (N_tiles_per_cell) :: count

end type regrid_map

contains
!
! ---------------------------------------------------------------------
!

  SUBROUTINE ESA2CLM_45 (nc, nr, gfile)

    implicit none

    integer  , intent (in) :: nc, nr
    character (*)          :: gfile
    
    integer  , parameter   :: N_lon_clm = 7200, N_lat_clm = 3600, lsmpft = 25
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
    
    REAL, ALLOCATABLE, DIMENSION (:,:,:) :: PCTPFT 
    integer, allocatable, dimension (:) :: density, loc_int
    real   , allocatable, dimension (:) :: loc_val
    logical, allocatable, dimension (:) :: unq_mask
    integer :: NBINS, NPLUS
    integer, allocatable, dimension (:,:) :: clm_veg
    integer :: esa_clm_veg (2)
    real    :: esa_clm_frac(2)

    ! These 2 values are assumed as same as they are in surfdata_0.23x0.31_simyr2000_c100406.nc

    EDGEW = -180.
    EDGES = -90. 

    ! Reading CLM pft data file
    !--------------------------

    ALLOCATE (PCTPFT      (1:N_lon_clm, 1:N_lat_clm, 1:lsmpft))
     
    status  = NF_OPEN ('data/CATCH/CLM45/mksrf_24pftNT_landuse_rc2000_c121207.nc', NF_NOWRITE, ncid)   
    status  = NF_INQ_VARID (ncid,'PCT_PFT',VarID) ; VERIFY_(STATUS)

    do k = 1, 25 ! Natural vegetation
       status  = NF_GET_VARA_REAL (ncid,VarID,(/1,1,k/),(/N_lon_clm, N_lat_clm, 1/),PCTPFT(:,:,k)) ; VERIFY_(STATUS)
    end do

    status = NF_CLOSE(ncid)

    ! CLM 4_5 description (25)                                CLM45-carbon description (27)                                    
    ! ------------------------                                ----------------------------- 

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
    ! 'WC4G'  15 	warm c4 grass                         16  warm c4 grass [moisture +  deciduous]
    ! 'WC4Gm'   	warm c4 grass                         17  warm c4 grass [moisture stress only]
    ! 'C3CROP' 16       c3_crop                               18
    ! 'C3IRR'  17       c3_irrigated                          19
    ! 'CORN'   18       corn                                  20
    ! 'ICORN'  19       irrigated corn                        21
    ! 'STCER'  20       spring temperate cereal               22
    ! 'ISTCER' 21       irrigated spring temperate cereal     23
    ! 'WTCER'  22       winter temperate cereal               24
    ! 'IWTCER' 23       irrigated winter temperate cereal     25
    ! 'SOYB'   24       soybean                               26
    ! 'ISOYB'  25       irrigated soybean                     27
    
!**    ! 'CROP'  16 	crop                                  18  crop [moisture +  deciduous]
!**    ! 'CROPm'   	crop                                  19  crop [moisture stress only]
!**    !         17        water

    dx_clm = 360./N_lon_clm
    dy_clm = 180./N_lat_clm

    do i = 1, N_lon_clm 
       x_min_clm (i) = (i-1)*dx_clm + EDGEW  
    end do

    do i = 1,  N_lat_clm
       y_min_clm (i) = (i-1)*dy_clm  + EDGES
    end do

    ! This data set is DE
    !PCTPFT (1:N_lon_clm/2             ,:,:) =  REAL (PCT_PFT_DBL(N_lon_clm/2 + 1: N_lon_clm,:,:))
    !PCTPFT (N_lon_clm/2 + 1: N_lon_clm,:,:) =  REAL (PCT_PFT_DBL(1:N_lon_clm/2             ,:,:))

    !DEALLOCATE (PCT_PFT_DBL)

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
                      
                      if ((esa_type ==  11).or. (esa_type ==  14).or.(esa_type ==  20).or. (esa_type == 190)) then

                         ! ESA type  11: Post-flooding or irrigated croplands 
                         ! ESA type  14: Rainfed croplands 
                         ! ESA type  20: Mosaic Cropland (50-70%) / Vegetation (grassland, shrubland, forest) (20-50%) 
                         ! ESA type 190:	Artificial surfaces and associated areas (urban areas >50%) 

                         if(sum(PCTPFT(ii,jj,16:25)) > 0.) then
                            do n = 16,25                                
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,16:25))
                            end do
                         else
                            clm_veg (tile_id(i), 16) = clm_veg (tile_id(i), 16) + 1.* density(k)   
                         endif
                      endif
                      
                      ! if (esa_type == 200)  clm_veg (tile_id(i), 11) = clm_veg (tile_id(i), 11) + 1.* density(k)   ! ESA type 200:	Bare areas
                      ! if (esa_type == 210)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ocean
                      ! if (esa_type == 220)  clm_veg (tile_id(i), 17) = clm_veg (tile_id(i), 17) + 1.* density(k)   ! ice  
                      ! gkw: bare soil excluded! only considering vegetated land                   
                      ! -----------------------------------------------------------------------------------------------------------------------------------------
                      
                      if (esa_type ==  30) then 
                         ! ESA type  30: Mosaic Vegetation (grassland, shrubland, forest) (50-70%) / Cropland (20-50%) 

                         if(sum(PCTPFT(ii,jj,16:25)) > 0.) then
                            do n = 16,25                                
                               clm_veg (tile_id(i), n) = clm_veg (tile_id(i), n) + 0.5*density(k)*(PCTPFT(ii,jj,n))/sum(PCTPFT(ii,jj,16:25))
                            end do
                         elseif(sum(PCTPFT(ii,jj,2:15)) > 0.) then
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

    open (10,file='clsm/CLM4.5_veg_typs_fracs',  &
         form='formatted',status='unknown')
    open (11, file = 'clsm/catchment.def', form = 'formatted', status = 'old', &
         action =  'read')
    
    read (11, *) maxcat   
 
    do k = 1, maxcat

       read (11,'(i8,i8,5(2x,f9.4))') tid,cid,minlon,maxlon,minlat,maxlat
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

! Now splitting CLM types for CLM-carbon model
! --------------------------------------------
 
! CLM types 2- 10,12,13 are not being splitted.
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

! CLM4_5 crop types - we don't split

       if ((esa_clm_veg (1) >= 16).and.(esa_clm_veg (1) <= 25)) then
          CPT1 = esa_clm_veg (1) + 2
          CPT2 = esa_clm_veg (1) + 2
          CPF1 = esa_clm_frac(1) 
          CPF2 = 0.
       endif

       if ((esa_clm_veg (2) >= 16).and.(esa_clm_veg (2) <= 25)) then
          CST1 = esa_clm_veg (2) + 2
          CST2 = esa_clm_veg (2) + 2
          CSF1 = esa_clm_frac(2) 
          CSF2 = 0.
       endif

! Now splitting (broadleaf deciduous temperate shrub )
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

! ............. (cool c3 grass)

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

! ............. warm c4 grass

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
! CLM_4.5 : we don't splot crop type anymore 16 has become 16-25 and they are now 18-27 in catchment-CN
!       if (esa_clm_veg (1) == 16) then
!          CPT1 = 18
!          CPT2 = 19
!          CPF1 = esa_clm_frac(1) * scale        
!          CPF2 = esa_clm_frac(1) * (1. - scale) 
!       endif
!
!       if (esa_clm_veg (2) == 16) then
!          CST1 = 18
!          CST2 = 19
!          CSF1 = esa_clm_frac(2) * scale        
!          CSF2 = esa_clm_frac(2) * (1. - scale) 
!       endif

       ! fractions must sum to 1
       ! -----------------------
       ftot = cpf1 + cpf2 + csf1 + csf2

       if(ftot /= 100.) then
          cpf1 = 100. * cpf1 / ftot  
          cpf2 = 100. * cpf2 / ftot 
          csf1 = 100. * csf1 / ftot 
          csf2 = 100. * csf2 / ftot 
       endif
    
       write (10,'(2I8,4I3,4f7.2,2I3,2f7.2)')     &
            tid,cid,cpt1, cpt2, cst1, cst2, cpf1, cpf2, csf1, csf2, &
            esa_clm_veg (1), esa_clm_veg (2), esa_clm_frac(1), esa_clm_frac(2)
    end do

    close (10, status = 'keep')
    close (11, status = 'keep')    

  END SUBROUTINE ESA2CLM_45

!
! ------------------------------------------------------------------------------------------------
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

    ! CLM description (17)                                       CLM-carbon description (19)                                    
    ! --------------------                                        -------------------------- 

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
 
    do k = 1, maxcat

       read (11,'(i8,i8,5(2x,f9.4))') tid,cid,minlon,maxlon,minlat,maxlat
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

! Now splitting CLM types for CLM-carbon model
! --------------------------------------------
 
! CLM types 2- 10,12,13 are not being splitted.
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
    
       write (10,'(2I8,4I3,4f7.2,2I3,2f7.2)')     &
            tid,cid,cpt1, cpt2, cst1, cst2, cpf1, cpf2, csf1, csf2, &
            esa_clm_veg (1), esa_clm_veg (2), esa_clm_frac(1), esa_clm_frac(2)
    end do

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
    logical :: jpl_height = .true.

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

       read (11,'(i8,i8,5(2x,f9.4))') tid,cid
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
          write (10,'(i8,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')     &
            tid,cid,mos1,mos2,100.*mfrac,100.*sfrac, z2(k), z0 (k)
          
       endif
    end do

    close (10,status='keep')    
    close (11,status='keep') 

    open (20,file='clsm/vegdyn.data',status='unknown',action='write',form='unformatted', &
         convert='little_endian')   
    write (20) real(ityp)
    write (20) z2 (:)
    write (20) z0 (:)
    
    close (20)    
    

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
   type (regrid_map), intent (inout), dimension (nc_data,nr_data) :: rmap
   integer :: i,j,n, i1,i2,j1,j2,ncatch, nbins, status, NPLUS
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
   rmap%NT    = 0

   if((nc_data  >= nc).and.(nr_data >= nr)) then

      allocate(iraster(nc_data,nr_data),stat=STATUS); VERIFY_(STATUS)
      call RegridRaster(tile_id,iraster)
  
      do j = 1,nr_data
         do i =  1,nc_data
            if((iraster (i,j) >=1).and.(iraster (i,j) <=ncatch)) then
               rmap(i,j)%NT = 1
               rmap(i,j)%TID  (rmap(i,j)%NT) = iraster (i,j)
               rmap(i,j)%count(rmap(i,j)%NT) = 1
            endif
         end do
      end do
      deallocate (iraster) ; VERIFY_(STATUS)

   else
      
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
                        rmap(i,j)%NT = rmap(i,j)%NT + 1
                        if(rmap(i,j)%NT > N_tiles_per_cell) then
                           print *,'N_tiles_per_cell exceeded :', rmap(i,j)%NT
                           print *,i,j,i1,i2,j1,j2
                           print *,'NT',rmap(i,j)%NT 
                           print *,rmap(i,j)%TID
                           print *,rmap(i,j)%count
                           stop
                        endif
                        rmap(i,j)%TID  (rmap(i,j)%NT) = NINT(loc_val(n))
                        rmap(i,j)%count(rmap(i,j)%NT) = density(n)
                     endif
                  END DO
                  deallocate (loc_val, density)
                  deallocate (loc_int, unq_mask)
               endif
               NULLIFY (subset)
            else
               if((tile_id (i1,j1) > 0).and.(tile_id(i1,j1).le.ncatch)) then
                  rmap(i,j)%NT       = 1
                  rmap(i,j)%TID(1)   = tile_id (i1,j1)
                  rmap(i,j)%COUNT(1) = 1
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
    integer, parameter :: i5_srtm = 267084 , i6_srtm =  291284   ! Australia

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
    REAL :: latt,lont,fr_gcm,fr_cat,tsteps,zth, slr,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    character*100 :: path,fname,fout,metpath
    character (*) :: gfile
    integer :: n,maxcat,ip
    integer :: yy,j,month
    integer, allocatable, dimension (:) :: vegcls 
    real, allocatable, dimension (:) :: &
         modisvf, modisnf,albvr,albnr,albvf,albnf,lat,lon, &
         green,lai,snw,lai_before,lai_after,grn_before,grn_after
    real, allocatable, dimension (:) :: &
         calbvr,calbnr,calbvf,calbnf
    character*300 :: ifile1,ifile2,ofile
    integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    integer :: day, hour, min ,secs_in_day,k
    real :: yr,mn,dy,yr1,mn1,dy1,dum, slice1,slice2
    logical :: save_sib = .false.

    fname='clsm/catchment.def'
    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)maxcat

    allocate (albvr    (1:maxcat))
    allocate (albvf    (1:maxcat))
    allocate (albnr    (1:maxcat))
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
    allocate (snw      (1:maxcat))
    close (10,status='keep')

    fname=trim(gfile)//'.til'

    open (10,file=fname,status='old',action='read',form='formatted')
    fname='clsm/mosaic_veg_typs_fracs'
    open (20,file=fname,status='old',action='read',form='formatted')

    read (10,*)ip
    read (10,*)j_dum
    read (10,'(a)')version
    read (10,*)nc_gcm
    read (10,*)nr_gcm
    read (10,'(a)')version
    read (10,*)nc_ocean
    read (10,*)nr_ocean

    do n = 1,ip
      if (ease_grid) then     
	 read(10,*,IOSTAT=ierr) typ,pfs,lon,lat,ig,jg,fr_gcm
      else
      read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lont,latt,ig,jg,fr_gcm,indx_dum,pfs,j_dum,fr_cat,j_dum
      endif
       if (typ == 100) then
          ip2 = n 
          read (20,'(i8,i8,2(2x,i3),2(2x,f6.4))')     &
            indr1,indr1,vegcls(ip2),indr1,fr_gcm,fr_gcm
       endif
       if(ierr /= 0)write (*,*)'Problem reading'
    end do
    close (10,status='keep')
    close (20,status='keep')

    albvr    =0.
    albvf    =0.
    albnr    =0.
    albnf    =0.
    calbvf   =0.
    calbnf   =0.
    modisvf  =0.
    modisnf  =0.
    snw      =0.

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
    albvr    =0.
    albnr    =0.
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
          
          call sibalb(                                    &
               albvr,albnr,albvf,albnf,                   &
               lai, green, 0.0, snw, vegcls, maxcat)  
          
          calbvf = calbvf + albvf
          calbnf = calbnf + albnf         
          tsteps = tsteps + 1.  
          call augment_date_time( 86400, date_time_new ) 

          if (datetime_le_refdatetime(date_time_new,af_modis_time)) then
             
          else
             bf_modis_time = af_modis_time
             calbvf = calbvf/tsteps
             calbnf = calbnf/tsteps

             modisvf = modisvf/calbvf
             modisnf = modisnf/calbnf
           
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
               albvr    =0.
               albnr    =0.
               albvf    =0.
               albnf    =0.
               tsteps   =0.   
            endif
         endif
      end do
  
      deallocate (modisvf,modisnf,albvr,albvf,albnf,albnr)
      deallocate (green,lai)
      deallocate (vegcls)
      deallocate (calbvf,calbnf)
      deallocate (lai_before,grn_before, lai_after,grn_after,snw)

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
  type (regrid_map), intent (in), dimension (nc_data,nr_data) :: rmap
  character*6 :: MA
  character(*)  :: gfiler
  integer :: n,maxcat,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj
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
                         if(net_data1(i-iLL +1 ,j - jLL +1) > 0) then
                            if(rmap(i,j)%nt > 0) then
                               do n = 1, rmap(i,j)%nt
                                  vec_AlbVis(rmap(i,j)%tid(n))  = vec_AlbVis(rmap(i,j)%tid(n)) +  &
                                       sf*net_data1(i-iLL +1 ,j - jLL +1)*rmap(i,j)%count(n) 
                                  count_AlbVis(rmap(i,j)%tid(n))= count_AlbVis(rmap(i,j)%tid(n)) + &
                                        1.*rmap(i,j)%count(n)
                               end do
                            endif
                         endif
                         if(net_data2(i-iLL +1 ,j - jLL +1) > 0) then
                            if(rmap(i,j)%nt > 0) then
                               do n = 1, rmap(i,j)%nt
                                  vec_AlbNir(rmap(i,j)%tid(n))  = vec_AlbNir(rmap(i,j)%tid(n)) +  &
                                       sf*net_data2(i-iLL +1 ,j - jLL +1)*rmap(i,j)%count(n) 
                                  count_AlbNir(rmap(i,j)%tid(n))= count_AlbNir(rmap(i,j)%tid(n)) + &
                                        1.*rmap(i,j)%count(n)
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
  type (regrid_map), intent (in), dimension (nc_data,nr_data) :: rmap
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
                           if(rmap(i,j)%nt > 0) then
                              do n = 1, rmap(i,j)%nt
                                 if(vec_lai(rmap(i,j)%tid(n)) == -9999.) vec_lai(rmap(i,j)%tid(n)) = 0.                                 
                                 vec_lai(rmap(i,j)%tid(n))   = vec_lai(rmap(i,j)%tid(n)) + &
                                      sf*net_data1(i-iLL +1 ,j - jLL +1)*rmap(i,j)%count(n)
                                 count_lai(rmap(i,j)%tid(n)) = &
                                      count_lai(rmap(i,j)%tid(n)) + 1.*rmap(i,j)%count(n)                                     
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
  type (regrid_map), intent (in), dimension (nc_data,nr_data) :: rmap
  real :: dum, gyr,gmn,gdy,gyr1,gmn1,gdy1, slice1,slice2
  character*100 :: fname,fout
  character*10 :: string
  character*2 :: VV,HH
  integer, allocatable, dimension (:,:) :: &
         net_data1
  REAL, ALLOCATABLE, dimension (:) :: vec_lai, count_lai
  character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
  logical :: regrid
  REAL :: sf
  logical :: first_entry = .true.
  type (date_time_type) :: date_time_new,bf_lai_time,   &
       af_lai_time
  integer, intent(in), optional :: merge

  if(trim(lai_name) == 'lai'  ) vid = 4
  if(trim(lai_name) == 'green') vid = 5

!
! Reading number of cathment-tiles from catchment.def file
! -------------------------------------------------------- 
      fname='clsm/catchment.def' 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
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
       !
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
                            if(rmap(i,j)%nt > 0) then
                               do n = 1, rmap(i,j)%nt
                                  if(vec_lai(rmap(i,j)%tid(n)) == -9999.) vec_lai(rmap(i,j)%tid(n)) = 0.                                 
                                  vec_lai(rmap(i,j)%tid(n))   = vec_lai(rmap(i,j)%tid(n)) + &
                                       sf*net_data1(i-iLL +1 ,j - jLL +1)*rmap(i,j)%count(n)
                                  count_lai(rmap(i,j)%tid(n)) = &
                                       count_lai(rmap(i,j)%tid(n)) + 1.*rmap(i,j)%count(n)                                     
                               end do
                            endif
                         endif
                      enddo
                   enddo
                   status = NF_CLOSE(ncid)
                endif
             end do
          end do
          
          write(31) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
          where (count_lai > 0.) vec_lai = vec_lai/count_lai
          where (count_lai == 0.)vec_lai = 0.0001
          write(31)  vec_lai(:)
       end do
       close(31,status='keep')

       deallocate (net_data1)
       deallocate (count_lai)
       deallocate (vec_lai)

  END SUBROUTINE hres_gswp2

!----------------------------------------------------------------------  

  SUBROUTINE soil_para_hwsd (nx,ny,gfiler)

!
! Processing NGDC-HWSD-STATSGO merged soil properties with Woesten Soil
! Parameters and produces tau_param.dat and soil_param.dat files
! 
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
      REAL, ALLOCATABLE, dimension (:) :: soildepth, grav_vec,soc_vec,poc_vec,&
             ncells_top,ncells_top_pro,ncells_sub_pro 
      integer(kind=2) , allocatable, dimension (:) :: ss_clay,    &
             ss_sand,ss_clay_all,ss_sand_all,ss_oc_all
      REAL, ALLOCATABLE :: count_soil(:)
      integer, allocatable, target, dimension (:,:) :: tile_id
      integer, pointer :: iRaster(:,:)
      integer :: tindex, pfafindex,fac,o_cl,o_clp,fac_surf,vtype
      real,dimension(4) :: cFamily
      real   ,dimension(5) :: cF_lim
      logical :: first_entry = .true.
      logical :: regrid,write_file
      INTEGER, allocatable, dimension (:) :: soil_class_top,soil_class_com
      REAL :: sf,factor,wp_wetness,fac_count

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
      cF_lim(1) =  0.
      cF_lim(2) =  0.4      ! 0.365    ! 0.3
      cF_lim(3) =  0.64     ! 0.585    ! 4.0	   
      cF_lim(4) =  15./1.72 ! 9.885    ! 8.5
      cF_lim(5) =  100.0

      nsoil_pcarbon(1) = 84 ! 84
      nsoil_pcarbon(2) = nsoil_pcarbon(1) + 84 ! 84
      nsoil_pcarbon(3) = nsoil_pcarbon(2) + 84 ! 57

      fname='clsm/catchment.def'
!
! Reading number of cathment-tiles from catchment.def file
! 
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      
      close (10,status='keep')

      fname =trim(c_data)//'SOIL-DATA/GSWP2_soildepth_H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_CLOSE(ncid); VERIFY_(STATUS)

      allocate(soildepth(1:maxcat))
      allocate(soil_high(1:i_highd,1:j_highd))  
      allocate(count_soil(1:maxcat))  
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
! reading soil depth data
!
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
      
! Interpolation or aggregation on to catchment-tiles

      soildepth =0.
      count_soil = 0.
 
      do j=1,ny_adj
         do i=1,nx_adj
            if((iRaster(i,j).gt.0).and.(iRaster(i,j).le.maxcat)) then
               if ((raster(i,j).gt.0)) then
                  soildepth(iRaster(i,j)) = &
                       soildepth(iRaster(i,j)) + sf*raster(i,j)
                  count_soil(iRaster(i,j)) = &
                       count_soil(iRaster(i,j)) + 1. 
               endif
            endif
         end do
      end do

      DO n =1,maxcat
	 	if(count_soil(n)/=0.) soildepth(n)=soildepth(n)/count_soil(n)	
	          soildepth(n) = max(soildepth(n),1334.)
!                  soildepth(n) = soildepth(n) + 2000.
!                  soildepth(n) = min(soildepth(n),8000.) 
       END DO

      deallocate (SOIL_HIGH)
      deallocate (count_soil)
      NULLIFY(Raster)
!
! Reading NGDC-HWSD-STATSGO merged Soil Properties
!
      fname =trim(c_data)//'SOIL-DATA/SoilProperties_H11V13.nc'
      status = NF_OPEN(trim(fname),NF_NOWRITE, ncid); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lon_global',i_highd); VERIFY_(STATUS)
      status = NF_GET_att_INT(ncid,NF_GLOBAL,'N_lat_global',j_highd); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,1,string, nc_10); VERIFY_(STATUS)
      status = NF_INQ_DIM (ncid,2,string, nr_10); VERIFY_(STATUS)
      status = NF_CLOSE(ncid)

      regrid = nx/=i_highd .or. ny/=j_highd
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

      sand_top = -9999.
      clay_top = -9999.
      oc_top   = -9999.
      sand_sub = -9999.
      clay_sub = -9999.
      oc_sub   = -9999.
      grav_grid= -9999.

      do jx = 1,18
      	 do ix = 1,36
	    write (vv,'(i2.2)')jx
	    write (hh,'(i2.2)')ix 
	    fname = trim(c_data)//'SOIL-DATA/SoilProperties_H'//hh//'V'//vv//'.nc'
            status = NF_OPEN(trim(fname),NF_NOWRITE, ncid)
	    if(status == 0) then
		status = NF_GET_att_INT  (ncid, NF_GLOBAL,'i_ind_offset_LL',iLL); VERIFY_(STATUS)
                status = NF_GET_att_INT  (ncid, NF_GLOBAL,'j_ind_offset_LL',jLL); VERIFY_(STATUS)
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

! now regridding

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

! Deallocate large arrays

      allocate(land_pixels(1:size(iRaster,1),1:size(iRaster,2)))
      land_pixels = (iRaster >=1).and.(iRaster<=maxcat)
      i1 = count(land_pixels)   
      deallocate (land_pixels) 

      allocate (tileid_vec(1:i1))
      allocate (data_vec1 (1:i1))
      allocate (data_vec2 (1:i1))
      allocate (data_vec3 (1:i1))
      allocate (data_vec4 (1:i1))
      allocate (data_vec5 (1:i1))
      allocate (data_vec6 (1:i1))
      allocate (grav_vec  (1:maxcat))
      allocate (soc_vec   (1:maxcat))
      allocate (poc_vec   (1:maxcat))
      allocate (ncells_top  (1:maxcat))
      allocate (ncells_top_pro  (1:maxcat))
      allocate (ncells_sub_pro  (1:maxcat))
      allocate(count_soil(1:maxcat))  
      count_soil = 0.
      grav_vec   = 0.
      soc_vec    = 0.
      poc_vec    = 0.
      ncells_top = 0.
      ncells_top_pro = 0.
      ncells_sub_pro = 0.

      n =1
      do j=1,ny_adj
         do i=1,nx_adj
            if((iRaster(i,j).ge.1).and.(iRaster(i,j).le.maxcat)) then

	       tileid_vec (n) =  iRaster(i,j)
	       data_vec1  (n) =  Raster1(i,j)
	       data_vec2  (n) =  Raster2(i,j)
	       data_vec3  (n) =  Raster3(i,j)
	       data_vec4  (n) =  Raster4(i,j)
	       data_vec5  (n) =  Raster5(i,j)
	       data_vec6  (n) =  Raster6(i,j)

               if ((raster(i,j).gt.0)) then 
                  grav_vec(iRaster(i,j)) = &
                       grav_vec(iRaster(i,j)) + sf*raster(i,j)
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

      deallocate (grav_grid)
      deallocate (count_soil)
      NULLIFY(Raster)

      NULLIFY(Raster1,Raster2,Raster3,Raster4,Raster5,Raster6)
      deallocate (clay_top,sand_top,oc_top,clay_sub,sand_sub,oc_sub) 
      deallocate (tile_id)

      allocate (arrayA    (1:i1))
      allocate (arrayB    (1:i1))

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
!
! Reading Woesten Soil Parameters and CLSM tau parameters
!
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

      fname = trim(c_data)//'SoilClasses-SoilHyd-TauParam.dat'
      table_map = 0
      open (11, file=trim(fname), form='formatted',status='old', &
           action = 'read')
      read (11,'(a)')fout
      do n =1,n_SoilClasses 
      	 read (11,'(4f7.3,4f8.4,e13.5,2f12.7,2f8.4,4f12.7)')a_sand(n),a_clay(n),a_silt(n),a_oc(n),a_bee(n),a_psis(n), &
              a_poros(n),a_wp(n),a_aksat(n),atau(n),btau(n),a_wpsurf(n),a_porosurf(n),atau_2cm(n),btau_2cm(n)

	 min_percs%clay_perc = a_clay(n)
	 min_percs%silt_perc = a_silt(n)
	 min_percs%sand_perc = a_sand(n)
	 if(n <= nsoil_pcarbon(1))                              table_map(soil_class (min_percs),1) = n  
	 if((n > nsoil_pcarbon(1)).and.(n <= nsoil_pcarbon(2))) table_map(soil_class (min_percs),2) = n  
         if((n > nsoil_pcarbon(2)).and.(n <= nsoil_pcarbon(3))) table_map(soil_class (min_percs),3) = n 

      end do
      close (11,status='keep') 
!
!  When Woesten Soil Parameters are not available for a particular Soil Class
!  ,as assumed by tiny triangles in HWSD soil triangle, Woesten Soil
!  parameters from the nearest available tiny triangle will be substituted.
!	     	  
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
!$OMP         soc_vec,poc_vec,ncells_top,ncells_top_pro,&
!$OMP         ncells_sub_pro)    &
!$OMP PRIVATE(n,i,j,k,icount,t_count,i1,i2,ss_clay,     &
!$OMP         ss_sand,ss_clay_all,ss_sand_all,          &
!$OMP         ss_oc_all,cFamily,factor,o_cl,o_clp,ktop, &
!$OMP         min_percs, fac_count, write_file)

    DO t_count = 1,n_threads
      DO n = low_ind(t_count),upp_ind(t_count)

	write_file = .false.

!	if (n==171010)  write_file = .true.

        if(n==low_ind(t_count)) then
             icount = 1
             do k=1,low_ind(t_count) - 1
                do while (tileid_vec(icount)== k)
                   icount = icount + 1
                end do
             end do
        endif
             
        i1 = icount 
        
	loop: do while (tileid_vec(icount)== n)
	   if(icount <= size(tileid_vec,1)) icount = icount + 1
           if(icount > size(tileid_vec,1)) exit loop
	end do loop 

	i2 = icount -1
 	i = i2 - i1 + 1

        allocate(ss_clay    (1:2*i))
        allocate(ss_sand    (1:2*i))
        allocate(ss_clay_all(1:2*i))
        allocate(ss_sand_all(1:2*i))
        allocate(ss_oc_all  (1:2*i))
	  
        ss_clay    = 0    
        ss_sand    = 0	
        ss_clay_all= 0
        ss_sand_all= 0
        ss_oc_all  = 0

	ss_clay_all (1:i)     = data_vec1(i1:i2)
	ss_sand_all (1:i)     = data_vec2(i1:i2)
	ss_oc_all   (1:i)     = data_vec3(i1:i2)	
	ss_clay_all (1+i:2*i) = data_vec4(i1:i2) 
	ss_sand_all (1+i:2*i) = data_vec5(i1:i2)
	ss_oc_all   (1+i:2*i) = data_vec6(i1:i2)	

	cFamily = 0.

	do j=1,i
	   if(j <= i) factor = 1.
	   if((ss_oc_all(j)*sf >=  cF_lim(1)).and. (ss_oc_all(j)*sf < cF_lim(2))) cFamily(1) = cFamily(1) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(2)).and. (ss_oc_all(j)*sf < cF_lim(3))) cFamily(2) = cFamily(2) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(3)).and. (ss_oc_all(j)*sf < cF_lim(4))) cFamily(3) = cFamily(3) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(4) ))                                   cFamily(4) = cFamily(4) + factor
	end do

	if (sum(cFamily) == 0.) o_cl  = 1
	if (sum(cFamily)  > 0.) o_cl  = maxloc(cFamily, dim = 1)

	cFamily = 0.

	do j=1,2*i
	   if(j <= i) factor = 1.
	   if(j  > i) factor = 2.33
	   if((ss_oc_all(j)*sf >=  cF_lim(1)).and. (ss_oc_all(j)*sf < cF_lim(2))) cFamily(1) = cFamily(1) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(2)).and. (ss_oc_all(j)*sf < cF_lim(3))) cFamily(2) = cFamily(2) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(3)).and. (ss_oc_all(j)*sf < cF_lim(4))) cFamily(3) = cFamily(3) + factor
	   if((ss_oc_all(j)*sf >=  cF_lim(4) ))                                   cFamily(4) = cFamily(4) + factor
	end do

	if (sum(cFamily) == 0.) o_clp = 1
	if (sum(cFamily)  > 0.) o_clp = maxloc(cFamily, dim = 1)

        if(o_cl == 4) then 
           soil_class_top(n) = n_SoilClasses
	   ktop = 0
	   do j=1,i
	     if(ss_oc_all(j)*sf >= cF_lim(4)) then
	        soc_vec (n) = soc_vec(n) + ss_oc_all(j)*sf
		ktop = ktop + 1
             endif
	   end do
	   if(ktop.ne.0) soc_vec (n)   = soc_vec(n)/ktop
	   ncells_top(n) = 100.*float(ktop)/float(i)
        else 
            k = 1
	    ktop = 1

	    do j=1,i
	        if((ss_oc_all(j)*sf >= cF_lim(o_cl)).and.(ss_oc_all(j)*sf < cF_lim(o_cl + 1))) then 
                   if((ss_clay_all(j)*sf >= 0.).and.(ss_sand_all(j)*sf >= 0.)) then   
    			   ss_clay (k) = ss_clay_all(j)
			   ss_sand (k) = ss_sand_all(j)
                           if((ss_clay (k) + ss_sand (k)) > 9999) then
                              if(ss_clay (k) >= ss_sand (k)) then
                                 ss_sand (k) = 10000 - ss_clay (k)
                               else
                                 ss_clay (k) = 10000 - ss_sand (k)
                               endif
                            endif
			   soc_vec (n) = soc_vec(n) + ss_oc_all(j)*sf
                           k = k + 1
                           ktop = ktop + 1
		   endif
                endif	
	    end do
	    
	    k = k - 1
	    ktop = ktop -1
	    if(ktop.ne.0) soc_vec (n) = soc_vec(n)/ktop
	    ncells_top(n) = 100.*float(ktop)/float(i)
            if (write_file) write(80+n,*)ktop,o_cl
            if(ktop > 0) then 
               if (write_file) write (80+n,*)ss_clay(1:ktop)
               if (write_file) write (80+n,*)ss_sand(1:ktop)
            endif
  	    j = center_pix_int0(sf, ktop,ktop, ss_clay(1:ktop),ss_sand(1:ktop))
            if (write_file) write(80+n,*)j

            if(j >=1) then 
               min_percs%clay_perc = ss_clay(j)*sf
               min_percs%sand_perc = ss_sand(j)*sf
               min_percs%silt_perc = 100. - ss_clay(j)*sf - ss_sand(j)*sf
               soil_class_top (n) = table_map(soil_class (min_percs),o_cl)   
            endif
        endif
            if (write_file) write(80+n,*)soil_class_top (n) 
        if(o_clp == 4) then 
           soil_class_com(n) = n_SoilClasses
	   fac_count = 0.
	   k =0
	   ktop =0
	   do j=1,2*i
	     if(ss_oc_all(j)*sf >= cF_lim(4)) then
                 if(j <= i) factor = 1.
	         if(j  > i) factor = 2.33
		 if(j  > i) k = k + 1
		 if(j <= i) ktop = ktop + 1
	         poc_vec (n) = poc_vec(n) + ss_oc_all(j)*sf*factor
		 fac_count = fac_count + factor
             endif
	   end do
	   if(fac_count.ne.0) poc_vec (n) = poc_vec (n)/fac_count
           ncells_sub_pro(n) = 100.*float(k)/float(i)
           ncells_top_pro(n) = 100.*float(ktop)/float(i)
        else
            k = 1
	    ktop = 1

            ss_clay=0
            ss_sand=0
	    fac_count = 0.

	    do j=1,2*i
	        if((ss_oc_all(j)*sf >=  cF_lim(o_clp)).and.(ss_oc_all(j)*sf < cF_lim(o_clp + 1))) then 
                   if((ss_clay_all(j)*sf >= 0.).and.(ss_sand_all(j)*sf >= 0.)) then 
  	              if(j <= i) factor = 1.
	              if(j  > i) factor = 2.33
	              poc_vec (n) = poc_vec(n) + ss_oc_all(j)*sf*factor
		      fac_count = fac_count + factor
		      if(j <= i) then
    			   ss_clay (k) = ss_clay_all(j)
			   ss_sand (k) = ss_sand_all(j)
                           if((ss_clay (k) + ss_sand (k)) > 9999) then
                              if(ss_clay (k) >= ss_sand (k)) then
                                 ss_sand (k) = 10000 - ss_clay (k)
                               else
                                 ss_clay (k) = 10000 - ss_sand (k)
                               endif
                            endif
                           k = k + 1
			   ktop = ktop + 1
		       else
    			   ss_clay (k) = ss_clay_all(j)
			   ss_sand (k) = ss_sand_all(j)
                           if((ss_clay (k) + ss_sand (k)) > 9999) then
                              if(ss_clay (k) >= ss_sand (k)) then
                                 ss_sand (k) = 10000 - ss_clay (k)
                               else
                                 ss_clay (k) = 10000 - ss_sand (k)
                               endif
                            endif
                           k = k + 1                         
		       endif   
		   endif
                endif	
	    end do
	    
	    k = k - 1
	    ktop = ktop -1
	    if(fac_count.ne.0) poc_vec (n) = poc_vec(n)/fac_count
	    ncells_top_pro(n) = 100.*float(ktop)/float(i)
	    ncells_sub_pro(n) = 100.*float(k-ktop)/float(i)

            if (write_file) write (80+n,*)ktop,k,o_cl
            if (write_file) write (80+n,*)ss_clay(1:k)
            if (write_file) write (80+n,*)ss_sand(1:k)
	    j = center_pix_int0 (sf, ktop,k, ss_clay(1:k),ss_sand(1:k))
            if (write_file) write(80+n,*) j
            if(j >=1) then 
               min_percs%clay_perc = ss_clay(j)*sf
               min_percs%sand_perc = ss_sand(j)*sf
               min_percs%silt_perc = 100. - ss_clay(j)*sf - ss_sand(j)*sf
               soil_class_com (n) = table_map(soil_class (min_percs),o_clp)  
            endif
            if (write_file) write(80+n,*) soil_class_com (n) 
            if (write_file) close(80+n)          
        endif
	deallocate (ss_clay,ss_sand,ss_clay_all,ss_sand_all,ss_oc_all)
      END DO
      END DO
!$OMP ENDPARALLELDO
            
      fname='clsm/catchment.def'
      open (10,file=fname,status='old',action='read',form='formatted')
      read(10,*) maxcat
      fname ='clsm/soil_param.first'
      open (11,file=trim(fname),form='formatted',status='unknown',action = 'write')

      fname ='clsm/tau_param.dat'
      open (12,file=trim(fname),form='formatted',status='unknown',action = 'write')

      fname ='clsm/mosaic_veg_typs_fracs'
      open (13,file=trim(fname),form='formatted',status='old',action = 'read')

      do n = 1, maxcat

      	 read (10,*) tindex,pfafindex
         read (13,*) tindex,pfafindex,vtype

         ! fill gaps from neighbor for rare missing values came from inconsistent masks
         if ((soil_class_top (n) == -9999).or.(soil_class_com (n) == -9999)) then

            ! if com-layer has data the issues is only with top-layer
            ! -------------------------------------------------------

            if(soil_class_com (n) >= 1) soil_class_top (n) = soil_class_com (n)

            ! if there is nothing look for the neighbor
            ! -----------------------------------------
            
            if (soil_class_com (n) == -9999) then
               do k = 1, maxcat
                  j  = 0
                  i1 = n - k
                  i2 = n + k
                  if((i1 >=     1).and.(soil_class_com (i1) >=1)) j = i1
                  if((i2 <=maxcat).and.(soil_class_com (i2) >=1)) j = i2

                  if (j > 0) then
                     soil_class_com (n) = soil_class_com (j)
                     soil_class_top (n) = soil_class_com (n)
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

         wp_wetness = a_wp(fac)     /a_poros(fac)
         write (11,'(i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)')tindex,pfafindex,      &
               soil_class_top(n),soil_class_com(n),a_bee(fac),a_psis(fac),a_poros(fac),&
               a_aksat(fac)/exp(-1.0*zks*gnu),wp_wetness,soildepth(n),                 &
               grav_vec(n),soc_vec(n),poc_vec(n), &
               a_sand(fac_surf),a_clay(fac_surf),a_sand(fac),a_clay(fac), &
	       a_wpsurf(fac_surf)/a_porosurf(fac_surf),a_porosurf(fac_surf)
	       	    
         write (12,'(i8,i8,4f10.7)')tindex,pfafindex, &
	       atau_2cm(fac_surf),btau_2cm(fac_surf),atau(fac_surf),btau(fac_surf)  
	    	    
      end do
      write (11,'(a)')'                    '
      write (11,'(a)')'FMT=i8,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4'
      write (11,'(a)')'TileIndex PfafID SoilClassTop SoilClassProfile BEE PSIS POROS Ks_at_SURF WPWET SoilDepth %Grav %OCTop %OCProf %Sand_top %Clay_top %Sand_prof %Clay_prof WPWET_SURF POROS_SURF'
      close (10, status = 'keep')	            
      close (11, status = 'keep')	            
      close (12, status = 'keep')	            
      close (13, status = 'keep')

      deallocate (data_vec1, data_vec2,data_vec3, data_vec4,data_vec5, data_vec6)
      deallocate (tileid_vec)
      deallocate (a_sand,a_clay,a_silt,a_oc,a_bee,a_psis,       &
            a_poros,a_wp,a_aksat,atau,btau,a_wpsurf,a_porosurf, &
            atau_2cm,btau_2cm)
      deallocate (soildepth, grav_vec,soc_vec,poc_vec,&
             ncells_top,ncells_top_pro,ncells_sub_pro,soil_class_top,soil_class_com)

  END SUBROUTINE soil_para_hwsd
!
! ====================================================================
!

INTEGER FUNCTION center_pix_int (sf,ktop, ktot, x,y,x0,y0,z0,ext_point)

implicit none

integer (kind =2), dimension (:), intent (in) :: x,y
integer, intent (in) :: ktop,ktot
real, intent (in) :: sf
real :: xi,xj,yi,yj,xx0,yy0,zz0
real, allocatable, dimension (:,:) :: length_m
real, allocatable, dimension (:) :: length
real, intent (inout) :: x0,y0,z0
integer :: i,j,npix
logical, intent(in) :: ext_point
real :: zi, zj

allocate (length_m (1:ktot,1:ktot))
allocate (length   (1:ktot))
length_m =0.
length   =0.

center_pix_int = -9999
if(ktot /= 0) then
   do i = 1,ktot
      xi = sf*x(i)
      yi = sf*y(i)
      zi = 100. - xi - yi
      if (.not. ext_point) then
         x0 = xi
         y0 = yi
         z0 = zi
      endif
      
      do j = 1,ktot
         xj = sf*x(j)
         yj = sf*y(j)
         zj = 100. - xj - yj
         xx0= xj - x0
         yy0= yj - y0
         zz0= zj - z0
 
         if(ktot > ktop) then 
            if(j <= ktop) then
               length_m (i,j) = (xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5
            else
               length_m (i,j) = 2.33*((xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5)
            endif
         else
            length_m (i,j) = (xx0*xx0 +  yy0*yy0 + zz0*zz0)**0.5
         endif
      end do
      length (i) = sum(length_m (i,:))
   end do
   
center_pix_int = minloc(length,dim=1)
endif

END FUNCTION center_pix_int

!
! ====================================================================
!

INTEGER FUNCTION center_pix_int0 (sf,ktop, ktot, x,y)

implicit none
! sf = 0.01 (integer to real scale factor), ktop = # of pixels in top layer
! ktot = total # of pixels, top + subsurface combined
! x (clay), y (sand_
integer (kind =2), dimension (:), intent (in) :: x,y
integer, intent (in) :: ktop,ktot
real, intent (in) :: sf
real :: xi,xj,yi,yj
real :: length

integer :: i,j,npix
real :: zi, zj, mindist,xc,yc,zc

length   =0.

center_pix_int0 = -9999

if(ktot /= 0) then
! There should be some data pixels
   if(ktot > ktop) then 
      ! Have both layers
      if(ktop > 0) then
         ! There are data in top layer
         xc = sf*0.3*sum(real(x(1:ktop)))/real(ktop) + sf*0.7*sum(real(x(ktop + 1 : ktot)))/real(ktot - ktop)  
         yc = sf*0.3*sum(real(y(1:ktop)))/real(ktop) + sf*0.7*sum(real(y(ktop + 1 : ktot)))/real(ktot - ktop)
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
   zc = 100. - xc - yc
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

!
! ====================================================================
!
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
    integer :: i, j, n, im, jm, lwi, idum, ntiles, nland, nv, ix, jx, itype, iband, isum, ntl, np, jalbx, ialbx

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

             READ (10,'(i8,i8,5(2x,f9.4), i4)')l,pfaf,mnx,mxx,mny,mxy

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
             
       write (10,'(2I8, i3, f8.4, f8.2, f10.2, f8.4)' ) tid, cid, abm_int, peatf_r, gdp_r, hdm_r, field_cap(sc_com)  

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

END MODULE  process_hres_data


! -------------------------------------------------------------------------------------------------------------------------------

#ifdef DO_CLM45
!gmake FOPT='-DDO_CLM45'
!mpiifort -o ../bin/mk_clm45 -DsysLinux -DESMA64 -DHAS_NETCDF4 -DHAS_NETCDF3 -DH5_HAVE_PARALLEL -DNETCDF_NEED_NF_MPIIO  -DHAVE_SHMEM  -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/esmf   -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/esmf  -DDO_CLM45     -fPIC -fpe0 -fp-model source -heap-arrays 32 -assume noold_maxminloc  -align dcommons  -DEIGHT_BYTE -I../mod/ -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/netcdf -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -L../lib/ -I../include/ -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/netcdf -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/hdf5 mod_clm45_routines.F90 ../lib/libraster.a -lz -L../lib -lraster -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lnetcdff -lnetcdf -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lnetcdf -ljpeg -lmfhdf -ldf -lhdf5_hl -lhdf5 -lm -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lmfhdf -ldf -lsz -ljpeg -lgpfs -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lcurl -lssl -lcrypto -lssl -lcrypto -ldl -lz -lz -lrt -lm -lm -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lcurl -lssl -lcrypto -lssl -lcrypto -ldl -lz -lz -lrt -lm -lirc -ldl -lc -lpthread -lrt

    PROGRAM mk_clm45
  
      use mod_clm45_routines
      implicit none
      character*400        :: arg(3), gfiler
      integer              :: n , iargc, nc, nr

      if(iargc() /= 3) then
         print *, "Wrong Number of arguments: ", iargc()
         print *, "Usage : ./mkclm45 NC NR rst/gfile"
         stop
      endif

      do n=1,3
         call getarg(n,arg(n))
      enddo

      read(arg(1),*) NC
      read(arg(2),*) NR
      read(arg(3),'(a)') gfiler

      call CLM45_fixed_parameters (nc, nr, gfiler)
      call CLM45_clim_parameters  (nc, nr, gfiler)

    END PROGRAM mk_clm45

#endif

! ------------------------------------------------------------------------------------------------------

#ifdef DO_GIMMS
!gmake FOPT='-DDO_GIMMS'
!mpiifort -o ../bin/mk_gimms -DsysLinux -DESMA64 -DHAS_NETCDF4 -DHAS_NETCDF3 -DH5_HAVE_PARALLEL -DNETCDF_NEED_NF_MPIIO  -DHAVE_SHMEM  -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/esmf   -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/esmf -DDO_GIMMS  -O3 -qopt-report0 -ftz -align all -fno-alias -qno-offload -traceback     -fPIC -fpe0 -fp-model source -heap-arrays 32 -assume noold_maxminloc  -align dcommons  -DEIGHT_BYTE -I../mod/ -I/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/include/netcdf -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -L../lib/ -I../include/ gimms_ndvi.F90 ../lib/libraster.a -lz -L../lib -lraster -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lnetcdff -lnetcdf -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lnetcdf -ljpeg -lmfhdf -ldf -lhdf5_hl -lhdf5 -lm -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lmfhdf -ldf -lsz -ljpeg -lgpfs -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lcurl -lssl -lcrypto -lssl -lcrypto -ldl -lz -lz -lrt -lm -lm -L/discover/nobackup/projects/gmao/share/gmao_ops/Baselibs/v4.0.6_build1/x86_64-unknown-linux-gnu/ifort_15.0.2.164-intelmpi_5.0.3.048/Linux/lib -lcurl -lssl -lcrypto -lssl -lcrypto -ldl -lz -lz -lrt -lm -lirc -ldl -lc -lpthread -lrt  -L/gpfsm/dnb32/mbhat/GCC/install/gcc-4.6.3/lib/gcc/x86_64-unknown-linux-gnu/4.6.3 -lstdc++

    PROGRAM mk_gimms
  
      use gimms_ndvi
      implicit none
      character*400        :: arg(4), gfiler, vegfile
      integer              :: n , iargc, nc, nr

      if(iargc() /= 4) then
         print *, "Wrong Number of arguments: ", iargc()
         print *, "Usage : ./mkclm45 NC NR rst/gfile vegfile"
         stop
      endif

      do n=1,4
         call getarg(n,arg(n))
      enddo

      read(arg(1),*) NC
      read(arg(2),*) NR
      read(arg(3),'(a)') gfiler
      read(arg(4),'(a)') vegfile

      call ascat_r0        (nc, nr, gfiler, vegfile)
      call gimms_clim_ndvi (nc, nr, gfiler)

    END PROGRAM mk_gimms

#endif
