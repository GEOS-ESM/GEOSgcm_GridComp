PROGRAM mkSMAPTilesPara
!     This program constructs land and lake tiles for the SMAP-EASE-M09 and M36 grids (just set MGRID) 
!         for CLSM implementation.
!     f90 -c create_smap_tiles.f90
!     f90 -c smapconv.f
!     f90 -o create_smap_tiles create_smap_tiles.o smapconv.o
!
      use easeV1_conv
      use rmTinyCatchParaMod
      use process_hres_data

      implicit none
      integer nc,nr,i,j,icount(21),ig,jg,i0,iop,n,d1,d2,j1
      parameter (nc=8640,nr=4320)
      integer,allocatable :: tileid_index(:,:),catid(:,:),catid_index(:,:)
      integer,allocatable, dimension (:) :: land_id,water_id,ice_id,pfaf_array
      real, allocatable, dimension (:) :: smap_grid_area,tile_area,lat_c,lon_c 
      integer l,imn,imx,jmn,jmx,mval,l_index,i_index,w_index,typ,pfaf,cindex
      character(3) :: easegrid
      integer*1,allocatable :: veg(:,:)
      real :: clat, clon, r_smap, s_smap, smap_convert, da!, smap_inverse, & 
!      	   ezlh_convert, ezlh_inverse, easeV1_convert, easeV1_inverse
      real :: fr_gcm
      integer :: ind_col, ind_row, status 
      REAL (kind=8), PARAMETER :: RADIUS=6371000.,pi=3.1415926535898
      character*100 :: veg_class (12)
      character*5 :: MGRID
      character*100 :: gfile,gtopo30,pathout
      integer :: nc_smap,nr_smap, N_args, iargc
      real :: EASE_grid_area, CELL_km
      real*4, dimension (:,:), allocatable ::    q0
      REAL, dimension (:), allocatable :: tile_ele, tile_area_land  
      REAL :: dx,dy,d2r,lats,mnx,mxx,mny,mxy,sum1,sum2,jgv
      character(40) :: arg
      character*200 :: tmpstring	      
      
      N_args = iargc()

      if(N_args < 1) then
        print *,'USAGE : bin/mkSMAPTiles -smap_grid MXX'
	print *,'Allowed SMAP grids are: M01 M03 M09 M25 M36'
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
      
   gtopo30 = 'data/CATCH/srtm30_withKMS_2.5x2.5min.data'

      allocate(q0(nc,nr)) 
      dx  = 360._8/nc
      dy  = 180._8/nr
      d2r = PI/180._8

      open (10,file=trim(gtopo30),form='unformatted',status='old',convert='little_endian')
      read (10) q0
      close (10,status='keep')
      pathout ='clsm/'

      if (trim(MGRID) == 'M09') then

         nc_smap = 3852
         nr_smap = 1632
         gfile = 'SMAP_EASE_'//trim(MGRID)//'_3852x1632'
         EASE_grid_area = 81.001812568020028

      elseif(trim(MGRID) == 'M36') then

         nc_smap = 963
         nr_smap = 408
         gfile = 'SMAP_EASE_'//trim(MGRID)//'_963x408'
         EASE_grid_area = 1296.029001087600 
 
      elseif(trim(MGRID) == 'M25') then
      
         nc_smap = 1383
         nr_smap = 586
         gfile = 'SMAP_EASE_M25_1383x586'
         EASE_grid_area = 628.3808529246343824
       
       else if (trim(MGRID) .eq. 'M03') then ! SMAP  3 km grid
          CELL_km = 3.00003356589     ! nominal cell size in kilometers
          nc_smap = 11556
          nr_smap = 4896
       	  gfile = 'SMAP_EASE_M03_11556x4896'
          EASE_grid_area = CELL_km*CELL_km

       else if (trim(MGRID) .eq. 'M01') then ! SMAP  1 km grid
          CELL_km = 1.00001118863     ! nominal cell size in kilometers
          nc_smap = 34668
          nr_smap = 14688
       	  gfile = 'SMAP_EASE_M01_34668x114688'
          EASE_grid_area = CELL_km*CELL_km
       
      else  !

         print *,'Unknown SMAP Grid stopping..'
         stop

      endif
            
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
    
      allocate(veg         (1:nc,1:nr))
      allocate(tileid_index(1:nc,1:nr))
      allocate(catid       (1:nc,1:nr))
      allocate(catid_index (1:nc,1:nr))
      allocate(land_id    (1:16330000))
      allocate(water_id   (1:16330000))
      allocate(ice_id     (1:16330000))
    
      da = radius*radius*pi*pi/24./24./180./180./1000000.    
    
      tileid_index=0
      land_id =0
      water_id=0
      ice_id=0
          
      ! reading SiB2 land cover classification data - the origin of the 
      ! 2.5'x2.5' vegetation raster file is global 1min IGBP data 
      ! (ftp://edcftp.cr.usgs.gov/pub/data/glcc/globe/latlon/sib22_0.leg)

      open (10,file='data/CATCH/sib22.5_v2.0.dat', &
           form='unformatted', &
	action='read', convert='big_endian',status='old')

      READ(10)veg

      close (10,status='keep')

!      ! reading 2.5'x2.5' global raster file of Pfafstetter Catchment IDs
!      ! In this version, the dateline has been overlaid over the catchments those straddle 
!      ! across. The numbers contain for
!      !  1 global ocean catchment                : Pfafstetter ID 0
!      !  36716 global land catchments            : Pfafstetter IDs 1000-5999900
!      !  1 global inland water (lakes) catchment : Pfafstetter ID 6190000
!      !  1 global ice catchment                  : Pfafstetter ID 6200000
! 
      open (10,file='data/CATCH/global.cat_id.catch.DL', form='formatted', &
	action='read', status='old')!

      do j=1,nr
         read(10,*)(catid(i,j),i=1,nc)
      end do

      close (10,status='keep')
      
      print *,'Read global.cat_id.catch.DL' 
      print *,'Min and Max of Pfafstetter IDs:', minval(catid),maxval(catid)

      ! reading the 2.5'x2.5' global raster file of tile indices for the 
      !  above Pfafstetter Catchments
      !  1 global ocean catchment                : tile_index 36719
      !  36716 global land catchments            : tile_index 1-36716
      !  1 global inland water (lakes) catchment : tile_index 36717
      !  1 global ice catchment                  : tile_index 36718

      open (10,file='data/CATCH/'  &
          //'PfafstatterDL.rst', form='unformatted',        &
          action='read',convert='little_endian', status='old')

      do j=1,nr
         read(10)(catid_index(i,j),i=1,nc)
      end do

      close (10,status='keep')

      print *,'Read PfafstatterDL.rst' 
      print *,'Min and Max of tile indices:',minval(catid_index),maxval(catid_index)

      ! Deriving SiB2 vegetation classification histogram

	icount=0

      do j=1,nr
         do i=1,nc
         icount(veg(i,j)+1)= icount(veg(i,j)+1) +1
         end do
      end do

      open (10,file='clsm/vegetation.hst2',form='formatted')

      do i=1,12
         write (10,'(i2,i10,1x,a100)')i-1,icount(i),veg_class(i)
      end do

      write (10,*)'-----------------------'
      write(10,*)'           ',sum(icount)
      close (10,status='keep')
      write (*,*)'End reading Sib2 vegetation classification'
!
! While looping through the 2.5'x2.5' grid, this section counts presence of 
! land, ice and water on the SMAP GIRD.
! Each SMAP grid cell is assigned with an ID = ind_row*10000 +  ind_col 
!         ind_col, ind_row are overlying SMAP grid cell indices
! This is just the prelimiminery assessment in the process of assigning separate  
!     tiles for land, water and ice fractions within the SMAP Grid cell
! The program checks all the underlying 2.5x2.5 cells for each SMAP EASE grid cell separately
! and counts the number of SMAP EASE water, land and ice grid cells were observed.
!
      do i = 1 ,nc

         clon = -180. + float(i-1)*2.5/60. + 1.25/60.
         
         do j =nr ,1 ,-1

            clat = -90. + float(j-1)*2.5/60. + 1.25/60.

	    call easeV1_convert(trim(MGRID), clat, clon, r_smap, s_smap)

!            if(trim(MGRID) == 'Ml') then
!	    	status = ezlh_convert(trim(MGRID), clat, clon, r_smap, s_smap)
!	    else                       
!            	status = smap_convert(trim(MGRID), clat, clon, r_smap, s_smap)
!	    endif

            ind_col = nint(r_smap) + 1 
            ind_row = nint(s_smap) + 1

            if((ind_row.ge.1).and.(veg(i,j).ge.1).and.(ind_row.le.1632)) then
               l=  ind_row*10000 +  ind_col

               if(veg(i,j)==10)  water_id(l) = 1
               if(veg(i,j)==11)  ice_id  (l) = 1
               if(veg(i,j).lt.10)land_id (l) = 1
            endif
           
         end do
      end do

      print *,'# of Land  pixels in SMAP: ',sum (land_id)
      print *,'# of water pixels in SMAP: ',sum (water_id)
      print *,'# of ice   pixels in SMAP: ',sum (ice_id)

      l_index=0
      w_index=sum (land_id)
      i_index=sum (land_id) + sum (water_id)

      land_id = 0
      water_id= 0
      ice_id  = 0

! While looping through the 2.5'x2.5' grid, this section derives land, ice and water tiles.
! Each SMAP grid cell is assigned with an ID = ind_row*10000 +  ind_col 
!         ind_col, ind_row are overlying SMAP grid cell indices
! Based on the above calculations: 
!         l_index Grid cells have land fractions (sum(land_id)) 
!         w_index SMAP Grid cells have inland water fractions (sum(water_id))
!         i_index SMAP Grid cells have ice fractions (sum(ice_id))
! hence, tile_index        1                     to l_index                     represent land tiles  
!         tile_index       l_index +1            to l_index + w_index           represent water (lakes) tiles  
!         tile_index       l_index + w_index +1  to l_index + w_index + i_index represent ice tiles
! global 2.5'x2.5' array of tileid_index(8640,4320) contains corresponding tile_index values which 
!        is derived in the below loop

      do i = 1 ,nc
 
         clon = -180. + float(i-1)*2.5/60. + 1.25/60.
         
         do j =nr ,1 ,-1

            clat = -90. + float(j-1)*2.5/60. + 1.25/60.

	    call easeV1_convert(trim(MGRID), clat, clon, r_smap, s_smap)
!            if(trim(MGRID) == 'Ml') then                       
!	        status = ezlh_convert(trim(MGRID), clat, clon, r_smap, s_smap)
!	    else
!                status = smap_convert(trim(MGRID), clat, clon, r_smap, s_smap)
!	    endif

            ind_col = nint(r_smap) + 1 
            ind_row = nint(s_smap) + 1

            if((ind_row.ge.1).and.(veg(i,j).ge.1).and.(ind_row.le.1632)) then
               
               l=  ind_row*10000 +  ind_col
               
               if(veg(i,j)==10) then
                  if(water_id(l)==0) then
                     w_index = w_index + 1
                     water_id(l) = w_index  
                     tileid_index(i,j)= water_id(l)
                  else
                     tileid_index(i,j)= water_id(l)
                  endif
               endif
               
               if(veg(i,j)==11) then
                  if(ice_id(l)==0) then
                     i_index = i_index + 1
                     ice_id  (l) = i_index
                     tileid_index(i,j)= ice_id  (l)  !i_index
                  else
                     tileid_index(i,j)= ice_id  (l)  !i_index
                  endif
               endif
               
               if(veg(i,j).lt.10) then
                  if(land_id(l)==0) then
                     l_index = l_index + 1     
                     land_id (l) = l_index        
                     tileid_index(i,j)= land_id (l) !1-l_index
                  else
                     tileid_index(i,j)= land_id (l) !1-l_index
                  endif
               endif
            endif
         end do
      end do

!      print *,l_index,w_index -l_index ,i_index - w_index
!      stop

      deallocate(land_id )
      deallocate(water_id)
      deallocate(ice_id  )
      allocate(smap_grid_area(1:16330000))
      allocate(land_id(1:i_index))
      allocate(lat_c(1:i_index))      
      allocate(lon_c(1:i_index))
      allocate(tile_area(1:i_index))
      allocate(water_id(1:i_index))
      allocate(ice_id  (1:i_index))
      allocate(pfaf_array  (1:i_index))

      land_id = 0
      water_id= 0
      ice_id  = 0
      smap_grid_area = 0.
      tile_area=0.
      lat_c = 0.
      lon_c = 0.

      do i = 1 ,nc    

         clon = -180. + float(i-1)*2.5/60. + 1.25/60.
         
         do j =nr ,1 ,-1

            clat = -90. + float(j-1)*2.5/60. + 1.25/60.

	    call easeV1_convert(trim(MGRID), clat, clon, r_smap, s_smap)

!            if(trim(MGRID) == 'Ml') then
!	        status = ezlh_convert(trim(MGRID), clat, clon, r_smap, s_smap)                                   
!	    else			   
!            	status = smap_convert(trim(MGRID), clat, clon, r_smap, s_smap)
!	    endif

            ind_col = nint(r_smap) + 1 
            ind_row = nint(s_smap) + 1

            if((ind_row.ge.1).and.(ind_row.le.1632)) then

               l=  ind_row*10000 +  ind_col
               if(veg(i,j).ge.1) then
               land_id(tileid_index(i,j)) = l
               ice_id (tileid_index(i,j)) = j*10000 + i  ! just  recycling the array ice_id, from this point ice_id contains infor to derive pfaf equivalent for SMAP cells
               tile_area(tileid_index(i,j))= tile_area(tileid_index(i,j)) + &
                    da*cos((-90.+float(j)/24. -1./48.)*pi/180.)
               endif
               smap_grid_area(l) = smap_grid_area(l) +             &
                    da*cos((-90.+float(j)/24. -1./48.)*pi/180.)
!               lat_c(tileid_index(i,j))=lat_c(tileid_index(i,j)) + & ! not being used
!                    (-90.+float(j)/24. -1./48.)
!               lon_c(tileid_index(i,j))=lon_c(tileid_index(i,j)) + & ! not being used
!                    (-180.+float(i)/24. -1./48.)              
!               water_id(tileid_index(i,j)) = water_id(tileid_index(i,j)) + 1 ! not being used
!               water_id(tileid_index(i,j)) = water_id(tileid_index(i,j)) + 1
               
            endif
         end do
      end do

 !     print *,minval(land_id),maxval(land_id)
      print *,'Creating ...', trim(gfile)//'rst'
      open (10, file ='rst/'//trim(gfile)//'.rst',form='unformatted',status='unknown',  &
           action='write')

      do j=1,nr
         write(10)(tileid_index(i,j),i=1,nc)
      end do

      close (10,status='keep')

      open (10, file ='til/'//trim(gfile)//'.til',form='formatted',status='unknown',action='write')
      write (10,*)i_index
      write (10,*)1
      write (10,*)'SMAP-EASE-'//trim(MGRID)
      write (10,*)nc_smap
      write (10,*)nr_smap
      write (10,*)'NO-OCEAN'
      write (10,*) -9999
      write (10,*) -9999      

      do l=1,i_index
         if (l <= l_index) typ = 100
         if ((l > l_index).and.(l <= w_index)) typ =19
         if (l > w_index) typ = 20

         pfaf  = catid(ice_id(l)-10000*(ice_id(l)/10000),ice_id(l)/10000)
         cindex= catid_index(ice_id(l)-10000*(ice_id(l)/10000),ice_id(l)/10000)
         ig    = land_id(l)-10000*(land_id(l)/10000)
         jg    = land_id(l)/10000
         pfaf_array(l) = pfaf

!         clat  = lat_c(l)/real(water_id(l))
!         clon  = lon_c(l)/real(water_id(l)) 

	  call easeV1_inverse (trim(MGRID), real(ig-1),  real(jg-1), clat, clon)

!	 if(trim(MGRID) == 'Ml') then
!             status = ezlh_inverse (trim(MGRID), real(ig-1),  real(jg-1), clat, clon)
!	 else
!             status = smap_inverse (trim(MGRID), real(ig-1),  real(jg-1), clat, clon) 
!         endif

         fr_gcm= tile_area(l)/smap_grid_area(jg*10000 +  ig)
         tile_area(l) = fr_gcm*EASE_grid_area
         write(10,'(i10,i9,2f10.4,2i5,f16.12,i10,f13.4,i8)') &
              typ,pfaf,clon,clat,ig-1,jg-1,fr_gcm ,cindex !,fr_gcm*EASE_grid_area
      end do

      close (10,status='keep')      
!      stop

! tile elevation

    allocate(tile_ele(1:l_index))
    allocate(tile_area_land(1:l_index))

    tile_ele = 0.
    tile_area_land = 0.    

    do j=1,nr

       lats = -90._8 + (j - 0.5_8)*dy

       do i=1,nc          
          if((tileid_index(i,j) > 0).and.(tileid_index(i,j) <= l_index))then
             tile_ele(tileid_index(i,j)) = tile_ele(tileid_index(i,j)) + q0(i,j)*   &	
                  (sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)
             tile_area_land(tileid_index(i,j)) = tile_area_land(tileid_index(i,j)) +  &	
                  (sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)		 
          endif
       enddo
    enddo

    tile_ele = tile_ele/tile_area_land

    ! adjustment Global Mean Topography to 614.649 (615.662 GTOPO 30) m
    ! ---------------------------
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

!
! Now catchment.def
!

    open (10,file='clsm/catchment.def',  &
         form='formatted',status='unknown')

    write (10,*)l_index

    do j=1,l_index

      ig    = land_id(j)-10000*(land_id(j)/10000)
      jg    = land_id(j)/10000

      call easeV1_inverse (trim(MGRID), real(ig-1),real(jg-1), clat, clon) 

!      if(trim(MGRID) == 'Ml') then       
!          status = ezlh_inverse (trim(MGRID), real(ig-1),real(jg-1), clat, clon) 
!      else
!          status = smap_inverse (trim(MGRID), real(ig-1),real(jg-1), clat, clon) 
!      endif

      mnx = clon - 180./real(nc_smap)
      mxx = clon + 180./real(nc_smap)
      
      jgv = real(jg-1) + 0.5
      
      call easeV1_inverse (trim(MGRID), real(ig-1),jgv, clat, clon)

!      if(trim(MGRID) == 'Ml') then 
!           status = ezlh_inverse (trim(MGRID), real(ig-1),jgv, clat, clon)
!      else
!           status = smap_inverse (trim(MGRID), real(ig-1),jgv, clat, clon) 
!      endif
      mny = clat
      
      jgv = real(jg-1) - 0.5

      call easeV1_inverse (trim(MGRID), real(ig-1),jgv, clat, clon) 

!      if(trim(MGRID) == 'Ml') then
!           status = ezlh_inverse (trim(MGRID), real(ig-1),jgv, clat, clon) 
!      else
!           status = smap_inverse (trim(MGRID), real(ig-1),jgv, clat, clon) 
!      endif

      mxy = clat

      write (10,'(i8,i8,5(2x,f9.4))')j,pfaf_array(j),mnx,mxx,mny,mxy,tile_ele(j)

   end do

      ! create Grid2Catch transfer file
      ! -------------------------------

      CALL CREATE_ROUT_PARA_FILE (NC, NR, trim(gfile), MGRID=MGRID)  

! now run mkCatchParam
   tmpstring = 'bin/mkCatchParam_openmp -e EASE -g '//trim(gfile)
   print *,trim(tmpstring)

   call system (tmpstring)   

   END PROGRAM mkSMAPTilesPara
