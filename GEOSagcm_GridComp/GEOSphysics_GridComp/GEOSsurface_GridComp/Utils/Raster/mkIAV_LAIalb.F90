#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

PROGRAM mkIAV_LAIalb

  use date_time_util  
  use leap_year
  use MAPL_ConstantsMod
  use module_sibalb, ONLY: sibalb
  use rmTinyCatchParaMod, only : regridrasterreal,Get_MidTime,Time_Interp_Fac,c_data

  implicit none
  include 'netcdf.inc'	
  character*300 :: GFILE = 'CF0090x6C_DE1440xPE0720-Pfafstetter'
  integer, parameter :: yearB = 1981, yearE = 2017, nc = 43200, nr = 21600
  integer     :: year, maxcat, tindex1,pfaf1,i,j, n
  REAL, ALLOCATABLE, dimension (:):: tile_lon, tile_lat
  real    :: minlat,maxlat,minlon,maxlon
  integer :: unit_l = 30, unit_v  =31, unit_n = 32, log_u = 99
  integer, dimension (:,:), allocatable, target :: tile_id
  integer, allocatable, dimension (:) :: vegcls
  character*8 :: YYYYMMDD1, YYYYMMDD2
  
  ! Reading number of cathment-tiles from catchment.def file
  !_________________________________________________________ 

 
  open (10,file='clsm/catchment.def',status='old',action='read',form='formatted')
  open (20,file='clsm/mosaic_veg_typs_fracs',status='old',action='read',form='formatted')
  
  read(10,*) maxcat
  allocate (tile_lon(1:maxcat)) 
  allocate (tile_lat(1:maxcat)) 
  allocate (vegcls  (1:maxcat))
  
  do n = 1, maxcat
     read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
     tile_lon(n) = (minlon + maxlon)/2.
     tile_lat(n) = (minlat + maxlat)/2.
     read (20,'(i8,i8,2(2x,i3),2(2x,f6.4))') tindex1,pfaf1,vegcls(n)
  end do

  close (10,status='keep')
  close (20,status='keep')
  
  ! Reading rst file
  !-----------------
  
  open (10,file='rst/'//trim(gfile)//'.rst',status='old',action='read',  &
       form='unformatted',convert='little_endian')
  allocate (tile_id    (1:nc,1:nr))           
  do j=1,nr
     read(10)tile_id(:,j)
  end do
  close (10,status='keep')
   
  ! writing GLASS LAI
  !
  open (log_u, file='clsm/IAV.log2', form='formatted',status='unknown', action = 'write')
  !write (log_u, *) 'Writing GLASS LAI'
  !
  !open (unit_l,file='clsm/lai.dat',  &
  !     form='unformatted',status='unknown',convert='little_endian')
  !do year = yearB, yearE
  !   call grid2tile_glass (year)  
  !end do
  !close(unit_l,status='keep')
      
  ! MODIS scale parameter files
  write (log_u, *) '  '
  write (log_u, *) 'MODIS scale parameters  '
  
  open (unit_l,file='clsm/lai.dat',status='old',action='read',form='unformatted', &
       convert='little_endian')      
  open (unit_v,file='clsm/visdf.dat',convert='little_endian', &
       action='write',status='unknown',form='unformatted')
  open (unit_n,file='clsm/nirdf.dat',convert='little_endian', &
       action='write',status='unknown',form='unformatted')
  
  call modis_scale_para_high 

  close(unit_l,status='keep')
  close(unit_v,status='keep')
  close(unit_n,status='keep')
  
contains

    SUBROUTINE grid2tile_glass (year)
      !
      ! Processing GLASS LAI (AVHRR or MODIS) and creating 8-day climatological data 
      !
      implicit none
      integer, intent (in) :: year
      integer, parameter :: N_lon_glass = 7200, N_lat_glass = 3600
      real, parameter :: dxy = 1.
      integer :: QSize
      integer :: n,i,j,k,ncid,i_highd,j_highd,nx_adj,ny_adj,ierr,nx,ny
      integer :: status,iLL,jLL,ix,jx,vid,nc_10,nr_10,n_tslices,d_undef,t,  &
           time_slice,time_slice_next,yr,mn,dd,yr1,mn1,dd1,i1,i2,tindex1,pfaf1
      character*100 :: fname,fout
      character*10 :: string
      character*2 :: VV,HH
      integer, allocatable, target,  dimension (:,:) :: net_data1
      real,  pointer, dimension (:,:) :: QSub
      real,  pointer, dimension (:,:) :: subset
      REAL, ALLOCATABLE, dimension (:):: vec_lai, count_lai, x, y !, distance
      real, allocatable, target, dimension (:,:) :: lai_grid, data_grid, data_grid2
      INTEGER ::imn,imx,jmn,jmx,mval,d1,d2,l, VarID
      character(len=4), dimension (:), allocatable :: MMDD, MMDD_next
      logical :: regrid
      REAL :: sf, dum,dist_save,tile_distance
      type (date_time_type) :: date_time_new,bf_lai_time,   &
           af_lai_time, date_time_this
      integer       ::  tileid_tile
      character*3   :: ddd
      integer, dimension (46)       :: &
           MODIS_DOYS = (/                              &
           1  ,  9, 17, 25, 33, 41, 49, 57, 65,         &
           73 , 81, 89, 97,105,113,121,129,137,         &
           145,153,161,169,177,185,193,201,209,         &
           217,225,233,241,249,257,265,273,281,         &
           289,297,305,313,321,329,337,345,353,361/)
      character*4 :: YYYY

      write (YYYY,'(i4.4)') year

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
      allocate (data_grid (1:NC,1:NR))
      allocate (data_grid2 (1 : N_lon_glass, 1 : N_lat_glass)) 

      do t =1,n_tslices
         
         time_slice = t
         yr = year
         yr1= year
         
         if(t == n_tslices) then 
            yr1 = year + 1
         endif
         
         read(mmdd(t),'(i2.2,i2.2)') mn,dd
         read(mmdd_next(t),'(i2.2,i2.2)') mn1,dd1
         
         date_time_this%year       = year
         date_time_this%month      = mn
         date_time_this%day        = dd
         date_time_this%hour   = 0            
         date_time_this%min    = 0            
         date_time_this%sec    = 0 
         call get_dofyr_pentad(date_time_this)                      

         write (ddd,'(i3.3)')  MODIS_DOYS(t) !date_time_this%dofyr
         !if (date_time_this%dofyr /= MODIS_DOYS (t)) then
         !   print *, date_time_this%dofyr, ' /= ', MODIS_DOYS(t)
         !   stop
         !end if
         ! Reading Interpolation or aggregation on to catchment-tiles
         
         vec_lai   = -9999.
         count_lai = 0.
         lai_grid  = -9999
         
         status  = NF_OPEN ('GLASSA/'//YYYY//'/GLASS01B02.V04.A'//YYYY//ddd//'.nc4', NF_NOWRITE, ncid) ; VERIFY_(STATUS)
         status  = NF_INQ_VARID (ncid,'LAI',VarID) ; VERIFY_(STATUS)
         status  = NF_GET_VARA_INT(ncid,VarID, (/1,1/),(/N_lon_glass, N_lat_glass/), net_data1) ; VERIFY_(STATUS)

         call RegridRasterReal(0.01*real(net_data1), data_grid)
         data_grid2 = 0.01*real(net_data1)

         status = NF_CLOSE(ncid)

         do j = 1,nr
            do i = 1, nc
               if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.MAXCAT)) then        
                  if((data_grid(i,j) >= 0.).and.(data_grid(i,j) <= 10.)) then 
                     if(vec_lai(tile_id(i,j)) == -9999.) vec_lai(tile_id(i,j)) = 0. 
                     vec_lai (tile_id(i,j)) = vec_lai (tile_id(i,j)) + data_grid(i,j)
                     count_lai (tile_id(i,j)) = count_lai (tile_id(i,j)) + 1.                     
                  endif
               endif
            end do
         end do

         write(unit_l) float((/yr,mn,dd,0,0,0,yr1,mn1,dd1,0,0,0,maxcat,1/))
         write(YYYYMMDD1,'(i8.8)')yr*10000 + mn*100 + dd
         write(YYYYMMDD2,'(i8.8)')yr1*10000 + mn1*100 + dd1
         write (log_u, *) '/GLASS01B02.V04.A'//YYYY//ddd
         
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
          write(unit_l)  vec_lai(:)
          write (log_u,*)YYYYMMDD1, ' ', YYYYMMDD2, minval(vec_lai), maxval(vec_lai)
       end do
 
       
       deallocate (net_data1, MMDD,MMDD_next,data_grid,data_grid2)
       deallocate (count_lai,x,y,lai_grid)
       deallocate (vec_lai) 

     END SUBROUTINE grid2tile_glass

!
!----------------------------------------------------------------------
!
  SUBROUTINE modis_scale_para_high 
    
    implicit none
    type (date_time_type) :: gf_green_time,af_green_time,end_time, &
              bf_lai_time,af_lai_time,date_time_new,bf_modis_time,   &
              af_modis_time, green_time_now, alb_time_now
    CHARACTER*20 :: version,resoln,continent
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: latt,lont,fr_gcm,fr_cat,tsteps,zth, slr,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    character*100 :: path,fname,fout,metpath
    integer :: n,ip
    integer :: yy,j,month
    real, allocatable, dimension (:) :: &
         modisvf, modisnf,albvf,albnf,lat,lon, &
         green,lai,lai_before,lai_after,grn_before,grn_after
    real, allocatable, dimension (:) :: &
         calbvf,calbnf
    character*300 :: ifile1,ifile2,ofile
    integer, dimension(12), parameter :: days_in_month_nonleap = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    integer :: day, hour, min ,secs_in_day,k
    real :: yr,mn,dy,yr1,mn1,dy1,dum, slice1,slice2
    logical :: save_sib = .false., update_modis = .false.

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

    ! MODIS Albedo files
    open (10,file='clsm/AlbMap.WS.8-day.tile.0.3_0.7.dat',&
         form='unformatted',convert='little_endian', &
         action='read',status='old')   
    open (11,file='clsm/AlbMap.WS.8-day.tile.0.7_5.0.dat',&
         form='unformatted',convert='little_endian', &
         action='read',status='old') 

    ! initialize LAI
    ! --------------
    
    read(unit_l) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(unit_l) lai_before
    lai_after = lai_before    
    call Get_MidTime(1980.-2001,12.,27.,yr-2001,mn,dy,bf_lai_time)
    call Get_MidTime(yr-2001,mn,dy,yr1-2001,mn1,dy1,af_lai_time)
    date_time_new     = af_lai_time
    date_time_new%day = 1
    end_time          = af_lai_time
    end_time%year     = yearE + 1
    end_time%month    = 1
    end_time%day      = 1

    write (log_u, *) 'START TIME : ', DT2STR(date_time_new)
    write (log_u, *) 'END TIME : ', DT2STR(end_time)
    write (log_u, *) 'LAI : ', DT2STR(bf_lai_time), ' ',DT2STR(af_lai_time) 

    ! Albedo files
    ! ------------
    
    write(unit_v) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
    write(unit_n) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.

    read(10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read (10) modisvf (:)
    read (11) modisnf (:)
     
    if (save_sib) then
       open (20,file='clsm/sib_visdf.dat',convert='little_endian', &
             action='write',status='unknown',form='unformatted')
       open (21,file='clsm/sib_nirdf.dat',convert='little_endian', &
             action='write',status='unknown',form='unformatted')  
       write(20) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
       write(21) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1. 
    endif
 
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,bf_modis_time)

    read(10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1

    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_modis_time)
    
    rewind (10)
    rewind (11)

    read(10) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read (10) modisvf (:)
    read (11) modisnf (:)
    
    alb_time_now     = af_modis_time
    alb_time_now%day = 1
    
    write (log_u, *) 'ALB : ', DT2STR(bf_modis_time), ' ',DT2STR(af_modis_time)
    
    open (41,file='clsm/green.dat',status='old',action='read',form='unformatted', &
         convert='little_endian')
    read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(41) grn_before
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,gf_green_time)
    
    read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
    read(41) grn_after
    call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)
    write (log_u, *) 'GREEN : ', DT2STR(gf_green_time), ' ',DT2STR(af_green_time)
    green_time_now = af_green_time
    green_time_now%day = 1

    calbvf   =0.
    calbnf   =0.
    albvf    =0.
    albnf    =0.
    tsteps   =0.    
    modisvf  =0.
    modisnf  =0.      
    
    do while (datetime_le_refdatetime(date_time_new,end_time))

       ! Check LAI
       if (datetime_le_refdatetime(date_time_new,af_lai_time)) then

       else
          update_modis = .true.
          ! LAI first
          read(unit_l,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
          if(ierr == 0) then 
             lai_before = lai_after
             read(unit_l) lai_after
             bf_lai_time = af_lai_time
             call Get_MidTime(yr-2001,mn,dy,yr1-2001,mn1,dy1,af_lai_time)
             write (log_u, *) 'UPDATING LAI FILE : ', DT2STR(date_time_new),' ', DT2STR(bf_lai_time), ' ',DT2STR(af_lai_time)
             
          else
             write (log_u, *) 'END OF LAI FILE : ',DT2STR(date_time_new)
             rewind(unit_l) 
             read(unit_l) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
             read(unit_l) lai_before
             call Get_MidTime(yr-2001,mn,dy,yr1-2001,mn1,dy1,bf_lai_time)
             read(unit_l) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
             read(unit_l) lai_after
             call Get_MidTime(yr-2001,mn,dy,yr1-2001,mn1,dy1,af_lai_time)

             ! if(date_time_new%dofyr < bf_lai_time%dofyr) then
             !    do while ((date_time_new%dofyr >  af_lai_time%dofyr))
             !       lai_before = lai_after
             !       bf_lai_time = af_lai_time
             !       read(unit_l) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
             !       read(unit_l) lai_after
             !       call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_lai_time)
             !    end do
             !endif
          endif          
       endif
       
       ! interpolate LAI         
       call Time_Interp_Fac (date_time_new, bf_lai_time, af_lai_time, slice1, slice2)
       lai    = (slice1*lai_before + slice2*lai_after)

       ! check Green date
         
       if (datetime_le_refdatetime(green_time_now,af_green_time)) then
          
       else
          
          read(41,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
          if(ierr == 0) then
             if(yr == 1) then
                grn_before = grn_after
                gf_green_time = af_green_time                  
                read(41) grn_after
                call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)
                write (log_u, *) 'UPDATING GREEN FILE : ', DT2STR(green_time_now),' ', DT2STR(gf_green_time), ' ',DT2STR(af_green_time)
             else
                rewind (41)
 !               read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
 !               read(41) grn_before
                read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read(41) grn_before
                call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,gf_green_time)
                
                read(41) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read(41) grn_after
                call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_green_time)
                green_time_now = af_green_time
                green_time_now%day =  date_time_new%day
                write (log_u, *) 'REWINDING GREEN FILE : ',DT2STR(green_time_now),' ',DT2STR(gf_green_time), ' ',DT2STR(af_green_time)
                write (log_u, *) 'Local time at GREEN rewind : ',DT2STR(date_time_new)
             endif
          else
             print *, 'End of Green file'
          endif
       endif
       
       ! interpolate green 
       call Time_Interp_Fac (green_time_now, gf_green_time, af_green_time, slice1, slice2)
       green  = (slice1*grn_before + slice2*grn_after)

       ! run SiB
       call sibalb (                  &
            MAXCAT,vegcls,lai,green,  &
            albvf, albnf)
       
       calbvf = calbvf + albvf
       calbnf = calbnf + albnf         
       tsteps = tsteps + 1.
       
       ! update time
       call augment_date_time( 86400, date_time_new ) 
       call augment_date_time( 86400, alb_time_now  )
       call augment_date_time( 86400, green_time_now)
       
       ! check LAI/alb time slice
       
       if (update_modis) then

          update_modis = .false.
          calbvf = calbvf/tsteps
          calbnf = calbnf/tsteps
          
          modisvf = modisvf/calbvf
          modisnf = modisnf/calbnf
          
          do n =1, maxcat
             if(modisvf(n).le.0.) modisvf(n) = 1.
             if(modisnf(n).le.0.) modisnf(n) = 1.
             if(modisvf(n).gt.100)modisvf(n)= 1.
             if(modisnf(n).gt.100)modisnf(n)= 1.
          enddo
          
          if (save_sib) then
             write (20) calbvf (:)
             write (21) calbnf (:)               
          endif
          
          write (unit_v) modisvf (:)
          write (unit_n) modisnf (:)
          
          ! Now update LAlb
          ! ------------------
          write(unit_v) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
          write(unit_n) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
             
          if (save_sib) then
             write(20) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1.
             write(21) yr,mn,dy,0.,0.,0.,yr1,mn1,dy1,0.,0.,0.,float(maxcat),1. 
          endif
          
          calbvf   =0.
          calbnf   =0.
          albvf    =0.
          albnf    =0.
          tsteps   =0.
          
          ! update alb file
 
          read(10,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1       
          
          if(ierr == 0) then
             if(yr == 1) then
                read(11) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read (10) modisvf (:)
                read (11) modisnf (:)                 
                bf_modis_time = af_modis_time
                call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_modis_time)
                
                write (log_u, *) 'UPDATING MODIS FILE : ', DT2STR(alb_time_now),' ', DT2STR(bf_modis_time), ' ',DT2STR(af_modis_time)
             else
                rewind (10)
                rewind (11)
!                read(10,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
!                read(11,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
!                read (10) modisvf (:)
!                read (11) modisnf (:)
                read(10,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read(11,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read (10) modisvf (:)
                read (11) modisnf (:)
                call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,bf_modis_time)
                read(10,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read(11,IOSTAT=ierr) yr,mn,dy,dum,dum,dum,yr1,mn1,dy1
                read (10) modisvf (:)
                read (11) modisnf (:)
                call Get_MidTime(yr,mn,dy,yr1,mn1,dy1,af_modis_time)
                alb_time_now = af_modis_time
                alb_time_now%day = date_time_new%day
                write (log_u, *) 'REWINDING MODIS FILE : ',DT2STR(alb_time_now),' ',DT2STR(bf_modis_time), ' ',DT2STR(af_modis_time)
                write (log_u, *) 'Local time at MODIS rewind : ',DT2STR(date_time_new)
             endif
          else
             print *, 'End of Alb file'
          endif
       endif
    end do
    
    deallocate (modisvf,modisnf,albvf,albnf)
    deallocate (green,lai)
    deallocate (calbvf,calbnf)
    deallocate (lai_before,grn_before, lai_after,grn_after)
    
    close (10, status='keep')
    close (11, status='keep')
    close (41, status='keep')
    close (unit_l, status='keep')
    close (unit_v, status='keep')
    close (unit_n, status='keep')
    if (save_sib) then
       close (20, status='keep')
       close (21, status='keep')
    endif
    
END SUBROUTINE modis_scale_para_high

CHARACTER*8 FUNCTION DT2STR (date_time)

  type (date_time_type) :: date_time
  write(DT2STR,'(i8.8)')date_time%year*10000 + date_time%month*100 + date_time%day
END FUNCTION DT2STR
 
END PROGRAM mkIAV_LAIalb
