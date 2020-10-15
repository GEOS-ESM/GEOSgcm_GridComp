#define VERIFY_(A) if(A /=0)then;print *,'ERROR code',A,'at',__LINE__;call exit(3);endif


! setenv LD_LIBRARY_PATH /ford1/local/EL6-64/lib
  ! source /ford1/local/intel/latest/bin/compilervars.csh intel64
  ! ifort -o mkLAIrst2tiles.x mkLAIrst2tiles.F90 -L/ford1/local/EL6-64/lib -I/ford1/local/EL6-64/include -lnetcdff -lnetcdf

!  USAGE:

! 1) Derive tile-spaced LAI/Albedo for using analysis ready MODIS raw data @ 15/30 arc-sec lat/lon
!    LAI
!    bin/mkMODISrst2tiles.x  -m MCD15A2H -y YYYYDOY  -r YYYYMMDD  -t time_step   -l minutes_since_20020704-0000
!    Albedo
!    bin/MODISrst2tiles.x  -m MCD43GF -y YYYYDOY  -r YYYYMMDD  -t time_step   -l minutes_since_20000101-0000
! 2) Smooth LAI using a 3 time step window
!    bin/MODISrst2tiles.x  -s Y -g GRID_NAME  -a lai_data.YYYYDOY_previous  -b lai_data.YYYYDOY  -c lai_data.YYYYDOY_next -t tstep

!#####################################################################################

PROGRAM mkMODISrst2tiles 

  use netcdf
  implicit none

  character*300, parameter ::                                                 &
       RAWDIR = '/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/',       &
       BCSDIR = '/discover/nobackup/smahanam/MERRA3/MODIS_DVG/python/bcs/',   &
       OUTDIR = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/MODIS_DVG/IAV/',  &
       SMDIR  = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/MODIS_DVG/IAV_smoothed/'
!       RAWDIR = '/l_data/model_parameters/LAI/MODIS6/',                &
!       BCSDIR = '/ldas/sarith/LAI/MODIS6/Python/bcs/' ,                &
!       OUTDIR = '/l_data/model_parameters/LAI/MODIS6/'
  
  real, parameter :: dxy = 0.5
  integer,     parameter :: NC15 = 86400, NR15 = 43200, NC30 = 43200, NR30 = 21600, NGRIDS = 7
! *** NOTE WARNING *** DE720 must be the first
  character*16, parameter,dimension (NGRIDS) :: GRIDS =(/ &
       'DE720           ', &
       'PE90x540-CF/    ', &
       'PE180x1080-CF/  ', &
       'PE360x2160-CF/  ', &
       'PE720x4320-CF/  ', &
       'SMAP-EASEv2-M09/', &
       'SMAP-EASEv2-M36/' /)

  character*20 :: MODIS_NAME, YYYYDOY, GRID_NAME, file_out, file1, file2, file3
  integer,dimension(1) :: tstep, time, refdate
  integer              :: n, iargc, NT, NXT
  character*256        :: arg
  logical              :: file_exists = .false.
  real,    allocatable, dimension (:)  :: tile_lat
  real,    allocatable, dimension (:)  :: tile_lon
  real,    allocatable, dimension (:)  :: vec_data, count_data
  integer, allocatable, dimension (:,:):: hr_data1, hr_data2
  character*1 :: opt
  logical     :: smooth
  
  ! Read arguments
  ! --------------
  
  smooth = .false.
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
     case ('m')
        MODIS_NAME = trim(arg)
     case ('y')
        YYYYDOY = trim(arg)
     case ('r')
        read (arg,* ) refdate(1)
     case ('t')
        read (arg,* ) tstep(1)
     case ('l')
        read (arg,* ) time(1)
     case ('s')
        smooth = .true.
     case ('a')
        file1 = trim(arg)
     case ('b')
        file2 = trim(arg)
        file_out = file2
     case ('c')
        file3 = trim(arg)
     case ('g')
        GRID_NAME = trim(arg)
     case default
        print *, " USAGE   : ./mkMODISrst2tiles.x -m MODIS_NAME -g GRID_NAME -y YYYYDOY -r REFDATE -t time_step -l lag_in_mins -s Y (to smooth) -a file1 -b file2 -c file3"
        print *, " PROCESS : ./mkMODISrst2tiles.x -m MODIS_NAME -y YYYYDOY -r REFDATE -t time_step -l lag_in_mins"
        print *, " SMOOTH  : ./mkMODISrst2tiles.x -g GRID_NAME -t time_step -s Y -a file1 -b file2 -c file3"
        call exit(1)
     end select

     nxt = nxt + 1
     call getarg(nxt,arg)

  end do

  MODIS_NAME = trim(MODIS_NAME)//'.006'
  
  if (.not. smooth) then
 
     call grid2tile (trim(MODIS_NAME), trim(YYYYDOY), tstep, time, refdate)
     
  else
     
     inquire(file=trim(SMDIR)//trim(GRID_NAME)//'/'//trim(file_out),exist=file_exists)
     if(file_exists) stop
     call smooth_data (trim(GRID_NAME), trim(file_out), trim(file1), trim(file2), trim(file3),tstep)
          
  endif

  STOP
  
contains

  SUBROUTINE smooth_data (GRID_NAME, file_out, file1, file2, file3, tstep)

    implicit none

    character(*), intent (in) :: GRID_NAME, file_out, file1, file2, file3
    integer, intent (in),dimension(:)   :: tstep
    integer, allocatable, dimension (:) :: ii,jj
    character*300             :: filename
    logical                   :: file_exists = .false.
    character*3               :: DOY, DOY1, DOY3
    character*4               :: YYYY
    integer                   :: n, nx, ny, NCID, STATUS, k, j
    real, allocatable, dimension   (:) :: yc, vec_lai
    real, allocatable, dimension (:,:) :: lai_grid

    open (20, file = trim(SMDIR)//trim(GRID_NAME)//'/'//trim(file_out), &
            form = 'unformatted', action = 'write', status = 'unknown')

    open (10,file = trim(BCSDIR)//trim(GRID_NAME)//'/catchment.def',status='old',action='read',form='formatted')
    read (10,*) NT
    close(10, status = 'keep')

    allocate (vec_lai (nt))
    allocate (vec_data(nt))
    allocate (yc      (nt))
    
    vec_lai = 0.0001
    vec_data= 0.
    yc      = 0.
    do n = 1,3
       if (n == 1) filename = trim(file1)
       if (n == 2) filename = trim(file2)
       if (n == 3) filename = trim(file3)
       open (10, file=trim(OUTDIR)//trim(GRID_NAME)//'/'//trim(filename),&
            form = 'unformatted', action = 'read', status = 'old')
       read (10) vec_data
       close(10, status = 'keep')
       do j = 1, NT
          if(vec_data (j) > 0.) then
             vec_lai(j) = vec_lai(j) + vec_data (j)
             yc(j)      = yc(j)  + 1
          endif
       end do
    end do
    where (yc > 0.) vec_lai  = vec_lai /yc
    write (20) vec_lai
    close(20, status = 'keep')

    if(trim(GRID_NAME) == 'DE720') then
       allocate (ii (1:NT))
       allocate (jj (1:NT))
       nx = nint (360./dxy)
       ny = nint (180./dxy)
       allocate (lai_grid (1 : nx, 1 : ny)) 
       call ReadTileFile_RealLatLon(trim(BCSDIR)//trim(GRID_NAME)//'/tilfile', NT, ii,jj)
       lai_grid =  1.e+15
       do n = 1,NT
          lai_grid(ii(n), jj(n)) =  vec_lai(n)
       end do
       status = NF90_OPEN(trim(SMDIR)//'/ExtData/MCD15A2H.006_LAI_ExtData.nc4',NF90_WRITE,    ncid); VERIFY_(STATUS)
       STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'MODIS_LAI'),     lai_grid,start = (/1,1,tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
       STATUS = NF90_CLOSE (NCID)
    endif

    ! derive climatological mean file
    
    DOY = file_out(14:16)
    inquire(file=trim(SMDIR)//trim(GRID_NAME)//'/lai_data.YYYY'//DOY,exist=file_exists)

    if(.not.file_exists) then
       DOY1 = file1(14:16)
       DOY3 = file3(14:16)
       vec_lai = 0.0001
       yc      = 0.
       do n = 2002, 2020
          do k = 1,3
             write (yyyy, '(i4.4)') n
             file_exists = .false.             
             if (k == 1) then
                if (doy == '001') write (yyyy, '(i4.4)') n -1
                inquire(file=trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY1,exist=file_exists)
                filename = trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY1
             endif
             if (k == 2) then
                inquire(file=trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY ,exist=file_exists)
                filename = trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY
             endif
             if (k == 3) then
                if (doy == '361') write (yyyy, '(i4.4)') n +1
                inquire(file=trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY3,exist=file_exists)
                filename = trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY3
             endif

             if(file_exists) then
                open (10, file=trim(filename), &
                     form = 'unformatted', action = 'read', status = 'old')
                read (10) vec_data
                close(10, status = 'keep')
                do j = 1, NT
                   if(vec_data (j) > 0.) then
                      vec_lai(j) = vec_lai(j) + vec_data (j)
                      yc(j)      = yc(j)  + 1
                   endif
                end do
             endif
          end do
       end do
       where (yc > 0.) vec_lai = vec_lai / yc
       open (20, file = trim(SMDIR)//trim(GRID_NAME)//'/lai_data.YYYY'//DOY,  &
            form = 'unformatted', action = 'write', status = 'unknown')
       write (20) vec_lai
       close(20, status = 'keep')
    endif
    
  END SUBROUTINE smooth_data

  ! -------------------------------------------------------------------
  
  SUBROUTINE grid2tile (MODIS_NAME, YYYYDOY, tstep, time, refdate)

    implicit none
    character(*), intent (in)                  :: MODIS_NAME, YYYYDOY
    integer, intent (in),dimension(:)          :: tstep, time, refdate
    integer, allocatable, dimension (:)        :: ii,jj
    integer*2, allocatable, dimension (:,:)    :: int2_data
    real, allocatable, target, dimension (:,:) :: data_grid, nir_grid
    real, allocatable, dimension (:,:)         :: data_nir, data_ext
    real,    pointer,     dimension (:,:)      :: subset
    real,    pointer,     dimension (:)        :: nirdf, nir_count
    real,    allocatable,    dimension(:)      :: x, y
    integer :: n,i,j,k, ncid, nx, ny,status, ng, val_min, val_max
    INTEGER :: imn,imx,jmn,jmx,mval,d1,d2,l, ix, jx
    integer         :: doy, doy_tstep,i_yyyydoy
    real            :: sf

    read (YYYYDOY, '(i7.7)') i_yyyydoy
    doy = i_yyyydoy - (i_yyyydoy/1000)*1000
    doy_tstep = (doy - 1)/8 + 1

    nx = nint (360./dxy)
    ny = nint (180./dxy)
    allocate (data_grid (1 : nx, 1 : ny)) 
    allocate (x(1:nx))
    allocate (y(1:ny))
    data_grid  = -9999
    FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
    FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy
    
    ! READ MODIS data
    ! ---------------
   
    if(trim(MODIS_NAME) == 'MCD15A2H.006') then
       
       allocate (hr_data1  (1:NC15, 1: NR15))
       allocate (hr_data2  (1:NC15, 1: NR15))
       allocate (int2_data (1:NC15, 1: NR15))
       status = NF90_OPEN(trim(RAWDIR)//trim(MODIS_NAME)//'/'//YYYYDOY(1:4)//'/'// &
            trim(MODIS_NAME)//'_LAI_'//trim(YYYYDOY)//'.nc4',NF90_NOWRITE, ncid); VERIFY_(STATUS)
       
       STATUS  = NF90_GET_VAR (NCID,VarID(NCID,'Lai_500m'  ), int2_data) ; VERIFY_(STATUS)
       hr_data1 = int2_data
       STATUS  = NF90_GET_VAR (NCID,VarID(NCID,'FparLai_QC'), int2_data) ; VERIFY_(STATUS)
       hr_data2 = int2_data
       deallocate (int2_data)
       STATUS  = NF90_CLOSE (NCID)
       val_min = 0
       val_max = 100
       sf      = 0.01

    elseif (trim(MODIS_NAME) == 'MCD43GF.006') then

       allocate (nir_grid  (1:NX,   1:   NY)) 
       allocate (hr_data1  (1:NC30, 1: NR30))
       allocate (hr_data2  (1:NC30, 1: NR30))
       status = NF90_OPEN(trim(RAWDIR)//trim(MODIS_NAME)//'/'//YYYYDOY(1:4)//'/'// &
            trim(MODIS_NAME)//'_WSA_'//trim(YYYYDOY)//'.nc4',NF90_NOWRITE, ncid); VERIFY_(STATUS)
       STATUS = NF90_GET_VAR (NCID,VarID(NCID,'Albedo_Map_0.3-0.7'), hr_data1) ; VERIFY_(STATUS)
       STATUS = NF90_GET_VAR (NCID,VarID(NCID,'Albedo_Map_0.7-5.0'), hr_data2) ; VERIFY_(STATUS)
       STATUS = NF90_CLOSE (NCID)
       val_min = 0
       val_max = 1000
       sf      = 0.001       
       nir_grid= -9999
       
    else

       print *, 'Unknown MODIS name ...: ', trim(MODIS_NAME)
       stop

    endif
    
    ! read BCs til/rst files and initialize
    ! -------------------------------------
    
    do ng = 1, ngrids

       if(trim(MODIS_NAME) == 'MCD15A2H.006') call process_grid (NG, val_min, val_max, sf)
       if(trim(MODIS_NAME) ==  'MCD43GF.006') call process_grid (NG, val_min, val_max, sf, &
            nirdf = nirdf, nir_count = nir_count)
       
       if(trim(GRIDS(ng)) == 'DE720') then
          allocate (ii (1:NT))
          allocate (jj (1:NT))
          call ReadTileFile_RealLatLon(trim(BCSDIR)//trim(GRIDS(NG))//'/tilfile', NT, ii,jj)
       endif

       if(trim(GRIDS(ng)) == 'DE720') then
          do k = 1, NT
             if(count_data(k) >= 1.)   data_grid(ii(k), jj(k)) =  vec_data(k)
             if(trim(MODIS_NAME) ==  'MCD43GF.006') then
                if(nir_count(k) >= 1.) nir_grid (ii(k), jj(k)) =  nirdf   (k)
             endif
          enddo
       end if
       
       do n = 1, NT
          if(count_data(n) == 0.) then
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
                   vec_data (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                   exit
                endif
                l = l + 1
                NULLIFY (subset)
             end do
          endif
       end do

       if(trim(MODIS_NAME) ==  'MCD43GF.006') then
          do n = 1, NT
             if(nir_count(n) == 0.) then
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
                   subset => nir_grid(imn: imx,jmn:jmx)
                   
                   if(maxval(subset) > 0.) then 
                      nirdf (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                      exit
                   endif
                   l = l + 1
                   NULLIFY (subset)
                end do
             endif
          end do          
       endif
       
       ! Write output
       ! ------------
           
       if(trim(MODIS_NAME) == 'MCD15A2H.006') open (10, file = trim(OUTDIR)//'/'//trim(GRIDS(ng))//'/lai_data.'//trim(YYYYDOY), &
            form = 'unformatted', action = 'write', status = 'unknown')
       if(trim(MODIS_NAME) == 'MCD43GF.006' ) open (10, file = trim(OUTDIR)//'/'//trim(GRIDS(ng))//'/alb_data.'//trim(YYYYDOY), &
            form = 'unformatted', action = 'write', status = 'unknown')
       write (10) vec_data
       close (10,status = 'keep')

       if(trim(GRIDS(ng)) == 'DE720') then
          allocate (data_ext (1 : nx, 1 : ny))
          data_ext = 1.e+15
          do k = 1, NT
             data_ext(ii(k), jj(k)) =  vec_data(k)
          enddo
          if(trim(MODIS_NAME) == 'MCD15A2H.006') then
             status = NF90_OPEN(trim(OUTDIR)//'/ExtData/'//trim(MODIS_NAME)//'_LAI_ExtData.nc4',NF90_WRITE,    ncid); VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'MODIS_LAI'),     data_ext, start = (/1,1,tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'REFERENCE_DATE'),refdate, start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'time'),          time,    start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
             STATUS = NF90_CLOSE (NCID)
          endif
          if(trim(MODIS_NAME) == 'MCD43GF.006') then
             allocate (data_nir (1 : nx, 1 : ny))
             data_nir = 1.e+15
             do k = 1, NT
                data_nir(ii(k), jj(k)) =  nirdf (k)
             enddo
             status = NF90_OPEN(trim(OUTDIR)//'/ExtData/'//trim(MODIS_NAME)//'_ALBEDO_ExtData.nc4',NF90_WRITE,    ncid); VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'MODIS_VISDF'),     data_ext, start = (/1,1,tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'MODIS_NIRDF'),     data_nir, start = (/1,1,tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'REFERENCE_DATE'),refdate, start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
             STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'time'),          time,    start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
             STATUS = NF90_CLOSE (NCID)
             deallocate (data_nir)
          endif
          deallocate (data_ext)
       endif
       deallocate (tile_lat, tile_lon, vec_data, count_data)
       if(trim(MODIS_NAME) ==  'MCD43GF.006') deallocate (nirdf, nir_count)
       NT = 0
    end do
       
    deallocate (hr_data1,hr_data2, data_grid)
    if(trim(MODIS_NAME) ==  'MCD43GF.006') deallocate (nir_grid)
    return
  END SUBROUTINE grid2tile
  
    ! -----------------------------------------------------------------------
    
  SUBROUTINE process_grid (ng, val_min, val_max, sf, nirdf, nir_count)
      
    implicit none
    integer, intent (in)                    :: NG, val_min, val_max
    real,    intent (in)                    :: sf
    real, dimension (:), pointer, optional  :: nirdf, nir_count
    integer, allocatable, dimension (:,:)   :: tile_id
    integer                                 :: i,j,k, tindex, pfaf, nc, nr
    real                                    :: minlon,maxlon,minlat,maxlat,dxm, dym
    character(8)                            :: qc_str
    
    nc = size (hr_data1, 1)
    nr = size (hr_data1, 2)
    allocate (tile_id (1:NC30, 1:NR30))    
    open (20,file=trim(BCSDIR)//trim(GRIDS(ng))//'/rstfile',status='old',action='read',form='unformatted',convert='little_endian')    
    open (10,file=trim(BCSDIR)//trim(GRIDS(ng))//'/catchment.def',status='old',action='read',form='formatted')
    read (10,*) NT
    allocate (tile_lat (1:NT))
    allocate (tile_lon (1:NT))
    allocate (vec_data  (1:NT))
    allocate (count_data(1:NT))

    if (present (nirdf)) then
       allocate (nirdf    (1:NT))
       allocate (nir_count(1:NT))
       nirdf     = -9999.
       nir_count = 0.
    endif
    
    vec_data   = -9999.
    count_data = 0.
    
    do k = 1, NT
       read (10,*) tindex,pfaf,minlon,maxlon,minlat,maxlat
       tile_lon(k) = (minlon + maxlon)/2.
       tile_lat(k) = (minlat + maxlat)/2.
    end do
    close(10, status = 'keep')
    
    do j=1,NR30
       read(20)tile_id(:,j)
    end do
    close(20, status = 'keep')
    
    ! process DATA
    dxm = real(nc) /real(NC30) 
    dym = real(nr) /real(NR30)
    
    do j = 1, nr
       do i = 1,nc 
          if((hr_data1(i,j) >=val_min).and.(hr_data1(i,j) <= val_max)) then
             ! Table 5 of https://lpdaac.usgs.gov/documents/2/mod15_user_guide.pdf
             ! (SCF_QC '000' OR '001') AND ( MODLAND_QC = 0)
             !if((FparLai_QC(i,j) <= 62).and.(MOD(FparLai_QC(i,j),2) ==0)) then
             ! Cloud free 
             !                qc_str = bit2str(INT(FparLai_QC(i,j)))
             !                if(qc_str(4:5) /= '01') then 
             k = tile_id (ceiling(i/dxm), ceiling (j/dym))
             if((k >= 1).and.(k <= NT)) then
                if(vec_data(k) == -9999.) vec_data(k) = 0.
                vec_data(k)   = vec_data(k) + hr_data1(i,j)*sf
                count_data(k) = count_data(k) + 1.
             endif             
          endif
          if (present (nirdf)) then
             if((hr_data2(i,j) >=val_min).and.(hr_data2(i,j) <= val_max)) then
                k = tile_id (ceiling(i/dxm), ceiling (j/dym))
                if((k >= 1).and.(k <= NT)) then
                   if(nirdf(k) == -9999.) nirdf(k) = 0.
                   nirdf(k)   = nirdf(k) + hr_data2(i,j)*sf
                   nir_count(k) = nir_count(k) + 1.
                endif
             endif
          endif
       end do
    end do
      
    where (count_data > 0.) vec_data = vec_data/count_data
    if (present (nirdf)) where (nir_count > 0.) nirdf = nirdf/nir_count  
  end SUBROUTINE process_grid
  
  ! -----------------------------------------------------------------------

  SUBROUTINE HANDLE_ERR(STATUS, Line)
    
    INTEGER,      INTENT (IN) :: STATUS
    CHARACTER(*), INTENT (IN) :: Line
    
    IF (STATUS .NE. NF90_NOERR) THEN
       PRINT *, trim(Line),': ',NF90_STRERROR(STATUS)
       STOP 'Stopped'
    ENDIF

  END SUBROUTINE HANDLE_ERR
  
  ! *****************************************************************************

  integer function VarID (NCFID, VNAME) 
     
    integer, intent (in)      :: NCFID
    character(*), intent (in) :: VNAME
    integer                   :: status
    
    STATUS = NF90_INQ_VARID (NCFID, trim(VNAME) ,VarID)
    IF (STATUS .NE. NF90_NOERR) &
         CALL HANDLE_ERR(STATUS, trim(VNAME))  
    
  end function VarID

  ! *****************************************************************************

  subroutine ReadTileFile_RealLatLon (InCNTileFile, ntiles, ii, jj,mask)
   
    ! read *.til tile definition file, return *real* lat/lon for slow but accurate processing
    
    implicit none
    character(*), intent (in) :: InCNTileFile
    integer , intent (in)     :: ntiles
    integer, dimension (:), intent(inout)    :: ii,jj
    integer, optional, intent(IN) :: mask
    integer :: n,icnt,ityp, ntl, umask, i
    real    :: xval,yval, pf, ik, jk
         
    if(present(mask)) then
       umask = mask
    else
       umask = 100
    endif
   
    open(11,file=InCNTileFile, &
         form='formatted',action='read',status='old')
    read (11,*, iostat=n) Ntl
    	  
    do n = 1,7 ! skip header
       read(11,*)
    end do
    
    icnt = 0
    
    do i=1,Ntl
       read(11,*) ityp,pf,xval,yval, ik,jk
       if(ityp == umask) then
          icnt = icnt + 1
          if (icnt > ntiles) then
             print *,trim(InCNTileFile)
             print *,icnt,ntiles
             stop
          endif
          ii(icnt) = ik
          jj(icnt) = jk
         endif
      end do
   
      close(11)
          
    end subroutine ReadTileFile_RealLatLon

    ! *********************************************************

    function bit2str(b8) result(str8)
      
      integer, intent(in)  :: b8
      character(8)         :: str8
      character(1)         :: xstr
      integer              :: j,x,a
      
      x = b8
      
      do j=8,1,-1 
         a = MOD(x,2)
         write(xstr,'(i1)'), a
         str8(8 -(j-1):8 -(j-1)) = xstr
         x = x/2
      end do
      
    end function bit2str
    
  END PROGRAM mkMODISrst2tiles
