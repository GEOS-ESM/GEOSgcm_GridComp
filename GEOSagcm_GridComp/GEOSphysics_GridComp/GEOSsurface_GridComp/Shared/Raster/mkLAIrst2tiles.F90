#define VERIFY_(A) if(A /=0)then;print *,'ERROR code',A,'at',__LINE__;call exit(3);endif
  ! setenv LD_LIBRARY_PATH /ford1/local/EL6-64/lib
  ! source /ford1/local/intel/latest/bin/compilervars.csh intel64
  ! ifort -o mkLAIrst2tiles.x mkLAIrst2tiles.F90 -L/ford1/local/EL6-64/lib -I/ford1/local/EL6-64/include -lnetcdff -lnetcdf

!  USAGE:

! 1) Derive tile-spaced LAI for using analysis ready MODIS raw data @ 15 arc-sec lat/lon
!    bin/mkLAIrst2tiles.x  YYYYDOY  YYYYMMDD  time_step   minutes_since_20020708-0000
! 2) Smmoth using 3 time step window
!    bin/mkLAIrst2tiles.x  GRID_NAME  lai_data.YYYYDOY  lai_data.YYYYDOY_previous  lai_data.YYYYDOY  lai_data.YYYYDOY_next

PROGRAM mkLAIrst2tiles 

  use netcdf
  implicit none

  character*300, parameter ::                                                 &
       RAWDIR = '/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/',       &
       BCSDIR = '/discover/nobackup/smahanam/MERRA3/MODIS_DVG/python/bcs/',   &
       OUTDIR = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/MCD15A2H.006/IAV/',  &
       SMOOTH = '/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/MCD15A2H.006/IAV_smoothed/'
!       RAWDIR = '/l_data/model_parameters/LAI/MODIS6/',                &
!       BCSDIR = '/ldas/sarith/LAI/MODIS6/Python/bcs/' ,                &
!       OUTDIR = '/l_data/model_parameters/LAI/MODIS6/'

  integer,     parameter :: NC = 86400, NR = 43200, IM = 43200, JM = 21600, NGRIDS = 7
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
  integer              :: n, iargc
  character*256        :: arg(5)
  logical              :: file_exists = .false.

  type :: regrid_catdef
     integer                             :: maxcat
     integer, allocatable, dimension (:) :: tile_lat
     integer, allocatable, dimension (:) :: tile_lon
     real,    allocatable, dimension (:) :: vec_lai, count_lai
  end type regrid_catdef
  
  if(iargc() ==4)  then

     do n = 1 ,4
        call getarg(n,arg(n))
     end do

     read (arg(1),'(a)') YYYYDOY
     read (arg(2),* ) refdate(1)
     read (arg(3),* )   tstep(1)
     read (arg(4),* )    time(1)
     
     print *, refdate, tstep, time
     
     MODIS_NAME = 'MCD15A2H.006'
     
     call grid2tile (MODIS_NAME, trim(YYYYDOY), tstep, time, refdate)
     !  call grid2tile_with_gaps (MODIS_NAME, YYYYDOY, tstep, time, refdate)
     
  elseif (iargc() ==5)  then
     
     do n = 1 ,5
        call getarg(n,arg(n))
     end do
     
     read (arg(1),'(a)') GRID_NAME
     read (arg(2),'(a)') file_out
     read (arg(3),'(a)') file1
     read (arg(4),'(a)') file2
     read (arg(5),'(a)') file3
     
     inquire(file=trim(SMOOTH)//trim(GRID_NAME)//'/'//trim(file_out),exist=file_exists)
     if(file_exists) stop
     call smooth_data (trim(GRID_NAME), trim(file_out), trim(file1), trim(file2), trim(file3))
     
  endif

contains

  SUBROUTINE smooth_data (GRID_NAME, file_out, file1, file2, file3)

    implicit none

    character(*), intent (in) :: GRID_NAME, file_out, file1, file2, file3
    character*100             :: filename
    logical                   :: file_exists = .false.
    character*3               :: DOY
    character*4               :: YYYY
    integer                   :: n, yc, nt
    real, allocatable         :: vec_lai(:), vec_data(:)

    open (20, file = trim(SMOOTH)//trim(GRID_NAME)//'/'//trim(file_out), &
            form = 'unformatted', action = 'write', status = 'unknown')

    open (10,file = trim(BCSDIR)//trim(GRID_NAME)//'/catchment.def',status='old',action='read',form='formatted')
    read (10,*) NT
    close(10, status = 'keep')

    allocate (vec_lai (nt))
    allocate (vec_data(nt))

    vec_lai = 0.
    vec_data= 0.

    do n = 1,3
       if (n == 1) filename = trim(file1)
       if (n == 2) filename = trim(file2)
       if (n == 3) filename = trim(file3)
       open (10, file=trim(OUTDIR)//trim(GRID_NAME)//'/'//trim(filename),&
            form = 'unformatted', action = 'read', status = 'old')
       read (10) vec_data
       close(10, status = 'keep')
       vec_lai = vec_lai + vec_data / 3.
    end do

    write (20) vec_lai
    close(20, status = 'keep')
    
    ! derive climatological mean file
    
    DOY = file_out(14:16)
    inquire(file=trim(SMOOTH)//trim(GRID_NAME)//'/lai_data.YYYY'//DOY,exist=file_exists)

    if(.not.file_exists) then
       vec_lai = 0.
       yc      = 0.
       do n = 2002, 2020
          write (yyyy, '(i4.4)') n
          file_exists = .false.
          inquire(file=trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY,exist=file_exists)
          if(file_exists) then
             open (10, file=trim(OUTDIR)//trim(GRID_NAME)//'/lai_data.'//YYYY//DOY, &
                  form = 'unformatted', action = 'read', status = 'old')
             read (10) vec_data
             close(10, status = 'keep')
             vec_lai = vec_lai + vec_data
             yc = yc + 1
          endif          
       end do
       vec_lai = vec_lai + vec_data / yc
       open (20, file = trim(SMOOTH)//trim(GRID_NAME)//'/lai_data.YYYY'//DOY,  &
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
    integer, allocatable, dimension (:,:,:)    :: tile_id
    integer, allocatable, dimension (:)        :: maxcat,ii,jj
    integer*2, allocatable, dimension (:,:)    :: Lai_500m, FparLai_QC
    type (regrid_catdef), dimension(:),pointer :: catdef
    real, allocatable, target, dimension (:,:) :: lai_grid
    real, allocatable, dimension (:,:)         :: lai_clim
    real,    pointer,     dimension (:,:)      :: subset
    real,    allocatable,    dimension(:)      :: x, y
    integer :: n,i,j,k, NT, tindex, pfaf, ncid, nx, ny,status,tileid_tile
    real    :: dxm, dym,minlon,maxlon,minlat,maxlat
    INTEGER :: imn,imx,jmn,jmx,mval,d1,d2,l, ix, jx
    real, parameter :: dxy = 0.5
    character(8)    :: qc_str
    integer         :: doy, doy_tstep,i_yyyydoy

    allocate (catdef(NGRIDS))
    read (YYYYDOY, '(i7.7)') i_yyyydoy
    doy = i_yyyydoy - (i_yyyydoy/1000)*1000
    doy_tstep = (doy - 1)/8 + 1
    print *,'DOY DT', doy_tstep
    allocate (tile_id (1:IM, 1:JM, 1:ngrids))
    
    ! read BCs til/rst files
    ! ----------------------

    do n = 1, ngrids

       open (20 +n,file=trim(BCSDIR)//trim(GRIDS(n))//'/rstfile',status='old',action='read',form='unformatted',convert='little_endian')    
       open (10,file=trim(BCSDIR)//trim(GRIDS(n))//'/catchment.def',status='old',action='read',form='formatted')
       read (10,*) NT

       catdef(n)%maxcat = NT
       allocate (catdef(n)%tile_lat (1:NT))
       allocate (catdef(n)%tile_lon (1:NT))
       allocate (catdef(n)%vec_lai  (1:NT))
       allocate (catdef(n)%count_lai(1:NT))
       catdef(n)%vec_lai   = -9999.
       catdef(n)%count_lai = 0.

       do k = 1, NT
          read (10,*) tindex,pfaf,minlon,maxlon,minlat,maxlat
          catdef(n)%tile_lon(k) = (minlon + maxlon)/2.
          catdef(n)%tile_lat(k) = (minlat + maxlat)/2.
       end do
       close(10, status = 'keep')
       if(trim(GRIDS(n)) == 'DE720') then
          allocate (ii (1:NT))
          allocate (jj (1:NT))
          call ReadTileFile_RealLatLon(trim(BCSDIR)//trim(GRIDS(n))//'/tilfile', NT, ii,jj)
       endif

    end do

    do j=1,JM
       do n = 1, ngrids
          read(20+n)tile_id(:,j,n)
       end do
    end do

    do n = 1, ngrids
       close(20+n, status = 'keep')
    end do

    dxm = real(nc) /real(im) 
    dym = real(nr) /real(jm)

    allocate (Lai_500m  (1:NC, 1: NR))
    allocate (FparLai_QC(1:NC, 1: NR))
    status = NF90_OPEN(trim(RAWDIR)//trim(MODIS_NAME)//'/'//trim(MODIS_NAME)//'_LAI_'//trim(YYYYDOY)//'.nc4',NF90_NOWRITE, ncid); VERIFY_(STATUS)
    STATUS = NF90_GET_VAR (NCID,VarID(NCID,'Lai_500m'  ),   Lai_500m) ; VERIFY_(STATUS)
    STATUS = NF90_GET_VAR (NCID,VarID(NCID,'FparLai_QC'), FparLai_QC) ; VERIFY_(STATUS)
    STATUS = NF90_CLOSE (NCID)
    print *, maxval(Lai_500m), minval (LAI_500m)
    nx = nint (360./dxy)
    ny = nint (180./dxy)
    allocate (lai_grid (1 : nx, 1 : ny)) 
    allocate (x(1:nx))
    allocate (y(1:ny))
    allocate (lai_clim (1 : nx, 1 : ny))

    status = NF90_OPEN('/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/'// &
         'V001/MCD15A2H.006/IAV/ExtData/MCD15A2H.006_LAIclim_wogapfilled_ExtData.nc4',NF90_NOWRITE, ncid); VERIFY_(STATUS)
    STATUS = NF90_GET_VAR (NCID,VarID(NCID,'MODIS_LAI'  ),lai_clim, start = (/1,1, doy_tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
    STATUS = NF90_CLOSE (NCID)
    print *,maxval(lai_clim), minval (lai_clim)
    FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
    FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy
    lai_grid  = -9999
 
    do j = 1, nr
       do i = 1, nc
          if((Lai_500m(i,j) >=0).and.(Lai_500m(i,j) <= 100)) then
             ! Table 5 of https://lpdaac.usgs.gov/documents/2/mod15_user_guide.pdf
             ! (SCF_QC '000' OR '001') AND ( MODLAND_QC = 0)
             if((FparLai_QC(i,j) <= 62).and.(MOD(FparLai_QC(i,j),2) ==0)) then
                ! Cloud free 
!                qc_str = bit2str(INT(FparLai_QC(i,j)))
!                if(qc_str(4:5) /= '01') then 
                tileid_tile = tile_id (ceiling(i/dxm), ceiling (j/dym),1)
                if((tileid_tile >= 1).and.(tileid_tile <= catdef(1)%maxcat)) then
                   do n = 1, ngrids
                      k = tile_id (ceiling(i/dxm), ceiling (j/dym),n)
                      if(catdef(n)%vec_lai(k) == -9999.) catdef(n)%vec_lai(k) = 0.
                      catdef(n)%vec_lai(k)   = catdef(n)%vec_lai(k) + Lai_500m(i,j)*0.1
                      catdef(n)%count_lai(k) = catdef(n)%count_lai(k) + 1. 
                   end do
                endif
!               endif
             endif
          endif
       end do
    end do

    do n = 1, ngrids
        where (catdef(n)%count_lai > 0.) catdef(n)%vec_lai = catdef(n)%vec_lai/catdef(n)%count_lai
        if(trim(GRIDS(n)) == 'DE720') then
           do k = 1, catdef(n)%maxcat
              if(lai_clim(ii(k), jj(k)) < 100.)then
                 lai_grid(ii(k), jj(k)) = amax1(lai_clim(ii(k), jj(k)),0.01)
 !                print *,ii(k), jj(k), lai_grid(ii(k), jj(k))
              endif
              if(catdef(n)%count_lai(k) >= 1.) lai_grid(ii(k), jj(k)) =  catdef(n)%vec_lai(k)
           enddo
        endif
    end do
    
    ! Filling gaps
    !-------------

    do k = 1, ngrids
       do n = 1, catdef(k)%maxcat
          if(catdef(k)%count_lai(n) == 0.) then
             DO i = 1,nx - 1
                if ((catdef(k)%tile_lon(n) >= x(i)).and.(catdef(k)%tile_lon(n) < x(i+1))) ix = i
             end do
             DO i = 1,ny -1
                if ((catdef(k)%tile_lat(n) >= y(i)).and.(catdef(k)%tile_lat(n) < y(i+1))) jx = i
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
                   catdef(k)%vec_lai (n) = sum(subset, subset>0.)/(max(1,count(subset>0.)))
                   exit
                endif
                l = l + 1
                NULLIFY (subset)
             end do
          endif
       end do
    end do    

    ! Write output
    ! ------------
    
    lai_grid =  1.e+15 
    print *, 'Writing'
    do n = 1, ngrids
       open (10, file = trim(OUTDIR)//'/'//trim(GRIDS(n))//'/lai_data.'//trim(YYYYDOY), &
            form = 'unformatted', action = 'write', status = 'unknown')
       write (10) catdef(n)%vec_lai
       close (10,status = 'keep')
       if(trim(GRIDS(n)) == 'DE720') then
           do k = 1, catdef(n)%maxcat
              lai_grid(ii(k), jj(k)) =  catdef(n)%vec_lai(k)
           enddo
        endif
    end do
    print *, 'Done calcs'
    status = NF90_OPEN(trim(OUTDIR)//'/ExtData/'//trim(MODIS_NAME)//'_LAI_ExtData.nc4',NF90_WRITE,    ncid); VERIFY_(STATUS)
    STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'MODIS_LAI'),     lai_grid,start = (/1,1,tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
    STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'REFERENCE_DATE'),refdate, start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
    STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'time'),          time,    start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
    STATUS = NF90_CLOSE (NCID)
    print *, 'The End'
    deallocate (catdef)
    deallocate (tile_id,Lai_500m, FparLai_QC,lai_grid,lai_clim)
    return
  END SUBROUTINE grid2tile

  ! -----------------------------------------------------------------------

  SUBROUTINE grid2tile_with_gaps (MODIS_NAME, YYYYDOY, tstep, time, refdate)

    implicit none
    integer, parameter                         ::  ngridsw = 1
    character(*), intent (in)                  :: MODIS_NAME, YYYYDOY
    integer, intent (in),dimension(:)          :: tstep, time, refdate
    integer, allocatable, dimension (:,:,:)    :: tile_id
    integer, allocatable, dimension (:)        :: maxcat,ii,jj
    integer*2, allocatable, dimension (:,:)    :: Lai_500m, FparLai_QC
    type (regrid_catdef), dimension(NGRIDSW)   :: catdef
    real, allocatable, target, dimension (:,:) :: lai_grid
    real,    pointer,     dimension (:,:)      :: subset
    real,    allocatable,    dimension(:)      :: x, y
    integer :: n,i,j,k, NT, tindex, pfaf, ncid, nx, ny,status,tileid_tile
    real    :: dxm, dym,minlon,maxlon,minlat,maxlat
    INTEGER :: imn,imx,jmn,jmx,mval,d1,d2,l, ix, jx
    real, parameter :: dxy = 0.5
    character(8)    :: qc_str

    allocate (tile_id (1:IM, 1:JM, 1:ngridsw))

    ! read BCs til/rst files
    ! ----------------------

    do n = 1, ngridsw

       open (20 +n,file=trim(BCSDIR)//trim(GRIDS(n))//'/rstfile',status='old',action='read',form='unformatted',convert='little_endian')    
       open (10,file=trim(BCSDIR)//trim(GRIDS(n))//'/catchment.def',status='old',action='read',form='formatted')
       read (10,*) NT

       catdef(n)%maxcat = NT
       allocate (catdef(n)%tile_lat (1:NT))
       allocate (catdef(n)%tile_lon (1:NT))
       allocate (catdef(n)%vec_lai  (1:NT))
       allocate (catdef(n)%count_lai(1:NT))
       catdef(n)%vec_lai   = -9999.
       catdef(n)%count_lai = 0.

       do k = 1, NT
          read (10,*) tindex,pfaf,minlon,maxlon,minlat,maxlat
          catdef(n)%tile_lon(k) = (minlon + maxlon)/2.
          catdef(n)%tile_lat(k) = (minlat + maxlat)/2.
       end do
       close(10, status = 'keep')
       if(trim(GRIDS(n)) == 'DE720') then
          allocate (ii (1:NT))
          allocate (jj (1:NT))
          call ReadTileFile_RealLatLon(trim(BCSDIR)//trim(GRIDS(n))//'/tilfile', NT, ii,jj)
       endif

    end do

    do j=1,JM
       do n = 1, ngridsw
          read(20+n)tile_id(:,j,n)
       end do
    end do

    do n = 1, ngridsw
       close(20+n, status = 'keep')
    end do

    dxm = real(nc) /real(im) 
    dym = real(nr) /real(jm)

    allocate (Lai_500m  (1:NC, 1: NR))
    allocate (FparLai_QC(1:NC, 1: NR))
    status = NF90_OPEN(trim(RAWDIR)//trim(MODIS_NAME)//'/'//trim(MODIS_NAME)//'_LAI_'//trim(YYYYDOY)//'.nc4',NF90_NOWRITE, ncid); VERIFY_(STATUS)
    STATUS = NF90_GET_VAR (NCID,VarID(NCID,'Lai_500m'  ),   Lai_500m) ; VERIFY_(STATUS)
    STATUS = NF90_GET_VAR (NCID,VarID(NCID,'FparLai_QC'), FparLai_QC) ; VERIFY_(STATUS)
    STATUS = NF90_CLOSE (NCID)
    print *, maxval(Lai_500m), minval (LAI_500m)
    nx = nint (360./dxy)
    ny = nint (180./dxy)
    allocate (lai_grid (1 : nx, 1 : ny)) 
    allocate (x(1:nx))
    allocate (y(1:ny))
    
    FORALL (i = 1:nx) x(i) =  -180. + dxy/2. + (i-1)*dxy
    FORALL (i = 1:ny) y(i) =   -90. + dxy/2. + (i-1)*dxy
    lai_grid  = 1.e+15
 
    do j = 1, nr
       do i = 1, nc
          if((Lai_500m(i,j) >=0).and.(Lai_500m(i,j) <= 100)) then
             ! Table 5 of https://lpdaac.usgs.gov/documents/2/mod15_user_guide.pdf
             ! (SCF_QC '000' OR '001') AND ( MODLAND_QC = 0)
             if((FparLai_QC(i,j) <= 62).and.(MOD(FparLai_QC(i,j),2) ==0)) then
                ! Cloud free 
!                qc_str = bit2str(INT(FparLai_QC(i,j)))
!                if(qc_str(4:5) /= '01') then 
                tileid_tile = tile_id (ceiling(i/dxm), ceiling (j/dym),1)
                if((tileid_tile >= 1).and.(tileid_tile <= catdef(1)%maxcat)) then
                   do n = 1, ngridsw
                      k = tile_id (ceiling(i/dxm), ceiling (j/dym),n)
                      if(catdef(n)%vec_lai(k) == -9999.) catdef(n)%vec_lai(k) = 0.
                      catdef(n)%vec_lai(k)   = catdef(n)%vec_lai(k) + Lai_500m(i,j)*0.1
                      catdef(n)%count_lai(k) = catdef(n)%count_lai(k) + 1. 
                   end do
                endif
!               endif
             endif
          endif
       end do
    end do

    do n = 1, ngridsw
        where (catdef(n)%count_lai > 0.) catdef(n)%vec_lai = catdef(n)%vec_lai/catdef(n)%count_lai
        if(trim(GRIDS(n)) == 'DE720') then
           do k = 1, catdef(n)%maxcat
              if(catdef(n)%count_lai(k) >= 1.) lai_grid(ii(k), jj(k)) =  catdef(n)%vec_lai(k)
           enddo
        endif
    end do
    

    ! Write output
    ! ------------
        
    status = NF90_OPEN(trim(OUTDIR)//trim(MODIS_NAME)//'/ExtData/'//trim(MODIS_NAME)//'_LAI_ExtData.nc4',NF90_WRITE,    ncid); VERIFY_(STATUS)
    STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'MODIS_LAI'),     lai_grid,start = (/1,1,tstep/),count =(/nx,ny,1/)) ; VERIFY_(STATUS)
    STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'REFERENCE_DATE'),refdate, start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
    STATUS = NF90_PUT_VAR (NCID,VarID(NCID,'time'),          time,    start = (/tstep/)    ,count =(/1/)      ) ; VERIFY_(STATUS)
    STATUS = NF90_CLOSE (NCID)

  END SUBROUTINE grid2tile_with_gaps
  
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
    integer :: n,icnt,ityp, nt, umask, i
    real    :: xval,yval, pf, ik, jk
         
    if(present(mask)) then
       umask = mask
    else
       umask = 100
    endif
   
    open(11,file=InCNTileFile, &
         form='formatted',action='read',status='old')
    read (11,*, iostat=n) Nt
    	  
    do n = 1,7 ! skip header
       read(11,*)
    end do
    
    icnt = 0
    
    do i=1,Nt
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
    
  END PROGRAM mkLAIrst2tiles
