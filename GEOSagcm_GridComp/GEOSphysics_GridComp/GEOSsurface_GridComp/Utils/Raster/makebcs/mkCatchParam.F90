#define I_AM_MAIN
#include "MAPL_ErrLog.h"

PROGRAM mkCatchParam

! !INTERFACE:
!
! !ARGUMENTS:
!
!  Usage = "mkCatchParam -x nx -y ny -g Gridname -b DL -v LBCSV "       
!     -x: Size of longitude dimension of input raster.            DEFAULT: 8640
!     -y: Size of latitude dimension of input raster.             DEFAULT: 4320
!     -b: Position of dateline w.r.t. first grid cell boundaries. DEFAULT: DC (dateline-on-center)
!     -g: Gridname (name of the .til or .rst file without file extension)  
!     -v: LBCSV : Land bcs version (F25, GM4, ICA, NL3, NL4, NL5, v06, v07, v08, v09)
!     
!
! This program is good to generate  
! model, vegetation, soil, and MODIS albedo parameter files for the 
! catchment model implementation
!  
! Sarith Mahanama - March 23, 2012 
! Email: sarith.p.mahanama@nasa.gov
  use MAPL
  use rmTinyCatchParaMod 
  use process_hres_data
  use LogRectRasterizeMod, only: SRTM_maxcat 

  !   use module_irrig_params, ONLY : create_irrig_params

  implicit none
  
  include 'netcdf.inc'	
  
  ! The default is NC=i_raster=8640, NR=j_raster=4320 via "use rmTinyCatchParaMod", but
  ! NC and NR are typically overwritten through command-line arguments "-x nx -y ny".
  !
  ! In make_bcs, nc=nx=43200 and nr=ny=21600 whenever the mask file is GEOS5_10arcsec_mask*.nc.

  integer              :: NC = i_raster, NR = j_raster    
  character*5          :: LBCSV = 'UNDEF'
  character*128        :: Gridname = ''
  character*128        :: withbcs  = ''
  character*128        :: ARG, MaskFile
  character*256        :: CMD
  character*1          :: opt
  character*7          :: PEATSOURCE   = 'GDLHWSD'
  character*3          :: VEGZSOURCE   = 'D&S'
  character*2          :: DL ='DC'    
  integer              :: II, JJ, Type
  integer              :: I, J, command_argument_count, nxt
  real*8               :: dx, dy, lon0
  logical              :: regrid
  character(len=128), dimension (7) ::  Usage 
  character*128        ::  Grid2
  character*2          :: poles
  character*128        :: fnameRst = ''        ! a.k.a. "gfile[r]" in mod_process_hres_data.F90
  character*128        :: fnameTil = ''
  logical              :: file_exists, file_exists2, file_exists3, file_exists4
  logical              :: F25Tag = .false.
  logical              :: ease_grid=.false., redo_modis=.false.
  character*40         :: lai_name 
  integer, parameter   :: log_file = 998
  type (regrid_map)    :: maparc30, mapgeoland2,maparc60
  character*200        :: tmpstring, tmpstring1, tmpstring2
  character*200        :: fname_tmp, fname_tmp2, fname_tmp3, fname_tmp4
  integer              :: n_land
  logical              :: process_snow_albedo = .false. 
  character(len=10)    :: nc_string, nr_string
  integer              :: nc_ease, nr_ease, unit, clock_rate, clock1, clock2
  real                 :: seconds
  integer, allocatable :: iTable(:,:), tile_pfs(:), tile_j_dum(:)
  integer, pointer     :: tile_id(:,:)
  real, allocatable    :: tile_lat(:), tile_lon(:), min_lon(:), max_lon(:), min_lat(:), max_lat(:)
  real                 :: minlon, minlat, maxlon, maxlat, elev
  integer              :: tindex1, pfaf1, n, status

! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
!
! NOTE: "!$" is for conditional compilation
!
logical :: running_omp = .false.
!
!$ integer :: omp_get_thread_num, omp_get_num_threads
!
integer :: n_threads=1

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


    USAGE(1) ="Usage: mkCatchParam -x nx -y ny -g Gridname -b DL -v LBCSV                       "
    USAGE(2) ="     -x: Size of longitude dimension of input raster. DEFAULT: 8640              "
    USAGE(3) ="     -y: Size of latitude dimension of input raster.  DEFAULT: 4320              "
    USAGE(4) ="     -g: Gridname  (name of the .til or .rst file *without* file extension)      "
    USAGE(5) ="     -b: Position of the dateline in the first grid box (DC or DE). DEFAULT: DC  "
    USAGE(6) ="     -v: Land bcs version (F25, GM4, ICA, NL3, NL4, NL5, v06, v07, v08 v09 )     "
    USAGE(7) ="     -p: if no, it creates catchment_def and nc4 tile files then exits           "

! Process Arguments                            
!------------------ 

    CALL system_clock(count_rate=clock_rate)
    CALL get_command (cmd)
    inquire(file='clsm/mkCatchParam.log', exist=file_exists)
    if(file_exists) then
       open (log_file, file ='clsm/mkCatchParam.log', status='old', position='append', form='formatted',action='write')
    else
       open (log_file, file ='clsm/mkCatchParam.log', status='unknown', form='formatted',action='write')
       write (log_file,'(a)')trim(cmd)
       write (log_file,'(a)')' '
    endif
    withbcs = 'yes'
    I = command_argument_count()
    if(I < 1 .or. I > 10) then
       write (log_file,'(a)') "Wrong Number of arguments: ", i
       do j = 1,size(usage)
          print "(sp,a100)", Usage(j)
       end do
       call exit(1)
    end if

    nxt = 1
    call get_command_argument(nxt,arg)
    do while(arg(1:1)=='-')
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zh')==0) then
             nxt = nxt + 1
             call get_command_argument(nxt,arg)
          endif
       else
          arg = arg(3:)
       end if
       select case (opt)
       case ('x')
          read(arg,'(i6)') nc
       case ('y')
          read(arg,'(i6)') nr
       case ('g')
          Gridname = trim(arg)
       case ('v')
          LBCSV = trim(arg)
          if (trim(arg).eq."F25") F25Tag = .true.
          call init_bcs_config (trim(LBCSV))       ! get bcs details from version string
       case ('b')
          DL = trim(arg)
       case ('p')
          withbcs = trim(arg)
       case default
          do j = 1,size(usage)
             print "(sp,a100)", Usage(j)
          end do
          call exit(1)
       end select
       nxt = nxt + 1
       call get_command_argument(nxt,arg)
    end do

   call get_environment_variable ("MASKFILE"        ,MaskFile        )
 
   if(trim(Gridname) == '') then
      write (log_file,'(a)')'Unable to create parameters without til/rst files.... !'
      stop
   endif
   
   regrid = nc/=i_raster .or. nr/=j_raster
  
   if (index(Gridname,'EASEv') /=0) then
      ! here Gridname has alias EASELabel
      call MAPL_ease_extent(Gridname, nc_ease, nr_ease, _RC)
      write(nc_string, '(i0)') nc_ease
      write(nr_string, '(i0)') nr_ease
      Gridname = trim(Gridname)//'_'//trim(nc_string)//'x'//trim(nr_string)
      ease_grid = .true.
   endif

   if(index(Gridname,'Pfaf.notiny')/=0) then 
      fnameRst='clsm/'//trim(Gridname)
      fnameTil='clsm/'//trim(Gridname)
   else
      fnameRst='rst/'//trim(Gridname)  
      fnameTil='til/'//trim(Gridname)  
    endif 

    if(use_PEATMAP)  PEATSOURCE   = 'PEATMAP'
    if(jpl_height)   VEGZSOURCE   = 'JPL'

    if (trim(SNOWALB)=='MODC061' .or. trim(SNOWALB) =='MODC061v2')  process_snow_albedo=.true.

!    if(n_threads == 1) then

       write (log_file,'(a)')trim(LAIBCS)
       write (log_file,'(a)')trim(MODALB)
       write (log_file,'(a)')trim(SOILBCS)   
       write (log_file,'(a)')trim(SNOWALB)
       write (log_file,'(a)')trim(MaskFile)
       write (log_file,'(a)')trim(PEATSOURCE)
       write (log_file,'(a)')trim(VEGZSOURCE)
       write (log_file,'(a)')'                              '
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'............ Begin CLSM parameter generation:...............'    
       write (log_file,'(a)')'                                               '
       write (log_file,'(a)')'CLSM parameters are being generated for the tile space :'
       write (log_file,'(a)')'     ',trim(fnameTil)
       write (log_file,'(a)')'                                               '
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'                                               '
       
       if(index(Gridname,'CF')/=0) then 
          DL = 'DE'
          write (log_file,'(a)')'Cube-Sphere Grid - assuming dateline-on-edge (DE)'
       endif
       
       ! ******************************************************************************
       !
       ! IMPORTANT: The top-level make_bcs script should not allow this program to
       !            run when ./clsm/ exists.  Consequently, across "Steps [xx]" below,
       !            the "inquire()" statements should be obsolete, and the case 
       !            "Using existing file" should never happen.
       !
       ! ******************************************************************************

       allocate(tile_id(nc, nr))
       fname_tmp = trim(fnameRst)//'.rst'
       open (newunit=unit,file=fname_tmp,status='old',action='read',form='unformatted',convert='little_endian', IOSTAT=status)
       if (status /=0) then
          write (log_file,'(a)')'         '//trim(fname_tmp) // 'cannot be opened, exit '
          call exit(1)
       endif
       do j = 1, nr
         read(unit)tile_id(:,j)
       end do        
       close(unit)
       ! Creating catchment.def 
       ! ----------------------
       
       tmpstring = 'Step 01: Supplemental tile attributes and nc4-formatted tile file'
       fname_tmp = 'clsm/catchment.def'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       if(.not.ease_grid) then  
          inquire(file=trim(fname_tmp), exist=file_exists)
          if (.not.file_exists) then
             write (log_file,'(a)')'         Creating catchment def and nc4 tile file...'
             call system_clock(clock1)
             call supplemental_tile_attributes(nc,nr,regrid,dl,fnameTil, tile_id) 
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
       else 
          write (log_file,'(a)')'Skipping step for EASE grid. '
          write (log_file,'(a)')'catchment.def file and tile file should already be created by mkEASETilesParam.x '
       endif
       write (log_file,'(a)')' '

       if (trim(withbcs) == 'no') then
         write (log_file,'(a)')'Skipping MOM BCs. BCs will be extracted from the corresponding BCs. '
         close (log_file,status='keep')
         call exit(0)
       endif
   
       call MAPL_ReadTilingNC4( trim(fnameTil)//".nc4", iTable = iTable) 

       N_land = count(iTable(:,0) == 100)          ! n_land = number of land tiles
       allocate(tile_j_dum, source = iTable(1:n_land,7)) ! possible used in cti_stats.dat
       deallocate (iTable)
   
       ! reading from catchment to preserve zero-diff
       open (newunit=unit,file='clsm/catchment.def',status='old',action='read',form='formatted', IOSTAT=status)
       if (status /=0) then
          write (log_file,'(a)')'         clsm/cathment.def cannot be opened, exit '
          call exit(1)
       endif
       read(unit,*) N
       if (n /= n_land) then
           write (log_file,'(a)')'n_land not consistent between tile file and catchment.def file, exit '
           write (log_file,*) n_land, n
           call exit(1) 
       endif

       allocate(min_lon(n_land), max_lon(n_land), min_lat(n_land), max_lat(n_land))
       allocate(tile_lat(n_land), tile_lon(n_land))
       allocate(tile_pfs(n_land))

       do n = 1, N_land
          read (unit,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat, elev
          min_lon(n) = minlon
          max_lon(n) = maxlon
          min_lat(n) = minlat
          max_lat(n) = maxlat
          tile_lon(n)= (minlon + maxlon)/2.0
          tile_lat(n)= (minlat + maxlat)/2.0
          tile_pfs(n)= pfaf1  
       end do
       close (unit,status='keep')

       inquire(file='clsm/catch_params.nc4', exist=file_exists)
       if (.not.file_exists) CALL open_landparam_nc4_files(N_land,process_snow_albedo)  

       ! Creating cti_stats.dat 
       ! ----------------------

       tmpstring = 'Step 02: Compound Topographic Index (CTI) stats'
       fname_tmp = 'clsm/cti_stats.dat'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file...'
          call system_clock(clock1)
          call cti_stat_file (MaskFile, n_land, tile_pfs, tile_j_dum)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif
       write (log_file,'(a)')' '
       
       ! Creating vegetation classification files
       !-----------------------------------------
       
       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then

          tmpstring = 'Step 03: Vegetation types using ESA land cover (MOSAIC/Catch)'
          fname_tmp = 'clsm/mosaic_veg_typs_fracs'
          write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
          inquire(file=trim(fname_tmp), exist=file_exists)
          if (.not.file_exists) then
             write (log_file,'(a)')'         Creating file...'
             call system_clock(clock1)
             call ESA2MOSAIC (nc,nr, n_land, tile_pfs, tile_id)
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
          write (log_file,'(a)')' '

          tmpstring = 'Step 04: Vegetation types using ESA land cover (CatchCNCLM40)'
          fname_tmp = 'clsm/CLM_veg_typs_fracs'
          write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
          inquire(file=trim(fname_tmp), exist=file_exists)
          if (.not.file_exists) then
             write (log_file,'(a)')'         Creating file...'
             call system_clock(clock1)
             call ESA2CLM (nc,nr, n_land, tile_lat, tile_pfs, tile_id)
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
          write (log_file,'(a)')' '

       else

          tmpstring = 'Step 03: Vegetation types using IGBP SiB2 land cover (MOSAIC/Catch)'
          fname_tmp = 'clsm/mosaic_veg_typs_fracs'
          write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
          inquire(file=trim(fname_tmp), exist=file_exists)
          if (.not.file_exists) then
             write (log_file,'(a)')'         Creating file...'
             call system_clock(clock1)
             call compute_mosaic_veg_types (nc, nr, regrid, n_land, tile_pfs, tile_id)
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
          write (log_file,'(a)')' '
          
          ! Per make_bcs, it looks like there are four possible mask files:
          !
          ! GEOS5_10arcsec_mask.nc
          ! global.cat_id.catch.DL
          ! global.cat_id.catch.GreatLakesCaspian_Updated.DL
          ! GEOS5_10arcsec_mask_freshwater-lakes.nc
          !
          ! If we are in this else block, we must be using one of the latter three masks.
          ! It looks like these latter masks only work for Catchment and not CatchCNCLM[xx]
          !
          ! - reichle, 11 Jan 2022
          
          write (log_file,'(a)')'NOTE: The selected mask works only for the Catchment model.'
          write (log_file,'(a)')'      Vegetation types *not* created for CatchCNCLM[xx].'
          write (log_file,'(a)')'      SKIPPING Step 04 and Step 05 !!!'
          write (log_file,'(a)')' '
          
       endif
       
       ! Processing Vegetation Climatology 
       ! ---------------------------------
       
       ! creating mapping arrays if necessary
                  
       tmpstring = 'Step 05: Vegetation climatologies'
       write (log_file,'(a,a,a)') trim(tmpstring),' ', trim(LAIBCS)
       
       if((trim(LAIBCS) == 'MODGEO').or.(trim(LAIBCS) == 'GEOLAND2')) then 
          fname_tmp = 'clsm/lai.GEOLAND2_10-DayClim'
          write (log_file,'(a,a)')'         --> ', trim(fname_tmp)
          inquire(file=trim(fname_tmp), exist=file_exists)          
          if (.not.file_exists) then
             write (log_file,'(a)')'         Creating file...'
             !allocate (mapgeoland2 (1:40320,1:20160))
             call system_clock(clock1)
             call create_mapping (nc,nr,40320,20160,mapgeoland2, n_land, tile_id )
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done create mapping mapgeoland2. Spent   ', seconds, "  seconds"
             lai_name = 'GEOLAND2_10-DayClim/geoland2_'
 
             write (log_file,'(a)')'         Creating '//lai_name
             call system_clock(clock1)
             if(trim(LAIBCS) == 'GEOLAND2') then
                call hres_lai_no_gswp (40320,20160,mapgeoland2, lai_name, n_land, tile_lon, tile_lat) 
             else
                call hres_lai_no_gswp (40320,20160,mapgeoland2, lai_name, n_land, tile_lon, tile_lat, merge=1) 
             endif
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
             ! if(allocated(mapgeoland2)) deallocate (mapgeoland2)
             deallocate (mapgeoland2%map)
             deallocate (mapgeoland2%ij_index)
             write (log_file,'(a)')'         Done.'
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
       endif
       
       if ((LAIBCS == 'MODGEO').or.(LAIBCS == 'MODIS').or.(MODALB == 'MODIS2')) then
          ! allocate (maparc30    (1:43200,1:21600))
          call system_clock(clock1)
          call create_mapping (nc,nr,43200,21600,maparc30, n_land,  tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done create mapping maparc30. Spent   ', seconds, "  seconds"
       endif
       
       fname_tmp = 'clsm/green.dat'
       write (log_file,'(a,a)')'         --> ', trim(fname_tmp)
       inquire(file=trim(fname_tmp), exist=file_exists)          
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file... (resolution will be added to file name later)'
          
          call system_clock(clock1)
          if (trim(LAIBCS) == 'GSWP2') then 
             call process_gswp2_veg (nc,nr,regrid,'grnFrac',n_land, tile_id)
          else
             if (size(maparc30%ij_index,1) /= 43200) then 
                ! allocate (maparc30    (1:43200,1:21600))
                call create_mapping (nc,nr,43200,21600,maparc30, n_land,  tile_id)
             endif
             call hres_gswp2 (43200,21600, maparc30, 'green', n_land, tile_lon, tile_lat) 
          endif
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent  ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif

       fname_tmp = 'clsm/lai.dat'
       write (log_file,'(a,a)')'         --> ', trim(fname_tmp)
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file... (resolution will be added to file name later)'
          redo_modis = .true.
          
          call system_clock(clock1)
          if (trim(LAIBCS) == 'GSWP2') call process_gswp2_veg (nc,nr,regrid,'LAI', n_land, tile_id) 
          if (trim(LAIBCS) == 'GSWPH') then
             if (size(maparc30%ij_index,1) /= 43200) then 
                ! allocate (maparc30    (1:43200,1:21600))
                call create_mapping (nc,nr,43200,21600,maparc30, n_land, tile_id)
             endif
             inquire(file='clsm/lai.MODIS_8-DayClim', exist=file_exists)
             if (.not.file_exists) call hres_gswp2 (43200,21600, maparc30, 'lai', n_land, tile_lon, tile_lat) 
          endif
          
          if (trim(LAIBCS) == 'MODIS') then
             lai_name = 'MODIS_8-DayClim/MODIS_'
             call hres_lai_no_gswp (43200,21600,maparc30,lai_name, n_land, tile_lon, tile_lon) 
          endif
          
          if (trim(LAIBCS) == 'MODGEO') then
             lai_name = 'MODIS_8-DayClim/MODIS_'
             inquire(file='clsm/lai.MODIS_8-DayClim', exist=file_exists)
             if (.not.file_exists)call hres_lai_no_gswp (43200,21600,maparc30,lai_name, n_land, tile_lon, tile_lat, merge=1)  
             call merge_lai_data (MaskFile, n_land, tile_pfs)
          endif
          
          if (trim(LAIBCS) == 'MODISV6') then
             lai_name = 'MCD15A2H.006/MODIS_'
             call grid2tile_modis6 (86400,43200,nc,nr,n_land, tile_lon, tile_lat, tile_id, lai_name)  
          endif

          if (trim(LAIBCS) == 'GLASSA') then
             lai_name = 'GLASS-LAI/AVHRR.v4/GLASS01B02.V04.AYYYY'
             call grid2tile_glass (nc,nr, tile_id,lai_name, n_land, tile_lon, tile_lat)  
          endif

          if (trim(LAIBCS) == 'GLASSM') then
             lai_name = 'GLASS-LAI/MODIS.v4/GLASS01B01.V04.AYYYY'
             call grid2tile_glass (nc,nr,tile_id,lai_name, n_land, tile_lon, tile_lat)  
          endif
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif

       fname_tmp = 'clsm/ndvi.dat'
       write (log_file,'(a,a)')'         --> ', trim(fname_tmp)
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file... (resolution will be added to file name later)'
          call system_clock(clock1)
          call gimms_clim_ndvi (nc,nr, n_land, tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif

       write (log_file,'(a)')' '
       
       ! -------------------------------------------------
       
       ! call modis_alb_on_tiles (nc,nr,ease_grid,regrid,fnameTil,fnameRst)
       ! call modis_scale_para (ease_grid,fnameTil)
       ! NOTE: modis_alb_on_tiles uses monthly climatological raster data on 8640x4320 to produce 
       ! MODIS albedo on tile space. The subroutine was replaced with "modis_alb_on_tiles_high" that process
       ! MODIS1 data on native grid and produces 8/16-day MODIS Albedo climatology


       tmpstring = 'Step 06: Albedo climatologies'
       write (log_file,'(a,a,a)') trim(tmpstring),' ', trim(MODALB)
       
       if(MODALB == 'MODIS1') then 
          fname_tmp = 'clsm/AlbMap.WS.16-day.tile.0.7_5.0.dat'
          write (log_file,'(a,a)')'         --> ', trim(fname_tmp)
          inquire(file=trim(fname_tmp), exist=file_exists)          
          if (.not.file_exists) then
             write (log_file,'(a)')'         Creating file...'
             call system_clock(clock1)
             if(F25Tag) then 
                call create_mapping (nc,nr,21600,10800,maparc60, n_land, tile_id)
                call modis_alb_on_tiles_high (21600,10800,maparc60,MODALB, n_land)
                deallocate (maparc60%map)
                deallocate (maparc60%ij_index)
             else
                !  This option is for legacy sets like Fortuna 2.1
                call modis_alb_on_tiles (nc,nr,regrid, n_land, tile_id)
             endif
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
       endif
       
       if(MODALB == 'MODIS2') then 
          fname_tmp  = 'clsm/AlbMap.WS.8-day.tile.0.3_0.7.dat'
          fname_tmp2 = 'clsm/AlbMap.WS.8-day.tile.0.7_5.0.dat'
          write (log_file,'(a,a,a,a)')'         --> ', trim(fname_tmp), ', ', trim(fname_tmp2)
          inquire(file=trim(fname_tmp ), exist=file_exists )          
          inquire(file=trim(fname_tmp2), exist=file_exists2)          
          if ((.not.file_exists).or.(.not.file_exists2)) then
             call system_clock(clock1)
             write (log_file,'(a)')'         Creating files...'
             call modis_alb_on_tiles_high (43200,21600,maparc30,MODALB, n_land)
             call system_clock(clock2)
             seconds = (clock2-clock1)/real(clock_rate)
             write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          else
             write (log_file,'(a)')'         Using existing file.'
          endif
       endif
       write (log_file,'(a)')' '
       
       ! ---------------------------------------------

       tmpstring = 'Step 07: Albedo scale factors'
       write (log_file,'(a,a,a)') trim(tmpstring),' ', trim(MODALB)
       
       ! NOTE: There are two files with albedo scale factors: "visdf.dat" and "nirdf.dat".
       !       Added check for "nirdf.dat", which was missing before.  - reichle, 13 Jan 2022

       fname_tmp  = 'clsm/visdf.dat'
       fname_tmp2 = 'clsm/nirdf.dat'
       write (log_file,'(a,a,a,a)')'         --> ', trim(fname_tmp), ', ', trim(fname_tmp2)
       inquire(file=trim(fname_tmp ), exist=file_exists )                        
       inquire(file=trim(fname_tmp2), exist=file_exists2)
       if ((redo_modis).or.(.not.file_exists).or.(.not.file_exists2)) then
          !   if(.not.F25Tag) then
          write (log_file,'(a)')'         Creating files... (resolution will be added to file name later)'
          call system_clock(clock1)
          call modis_scale_para_high (MODALB, n_land)
          !  else
          !     This option is for legacy sets like Fortuna 2.1
          !     inquire(file='clsm/modis_scale_factor.albvf.clim', exist=file_exists)
          !     if ((redo_modis).or.(.not.file_exists)) then
          !        call modis_scale_para (ease_grid,fnameTil)
          !        call REFORMAT_VEGFILES
          !     endif
          !  endif
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing files.'
       endif
       write (log_file,'(a)')' '
       
       !       tmpstring1 = '-e EASE -g '//trim(gfile) 
       !       write(tmpstring2,'(2(a2,x,i5,x))')'-x',nc,'-y',nr
       !       tmpstring = 'bin/mkCatchParam_openmp '//trim(tmpstring2)//' '//trim(tmpstring1)
       
!    else      
       
       ! this block is for n_threads>1
       !==============================
       
       if(trim(SOILBCS)=='NGDC') then
          write (log_file,'(a)')'Creating (intermediate) NGDC soil types file...'
          call system_clock(clock1)
          call create_soil_types_files (nc,nr, n_land, tile_pfs, tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
          write (log_file,'(a)')' '
       endif
       
       ! Creating soil_param.first and tau_param.dat files that has 2 options: 
       !  1) NGDC soil properties, 2) HWSD-STATSGO2 Soil Properties
       ! ---------------------------------------------------------------------
       
       tmpstring = 'Step 08: Soil parameters ' // trim(SOILBCS) 
       fname_tmp = 'clsm/soil_param.first'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file...'
          call system_clock(clock1)
          if(trim(SOILBCS)=='NGDC')  then 
             if(     F25Tag) call soil_para_high (nc,nr,regrid, n_land, tile_id,F25Tag=F25Tag)
             if(.not.F25Tag) call soil_para_high (nc,nr,regrid, n_land, tile_id)
          endif
          if(SOILBCS(1:4)=='HWSD') call soil_para_hwsd (nc,nr, n_land, tile_pfs, tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a,a)')'         Using existing file.'
       endif
       write (log_file,'(a)')' '
              
       tmpstring  = 'Step 09: CLSM model parameters ' // trim(SOILBCS) 
       fname_tmp  = 'clsm/ar.new'
       fname_tmp2 = 'clsm/bf.dat'
       fname_tmp3 = 'clsm/ts.dat'
       fname_tmp4 = 'clsm/soil_param.dat'
       tmpstring1 = trim(fname_tmp) // ', ' // trim(fname_tmp2) // ', ' // trim(fname_tmp3) // ', ' // trim(fname_tmp4) 
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(tmpstring1), ')'
       inquire(file=trim(fname_tmp ), exist=file_exists )
       inquire(file=trim(fname_tmp2), exist=file_exists2)
       inquire(file=trim(fname_tmp3), exist=file_exists3)
       inquire(file=trim(fname_tmp4), exist=file_exists4)
       if ((.not.file_exists).or.(.not.file_exists2).or.(.not.file_exists3).or.(.not.file_exists4)) then
          write (log_file,'(a)')'         Creating files...'
          call system_clock(clock1)
          if(trim(SOILBCS)=='NGDC') call create_model_para(        MaskFile, n_land, tile_lon, tile_lat, tile_pfs)
          if(SOILBCS(1:4) =='HWSD') call create_model_para_woesten(MaskFile, n_land, tile_lon, tile_lat, tile_pfs)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) '         Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a,a)')'         Using existing files.'
       endif
       write (log_file,'(a)')' '
       
       ! Commented out this call because 7.5-minute raster file is only used
       ! for plotting purposes
       !  call make_75 (nc,nr,regrid,c_data,fnameRst)
       !    write (log_file,'(a)')'Done creating 7.5 minute raster file ......................'
       write (log_file,'(a)')'NOTE: 7.5 minute raster file not created (only needed for diagnostic plotting).'
       write (log_file,'(a)')'      Uncomment associated lines in source to generate 7.5 minute raster file.'
       write (log_file,'(a)')' '
       
       tmpstring = 'Step 10: CatchCNCLM40 NDep T2m SoilAlb parameters'
       fname_tmp = 'clsm/CLM_NDep_SoilAlb_T2m'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       ! create this file only if matching veg types file already exists
       inquire(file='clsm/CLM_veg_typs_fracs', exist=file_exists)
       if (file_exists) then
          write (log_file,'(a)')'         Creating file...'
          call system_clock(clock1)
          call grid2tile_ndep_t2m_alb (nc,nr, n_land,tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) 'Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'Skipping step for lack of matching veg types file.'
       endif
       write (log_file,'(a)')' '       
       
       tmpstring = 'Step 11: CatchCNCLM45 abm peatf gdp hdm fc parameters'
       fname_tmp = 'clsm/CLM4.5_abm_peatf_gdp_hdm_fc'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file...'
          call system_clock(clock1)
          call CLM45_fixed_parameters (nc,nr, n_land, tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) 'Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif
       write (log_file,'(a)')' '
       
       tmpstring = 'Step 12: CatchCNCLM45 lightning frequency'
       fname_tmp = 'clsm/lnfm.dat'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file... (resolution will be added to file name later)'
          call system_clock(clock1)
          call CLM45_clim_parameters (nc,nr,n_land,tile_id)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) 'Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif
       write (log_file,'(a)')' '

       tmpstring = 'Step 13: Country and state codes'
       fname_tmp = 'clsm/country_and_state_code.data'
       write (log_file,'(a,a,a,a)') trim(tmpstring), ' (', trim(fname_tmp), ')'
       inquire(file=trim(fname_tmp), exist=file_exists)
       if (.not.file_exists) then
          write (log_file,'(a)')'         Creating file...'
          call system_clock(clock1)
          call map_country_codes (nc,nr,n_land, tile_lon, tile_lat)
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) 'Done. Spent   ', seconds, "  seconds"
       else
          write (log_file,'(a)')'         Using existing file.'
       endif
       write (log_file,'(a)')' '

       if(process_snow_albedo)then
          tmpstring = 'Step 14: Static snow albedo from MODIS' 
          write (log_file,'(a)') trim(tmpstring)
          write (log_file,'(a)')'         Creating file...'
          call system_clock(clock1)
          if (trim(SNOWALB)=='MODC061') then 
             call MODIS_snow_alb (n_land, min_lon, max_lon, min_lat, max_lat)
          elseif (trim(SNOWALB)=='MODC061v2') then
             if (size(maparc30%ij_index,1) /= 43200) then 
                call create_mapping (nc,nr,43200,21600,maparc30, n_land, tile_id)
             end if
             call MODIS_snow_alb_v2(43200,21600,maparc30, n_land) 
          else
             write (log_file,'(a)')'Unknown SNOWALB... stopping!'
             stop
          endif
          call system_clock(clock2)
          seconds = (clock2-clock1)/real(clock_rate)
          write (log_file, *) 'Done. Spent   ', seconds, "  seconds"
          write (log_file,'(a)')' '
       endif

       !      inquire(file='clsm/irrig.dat', exist=file_exists)
       !      if (.not.file_exists) call create_irrig_params (nc,nr,fnameRst)
       !      write (log_file,'(a)')'Done computing irrigation model parameters ...............13'
       
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'DONE creating CLSM data files...............................'
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')' '
       
       !       call execute_command_line ('chmod 755 bin/create_README.csh ; bin/create_README.csh')
!    endif

    close (log_file,status='keep') 

contains
    subroutine nearest_index_vector(candidates, targets, idx)
    ! For each targets(k), find argmin_j |candidates(j) - targets(k)|
    real, intent(in)  :: targets(:)
    real, intent(in)  :: candidates(:)
    integer, intent(out) :: idx(size(candidates))   ! 1-based indices
    integer :: k, j, nT, nC, best_j
    real :: best_d, d

    nT = size(targets)
    nC = size(candidates)

    do k = 1, nC
       best_d = huge(1.0d0)
       best_j = 1
       do j = 1, nT
          d = abs(candidates(k) - targets(j))
          if (d < best_d) then
             best_d = d
             best_j = j
          end if
       end do
       idx(k) = best_j    ! already 1-based
    end do

  end subroutine nearest_index_vector

  subroutine create_catchmapping()
    integer,parameter :: nlon=21600         ! Number of longitude grid points in the original grid
    integer,parameter :: nlat=10800         ! Number of latitude grid points in the original grid
    
    
    ! Variable declarations:
    integer              :: id, xi, yi, i,j, flag, subi,nmax,nmax2,ntot,it,ns_tot
    integer,allocatable  :: nsub(:)                   ! Array to store the number of sub-areas for each catchment
    integer, allocatable :: xsub(:,:), ysub(:,:), isub(:,:)
    ! xsub and ysub: 2D arrays to store mapped x and y indices for sub-catchments (not using subi_global in this code)
    real, allocatable    :: asub(:,:)        ! 2D array to store aggregated area for each sub-catchment
    
    real, allocatable  :: lon(:), lat(:)  ! Arrays to hold longitude and latitude values from the NetCDF file
    real, allocatable  :: lon36(:), lat36(:), tmp  ! Arrays to hold longitude and latitude values from the NetCDF file
    integer, allocatable :: loni(:), lati(:)
    integer, allocatable :: loni2(:), lati2(:)
    ! loni and lati: Arrays holding mapping indices from 1-minute resolution data files
    integer, allocatable :: catchind(:,:) ! 2D array holding catchment indices for each grid cell
    real, allocatable    :: cellarea(:,:)    ! 2D array containing the area of each grid cell
    
    ! Declare allocatable arrays for catchment ID, parent catchment ID, and boundary coordinates
    integer, allocatable, dimension(:) :: tid, catid, lati_tile, loni_tile
    real, allocatable, dimension(:)    :: lon_left, lon_right, lat_bottom, lat_top, latc,lonc
    integer,allocatable,dimension(:,:) :: map_tile
    
    real,allocatable,dimension(:) :: area_cat, pfaf_frac
    real,allocatable,dimension(:,:) :: frac,frac_tile
    integer,allocatable,dimension(:) ::  nsub_tile
    integer,allocatable ::csub(:,:),flag2(:,:), tile_id(:), pfaf_index(:)
    type(NetCDF4_FileFormatter) :: formatter
    type (FileMetadata) :: meta

    ! Define file path for input routing data:
    character(len=128):: pfaf_file ="/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/SRTM_PfafData.nc"
    
    ! Allocate arrays with the specified dimensions:
    allocate(catchind(nlon, nlat), cellarea(nlon, nlat))
    allocate(lon(nlon), lat(nlat))
    allocate(loni(nlon), lati(nlat),loni2(nlon), lati2(nlat))
    allocate(lon36(nc_ease),lat36(nr_ease))
    allocate(nsub(SRTM_maxcat))
    
    ! Read grid information from the NetCDF file "CatchIndex.nc":
    call formatter%open(trim(pfaf_file), PFIO_READ)
    call formatter%get_var("latitude", lat)
    call formatter%get_var("longitude", lon)
    call formatter%get_var("CatchIndex", catchind)
    call formatter%close()

    do i = 1, nc_ease
       call MAPL_ease_inverse( trim(GridName), real(i-1), 0.0, tmp, lon36(i))
    enddo
    do j = 1, nr_ease
       call MAPL_ease_inverse( trim(GridName), 0.0, real(j-1), lat36(j), tmp)
    enddo

    call nearest_index_vector(lat, lat36, lati)
    call nearest_index_vector(lon, lon36, loni)
    
    call formatter%open("/discover/nobackup/yzeng3/data/river_preproc_input/cellarea.nc", PFIO_READ)
    call formatter%get_var("data", cellarea)
    call formatter%close()
    
    cellarea = cellarea / 1.e6  ! Convert cell area (e.g., from m^2 to km^2)
    ! Initialize aggregation arrays to zero:
    allocate(xsub(9999, SRTM_maxcat), ysub(9999, SRTM_maxcat))
    nsub = 0
    ! Loop over all grid cells to aggregate cell areas by catchment and sub-area:
    do xi = 1, nlon
      do yi = 1, nlat
        if (catchind(xi, yi) >= 1) then
          ! The grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0              ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub(i, id) .and. lati(yi) == ysub(i, id)) then
                flag = 1
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(id) = nsub(id) + 1
            xsub(nsub(id), id) = loni(xi)
            ysub(nsub(id), id) = lati(yi)
          endif
        endif
      end do
    end do
    nmax = maxval(nsub)
    !print *,nmax
    deallocate(xsub,ysub)
    allocate(xsub(nmax, SRTM_maxcat), ysub(nmax, SRTM_maxcat), asub(nmax, SRTM_maxcat))
    ! Initialize aggregation arrays to zero:
    nsub = 0;xsub = 0;ysub = 0;asub = 0.
    ! Loop over all grid cells to aggregate cell areas by catchment and sub-area:
    do xi = 1, nlon
      do yi = 1, nlat
        if (catchind(xi, yi) >= 1) then
          ! The grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0              ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub(i, id) .and. lati(yi) == ysub(i, id)) then
                flag = 1
                ! If a match is found, accumulate the cell area into the existing sub-area:
                asub(i, id) = asub(i, id) + cellarea(xi, yi)
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(id) = nsub(id) + 1
            xsub(nsub(id), id) = loni(xi)
            ysub(nsub(id), id) = lati(yi)
            asub(nsub(id), id) = cellarea(xi, yi)
          endif
        endif
      end do
    end do
    
    ! Open the catchment definition file for the M36 grid and read the total number of catchments (header)
    open(77, file="clsm/catchment.def");read(77, *) ntot
    ! Allocate arrays with size nt
    allocate(tid(ntot), catid(ntot), lon_left(ntot), lon_right(ntot), lat_bottom(ntot), lat_top(ntot),latc(ntot),lonc(ntot))
    allocate(lati_tile(ntot),loni_tile(ntot),map_tile(nc_ease,nr_ease))
    ! Loop over each catchment and read: id, catchment id, left/right longitudes, bottom/top latitudes
    do i = 1, ntot
      read(77, *) tid(i), catid(i), lon_left(i), lon_right(i), lat_bottom(i), lat_top(i)
    end do
    latc = (lat_bottom + lat_top) / 2.
    lonc = (lon_left + lon_right) / 2.
    call nearest_index_vector(latc,lat36,lati_tile)
    call nearest_index_vector(lonc,lon36,loni_tile)
    map_tile=-9999
    do i =1, ntot
      map_tile(loni_tile(i),lati_tile(i)) = i
    enddo
    
    allocate(isub(nmax, SRTM_maxcat))
    isub=0
    ! Loop over each catchment
    do i = 1, SRTM_maxcat
      ! Loop over each potential sub-area within the current catchment
      do j = 1, nmax
        xi = xsub(j, i)  ! Retrieve the x-coordinate for the sub-area
        yi = ysub(j, i)  ! Retrieve the y-coordinate for the sub-area
        if (xi /= 0) then  ! Check if a valid sub-area exists (non-zero x-coordinate)
          isub(j, i) = map_tile(xi, yi)  ! Map the sub-area coordinates to a global tile index using map_tile
        endif
      enddo
    enddo
    
    allocate(area_cat(SRTM_maxcat),frac(nmax,SRTM_maxcat),nsub_tile(ntot),flag2(nmax,SRTM_maxcat))
    where(isub==-9999) asub=0.
    area_cat=0.
    do i=1, SRTM_maxcat
      area_cat(i)=sum(asub(:,i))
    enddo
    nsub_tile=0
    frac=0.
    do i=1, SRTM_maxcat
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1
          frac(j,i)=asub(j,i)/area_cat(i)
        endif
        if(isub(j,i)==0)exit
      enddo
    enddo
    nmax2=maxval(nsub_tile)
    allocate(csub(nmax2,ntot),frac_tile(nmax2,ntot))
    csub=0
    nsub_tile=0
    frac_tile=0.
    do i=1, SRTM_maxcat
      do j=1,nmax
        if(isub(j,i)>0)then
          it=isub(j,i)
          nsub_tile(it)=nsub_tile(it)+1
          csub(nsub_tile(it),it)=i
          frac_tile(nsub_tile(it),it)=frac(j,i)
        endif
        if(isub(j,i)==0)exit
      enddo
    enddo
    ns_tot=sum(nsub_tile)
  
    allocate(tile_id(ns_tot), pfaf_index(ns_tot), pfaf_frac(ns_tot))
    id = 0
    do i=1,ntot
      do j=1,nmax2
        if(csub(j,i)>0)then
          id=id+1
          tile_id(id)    = i
          pfaf_index(id) = csub(j,i)
          pfaf_frac(id)  = frac(j,i)
        endif
        if(csub(j,i)==0)exit
      enddo
    enddo

    call meta%add_dimension('tile', ns_tot)
    call meta%add_variable('tile_id', Variable(type=pFio_Int32, dimensions='tile'))
    call meta%add_variable('pfaf_index', Variable(type=pFio_Int32, dimensions='tile'))
    call meta%add_variable('pfaf_frac', Variable(type=pFio_REAL32, dimensions='tile'))
    call formatter%create('clsm/EASEv2_M36_tile_pfaf.nc4', mode=PFIO_CLOBBER)
    call formatter%write(meta)
    call formatter%put_var('tile_id', tile_id)
    call formatter%put_var('pfaf_index', pfaf_index)
    call formatter%put_var('pfaf_frac',  pfaf_frac)
    call formatter%close()
         
  end subroutine


END PROGRAM mkCatchParam
