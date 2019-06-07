#include "Raster.h"

 PROGRAM mkCatchParam

! !INTERFACE:
!
! !ARGUMENTS:
!
!  Usage = "mkCatchParam -x nx -y ny -g Gridname -b DL -m MA -l LD -s SD -e EASE"       
!     -x: Size of longitude dimension of input raster. DEFAULT: 8640
!     -y: Size of latitude dimension of input raster.  DEFAULT: 4320
!     -b: position of the dateline in the first box. DEFAULT: DC 
!     -g: Gridname  (name of the .til or .rst file without file extension)  
!     -l: Choice of LAI data set. DEFAULT : MODISV6
!         MODISV6 : 8-day climatology from the period 2002.01-2016.10 on  86400x43200 grid
!         MODGEO  : MODIS with GEOLAND2 overlaid on South America, Afirca and Australia
!         GEOLAND2: 10-day climatology from the period 1999-2011 on 40320x20160 grid               
!         GSWP2   : Monthly climatology from the period 1982-1998 on 360x180 grid                  
!         MODIS   : 8-day climatology from the period 2000-2013 on  43200x21600 grid
!         GSWPH   : Monthly climatology from the period 1982-1998 on 43200x21600 grid              "
!     -m: Choice of MODIS Albedo data. DEFAULT : MODIS2                                            
!         MODIS1 : 16-day Climatology from 1'x1 (21600x10800) MODIS data from the period 2000-2004 
!         MODIS2 : 8-day Climatology from 30"x30"(43200x21600) MODIS data from the period 2001-2011 
!     -s: Choice of soil data. DEFAULT :HWSD                                                       
!         HWSD : Merged HWSD-STATSGO2 soil properties on 43200x21600 with Woesten et al. (1999) Parameters   
!         NGDC : Reynolds soil texture clsses on 4320x2160 with GSWP2 soil hydraulic  parameters                   
!     -e: EASE : This is optional if catchment.def file is available already or                    
!         the til file format is pre-Fortuna-2.                                                    
!     
!
! This program is good to generate  
! model, vegetation, soil, and MODIS albedo parameter files for the 
! catchment model implementation
!  
! Sarith Mahanama - March 23, 2012 
! Email: sarith.p.mahanama@nasa.gov

   use rmTinyCatchParaMod
   use process_hres_data

  implicit none
    integer              :: NC = i_raster, NR = j_raster
    character*128        :: GridName = ''
    character*128        :: ARG, MaskFile
    character*256        :: CMD
    character*1          :: opt
    character*8          :: LD = 'MODGEO'
    character*4          :: SD = 'HWSD'
    character*4          :: EASE ='    '
    character*2          :: DL ='DC'
    character*6          :: MA = 'MODIS2'
    integer              :: II, JJ, Type
    integer              :: I, J, iargc, nxt
    real*8               :: dx, dy, lon0
    logical :: regrid
    character(len=400), dimension (20) ::  Usage 
    character*128        ::  Grid2
    character*2 :: poles
    CHARACTER*100 :: gfile,fname,pdir,rstdir
    character*128        :: GridNameR = ''
    character*128        :: GridNameT = ''
    logical :: file_exists
    logical, parameter :: F25Tag = .false.
    logical :: ease_grid=.false., redo_modis=.false.
    character*40       :: lai_name 
    integer, parameter :: log_file = 998
    include 'netcdf.inc'	
    type (regrid_map), allocatable, dimension (:,:) :: maparc30, mapgeoland2,maparc60
    logical :: running_omp = .false.
    character*200 :: tmpstring, tmpstring1, tmpstring2

    ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
    !
    ! FIND OUT WHETHER -omp FLAG HAS BEEN SET DURING COMPILATION

!$ running_omp = .true.         ! conditional compilation


!   call system('cd data/ ; ln -s /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/ CATCH')
!   call system('cd ..')

    USAGE(1) ="Usage: mkCatchParam -x nx -y ny -g Gridname -b DL -m MA -l LD -s SD -e EASE                       "
    USAGE(2) ="     -x: Size of longitude dimension of input raster. DEFAULT: 8640                               "
    USAGE(3) ="     -y: Size of latitude dimension of input raster.  DEFAULT: 4320                               "
    USAGE(4) ="     -g: Gridname  (name of the .til or .rst file without file extension)                         "
    USAGE(5) ="     -b: Position of the dateline in the first grid box (DC or DE). DEFAULT: DC                   "
    USAGE(6) ="     -l: Choice of LAI data set. DEFAULT : MODIS                                                  "
    USAGE(7) ="         MODISV6 : 8-day climatology from the period 2002.01-2016.10 on 86400x43200 grid          "
    USAGE(8) ="         MODGEO  : MODIS with GEOLAND2 overlaid on South America, Africa and Australia            "
    USAGE(9) ="         GEOLAND2: 10-day climatology from the period 1999-2011 on 40320x20160 grid               "
    USAGE(10)="         GSWP2   : Monthly climatology from the period 1982-1998 on 360x180 grid                  "
    USAGE(11)="         GSWPH   : Monthly climatology from the period 1982-1998 on 43200x21600 grid              "
    USAGE(12)="         MODIS   : 8-day climatology from the period 2000-2013  on 43200x21600 grid               "
    USAGE(13)="     -s: Choice of soil data. DEFAULT :HWSD                                                       "
    USAGE(14)="         HWSD : Merged HWSD-STATSGO2 soil properties on 43200x21600 with Woesten (1999) Parameters"
    USAGE(15)="         NGDC : Reynolds soil texture classes on 4320x2160 with GSWP2  soil hydraulic parameters  "
    USAGE(16)="     -m: Choice of MODIS Albedo data. DEFAULT : MODIS2                                            "
    USAGE(17)="         MODIS1: 16-day Climatology from 1'x1'(21600x10800) MODIS data from the period 2000-2004  "
    USAGE(18)="         MODIS2: 8-day Climatology from 0.5'x0.5'(43200x21600) MODIS data from the period 2001-2011"
    USAGE(19)="     -e: EASE : This is optional if catchment.def file is available already or                    "          
    USAGE(20)="                    the til file format is pre-Fortuna-2.                                         "           

! Process Arguments                            
!------------------ 

    CALL get_command (cmd)
    inquire(file='clsm/mkCatchParam.log', exist=file_exists)
    if(file_exists) then
       open (log_file, file ='clsm/mkCatchParam.log', status='old', position='append', form='formatted',action='write')
    else
       open (log_file, file ='clsm/mkCatchParam.log', status='unknown', form='formatted',action='write')
       write (log_file,'(a)')trim(cmd)
       write (log_file,'(a6,a3,L1)')'F25Tag',' : ' ,F25Tag
    endif

    I = iargc()

    if(I < 1 .or. I > 16) then
       write (log_file,'(a)') "Wrong Number of arguments: ", i
       do j = 1,size(usage)
          print "(sp,a100)", Usage(j)
       end do
       call exit(1)
    end if

    nxt = 1
    call getarg(nxt,arg)
    do while(arg(1:1)=='-')
       opt=arg(2:2)
       if(len(trim(arg))==2) then
          if(scan(opt,'zvh')==0) then
             nxt = nxt + 1
             call getarg(nxt,arg)
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
          GridName = trim(arg)
       case ('l')
          LD = trim(arg)
       case ('s')
          SD = trim(arg)
       case ('b')
          DL = trim(arg)
       case ('m')
          MA = trim(arg)
       case ('e')
          EASE = trim(arg)
	  if(EASE=='EASE') ease_grid=.true.
       case default
          do j = 1,size(usage)
             print "(sp,a100)", Usage(j)
          end do
          call exit(1)
       end select
       nxt = nxt + 1
       call getarg(nxt,arg)
    end do

    if(F25Tag) then
!
! Going back to f2.5 tag
!  
      LD ='GSWP2'
      MA ='MODIS1'
      SD ='HWSD'

    endif

   call getenv ("MASKFILE"        ,MaskFile        )
 
  if(trim(Gridname) == '') then
      write (log_file,'(a)')'Unable to create parameters without til/rst files.... !'
      stop
   endif
   
   regrid = nc/=i_raster .or. nr/=j_raster

   if(index(Gridname,'Pfaf.notiny')/=0) then 
      GridnameR='clsm/'//trim(Gridname)
      GridnameT='clsm/'//trim(Gridname)
      else
      GridnameR='rst/'//trim(Gridname)  
      GridnameT='til/'//trim(Gridname)  
    endif 

    if(running_omp .eq. .false.) then

       write (log_file,'(a)')trim(LD)
       write (log_file,'(a)')trim(MA)
       write (log_file,'(a)')trim(SD)   
       write (log_file,'(a)')trim(MaskFile)
       write (log_file,'(a)')'                              '
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'............ Begin CLSM parameter generation:...............'    
       write (log_file,'(a)')'                                               '
       write (log_file,'(a)')'CLSM parameters are being generated for the tile space :'
       write (log_file,'(a)')'     ',trim(GridnameT)
       write (log_file,'(a)')'                                               '
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'                                               '
       
       if(index(Gridname,'CF')/=0) then 
          DL = 'DE'
          write (log_file,'(a)')'Cube Grid - assuming DE'
       endif
       
       ! Creating catchment.def 
       ! ----------------------
       
       if(.not.ease_grid) then  
          inquire(file='clsm/catchment.def', exist=file_exists)
          if (.not.file_exists) call catchment_def (nc,nr,regrid,dl,gridnamet,gridnamer) 
          write (log_file,'(a)')'Done creating catchment.def file ..........................1'
       endif
       
       ! Creating cti_stats.dat 
       ! ----------------------
       
       inquire(file='clsm/cti_stats.dat', exist=file_exists)
       if (.not.file_exists) call cti_stat_file (ease_grid,gridnamet, MaskFile)
       write (log_file,'(a)')'Done creating CTI stat file ...............................2'	 
       
       ! Creating vegetation classification files
       !-----------------------------------------
       
       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
          
          inquire(file='clsm/mosaic_veg_typs_fracs', exist=file_exists)      
          if (.not.file_exists) call ESA2MOSAIC (nc,nr,gridnamer)
          inquire(file='clsm/CLM_veg_typs_fracs', exist=file_exists) 
          if (.not.file_exists) call ESA2CLM (nc,nr,gridnamer)    
          inquire(file='clsm/CLM4.5_veg_typs_fracs', exist=file_exists) 
          if (.not.file_exists) call ESA2CLM_45 (nc,nr,gridnamer)           
          write (log_file,'(a)')'Done creating vegetation types using ESA land cover........3'
          
       else
          
          inquire(file='clsm/mosaic_veg_typs_fracs', exist=file_exists)
          call compute_mosaic_veg_types (nc,nr,ease_grid,regrid,gridnamet,gridnamer)
          
          write (log_file,'(a)')'Done creating vegetation types using IGBP SiB2 land cover..3'
          
       endif
       
       ! Processing Vegetation Climatology 
       ! ---------------------------------
       
       ! creating mapping arrays if necessary
       
       if((trim(LD) == 'MODGEO').or.(trim(LD) == 'GEOLAND2')) then 
          inquire(file='clsm/lai.GEOLAND2_10-DayClim', exist=file_exists)
          if (.not.file_exists) then
             allocate (mapgeoland2 (1:40320,1:20160))
             call create_mapping (nc,nr,40320,20160,mapgeoland2, gridnamer)         
             lai_name = 'GEOLAND2_10-DayClim/geoland2_' 
             if(trim(LD) == 'GEOLAND2') then
                call hres_lai_no_gswp (40320,20160,mapgeoland2,gridnamer, lai_name) 
             else
                call hres_lai_no_gswp (40320,20160,mapgeoland2,gridnamer, lai_name, merge=1) 
             endif
             if(allocated(mapgeoland2)) deallocate (mapgeoland2)
          endif
       endif
       
       if ((LD == 'MODGEO').or.(LD == 'MODIS').or.(MA == 'MODIS2')) then
          allocate (maparc30    (1:43200,1:21600))
          call create_mapping (nc,nr,43200,21600,maparc30,    gridnamer)
       endif
       
       inquire(file='clsm/green.dat', exist=file_exists)
       
       if (.not.file_exists) then
          if (trim(LD) == 'GSWP2') then 
             call process_gswp2_veg (nc,nr,regrid,'grnFrac',gridnamer)
          else
             if (.not. allocated(maparc30)) then 
                allocate (maparc30    (1:43200,1:21600))
                call create_mapping (nc,nr,43200,21600,maparc30,    gridnamer)
             endif
             call hres_gswp2 (43200,21600, maparc30, gridnamer,'green') 
          endif
       endif
       
       inquire(file='clsm/lai.dat', exist=file_exists)
       
       if (.not.file_exists) then
          redo_modis = .true.
          
          if (trim(LD) == 'GSWP2') call process_gswp2_veg (nc,nr,regrid,'LAI',gridnamer) 
          if (trim(LD) == 'GSWPH') then
             if (.not. allocated(maparc30)) then 
                allocate (maparc30    (1:43200,1:21600))
                call create_mapping (nc,nr,43200,21600,maparc30,    gridnamer)
             endif
             inquire(file='clsm/lai.MODIS_8-DayClim', exist=file_exists)
             if (.not.file_exists) call hres_gswp2 (43200,21600, maparc30, gridnamer,'lai') 
          endif
          
          if (trim(LD) == 'MODIS') then
             lai_name = 'MODIS_8-DayClim/MODIS_'
             call hres_lai_no_gswp (43200,21600,maparc30,gridnamer,lai_name) 
          endif
          
          if (trim(LD) == 'MODGEO') then
             lai_name = 'MODIS_8-DayClim/MODIS_'
             inquire(file='clsm/lai.MODIS_8-DayClim', exist=file_exists)
             if (.not.file_exists)call hres_lai_no_gswp (43200,21600,maparc30,gridnamer,lai_name, merge=1)  
             call merge_lai_data (MaskFile)
          endif
          
          if (trim(LD) == 'MODISV6') then
             lai_name = 'MCD15A2H.006/MODIS_'
             call grid2tile_modis6 (86400,43200,nc,nr,gridnamer,lai_name)  
          endif
          
       endif

       inquire(file='clsm/ndvi.dat', exist=file_exists)
       if (.not.file_exists)  call gimms_clim_ndvi (nc,nr,gridnamer)

       write (log_file,'(a,a,a)')'Done computing ', trim(LD),' vegetation climatologies ............4'
  
       ! call modis_alb_on_tiles (nc,nr,ease_grid,regrid,gridnamet,gridnamer)
       ! call modis_scale_para (ease_grid,gridnamet)
       ! NOTE: modis_alb_on_tiles uses monthly climatological raster data on 8640x4320 to produce 
       ! MODIS albedo on tile space. The subroutine was replaced with "modis_alb_on_tiles_high" that process
       ! MODIS1 data on native grid and produces 8/16-day MODIS Albedo climatology
       
       if(MA == 'MODIS1') then 
          inquire(file='clsm/AlbMap.WS.16-day.tile.0.7_5.0.dat', exist=file_exists)
          if (.not.file_exists) then
             if(.not.F25Tag) then 
                allocate (maparc60    (1:21600,1:10800))
                call create_mapping (nc,nr,21600,10800,maparc60,    gridnamer)
                call modis_alb_on_tiles_high (21600,10800,maparc60,MA,gridnamer)
                if(allocated (maparc30)) deallocate (maparc60)
             else
                call modis_alb_on_tiles (nc,nr,ease_grid,regrid,gridnamet,gridnamer)
             endif
          endif
       endif
       
       if(MA == 'MODIS2') then 
          inquire(file='clsm/AlbMap.WS.8-day.tile.0.7_5.0.dat', exist=file_exists)
          if (.not.file_exists) call modis_alb_on_tiles_high (43200,21600,maparc30,MA,gridnamer)
       endif
       write (log_file,'(a,a,a)')'Done putting ',trim(MA), ' Albedo on the tile space  .............5'
       
       if(allocated (maparc30)) deallocate (maparc30)
       
       inquire(file='clsm/visdf.dat', exist=file_exists)
       if ((redo_modis).or.(.not.file_exists)) then
          if(.not.F25Tag) then
             call modis_scale_para_high (ease_grid,MA,gridnamet)
          else
             inquire(file='clsm/modis_scale_factor.albvf.clim', exist=file_exists)
             if ((redo_modis).or.(.not.file_exists)) then
                call modis_scale_para (ease_grid,gridnamet)
                call REFORMAT_VEGFILES
             endif
          endif
       endif
       
       write (log_file,'(a,a,a)')'Done computing ',trim(MA), ' scale factors .......................6'
       tmpstring1 = '-e EASE -g '//trim(gfile) 
       write(tmpstring2,'(2(a2,x,i5,x))')'-x',nc,'-y',nr
       tmpstring = 'bin/mkCatchParam_openmp '//trim(tmpstring2)//' '//trim(tmpstring1)

    else      
       if(SD=='NGDC') call create_soil_types_files (nc,nr,ease_grid,gridnamet,gridnamer)    
       if(SD=='NGDC') write (log_file,'(a)')'Done creating NGDC soil types file .......................7a'	   
       
       ! Creating soil_param.first and tau_param.dat files that has 2 options: 
       !  1) NGDC soil properties, 2) HWSD-STATSGO2 Soil Properties
       ! ---------------------------------------------------------------------
       
       inquire(file='clsm/soil_param.first', exist=file_exists)
       if (.not.file_exists) then
          if(SD=='NGDC')  then 
             if(F25Tag) call soil_para_high (nc,nr,regrid,gridnamer,F25Tag=F25Tag)
             if(.not.F25Tag) call soil_para_high (nc,nr,regrid,gridnamer)
          endif
          
          if(SD=='HWSD')  call soil_para_hwsd (nc,nr,gridnamer) 
       endif
       write (log_file,'(a,a,a)')'Done computing ',trim(SD),' soil parameters .......................7'
       
       
       inquire(file='clsm/ts.dat', exist=file_exists)
       if (.not.file_exists) then
          if(SD=='NGDC') call create_model_para (MaskFile)
          if(SD=='HWSD') call create_model_para_woesten (MaskFile) 
       endif
       write (log_file,'(a,a,a)')'Done computing CLSM model parameters based on ',trim(SD),'.........8'

       ! Commented out this call because 7.5-minute raster file is only used
       ! for plotting purposes
       !  call make_75 (nc,nr,regrid,c_data,gridnamer)
       !    write (log_file,'(a)')'Done creating 7.5 minute raster file ......................'
       !    write (log_file,'(a)')'Not created 7.5 minute raster file   ......................'
       
       inquire(file='clsm/CLM_veg_typs_fracs', exist=file_exists)
       if (file_exists) then
          
          call grid2tile_ndep_t2m_alb (nc,nr,gridnamer)  
          write (log_file,'(a)')'Done computing CLSM-CN NDep T2m SoilAlb ...................9'
          
       endif

       inquire(file='clsm/CLM4.5_abm_peatf_gdp_hdm_fc', exist=file_exists) 
       if (.not.file_exists) call CLM45_fixed_parameters (nc,nr,gridnamer)           
       write (log_file,'(a)')'Done creating CLM4.5_abm_peatf_gdp_hdm_fc ................10'

       inquire(file='clsm/lnfm.dat', exist=file_exists)
       if (.not.file_exists) call CLM45_clim_parameters (nc,nr,gridnamer)   
       write (log_file,'(a)')'Done creating CLM4.5 lightening frequency clim ...........11'
       
       
       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'DONE creating CLSM data files...............................'
       write (log_file,'(a)')'============================================================'
              
       call system ('chmod 755 src/create_README.csh ; src/create_README.csh')
    endif

    close (log_file,status='keep') 

END PROGRAM mkCatchParam
