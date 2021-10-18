 PROGRAM mkCatchParam

! !INTERFACE:
!
! !ARGUMENTS:
!
!  Usage = "mkCatchParam -x nx -y ny -g Gridname -b DL -v LBSV -e EASE"       
!     -x: Size of longitude dimension of input raster. DEFAULT: 8640
!     -y: Size of latitude dimension of input raster.  DEFAULT: 4320
!     -b: position of the dateline in the first box. DEFAULT: DC 
!     -g: Gridname  (name of the .til or .rst file without file extension)  
!     -v: LBCSV : use a configuration from GEOS5 bcs directory ICA, NL3, NL4, or NL4p              
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
   use comp_CATCHCN_AlbScale_parameters, ONLY : albedo4catchcn
!   use module_irrig_params, ONLY : create_irrig_params

  implicit none
  integer                :: NC = i_raster, NR = j_raster
    character*4          :: LBSV = 'DEF'
    character*128        :: GridName = ''
    character*128        :: ARG, MaskFile
    character*256        :: CMD
    character*1          :: opt
    character*7          :: PEATSOURCE   = 'GDLHWSD'
    character*3          :: VEGZSOURCE   = 'D&S'
    character*4          :: EASE ='    '
    character*2          :: DL ='DC'    
    integer              :: II, JJ, Type
    integer              :: I, J, iargc, nxt
    real*8               :: dx, dy, lon0
    logical :: regrid
    character(len=400), dimension (8) ::  Usage 
    character*128        ::  Grid2
    character*2 :: poles
    CHARACTER*100 :: gfile,fname,pdir,rstdir
    character*128        :: GridNameR = ''
    character*128        :: GridNameT = ''
    logical :: file_exists
    logical             :: F25Tag = .false.
    logical :: ease_grid=.false., redo_modis=.false.
    character*40       :: lai_name 
    integer, parameter :: log_file = 998
    include 'netcdf.inc'	
    type (regrid_map) :: maparc30, mapgeoland2,maparc60
    character*200 :: tmpstring, tmpstring1, tmpstring2
    

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

!   call system('cd data/ ; ln -s /discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/ CATCH')
!   call system('cd ..')

    USAGE(1) ="Usage: mkCatchParam -x nx -y ny -g Gridname -b DL -v LBCSV -e EASE                                "
    USAGE(2) ="     -x: Size of longitude dimension of input raster. DEFAULT: 8640                               "
    USAGE(3) ="     -y: Size of latitude dimension of input raster.  DEFAULT: 4320                               "
    USAGE(4) ="     -g: Gridname  (name of the .til or .rst file without file extension)                         "
    USAGE(5) ="     -b: Position of the dateline in the first grid box (DC or DE). DEFAULT: DC                   "
    USAGE(6) ="     -e: EASE : This is optional if catchment.def file is available already or                    "          
    USAGE(7) ="                the til file format is pre-Fortuna-2.                                             "
    USAGE(8) ="     -v  LBCSV : use a configuration from GEOS5 bcs directory F25, GM4, ICA, NL3, NL4, or NL4p    "

! Process Arguments                            
!------------------ 

    CALL get_command (cmd)
    inquire(file='clsm/mkCatchParam.log', exist=file_exists)
    if(file_exists) then
       open (log_file, file ='clsm/mkCatchParam.log', status='old', position='append', form='formatted',action='write')
    else
       open (log_file, file ='clsm/mkCatchParam.log', status='unknown', form='formatted',action='write')
       write (log_file,'(a)')trim(cmd)
       write (log_file,'(a)')' '
    endif

    I = iargc()
    if(I < 1 .or. I > 10) then
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
          if(scan(opt,'zh')==0) then
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
       case ('v')
          LBSV = trim(arg)
          if (trim(arg).eq."F25") F25Tag = .true.
          call init_bcs_config (trim(LBSV))
       case ('b')
          DL = trim(arg)
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

    if(process_peat) PEATSOURCE   = 'PEATMAP'
    if(jpl_height)   VEGZSOURCE   = 'JPL'

    if(n_threads == 1) then

       write (log_file,'(a)')trim(LAIBCS)
       write (log_file,'(a)')trim(MODALB)
       write (log_file,'(a)')trim(SOILBCS)   
       write (log_file,'(a)')trim(MaskFile)
       write (log_file,'(a)')trim(PEATSOURCE)
       write (log_file,'(a)')trim(VEGZSOURCE)
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
       
       inquire(file='clsm/catch_params.nc4', exist=file_exists)
       if (.not.file_exists) CALL open_landparam_nc4_files 

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
       
       if((trim(LAIBCS) == 'MODGEO').or.(trim(LAIBCS) == 'GEOLAND2')) then 
          inquire(file='clsm/lai.GEOLAND2_10-DayClim', exist=file_exists)
          if (.not.file_exists) then
             !allocate (mapgeoland2 (1:40320,1:20160))
             call create_mapping (nc,nr,40320,20160,mapgeoland2, gridnamer)         
             lai_name = 'GEOLAND2_10-DayClim/geoland2_' 
             if(trim(LAIBCS) == 'GEOLAND2') then
                call hres_lai_no_gswp (40320,20160,mapgeoland2,gridnamer, lai_name) 
             else
                call hres_lai_no_gswp (40320,20160,mapgeoland2,gridnamer, lai_name, merge=1) 
             endif
             ! if(allocated(mapgeoland2)) deallocate (mapgeoland2)
             deallocate (mapgeoland2%map)
             deallocate (mapgeoland2%ij_index)
          endif
       endif
       
       if ((LAIBCS == 'MODGEO').or.(LAIBCS == 'MODIS').or.(MODALB == 'MODIS2')) then
          ! allocate (maparc30    (1:43200,1:21600))
          call create_mapping (nc,nr,43200,21600,maparc30,    gridnamer)
       endif
       
       inquire(file='clsm/green.dat', exist=file_exists)
       
       if (.not.file_exists) then
          if (trim(LAIBCS) == 'GSWP2') then 
             call process_gswp2_veg (nc,nr,regrid,'grnFrac',gridnamer)
          else
             if (size(maparc30%ij_index,1) /= 43200) then 
                ! allocate (maparc30    (1:43200,1:21600))
                call create_mapping (nc,nr,43200,21600,maparc30,    gridnamer)
             endif
             call hres_gswp2 (43200,21600, maparc30, gridnamer,'green') 
          endif
       endif
       
       inquire(file='clsm/lai.dat', exist=file_exists)
       
       if (.not.file_exists) then
          redo_modis = .true.
          
          if (trim(LAIBCS) == 'GSWP2') call process_gswp2_veg (nc,nr,regrid,'LAI',gridnamer) 
          if (trim(LAIBCS) == 'GSWPH') then
             if (size(maparc30%ij_index,1) /= 43200) then 
                ! allocate (maparc30    (1:43200,1:21600))
                call create_mapping (nc,nr,43200,21600,maparc30,    gridnamer)
             endif
             inquire(file='clsm/lai.MODIS_8-DayClim', exist=file_exists)
             if (.not.file_exists) call hres_gswp2 (43200,21600, maparc30, gridnamer,'lai') 
          endif
          
          if (trim(LAIBCS) == 'MODIS') then
             lai_name = 'MODIS_8-DayClim/MODIS_'
             call hres_lai_no_gswp (43200,21600,maparc30,gridnamer,lai_name) 
          endif
          
          if (trim(LAIBCS) == 'MODGEO') then
             lai_name = 'MODIS_8-DayClim/MODIS_'
             inquire(file='clsm/lai.MODIS_8-DayClim', exist=file_exists)
             if (.not.file_exists)call hres_lai_no_gswp (43200,21600,maparc30,gridnamer,lai_name, merge=1)  
             call merge_lai_data (MaskFile)
          endif
          
          if (trim(LAIBCS) == 'MODISV6') then
             lai_name = 'MCD15A2H.006/MODIS_'
             call grid2tile_modis6 (86400,43200,nc,nr,gridnamer,lai_name)  
          endif

          if (trim(LAIBCS) == 'GLASSA') then
             lai_name = 'GLASS-LAI/AVHRR.v4/GLASS01B02.V04.AYYYY'
             call grid2tile_glass (nc,nr,gridnamer,lai_name)  
          endif

          if (trim(LAIBCS) == 'GLASSM') then
             lai_name = 'GLASS-LAI/MODIS.v4/GLASS01B01.V04.AYYYY'
             call grid2tile_glass (nc,nr,gridnamer,lai_name)  
          endif
          
       endif

       inquire(file='clsm/ndvi.dat', exist=file_exists)
       if (.not.file_exists)  call gimms_clim_ndvi (nc,nr,gridnamer)

       write (log_file,'(a,a,a)')'Done computing ', trim(LAIBCS),' vegetation climatologies .............4'
  
       ! call modis_alb_on_tiles (nc,nr,ease_grid,regrid,gridnamet,gridnamer)
       ! call modis_scale_para (ease_grid,gridnamet)
       ! NOTE: modis_alb_on_tiles uses monthly climatological raster data on 8640x4320 to produce 
       ! MODIS albedo on tile space. The subroutine was replaced with "modis_alb_on_tiles_high" that process
       ! MODIS1 data on native grid and produces 8/16-day MODIS Albedo climatology
       
       if(MODALB == 'MODIS1') then 
          inquire(file='clsm/AlbMap.WS.16-day.tile.0.7_5.0.dat', exist=file_exists)
          if (.not.file_exists) then
             if(F25Tag) then 
                call create_mapping (nc,nr,21600,10800,maparc60,    gridnamer)
                call modis_alb_on_tiles_high (21600,10800,maparc60,MODALB,gridnamer)
                deallocate (maparc60%map)
                deallocate (maparc60%ij_index)
             else
              !  This option is for legacy sets like Fortuna 2.1
                call modis_alb_on_tiles (nc,nr,ease_grid,regrid,gridnamet,gridnamer)
             endif
          endif
       endif
       
       if(MODALB == 'MODIS2') then 
          inquire(file='clsm/AlbMap.WS.8-day.tile.0.7_5.0.dat', exist=file_exists)
          if (.not.file_exists) call modis_alb_on_tiles_high (43200,21600,maparc30,MODALB,gridnamer)
       endif
       write (log_file,'(a,a,a)')'Done putting ',trim(MODALB), ' Albedo on the tile space  .............5'
       
       if(.not.F25Tag) then 
          deallocate (maparc30%map)
          deallocate (maparc30%ij_index)
       endif

       inquire(file='clsm/visdf.dat', exist=file_exists)
       if ((redo_modis).or.(.not.file_exists)) then
       !   if(.not.F25Tag) then
             call modis_scale_para_high (ease_grid,MODALB,gridnamet)
        !  else
        !     This option is for legacy sets like Fortuna 2.1
        !     inquire(file='clsm/modis_scale_factor.albvf.clim', exist=file_exists)
        !     if ((redo_modis).or.(.not.file_exists)) then
        !        call modis_scale_para (ease_grid,gridnamet)
        !        call REFORMAT_VEGFILES
        !     endif
        !  endif
       endif
       
       write (log_file,'(a,a,a)')'Done computing ',trim(MODALB), ' scale factors .......................6'
!       tmpstring1 = '-e EASE -g '//trim(gfile) 
!       write(tmpstring2,'(2(a2,x,i5,x))')'-x',nc,'-y',nr
!       tmpstring = 'bin/mkCatchParam_openmp '//trim(tmpstring2)//' '//trim(tmpstring1)

    else      
 
       if(SOILBCS=='NGDC') call create_soil_types_files (nc,nr,ease_grid,gridnamet,gridnamer)    
       if(SOILBCS=='NGDC') write (log_file,'(a)')'Done creating NGDC soil types file .......................7a'	   
       
       ! Creating soil_param.first and tau_param.dat files that has 2 options: 
       !  1) NGDC soil properties, 2) HWSD-STATSGO2 Soil Properties
       ! ---------------------------------------------------------------------
       
       inquire(file='clsm/soil_param.first', exist=file_exists)
       if (.not.file_exists) then
          if(SOILBCS=='NGDC')  then 
             if(F25Tag) call soil_para_high (nc,nr,regrid,gridnamer,F25Tag=F25Tag)
             if(.not.F25Tag) call soil_para_high (nc,nr,regrid,gridnamer)
          endif
          
          if(SOILBCS=='HWSD')  call soil_para_hwsd (nc,nr,gridnamer) 
       endif
       write (log_file,'(a,a,a)')'Done computing ',trim(SOILBCS),' soil parameters .......................7'
       
       
       inquire(file='clsm/ts.dat', exist=file_exists)
       if (.not.file_exists) then
          if(SOILBCS=='NGDC') call create_model_para (MaskFile)
          if(SOILBCS=='HWSD') call create_model_para_woesten (MaskFile) 
       endif
       write (log_file,'(a,a,a)')'Done computing CLSM model parameters based on ',trim(SOILBCS),'.........8'

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

       inquire(file='clsm/country_and_state_code.data', exist=file_exists)
       if (.not.file_exists) call map_country_codes (nc,nr,gridnamer)
       write (log_file,'(a)')'Done mapping country and state codes .....................12'

 !      inquire(file='clsm/irrig.dat', exist=file_exists)
 !      if (.not.file_exists) call create_irrig_params (nc,nr,gridnamer)
 !      write (log_file,'(a)')'Done computing irrigation model parameters ...............13'

     !   call albedo4catchcn (gridnamet)

       write (log_file,'(a)')'============================================================'
       write (log_file,'(a)')'DONE creating CLSM data files...............................'
       write (log_file,'(a)')'============================================================'
              
!       call system ('chmod 755 bin/create_README.csh ; bin/create_README.csh')
    endif

    close (log_file,status='keep') 

END PROGRAM mkCatchParam
