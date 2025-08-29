#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif
!
! A Collection subroutines needed by mkCatchParam.F90
!   Contact: Sarith Mahanama  sarith.p.mahanama@nasa.gov
!   Email  : sarith.p.mahanama@nasa.gov

module rmTinyCatchParaMod

  use LDAS_DateTimeMod
  use MAPL_ConstantsMod
  use MAPL_Base,           ONLY: MAPL_UNDEF
  use MAPL,                only: MAPL_WriteTilingNC4
  use lsm_routines,        ONLY: sibalb
  use LogRectRasterizeMod, ONLY: SRTM_maxcat, MAPL_UNDEF_R8 
  use, intrinsic :: iso_fortran_env, only: REAL64 
  implicit none
  
  logical, parameter :: error_file=.true.
  integer, parameter :: n_SoilClasses = 253
  real,    parameter :: zks = 2.0
  integer, parameter :: i_raster = 8640, j_raster=4320
  integer, parameter :: ncat_gswp2 = 15238
  REAL,    PARAMETER :: undef = 1.e+20
  integer, parameter :: arr_len = 1734915
  integer, parameter :: ip1 =0                       ! index offset for land tiles within all vector of all tiles (ip1=0 => land tiles are first) 
  real,    parameter :: dx_gswp2 =1.,dy_gswp2=1. 
  integer, parameter :: MAX_NOF_GRID = ncat_gswp2
  integer, PARAMETER :: nbdep=150, NAR=1000,nwt=81,nrz=41
  real,    parameter :: slice=0.1, lim =5.,grzdep =1.1
  logical, parameter :: bug =.false.

  include 'netcdf.inc'

  logical :: preserve_soiltype = .false.
  
  ! Bugfix for Target_mean_land_elev: 
  ! Previously, the hardcoded value 614.649 m was used as the target mean land elevation.
  ! This was incorrect because it did not account for the proper cosine-lat-weighted mean over land.
  ! The correct value, 656.83 m, is derived from NCAR GMTED TOPO 30arcsec dataset.
  ! This ensures that land elevation adjustment is based on the correct reference mean land elevation.

  real*8, parameter :: Target_mean_land_elev = 656.83D0  ! cosine-lat-weighted mean land elev from NCAR GMTED TOPO 30arcsec 
  
  private
  
  public Target_mean_land_elev
  public modis_alb_on_tiles
  public supplemental_tile_attributes, soil_para_high
  public create_soil_types_files, compute_mosaic_veg_types
  public cti_stat_file, create_model_para_woesten
  public create_model_para, regridraster, regridrasterreal
  public i_raster, j_raster, regridraster1, regridraster2, n_SoilClasses, zks
  public mineral_perc, process_gswp2_veg, center_pix, soil_class
  public REFORMAT_VEGFILES
  public Get_MidTime, Time_Interp_Fac
  public ascat_r0, jpl_canoph, NC_VarID, init_bcs_config  
  
  ! The following variables define the details of the BCS version (data sources).
  ! Initialize to dummy values here and set to desired values in init_bcs_config().
  
  logical,      public, save :: use_PEATMAP = .false.
  logical,      public, save :: jpl_height  = .false.
  character*8,  public, save :: LAIBCS      = 'UNDEF'
  character*6,  public, save :: SOILBCS     = 'UNDEF'
  character*6,  public, save :: MODALB      = 'UNDEF'
  character*10, public, save :: SNOWALB     = 'UNDEF'
  character*5,  public, save :: OUTLETV     = 'UNDEF'  
  REAL,         public, save :: GNU         = MAPL_UNDEF
  
  character*512              :: MAKE_BCS_INPUT_DIR
  
  type :: mineral_perc
     real :: clay_perc
     real :: silt_perc
     real :: sand_perc
  end type mineral_perc

contains
  
  SUBROUTINE init_bcs_config(LBCSV)
    
    ! determine BCs details from land BCs version string (LBCSV)
    !
    ! LAIBCS:  Leaf-Area-Index data set.        DEFAULT : MODGEO
    !   GLASSA    : 8-day AVHRR clim, 1981-2017,        7200x3600  grid
    !   GLASSM    : 8-day MODIS clim, 2000-2017,        7200x3600  grid
    !   MODISV6   : 8-day       clim, 2002.01-2016.10, 86400x43200 grid
    !   MODGEO    : MODIS with GEOLAND2 overlaid on South America, Africa, and Australia
    !   GEOLAND2  : 10-day  clim,     1999-2011,       40320x20160 grid               
    !   GSWP2     : Monthly clim,     1982-1998,         360x180   grid                  
    !   MODIS     : 8-day   clim,     2000-2013,       43200x21600 grid
    !   GSWPH     : Monthly clim,     1982-1998,       43200x21600 grid           
    !
    ! MODALB:  MODIS Albedo data (snow-free).   DEFAULT : MODIS2                                            
    !   MODIS1    : 16-day clim,  1'x1' (21600x10800) MODIS data, 2000-2004 
    !   MODIS2    :  8-day clim, 30"x30"(43200x21600) MODIS data, 2001-2011 
    !
    ! SNOWALB: Snow albedo data.                DEFAULT : LUT
    !   LUT       : Parameterization based on look-up table values. 
    !   MODC061   : Static snow albedo derived from MODIS Collection 6.1 data where available, fill value of 0.56 elsewhere. 
    !   MODC061v2 : Same as MODC061 but using tile ID instead of tile bounding box for mapping from raster to tile.
    !
    ! SOILBCS: Soil parameter data.             DEFAULT : HWSD    
    !   NGDC      : Soil parameters from Reynolds et al. 2000, doi:10.1029/2000WR900130 (MERRA-2, Fortuna, Ganymed, Icarus)
    !   HWSD      : Merged HWSDv1.21-STATSGO2 soil properties on 43200x21600 with Woesten et al. (1999) parameters   
    !   HWSD_b    : As in HWSD but with surgical fix of Argentina peatland issue (38S,60W)
    !
    ! OUTLETV: Definition of outlet locations.  DEFAULT : N/A
    !   N/A       : No information (do not create routing "TRN" files).
    !   v1        : Outlet locations file produced manually by Randy Koster.
    !   v2        : Outlet locations file produced by run_routing_raster.py using routing information encoded 
    !               in SRTM-based Pfafstetter catchments and Greenland outlets info provided by Lauren Andrews.

    character(*), intent (in) :: LBCSV     ! land BCs version 

    select case (trim(LBCSV))
       
    case ("F25")
       LAIBCS  = 'GSWP2'
       SOILBCS = 'NGDC'
       MODALB  = 'MODIS1'
       SNOWALB = 'LUT'
       OUTLETV = "N/A"
       GNU     = 2.17
       use_PEATMAP = .false.
       jpl_height  = .false.
       
    case ("GM4", "ICA")
       LAIBCS  = 'GSWP2'
       SOILBCS = 'NGDC'
       MODALB  = 'MODIS2'
       SNOWALB = 'LUT'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .false.
       jpl_height  = .false.

    case ("NL3")
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'LUT'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .false.
       jpl_height  = .false.

    case ("NL4")
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'      
       SNOWALB = 'LUT'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .false.
       jpl_height  = .true.

    case ("NL5")
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'LUT'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .true.

    case ("v06")   
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'MODC061'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .true.

    case ("v07")   
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'LUT'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .false.
       
    case ("v08")   
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'MODC061'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .false.
       jpl_height  = .false.
       
    case ("v09")   
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'MODC061'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .false.

    case ("v10")   
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'MODC061v2'
       OUTLETV = "N/A"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .false.

    case ("v11")   
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD'
       MODALB  = 'MODIS2'
       SNOWALB = 'MODC061v2'
       OUTLETV = "v1"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .true.

    case ("v12","v13","v14")  

       ! "v12", "v13", and "v14" are identical except for:
       ! - topography used for the atm (processed outside of make_bcs)
       ! - bug fix for land elevation in catchment.def file
       ! - generation of nc4-formatted tile file
       ! - v14 is used for coupled atm-ocean-seaice with MOM6/v2 (OM4) ocean bathymetry
 
       LAIBCS  = 'MODGEO'
       SOILBCS = 'HWSD_b'
       MODALB  = 'MODIS2'
       SNOWALB = 'MODC061v2'
       OUTLETV = "v2"
       GNU     = 1.0
       use_PEATMAP = .true.
       jpl_height  = .true.

    case default
       
       print *,'init_bcs_config(): unknown land boundary conditions version (LBCSV)'
       stop
       
    end select
    
  END SUBROUTINE init_bcs_config
  
  ! --------------------------------------------------------------------------------------------
  
  SUBROUTINE Get_MidTime (                           &
                          yr1,mn1,dy1,yr2,mn2,dy2,   &
                          MIDT)
    
    real,                 intent(in)  :: yr1,mn1,dy1,yr2,mn2,dy2
    type(date_time_type), intent(out) :: MIDT
    
    ! ------------
    
    type(date_time_type)              :: TIME1, TIME2
    integer                           :: TIMEDIF
    
    TIME1%year  = NINT(yr1) + 2001
    TIME1%month = NINT(mn1)
    TIME1%day   = NINT(dy1)
    TIME1%hour  = 0
    TIME1%min   = 0
    TIME1%sec   = 0
    
    call get_dofyr_pentad(TIME1)
    MIDT = TIME1

    TIME2%year  = NINT(yr2) + 2001
    TIME2%month = NINT(mn2)
    TIME2%day   = NINT(dy2)
    TIME2%hour  = 23
    TIME2%min   = 59
    TIME2%sec   = 59
    call get_dofyr_pentad(TIME2)
    
    TIMEDIF = datetime2_minus_datetime1(TIME1,TIME2)
    TIMEDIF = TIMEDIF/2
    
    call augment_date_time(TIMEDIF, MIDT) 
    
    !    print *,'MIDTIME'
    !    print *,'TIME1:', time1
    !    print *,'MIDT :', midt
    !    print *,'TIME2:', time2
    
  END SUBROUTINE Get_MidTime
  
  ! ---------------------------------------------------------------------------------------------
  
  SUBROUTINE Time_Interp_Fac (TIME0, TIME1, TIME2, FAC1, FAC2)

    !  PURPOSE:
    !  ========
    !
    !    Compute interpolation factors, fac, to be used 
    !    in the calculation of the instantaneous boundary 
    !    conditions, ie:
    !
    !     q(i,j) = fac1*q1(i,j) + (1.-fac1)*q2(i,j)
    !
    !    where:
    !     q(i,j)  => Boundary Data valid    at time0
    !     q1(i,j) => Boundary Data centered at time1
    !     q2(i,j) => Boundary Data centered at time2
    
    !  INPUT:
    !  ======
    !    time0    : Time of current timestep
    !    time1    : Time of boundary data 1 
    !    time2    : Time of boundary data 2 
    
    !  OUTPUT:
    !  =======
    !     fac1    : Interpolation factor for Boundary Data 1
    !
    
    type(date_time_type),   intent(in ) :: TIME0, TIME1, TIME2
    real,                   intent(out) :: FAC1
    real,                   intent(out) :: FAC2

    real        :: TimeDif1
    real        :: TimeDif

    !    print *,'Interpolation'
    !    print *,'TIME1:', time1
    !    print *,'TIME0:', time0
    !    print *,'TIME2:', time2
    
    TimeDif1 = real(datetime2_minus_datetime1(TIME0,TIME2))
    TimeDif  = real(datetime2_minus_datetime1(TIME1,TIME2))
    
    FAC1 = TimeDif1/TimeDif
    
    FAC2 = 1.-FAC1
    
    !    print *,fac1,fac2
    
  END SUBROUTINE Time_Interp_Fac
  
  ! ---------------------------------------------------------------------
  
  SUBROUTINE process_gswp2_veg (nc,nr,regrid,vname, n_land, tile_id,merge)
    
    integer,      intent(in)           :: nc, nr
    logical,      intent(in)           :: regrid
    integer,      intent(in)           :: n_land
    integer,      intent(in)           :: tile_id(:,:)
    character(*), intent(in)           :: vname

    integer,      intent(in), optional :: merge
    
    ! -------------------------------------------------------------    

    integer, parameter                           :: MAX_NOF_GRID = 15238
    REAL,    PARAMETER                           :: undef        = 1.e+20
    REAL,    PARAMETER                           :: UNDEF_GSWP2  = -9999.	
    
    real,    allocatable,         dimension(:)   :: catforc,vecforc,catcount	  
    integer, allocatable, target, dimension(:,:) :: gswp2_mask
    REAL,    ALLOCATABLE                         :: mon_climate(:,:)	
    integer                                      :: ierr, ncid,iret
    integer                                      :: i1,k1,n,iv,year,smon,imon,mon,i,j,status
    integer                                      :: k,ncatch
    integer                                      :: yr,mn,yr1,mn1
    integer, pointer                             :: Raster(:,:)
    character*512                                :: fname

    ! -----------------------------------------------------------------
 
    
    allocate (gswp2_mask (1:i_raster,1:j_raster))
    
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/shared/mask/mapping_2.5_grid_to_gswp2_tile_index.rst',&
         form='unformatted',status='old',action='read',convert='little_endian')
    
    do j =1,j_raster
       read (10) gswp2_mask(:,j)
    end do
    close (10,status='keep')
    
    if(regrid) then
       allocate(raster(nc,nr),stat=STATUS); VERIFY_(STATUS)
    else
       raster => gswp2_mask
    end if
    
    if(regrid) then
       call RegridRaster(gswp2_mask,raster)
    endif
   
    allocate(vecforc(1:MAX_NOF_GRID))
    allocate(catforc(n_land))
    allocate(catcount(n_land))
    allocate(mon_climate(1:n_land,1:12))
    mon_climate(:,:)=0.
    catforc=0.
    
    mon_climate(:,:)=0.
    
    iret = NF_OPEN(trim(MAKE_BCS_INPUT_DIR)//'/land/veg/lai_grn/'//trim(vname)//'_uk.nc',NF_NOWRITE, ncid)
    
    ASSERT_(iret==NF_NOERR)
    
    if (present (merge)) then
       open (31,file='clsm/lai.dat.gswp2',  &
            form='unformatted',status='unknown',convert='little_endian')
    else
       if(trim(vname) == 'LAI') open (31,file='clsm/lai.dat',  &
            form='unformatted',status='unknown',convert='little_endian')
       if(trim(vname) == 'grnFrac') open (31,file='clsm/green.dat',  &
            form='unformatted',status='unknown',convert='little_endian')
    endif
    
    do year=82,98
       
       smon=(year-82)*12
       imon=0
       do mon=smon+1,smon+12
          imon=imon+1
	  
	  iret = NF_GET_VARA_REAL(ncid, 6,(/1,mon/),(/MAX_NOF_GRID,1/),vecforc)
          ASSERT_(iret==NF_NOERR)		          
          catforc =1.e-20
          catcount=0
          DO j =1,nr
             DO I = 1,nc
                if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.n_land)) then
                   if ((Raster(i,j).ge.1).and.(Raster(i,j).le.MAX_NOF_GRID)) then
                      catforc(tile_id(i,j)) = catforc(tile_id(i,j)) +   &
                           vecforc(Raster(i,j))  
                      catcount(tile_id(i,j)) = catcount(tile_id(i,j)) + 1. 
                   endif
                endif
             END DO
          END DO
	  
	  do i = 1, n_land
             if(catcount(i).gt.0.) catforc(i) = catforc (i) /catcount(i)
          end do
          mon_climate(:,imon)=mon_climate(:,imon)+catforc(:)/17.
          
       END DO
    END DO ! Year
    iret = NF_CLOSE(ncid)
    ASSERT_(iret==NF_NOERR)

    ncatch = n_land
    
    do K=0,13
       yr = (k+11)/12
       mn = mod(k+11,12)+1
       yr1= (k+12)/12
       mn1= mod(k+12,12)+1
       write(31) float((/yr,mn,1,0,0,0,yr1,mn1,1,0,0,0,ncatch,1/))
       write(31) mon_climate(:,mod(k+11,12)+1)
    end do
    close(31,status='keep')
    
    deallocate(catforc)
    deallocate(mon_climate) 
    deallocate(vecforc)
    deallocate(catcount)  
    deallocate(gswp2_mask)
    if(regrid) then
       deallocate(raster)
    endif 
    
  END SUBROUTINE process_gswp2_veg

  !----------------------------------------------------------------------
  !
  ! SUBROUTINE modis_lai (nx,ny,regrid,gfile)
  !
  ! apparently not used; removed by reichle, 24 Dec 2024 
  !
  !----------------------------------------------------------------------  
  
  SUBROUTINE soil_para_high (nx,ny,regrid, n_land, tile_id, F25Tag)
    
    integer,      intent(in)            :: nx, ny
    logical,      intent(in)            :: regrid
    integer,      intent(in)            :: n_land
    integer,      intent(in)            :: tile_id(:,:)
    logical,      intent (in), optional :: F25Tag

    ! -----------------------------------------------------------
    
    real, dimension(12) :: lbee,lpsis,lporo,lcond,lwpwet, &
         atau2,btau2,atau5,btau5
    REAL, ALLOCATABLE :: soildepth (:)
    INTEGER :: soil_class_top,soil_class_com,soil_gswp,swit
    REAL :: BEE, PSIS, POROS,COND,WPWET
    integer :: n,count,k1,i1,i,j
    character*512 :: path,fname,fout,metpath

    CHARACTER*512 :: version,resoln,continent
    integer :: iret,ncid,ncid1
    real, allocatable, target, dimension (:,:) :: SOIL_HIGH
    REAL, ALLOCATABLE :: count_soil(:)
    integer :: tindex, pfafindex,i_sib,j_sib
    integer :: status
    real, allocatable, dimension(:) :: soildepth_gswp2 
    integer, allocatable, dimension (:) :: land_gswp2

    real, pointer :: Raster(:,:)

    logical                            :: file_exists
    real, allocatable, dimension (:,:) :: parms4file
    
    
    ! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
    !
    ! NOTE: "!$" is for conditional compilation
    !
    logical :: running_omp = .false.
    !
    !$ integer :: omp_get_thread_num, omp_get_num_threads
    !
    integer :: n_threads=1
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
    ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
    
    data lbee /3.30, 3.80, 4.34, 5.25, 3.63, 5.96, 7.32,       &
               8.41, 8.34, 9.70, 10.78, 12.93/
    data lpsis /-0.05, -0.07, -0.16, -0.65, -0.84, -0.24,      &
                -0.12, -0.63, -0.28, -0.12, -0.58, -0.27/
    data lporo /0.373, 0.386, 0.419, 0.476, 0.471, 0.437,      &
                0.412, 0.478, 0.447, 0.415, 0.478, 0.450/
    data lcond /2.45e-05, 1.75e-05, 8.35e-06, 2.36e-06,        &
                1.1e-06, 4.66e-06, 6.31e-06, 1.44e-06,         &
                2.72e-06, 4.25e-06, 1.02e-06, 1.33e-06/
    data lwpwet /0.033,0.051,0.086,0.169,0.045,0.148,0.156,    &
                 0.249,0.211,0.199,0.286,0.276/
    
    data atau2/0.0030065,0.0276075,0.0200614,0.0165152,   &
               0.0165152,0.0168748,0.0308809,0.0329365,    &
               0.0437085,0.0466403,0.0956670,0.1257360/
    
    data btau2/0.0307900,0.0196558,0.0299702,0.0443406,   &
               0.0443406,0.0359961,0.0234851,0.0370919,    &
               0.0312746,0.0249973,0.0222786,0.0193874/
    
    data atau5/0.0067424,0.0766189,0.0540989,0.0439714,   &
               0.0439714,0.0457011,0.0589881,0.0885157,    &
               0.1175960,0.0692305,0.1348880,0.1535540/
    
    data btau5/0.0569718,0.0492634,0.0678898,0.0786387,   &
               0.0786387,0.0737872,0.0713841,0.0742609,    &
               0.0693533,0.0745496,0.0732726,0.0718882/

    i_sib = i_raster
    j_sib = j_raster
   
    allocate(soildepth(n_land))
    allocate(soil_high(1:i_raster,1:j_raster))  
    allocate(count_soil(1:n_land))  
    
    inquire(file='clsm/catch_params.nc4', exist=file_exists)
    
    if(file_exists) then
       status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid) ; VERIFY_(STATUS)
       allocate (parms4file (1:n_land, 1:10))
    endif
    
    soil_high =-9999.
    
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    
    if (present(F25Tag)) then 
       
       iret = NF_OPEN(trim(MAKE_BCS_INPUT_DIR)//'/land/soil/SOIL-DATA/soil_depth/v1/SoilDepth.nc',NF_NOWRITE, ncid1)
       ASSERT_(iret==NF_NOERR)
       allocate (soildepth_gswp2(1: ncat_gswp2))
       allocate (land_gswp2     (1: ncat_gswp2)) 
       iret = NF_GET_VARA_INT (ncid1, 3,(/1/),(/ncat_gswp2/),land_gswp2)
       ASSERT_(iret==NF_NOERR)	
       iret = NF_GET_VARA_REAL(ncid1, 4,(/1/),(/ncat_gswp2/),soildepth_gswp2)
       ASSERT_(iret==NF_NOERR)		          
       iret = NF_CLOSE(ncid1)
       ASSERT_(iret==NF_NOERR)
       
       k1 = i_raster/360
       
       do n = 1,ncat_gswp2
          
          j = (land_gswp2(n)-1)/360  + 1
          i = land_gswp2(n) - (j - 1)*360
          j = 181 - j
          soil_high((i-1)*k1+1:i*k1,(j-1)*k1+1:j*k1) = soildepth_gswp2(n)
          
       end do
       deallocate (soildepth_gswp2,land_gswp2)
    else
       
       open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/soil/SOIL-DATA/soil_depth/v1/soil_depth_2.5.rst',&
            form='unformatted',status='old',action='read',convert='little_endian')
       
       do j =1,j_raster
          read (10) soil_high(:,j)
       end do
       close (10,status='keep')
       
    endif
    
    if(regrid) then
       allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
    else
       raster => soil_high
    end if
    
    if(regrid) then
       call RegridRasterReal(soil_high,raster)
    endif
    
    soildepth =0.
    count_soil = 0.
    
    do j=1,ny
       do i=1,nx
          if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.n_land)) then
             if(raster(i,j).eq.-9999.)   then
                !		   write (*,*)'soil_high UNDEF',i,j,tile_id(i,j),raster(i,j) 
                !                           stop
             endif
             if (raster(i,j).gt.0.) then
                
                soildepth(tile_id(i,j)) = &
                     soildepth(tile_id(i,j)) + raster(i,j)
                count_soil(tile_id(i,j)) = &
                     count_soil(tile_id(i,j)) + 1. 
             endif
          endif
       end do
    end do
    
    DO n =1,n_land
       if(count_soil(n)/=0.) soildepth(n)=soildepth(n)/count_soil(n)	
       if (present(F25Tag)) then
          soildepth(n) = max(soildepth(n),1.)
       else
          soildepth(n) = max(soildepth(n),1.334)
       endif
    END DO
    
    soildepth = soildepth*1000. 
    
    !     Openning files
    
    fname='clsm/soil_text.top'
    open (10,file=fname,status='old',action='read',form='formatted')
    fname='clsm/soil_text.com'
    open (11,file=fname,status='old',action='read',form='formatted')
    fout='clsm/soil_param.first'          
    open (21,file=fout,status='unknown',action='write',form='formatted')
    fout='clsm/tau_param.dat'          
    open (22,file=fout,status='unknown',action='write',form='formatted')
    
    swit =0
    DO n=1 , n_land
       read (10,*) tindex,pfafindex, soil_class_top
       write (22,'(i10,i8,4f10.7)')tindex,pfafindex,atau2(soil_class_top), &
            btau2(soil_class_top),atau5(soil_class_top),btau5(soil_class_top)
       read (11,*) tindex,pfafindex, soil_class_com
       
       !if (soil_class_com.eq.4) then
       !   soil_gswp = 5
       !elseif (soil_class_com.eq.5) then
        !   soil_gswp = 6
       !elseif (soil_class_com.eq.6) then
        !   soil_gswp = 4
       !elseif (soil_class_com.eq.8) then
       !   soil_gswp = 9
       !elseif (soil_class_com.eq.9) then
       !   soil_gswp = 8
       !else
       !   soil_gswp = soil_class_com
       !endif
       
       soil_gswp = soil_class_com
       
       cond=lcond(soil_gswp)/exp(-1.*zks*gnu)
       wpwet=lwpwet(soil_gswp)/lporo(soil_gswp) 
       write (21,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)')tindex,pfafindex,   &
            soil_class_top,soil_class_com,lBEE(soil_gswp), lPSIS(soil_gswp),          &
            lPORO(soil_gswp),COND,WPWET,soildepth(n)
       
       if (allocated (parms4file)) then
          
          parms4file (n, 1) = lBEE(soil_gswp)
          parms4file (n, 2) = COND
          parms4file (n, 3) = lPORO(soil_gswp)
          parms4file (n, 4) = lPSIS(soil_gswp)
          parms4file (n, 5) = WPWET         
          parms4file (n, 6) = soildepth(n)
          parms4file (n, 7) = atau2(soil_class_top)
          parms4file (n, 8) = btau2(soil_class_top)
          parms4file (n, 9) = atau5(soil_class_top)
          parms4file (n,10) = btau5(soil_class_top)
          
       endif
       
    END DO
    close (10,status='delete')
    close (11,status='delete')
    close (21,status='keep')
    close (22,status='keep')
    
    if(file_exists) then
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ) ,(/1/),(/n_land/), parms4file (:, 1)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'COND' ) ,(/1/),(/n_land/), parms4file (:, 2)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'POROS') ,(/1/),(/n_land/), parms4file (:, 3)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ) ,(/1/),(/n_land/), parms4file (:, 4)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'WPWET') ,(/1/),(/n_land/), parms4file (:, 5)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR') ,(/1/),(/n_land/), parms4file (:, 6)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ATAU2') ,(/1/),(/n_land/), parms4file (:, 7)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BTAU2') ,(/1/),(/n_land/), parms4file (:, 8)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ATAU5') ,(/1/),(/n_land/), parms4file (:, 9)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BTAU5') ,(/1/),(/n_land/), parms4file (:,10)) ; VERIFY_(STATUS)
       STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
       DEALLOCATE (parms4file)
    endif
    
    deallocate (soildepth, soil_high)
    if(regrid) then
       deallocate(raster)
    endif
  END SUBROUTINE soil_para_high

  ! ====================================================================
  !
  !  SUBROUTINE remove_tiny_tiles (                         &
  !       dateline,poles,gout) 
  !
  !  ***** subroutine not used as of Dec 2024; removed by rreichle, 20 Dec 2024 *****
  !
  !  END SUBROUTINE remove_tiny_tiles
  !
  ! ---------------------------------------------------------------------
  ! ---------------------------------------------------------------------
  ! ---------------------------------------------------------------------
  
  SUBROUTINE modis_alb_on_tiles (nx,ny,regrid, n_land, tile_id)
    
    integer,      intent(in) :: nx, ny
    logical,      intent(in) :: regrid
    integer,      intent(in) :: n_land
    integer,      intent(in) :: tile_id(:,:)
    
    ! -----------------------------------------------

    CHARACTER*512 :: version,resoln,continent
    character*512 :: path,fname,fout,metpath
    integer :: n,count,k1,i1,i
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3,ip2
    INTEGER :: laiid,year,mon,smon,imon,iret  
    integer :: ialbt,ialbs,yy,j,month
    character*2 :: bw
    character*5 :: cyy
    character*512 :: albtype, albspec
    real, allocatable, target, dimension (:,:) :: alb_in
    real, allocatable, dimension (:) :: alb_count,alb_out
    character*512 :: ifile,ofile
    integer :: status
    real,pointer :: raster (:,:)

    ip2 = ip1 + n_land

    allocate(alb_in(1:i_raster,1:j_raster))  
    allocate(alb_out(1:n_land)) 
    allocate(alb_count(1:n_land)) 
    if(regrid) then
       allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
    else
       raster => alb_in
    end if
    
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    do  ialbt = 2,2
       do ialbs = 1,2
          do yy = 2005,2005
             if (yy.eq.2005)cyy='00-04'
             if(ialbt.eq.1)albtype='BlackSky/'
             if(ialbt.eq.2)albtype='WhiteSky/'
             if(ialbt.eq.1)bw='BS'
             if(ialbt.eq.2)bw='WS'
             
             if(ialbs.eq.1)albspec='0.3_0.7/'
             if(ialbs.eq.2)albspec='0.7_5.0/'           
             ifile=trim(MAKE_BCS_INPUT_DIR)//'/land/albedo/AlbMap.'//bw//'.2x5.'//trim(cyy)//        &
                  '.monthly.'//albspec(1:index(albspec,'/')-1)//'.dat' 
             ofile='clsm/AlbMap.'//bw//'.2x5.'//trim(cyy)//'.monthly-tile.' &
                  //albspec(1:index(albspec,'/')-1)//'.dat'

             open (20,file=trim(ifile),form='unformatted',&
                  convert='big_endian', &
                  action='read',status='old')
             open (30,file=trim(ofile),form='unformatted', &
                  convert='big_endian', &
                  action='write',status='unknown')
             
             do month =1,12
                read (20) alb_in
                if(regrid) then
                   call RegridRasterReal(alb_in,raster)
                else
                   raster = alb_in
                endif
                
                alb_out = 0.
                alb_count = 0.
                do j=1,ny
                   do i=1,nx
                      if((tile_id(i,j).gt.ip1).and.(tile_id(i,j).le.ip2)) then
                         if(raster(i,j).eq.undef)   then
                            !                           write (*,*)'raster UNDEF',i,j,month,albtype,albspec
                            !                           stop
                         endif
                         if ((raster(i,j).gt.0.).and.(raster(i,j).ne.undef)) then
                            alb_out(tile_id(i,j)-ip1) = &
                                 alb_out(tile_id(i,j)-ip1) + raster(i,j)
                            alb_count(tile_id(i,j)-ip1) = &
                                 alb_count(tile_id(i,j)-ip1) + 1. 
                         endif
                      endif
                   end do
                end do
                
                do n = 1,n_land
                   if (alb_count(n).gt.0)then
                      alb_out(n) = alb_out(n)/alb_count(n)
                   else
                      !                     print *,'No albedo for the tile :',n
                      alb_out(n) = alb_out(n-1)
                   endif
                end do
                write (30) alb_out
             end do
             close (20,status='keep')
             close (30,status='keep')           
             
          end do
       end do
    end do
    
    deallocate (alb_in,alb_out,alb_count)
    if(regrid) then
       deallocate(raster)
    endif
  END SUBROUTINE modis_alb_on_tiles
  
  !----------------------------------------------------------------------
  !
  ! The following subroutines were already commented out as of 24 Dec 2024.
  ! Removed by rreichle, 24 Dec 2024.
  !
  ! SUBROUTINE modis_scale_para (ease_grid,gfile)
  !
  ! SUBROUTINE make_75 (nx,ny,regrid,path,gfile)
  ! 
  ! subroutine pick_cat(sam,clr)
  !
  !----------------------------------------------------------------------
  
  SUBROUTINE supplemental_tile_attributes(nx,ny,regrid,dateline,fnameTil, Rst_id)
    
    ! 1) get supplemental tile attributes not provided in MAPL-generated (ASCII) tile file,
    !    incl. min/max lat/lon of each tile and tile elevation
    ! 2) write nc4-formatted til file (incl. supplemental tile attributes)
    
    integer,      intent(in) :: nx, ny

    logical,      intent(in) :: regrid    

    character(*), intent(in) :: dateline

    character(*), intent(in) :: fnameTil   ! file name (w/o extension) of tile file
    integer, intent(in)      :: Rst_id(:,:)   
 
    ! ---------------------------------------------------------

    INTEGER, allocatable, dimension(:) :: CATID  
    integer                            :: n, ip, n_land, i, j, i_sib, j_sib, status
    INTEGER, allocatable, dimension(:) :: id, I_INDEX, J_INDEX 
    integer                            :: nc_gcm, nr_gcm, nc_ocean, nr_ocean
    REAL                               :: lat, lon, fr_gcm, fr_cat, tarea
    INTEGER                            :: typ, pfs, ig, jg, j_dum, i_dum, ierr, indx_dum, ip2, n_grid

    REAL (REAL64), PARAMETER           :: RADIUS=MAPL_RADIUS, pi= MAPL_PI 

    character*512                      :: fname
    character*512                      :: gtopo30
    CHARACTER*512                      :: version

    REAL, allocatable                  :: limits(:,:)

    REAL :: mnx,mxx,mny,mxy,dx,dy,d2r,lats,sum1,dx_gcm,area_rst

    REAL,    allocatable, dimension(:) :: tile_ele, tile_area,tile_area_rst  
    integer                            :: IM(2), JM(2)

    real,    pointer                   :: Raster(:,:)
    real                               :: mean_land_elev

    real*4,            allocatable, target :: q0 (:,:)
    real(REAL64),      allocatable         :: rTable(:,:)
    integer,           allocatable         :: iTable(:,:)
    character(len=128)                     :: gName(2)
    logical,           allocatable         :: IsOcean(:)

    ! -----------------------------------------------------
    !
    ! get elevation (q0) from "gtopo30" raster file ("srtm30_withKMS_2.5x2.5min.data")
    
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    gtopo30   = trim(MAKE_BCS_INPUT_DIR)//'/land/topo/v1/srtm30_withKMS_2.5x2.5min.data'
    allocate (q0(1:i_raster,1:j_raster))

    i_sib = nx
    j_sib = ny
    
    dx  = 360._8/i_sib
    dy  = 180._8/j_sib
    d2r = PI/180._8

    open (10,file=trim(gtopo30),form='unformatted',status='old')
    read (10) q0
    close (10,status='keep')
    
    if(regrid) then
       allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
    else
       raster => q0
    end if
    
    if(regrid) then
       call RegridRasterReal(q0,raster)
    endif

    ! -----------------------------------------------------------
    !
    ! read ASCII-formatted tile file (*.til)
    !
    ! ip  = number of tiles in global domain (all types, incl. land, landice, lake, & ocean)
    ! ip1 = index offset for land tiles in *.til files  (ip1=0 implies that land tiles first in *.til file)
    ! ip2 = ip1 + n_land = end index of land tiles (where n_land is number of land tiles in global *.til file)
    
    allocate (catid(1:i_sib))
    catid=0
    fname=trim(fnameTil)//'.til'

    open (10,file=fname,status='old',action='read',form='formatted')
    read (10,*)ip

    allocate(id(       ip))
    allocate(i_index(  ip))
    allocate(j_index(  ip))
    allocate(tile_area(ip))
  
    id=0
    read (10,*) n_grid
    IM = 0
    JM = 0
    gName = ['','']
    do n = 1, n_grid
       read (10,'(a)')version
       read (10,*)nc_gcm
       read (10,*)nr_gcm
       gName(n) = trim(adjustl(version))
       IM(n)    = nc_gcm
       JM(n)    = nr_gcm
    end do
    
!    dx_gcm = 360./float(nc_gcm)
    
    allocate(iTable(ip,0:7))
    allocate(rTable(ip,10))    
    rTable = MAPL_UNDEF_r8
    
    allocate(IsOcean(ip))
    IsOcean = .false.
    
    do n = 1,ip
       
       read(10,'(I10,3E20.12,9(2I10,E20.12,I10))',IOSTAT=ierr)     &    
            typ,tarea,lon,lat,ig,jg,fr_gcm,indx_dum,pfs,i_dum,fr_cat,j_dum
       
       if(ierr /= 0) write (*,*)'Problem reading ' // trim(fname)
       
       tile_area(n) = tarea
       id(n)        = pfs
       i_index(n)   = ig 
       j_index(n)   = jg 

       if (typ == 100) ip2 = n
       if (typ == 0  ) IsOcean(n) = .true.

       iTable(n,0) = typ
       rTable(n,3) = tarea
       rTable(n,1) = lon
       rTable(n,2) = lat
       iTable(n,2) = ig
       iTable(n,3) = jg
       rTable(n,4) = fr_gcm
       iTable(n,6) = indx_dum
       iTable(n,4) = pfs
       iTable(n,5) = i_dum
       rTable(n,5) = fr_cat
       iTable(n,7) = j_dum           
    end do
    close (10,status='keep')
    
    n_land=ip2-ip1                 ! = number of land tiles
    
    ! ---------------------------------------------------------------
    !
    ! compute supplemental tile info: mean elevation and min/max lat/lon of each tile
    
    allocate(tile_ele(     1:ip))
    allocate(tile_area_rst(1:ip))

    tile_ele       =    0.
    tile_area_rst  =    0.   ! total area of raster grid cells contributing to each tile
    
    allocate(limits(       1:ip,1:4))

    limits(:,1)    =  360.
    limits(:,2)    = -360.
    limits(:,3)    =   90.
    limits(:,4)    =  -90.

    ! read raster file with tile IDs

    do j=1,j_sib
       
       ! latitude and area of raster grid cells associated with lat index j
       
       lats     = -90._8 + (j - 0.5_8)*dy
        
       ! preserve zero-diff
       !area_rst = (sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)

       ! read tile IDs for lat index j
       catid(:) = rst_id(:,j) 

       ! compute average elevation weighted by area of contributing raster grid cells
       
       do i=1,i_sib          
          if (.not. IsOcean(catid(i)-ip1)) then
             
             tile_ele(     catid(i)-ip1) = tile_ele(     catid(i)-ip1) + raster(i,j)* &
                                           (sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)
             
             tile_area_rst(catid(i)-ip1) = tile_area_rst(catid(i)-ip1) +             &
                                            (sin(d2r*(lats+0.5*dy)) -sin(d2r*(lats-0.5*dy)))*(dx*d2r)
             
          endif
       enddo
      
       mny=-90. + float(j-1)*180./float(j_sib)
       mxy=-90. + float(j)  *180./float(j_sib)
       
       if (index(dateline,'DE')/=0) then
          do i=1,i_sib          
             if( .not. IsOcean(catid(i)- ip1))then
                mnx =-180. + float(i-1)*360./float(i_sib)
                mxx =-180. + float(i)  *360./float(i_sib)
                if(mnx .lt.limits(catid(i)-ip1,1))limits(catid(i)-ip1,1)=mnx 
                if(mxx .gt.limits(catid(i)-ip1,2))limits(catid(i)-ip1,2)=mxx 
                if(mny .lt.limits(catid(i)-ip1,3))limits(catid(i)-ip1,3)=mny 
                if(mxy .gt.limits(catid(i)-ip1,4))limits(catid(i)-ip1,4)=mxy 
             endif
          end do
       else
          do i=1,i_sib- i_sib/nc_gcm/2         
             if( .not. IsOcean(catid(i) - ip1)) then
                mnx =-180. + float(i-1)*360./float(i_sib)
                mxx =-180. + float(i)  *360./float(i_sib)
                if(mnx .lt.limits(catid(i)-ip1,1))limits(catid(i)-ip1,1)=mnx 
                if(mxx .gt.limits(catid(i)-ip1,2))limits(catid(i)-ip1,2)=mxx 
                if(mny .lt.limits(catid(i)-ip1,3))limits(catid(i)-ip1,3)=mny 
                if(mxy .gt.limits(catid(i)-ip1,4))limits(catid(i)-ip1,4)=mxy 
             endif
          end do
          do i=i_sib- i_sib/nc_gcm/2  +1,i_sib       
             if( .not. IsOcean(catid(i) - ip1)) then
                mnx =-360. -180. + float(i-1)*360./float(i_sib)
                mxx =-360. -180. + float(i)  *360./float(i_sib)
                if(mnx <  -180.) mnx = mnx + 360.
                if(mxx <= -180.) mxx = mxx + 360.
                if(mnx .lt.limits(catid(i)-ip1,1))limits(catid(i)-ip1,1)=mnx 
                if(mxx .gt.limits(catid(i)-ip1,2))limits(catid(i)-ip1,2)=mxx 
                if(mny .lt.limits(catid(i)-ip1,3))limits(catid(i)-ip1,3)=mny 
                if(mxy .gt.limits(catid(i)-ip1,4))limits(catid(i)-ip1,4)=mxy 
             endif
          end do
       endif
    enddo
    
    ! finalize min/max lat/lon
    
    where (limits(:,1).lt.-180.) limits(:,1) = limits(:,1) + 360.0
    where (limits(:,2).le.-180.) limits(:,2) = limits(:,2) + 360.0
    
    ! finalize elevation
    
    where ( .not. IsOcean)  tile_ele = tile_ele/tile_area_rst
    
    ! adjust global mean (land) topography to 614.649 (615.662 GTOPO 30) m
    
    sum1=0.
    
    do j=1,n_land
       sum1 = sum1 + tile_ele(j)*tile_area(j)
    enddo
    
    mean_land_elev = sum1/sum(tile_area(1:n_land))
   
    if ( mean_land_elev .ne. Target_mean_land_elev ) then
       
       print *, 'Global mean land elevation before adjustment     [m]: ', mean_land_elev
       
       tile_ele(1:n_land) = tile_ele(1:n_land)*(Target_mean_land_elev / mean_land_elev) 
       
       ! verify adjustment
       
       sum1=0.
       
       do j=1,n_land
          sum1 = sum1 + tile_ele(j)*tile_area(j)
       enddo
       
       print *, 'Global mean land elevation after scaling to SRTM [m]: ', sum1/sum(tile_area(1:n_land))
       
    endif
    
    ! --------------------------------------------------------------------------
    !
    ! write (ASCII) catchment.def file (land tiles only!)
    
    open (10,file='clsm//catchment.def',  &
         form='formatted',status='unknown')
    write (10,*) n_land
    
    do j=1,n_land
 !      if(trim(dateline)=='DC')then
 !         limits(j,1) = max(limits(j,1),(i_index(j)-1)*dx_gcm -180. - dx_gcm/2.)       
 !         limits(j,2) = min(limits(j,2),(i_index(j)-1)*dx_gcm -180. + dx_gcm/2.)  
 !      endif
       write (10,'(i10,i8,5(2x,f9.4))')j+ip1,id(j+ip1),limits(j,1),   &
            limits(j,2),limits(j,3),limits(j,4),tile_ele(j)       
    end do
    close(10,status='keep')
    
    ! --------------------------------------------------------------------------
    !
    ! write nc4-formatted tile file (all tile types)
    
    rTable(1:ip,6:9) = limits
    rTable(1:ip, 10) = tile_ele(1:ip)
    ! re-define rTable(:,4) and rTable(:,5).
    ! fr will be re-created in WriteTilingNC4
    where (rTable(:,4) /=0.0)
       rTable(:,4) = rTable(:,3)/rTable(:,4)
    endwhere
    where (rTable(:,5) /=0.0)
       rTable(:,5) = rTable(:,3)/rTable(:,5)
    endwhere
    
    fname=trim(fnameTil)//'.nc4'
    call MAPL_WriteTilingNC4(fname,  gName(1:n_grid), im(1:n_grid), jm(1:n_grid), nx, ny, iTable, rTable, N_PfafCat=SRTM_maxcat, rc=status)
    
    deallocate (rTable, iTable)
    deallocate (limits)
    deallocate (catid)
    deallocate (q0)
    if(regrid) then
       deallocate(raster)
    endif
    
  END SUBROUTINE supplemental_tile_attributes
 
  !----------------------------------------------------------------------
  
  SUBROUTINE create_soil_types_files( nx, ny, n_land, tile_pfs, catid )
    
    integer,        intent(in) :: nx, ny
    integer,        intent(in) :: n_land    
    integer,        intent(in) :: tile_pfs(:)    
    INTEGER, target,intent(in) :: CATID(:,:)

    
    !     This program reads global 5'x5' soil texture classification,
    !     then find the dominant Soil Classes for the GCM
    !     http://www.ngdc.noaa.gov/seg/eco/cdroms/reynolds/reynolds/reynolds.htm
    !     http://www.ngdc.noaa.gov/ecosys/cdroms/reynolds/reynolds/reynolds.htm
    !     Soil texture classification follows USDA classification (13 classes)
    !     min. value  : 0
    !     max. value  : 12
    !     legend cats : 13
    !     category 0  : Ocean/No Data
    !     category 1  : Sand
    !     category 2  : Loamy Sand
    !     category 3  : Sandy Loam
    !     category 4  : Silt Loam
    !     category 5  : Silt
    !     category 6  : Loam
    !     category 7  : Sandy Clay Loam
    !     category 8  : Silty Clay Loam
    !     category 9  : Clay Loam
    !     category 10 : Sandy Clay
    !     category 11 : Silty Clay
    !     category 12 : Clay
    !     

    INTEGER col,row,i,j,k,ii,n,ip
    INTEGER colsib,rowsib,isol,jsol,cls,clr1(1),clr2(1)
    PARAMETER(col=4320,row=2160)
    
    INTEGER, allocatable :: SIB_LAY(:,:)
    INTEGER, allocatable :: SOIL1(:,:)
    INTEGER, allocatable :: SOIL2(:,:)
    INTEGER tem1 (13),tem2(13),tem3(13)
    INTEGER, ALLOCATABLE :: TOP(:,:),COM(:,:)
    INTEGER IDVAL,STEX
    INTEGER (kind=1), allocatable :: gtext(:,:)
    INTEGER irrecs, c1,c2,r1,r2
    CHARACTER*512 ifile,ifile2,ofile1,ofile2,fname
    CHARACTER*512 :: version,resoln    
    INTEGER, allocatable, dimension (:) :: id !indx,id,indx_old
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: typ,pfs,ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3 ,ip2
    integer :: status

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
    ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------
    
    colsib = nx
    rowsib = ny
    !
    !	Compute the number of input records per row.
    irrecs = nint (col / 4.0)      
    !
    
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    

    ifile=trim(MAKE_BCS_INPUT_DIR)//'/land/soil/SOIL-DATA/soil_properties/v1/'//'dtex_tp1.bin'
    ifile2=trim(MAKE_BCS_INPUT_DIR)//'/land/soil/SOIL-DATA/soil_properties/v1/'//'dtex_sb1.bin'
    ofile1='clsm/soil_text.top'
    ofile2='clsm/soil_text.com'     

    ip = size(tile_pfs,1)
    allocate(id(1:ip), source = tile_pfs)

    ip2 = ip1 + n_land 

    ! write(*,*)'Finished reading CAT_IDs'
    
    !     Top layer soil classification 0-30cm
    !
    open (unit=11, file=ifile, form='unformatted', status='old',  &
         convert = 'big_endian')
    
    !
    allocate(gtext(1:col,1:row))
    allocate(SIB_LAY(1:nx,1:ny))
    gtext(:,:)=0
    SIB_LAY(:,:)=0
    k=0
    do j=row,1,-1
       !           do i=1,irrecs
       !              k=k+1
       !              c1 = (4*i)-3
       !              c2 = (4*i)
       !              read (unit=11, rec=k) (gtext(ii,j), ii=c1,c2)
       read (unit=11) (gtext(i,j), i=1,col)
       !           end do
    end do
    
    close (11,status='keep')
    !
    do j=1,rowsib
       jsol=CEILING(j/(ny/real(row)))
       do i=1,colsib
          isol=CEILING(i/(nx/real(col)))
          SIB_LAY(i,j)=gtext(isol,jsol)
       end do
    end do
    
    deallocate(gtext)
    !
    !     Top layer on 2x2.5
    allocate(soil1(n_land,1:13))
    soil1(:,:)=0
    do j=1,rowsib
       do i=1,colsib
          if((catid(i,j) > ip1).and.(catid(i,j) <= ip2))then
             IDVAL=catid(i,j) 
             STEX=SIB_LAY(i,j)
             SOIL1(IDVAL,STEX+1)=SOIL1(IDVAL,STEX+1)+1
          end if
       end do
    end do
    !
    !  write(*,*)'Finished reading top layer'
    deallocate(sib_lay)
    !
    !     Bottom layer soil classification 30-100cm
    !
    ! 	open (unit=11, file=ifile2, form='unformatted', status='old',access='direct',recl=1,  &
    !             convert = 'big_endian')
    open (unit=11, file=ifile2, form='unformatted', status='old',  &
         convert = 'big_endian')
    
    !
    allocate(gtext(1:col,1:row))
    allocate(SIB_LAY(1:colsib,1:rowsib))
    gtext(:,:)=0
    SIB_LAY(:,:)=0
    k=0
    do j=row,1,-1
       !           do i=1,irrecs
       !              k=k+1
       !              c1 = (4*i)-3
       !              c2 = (4*i)
       !              read (unit=11, rec=k) (gtext(ii,j), ii=c1,c2)
       read (unit=11) (gtext(i,j), i=1,col)
       !           end do
    end do
    !
    close (11,status='keep')
    
    do j=1,rowsib
       jsol=CEILING(j/(ny/real(row)))
       do i=1,colsib
          isol=CEILING(i/(nx/real(col)))
          SIB_LAY(i,j)=gtext(isol,jsol)
       end do
    end do
    deallocate(gtext)
    ! write(*,*)'Finished reading bottom layer'
    !
    !     Bottom layer on 2x2.5
    allocate(soil2(n_land,1:13))
    soil2(:,:)=0
    do j=1,rowsib
       do i=1,colsib
          if((catid(i,j) > ip1).and.(catid(i,j) <= ip2))then
             IDVAL=catid(i,j)  
             
             STEX=SIB_LAY(i,j)
             SOIL2(IDVAL,STEX+1)=SOIL2(IDVAL,STEX+1)+1
          endif
       end do
    end do
    deallocate(sib_lay)
    !
    !         write(*,*)'Finished counting pixels for each catchment'
    k=0
    allocate(top(n_land,2))
    allocate(com(n_land,2))
    top=0
    com=0
    do j=1,n_land
       tem1(1:13)=SOIL1(j,1:13)
       tem2(1:13)=SOIL2(j,1:13)
       
       tem3(:)=3*tem1(:)+7*tem2(:)
       if((sum(tem3).gt.0).and.(sum(tem1).eq.0))then
          tem1(:)=tem3(:)
          write(*,*)'Filled from the bottom layer',j
       end if
       if(sum(tem1).gt.0)then
          !              k=k+1
          !              ! 
          !             clr1=maxloc(tem1)
          !             clr2=maxloc(tem3)
          !             top(k,1)=j 
          !             top(k,2)=clr1(1)-1
          !             com(k,1)=j 
          !             com(k,2)=clr2(1)-1             
          k=k+1
          ! 
          clr1=maxloc(tem1)
          clr2=maxloc(tem3)
          top(j,1)=j 
          top(j,2)=clr1(1)-1
          com(j,1)=j 
          com(j,2)=clr2(1)-1      
       end if
    end do
    !
    open (unit=11, file=ofile1, form='formatted', status='unknown')
    open (unit=12, file=ofile2, form='formatted', status='unknown')
    
    !
    if(top(1,2).eq.0)top(1,2)= 3
    if(com(1,2).eq.0)com(1,2)= 9
    
    do j=1,n_land
       
       if(top(j,2).eq.0)top(j,2)=top(j-1,2)
       if(com(j,2).eq.0)com(j,2)=com(j-1,2)
       
       !         if(com(j,1).gt.0)then
       !            if(j.gt.1)then
       !               if(top(j,2).eq.0)top(j,2)=top(j-1,2)
       !               if(com(j,2).eq.0)com(j,2)=com(j-1,2)
       !            end if
       !
       write(11,*)j,id(j),top(j,2)
       write(12,*)j,id(j),com(j,2)
       
    end do
    close(11)
    close(12) 
    deallocate (soil1,soil2,top,com,id)
    
  END SUBROUTINE create_soil_types_files
  
  !----------------------------------------------------------------------
  
  SUBROUTINE compute_mosaic_veg_types( nx, ny, regrid, n_land, tile_pfs, Rst_id)
    
    integer,      intent(in) :: nx, ny
    
    logical,      intent(in) :: regrid
    
    integer,      intent(in) :: n_land
    integer,      intent(in) :: tile_pfs(:)
    integer,      intent(in) :: Rst_id(:,:)   
 
    ! -----------------------------
    
    integer*1, allocatable , dimension (:,:)  :: sib_veg2
    integer,   allocatable , target , dimension (:,:)  :: sib_veg   
    integer, allocatable ::  mos_veg(:,:)
    real, allocatable :: veg_frac(:,:)
    integer :: j,k,mcls
    INTEGER CATID(nx),cls,cls2,bcls
    REAL, allocatable :: veg(:,:),bare_frac(:),zdep2_g(:,:)
    REAL :: fmax0,dummy,tem(6),mfrac,sfrac,bfrac
    
    integer :: n,ip,count,k1,i1,i
    INTEGER, allocatable, dimension (:) :: id ! indx,id,indx_old
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    INTEGER :: ig,jg,i_dum,j_dum,ierr,indx_dum,indr1,indr2,indr3,ip2
    character*512 :: fname,fout
    CHARACTER*512 :: version,resoln,continent
    character*2 :: chyear
    integer :: mon,smon,imon,year
    integer :: status
    integer, pointer :: Raster(:,:)
    real, pointer, dimension (:)  :: z2, z0
    real, dimension (6) :: VGZ2 = (/35.0, 20.0, 17.0, 0.6, 0.5, 0.6/) ! Dorman and Sellers (1989)
    logical                            :: file_exists
    integer                            :: ncid

    ip = size(tile_pfs,1)

    ip2 = ip1 + n_land

    allocate(id(1:ip), source = tile_pfs)
    
    allocate(sib_veg2(1:i_raster,1:j_raster))
    allocate(sib_veg (1:i_raster,1:j_raster))
    
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/veg/pft/v1/sib22.5_v2.0.dat',form='unformatted',      &
         status='old',action='read',convert='big_endian')
    READ(10)sib_veg2
    
    close (10,status='keep')
    sib_veg = sib_veg2
    if(regrid) then
       allocate(raster(nx,ny),stat=STATUS); VERIFY_(STATUS)
    else
       raster => sib_veg
    end if
    
    if(regrid) then
       call RegridRaster(sib_veg,raster)
    endif
    
    allocate(veg(1:n_land,1:6))
    allocate(zdep2_g(1:n_land,1:1))
    
    veg=0.
    zdep2_g=0.
    
    n=1

    do j=1,ny
       
       catid(:) = Rst_id(:,j)
       
       do i=1,nx
          
          if((catid(i) > ip1).and.(catid(i) <= ip2))then
             zdep2_g(catid(i)-ip1,1)=zdep2_g(catid(i)-ip1,1)+1.
             if(raster(i,j).eq.0) then
                !                write (*,*)'Warning : SiB2 =0, an ocean pixel found !'
             elseif (raster(i,j).eq.1) then
                veg(catid(i)-ip1,1)=veg(catid(i)-ip1,1) + 1.
             elseif (raster(i,j).eq.2) then
                veg(catid(i)-ip1,2)=veg(catid(i)-ip1,2) + 1.
             elseif (raster(i,j).eq.3) then
                veg(catid(i)-ip1,2)=veg(catid(i)-ip1,2) + 0.5
                veg(catid(i)-ip1,3)=veg(catid(i)-ip1,3) + 0.5                  
             elseif (raster(i,j).eq.4) then
                veg(catid(i)-ip1,3)=veg(catid(i)-ip1,3) + 1.
             elseif (raster(i,j).eq.5) then
                veg(catid(i)-ip1,3)=veg(catid(i)-ip1,3) + 1.
             elseif (raster(i,j).eq.6) then
                veg(catid(i)-ip1,4)=veg(catid(i)-ip1,4) + 1.
             elseif (raster(i,j).eq.7) then
                veg(catid(i)-ip1,5)=veg(catid(i)-ip1,5) + 1.
             elseif (raster(i,j).eq.8) then
                !               if (j >= NINT(float(ny)*(140./180.))) then 
                !                   veg(catid(i)-ip1,6)=veg(catid(i)-ip1,6) + 1.
                !                else
                !                   veg(catid(i)-ip1,5)=veg(catid(i)-ip1,5) + 1.   
                !                endif
                if ((j > NINT(float(ny)*(40./180.))).and.(j < NINT(float(ny)*(140./180.)))) then 
                   veg(catid(i)-ip1,5)=veg(catid(i)-ip1,5) + 1.
                else
                   veg(catid(i)-ip1,6)=veg(catid(i)-ip1,6) + 1.   
                endif
             elseif (raster(i,j).eq.9) then
                veg(catid(i)-ip1,4)=veg(catid(i)-ip1,4) + 1.
             elseif (raster(i,j).eq.10) then
                !                write (*,*)'Warning : SiB2 =10, a water pixel found !'
             elseif (raster(i,j).eq.11) then
                !                write (*,*)'Warning : SiB2 =11, an ice pixel found !'
             elseif (raster(i,j).eq.100) then
                !                write (*,*)'Warning : SiB2 =100, NODATA pixel found !'
             endif
          endif
       enddo
    enddo
    
    allocate(mos_veg(1:n_land,1:2))
    allocate(veg_frac(1:n_land,1:3))
    mos_veg=0
    veg_frac=0.
    
    k=0
    do j=1,n_land
       tem(1:6)=veg(j,1:6)
       
       if(sum(tem).le.0.)write(*,*) 'Warning no veg types',j
       !       if(sum(tem).le.0.) stop
       if(sum(tem).gt.0)then
          
          k=k+1
          
          mfrac=-10.
          sfrac=-10.
          bfrac=0.
          cls=100
          cls2=100
          do n=1,6
             
             if(mfrac.le.tem(n))then
                sfrac=mfrac
                cls2=cls
                mfrac=tem(n)
                cls=n
             elseif(sfrac.le.tem(n)) then
                if(tem(n).lt.mfrac)then
                   sfrac=tem(n)
                   cls2=n
                endif
             endif
          enddo
          
          mos_veg(k,1)=cls
          mos_veg(k,2)=cls2
          veg_frac(k,1)=mfrac/zdep2_g(k,1)
          veg_frac(k,2)=sfrac/zdep2_g(k,1)
          veg_frac(k,3)=1.-mfrac-sfrac
       else
          k = k + 1
          mos_veg(k,1) = mos_veg(k-1,1)
          mos_veg(k,2) = mos_veg(k-1,2)
          veg_frac(k,1)= veg_frac(k-1,1)
          veg_frac(k,2)= veg_frac(k-1,2)
          veg_frac(k,3)= veg_frac(k-1,3)
       endif
       if(veg_frac(k,1).eq.0.)then
          write(*,*)'Checking for 100% desert soil tiles',k,mos_veg(k,1)&
               ,mos_veg(k,2),veg_frac(k,1),veg_frac(k,2),veg_frac(k,3)
          if(veg_frac(k,3).ne.1.)then
             write(*,*)'it is not 100% desert soil either'
          endif
          
       endif
       if(mos_veg(k,1).eq.100) then
          write(*,*)'Checking for 100% desert soil tiles',k,mos_veg(k,1)&
               ,mos_veg(k,2),veg_frac(k,1),veg_frac(k,2),veg_frac(k,3)
          write(*,*) 'Prob1'
       endif
       if(mos_veg(k,2).eq.100) then
          write(*,*) 'Prob1'
          mos_veg(k,2)=7
          veg_frac(k,2)=0.
          write(*,*)k,tem
          write(*,*)mos_veg(j,1),mos_veg(j,2),veg_frac(j,1),veg_frac(j,2),veg_frac(j,3)
       endif
    end do
    deallocate(veg)
    
    ! Canopy height and ASCAT roughness length
    
    call ascat_r0 (nx,ny, n_land, Rst_id, z0)
    
    if(jpl_height) then
       call jpl_canoph (nx,ny, n_land, Rst_id, z2)
    else
       allocate (z2(1:n_land))       
    endif
    
    open (10,file='clsm/mosaic_veg_typs_fracs',  &
         form='formatted',status='unknown')
    do j=1,n_land
       if (mos_veg(j,1) == 0) then
          if(.not.jpl_height) z2(j) = VGZ2(mos_veg(j,1))
          mos_veg(j,1) = mos_veg(j-1,1)
          mos_veg(j,2) = mos_veg(j-1,2)
          write (10,'(i10,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')     &
               j+ip1,id(j+ip1),mos_veg(j-1,1),mos_veg(j-1,2),veg_frac(j,1),veg_frac(j,2),z2(j), z0 (j)           
       else
          if(.not.jpl_height) z2(j) = VGZ2(mos_veg(j,1))
          write (10,'(i10,i8,2(2x,i3),2(2x,f6.2),2x,f6.3,2x,f10.7)')     &
               j+ip1,id(j+ip1),mos_veg(j,1),mos_veg(j,2),veg_frac(j,1),veg_frac(j,2),z2(j), z0 (j)
       endif
    end do
    close(10,status='keep')
    
    inquire(file='clsm/catch_params.nc4', exist=file_exists)
    
    if(file_exists) then
       status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid                                  ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'OLD_ITY'),(/1/),(/n_land/),real(mos_veg(:,1))) ; VERIFY_(STATUS)
       STATUS = NF_CLOSE (NCID) ; VERIFY_(STATUS)
    endif
    
    inquire(file='clsm/vegdyn.data', exist=file_exists)
    
    if(file_exists) then
       status = NF_OPEN ('clsm/vegdyn.data', NF_WRITE, ncid                                       ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ITY'    ),(/1/),(/n_land/),real(mos_veg(:,1))) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'Z2CH'   ),(/1/),(/n_land/),z2                ) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ASCATZ0'),(/1/),(/n_land/),Z0                ) ; VERIFY_(STATUS)
       STATUS = NF_CLOSE (NCID) ; VERIFY_(STATUS)
    else
       open (20,file='clsm/vegdyn.data',status='unknown',action='write',form='unformatted', &
            convert='little_endian')      
       
       write (20) real(mos_veg(:,1))
       write (20) z2 (:)
       write (20) z0 (:)    
       close (20)   
    endif
    
    deallocate (sib_veg2,sib_veg,mos_veg,veg_frac,zdep2_g,id,  z0, z2)
    if(regrid) then
       deallocate(raster)
    endif
    
  END SUBROUTINE compute_mosaic_veg_types
  
  !----------------------------------------------------------------------
  
  SUBROUTINE cti_stat_file ( MaskFile, n_land, tile_pfs, til_j_dum)
    character(*), intent(in) :: MaskFile
    integer, intent(in)      :: n_land, tile_pfs(:), til_j_dum(:)

    ! ----------------------------------------------

    INTEGER, PARAMETER :: nbcat=36716,nofvar=6
    INTEGER :: n,i,ip, itext(SRTM_maxcat,2),ix, jx,ip2
    INTEGER :: pfs, ig,jg,j_dum,ierr,indx_dum,indr1,indr2,indr3
    INTEGER*8 :: idum8
    INTEGER :: ncat,i_dum
    INTEGER, dimension(:), allocatable :: colin2cat 
    INTEGER, allocatable, dimension (:) :: id,indx_old
    integer :: idum1,idum2
    real :: dum1,dum2,dum3,dum4,dum5,dum6
    integer :: nc_gcm,nr_gcm,nc_ocean,nr_ocean
    REAL :: lat,lon,fr_gcm,fr_cat,tarea
    REAL :: fr
    REAL, allocatable, dimension (:,:) :: var
    REAL, allocatable, dimension (:) :: dummy
    CHARACTER*512 :: version
    character*512 :: fname

    ip = size(tile_pfs,1)
    allocate(indx_old(ip))
    allocate(id(ip))

    ip2      = ip1 + n_land
    id       = tile_pfs
    indx_old = til_j_dum
    if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) indx_old = tile_pfs

    allocate(colin2cat(1:6000000))
    colin2cat=0
    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
       open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/topo/v1/SRTM-TopoData/Pfafcatch-routing.dat',   &
            form='formatted', status='old',action='read')       
    else
       open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/misc/old_land/catchment.def',   &
            form='formatted', status='old',action='read')
    endif

    read (10,*) ncat
    do n=1,ncat
       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
          read (10,*)indx_dum,idum8
       else
          read (10,*)j_dum,indx_dum
       endif
       colin2cat(indx_dum)=n
    end do
    close (10,status='keep')

    if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
       open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/topo/v1/SRTM-TopoData/SRTM_cti_stats.dat',       &
            form='formatted', status='old',action='read')
    else
       open (10,file=trim(MAKE_BCS_INPUT_DIR)//'/land/misc/old_land/cti_stats.dat',       &
            form='formatted', status='old',action='read')
    endif

    fname='clsm/cti_stats.dat'
    open (20,file=fname,form='formatted', status='unknown')
    write (20,*) n_land

    read (10,*)ncat
    allocate(var(1:ncat,1:nofvar))
    var=0.
    pfs=1
    do i=1,ncat
       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
          read(10,'(i8,i15,5(1x,f8.4),i5,e18.3)')n,idum8,dum1,dum2,dum3,dum4,dum5,idum2,dum6
          idum1 = n
       else
          read(10,'(i8,i8,5(1x,f8.4),i5,e18.3)')n,idum1,dum1,dum2,dum3,dum4,dum5,idum2,dum6
       endif
       if(colin2cat(idum1).gt.0)then
          itext(pfs,1)=idum1
          var(pfs,1)=dum1
          var(pfs,2)=dum2
          var(pfs,3)=dum3
          var(pfs,4)=dum4
          var(pfs,5)=dum5
          itext(pfs,2)=idum2
          var(pfs,6)=dum6
          pfs=pfs+1
       endif
    end do

    close (10,status='keep')
    !
    do i=ip1+1,ip2
       if(((id(i).ge.5000142).and.(id(i).le.5025829)))then
          write(20,'(i10,i8,5(1x,f8.4),i5,e18.3)')i,id(i),var(indx_old(i),1)*11.1/9.1,var(indx_old(i),2), &
               var(indx_old(i),3),var(indx_old(i),4),var(indx_old(i),5),itext(indx_old(i),2),            &
               var(indx_old(i),6)
       else

          write(20,'(i10,i8,5(1x,f8.4),i5,e18.3)')i,id(i),var(indx_old(i),1),var(indx_old(i),2), & 
               var(indx_old(i),3),var(indx_old(i),4),var(indx_old(i),5),itext(indx_old(i),2),   &
               var(indx_old(i),6)
       endif
    end do
    close (20,status='keep')
    deallocate (colin2cat,var,id,indx_old)

  END SUBROUTINE cti_stat_file
  
  !---------------------------------------------------------------------

  SUBROUTINE create_model_para (MaskFile, nbcatch, tile_lon, tile_lat, tile_pfs)

    character(*), intent(in) :: MaskFile
    integer,      intent(in) :: nbcatch
    real,         intent(in) :: tile_lon(:), tile_lat(:)
    integer,      intent(in) :: tile_pfs(:)

    ! --------------------------------------------
    
    integer i,n,k, tindex1,pfaf1
    integer soil_gswp
    real meanlu,stdev,minlu,maxlu,coesk,rzdep
    real minlat,maxlat,minlon,maxlon
    real,allocatable, dimension (:) ::   &
         BEE, PSIS,POROS,COND,WPWET,soildepth
    REAL, allocatable, dimension(:) :: TOPMEAN, TOPVAR, TOPSKEW
    REAL ST(NAR), AC(NAR),COESKEW
    REAL, allocatable, dimension (:) ::   &
         ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4, &
         ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,   &
         tsa1, tsa2,tsb1, tsb2,               &
         taberr1,taberr2,normerr1,normerr2,   &
         taberr3,taberr4,normerr3,normerr4
    integer, dimension(12) :: tile_pick
    integer, allocatable, dimension (:) :: soil_class_top,soil_class_com,tindex2,pfaf2
    real watdep(nwt,nrz),wan(nwt,nrz),rzexcn(nwt,nrz),frc(nwt,nrz)
    real, allocatable, dimension  (:,:,:,:) :: &
         gwatdep,gwan,grzexcn,gfrc
    real :: wtdep,wanom,rzaact,fracl,profdep,dist_save,tile_distance 
    character*512 :: pathout,fname,fout,losfile
    CHARACTER*512 :: version,resoln,continent
    character*6 :: rdep,ext
    integer :: iwt,irz,group
    logical :: picked

    integer :: ncid, status
    logical :: file_exists
    real, allocatable, dimension (:,:) :: parms4file

    ! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
    !
    ! NOTE: "!$" is for conditional compilation
    !
    logical :: running_omp = .false.
    !
    !$ integer :: omp_get_thread_num, omp_get_num_threads
    !
    integer :: n_threads=1, li, ui
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
    ! ----------- OpenMP PARALLEL ENVIRONMENT ----------------------------

    !c-------------------------------------------------------------------------

    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    losfile =trim(MAKE_BCS_INPUT_DIR)//'/land/soil/soil_water_loss/v1/loss_perday'
    !c     opening files


    allocate (gwatdep(1:nwt,1:nrz,1:12,1:2))
    allocate (gwan   (1:nwt,1:nrz,1:12,1:2))
    allocate (grzexcn(1:nwt,1:nrz,1:12,1:2))
    allocate (gfrc   (1:nwt,1:nrz,1:12,1:2))

    do n =1,12
       if(n.lt.10)write(ext,'(i1.1)')n
       if(n.ge.10)write(ext,'(i2.2)')n
       do i =1,2
          if (i==1) rdep='.rz75.'
          if (i==2) rdep='.rz1.'
          open (120,file=trim(losfile)//trim(rdep)//trim(ext),  &
               form='formatted',status='old')

          do iwt=1,nwt
             do irz=1,nrz
                read(120,2000) wtdep,wanom,rzaact,fracl
2000            format(1x,4e16.8)
                gwatdep(iwt,irz,n,i)=wtdep
                gwan(iwt,irz,n,i)=wanom
                grzexcn(iwt,irz,n,i)=rzaact
                gfrc(iwt,irz,n,i)=amin1(fracl,1.)
             enddo
	  enddo
          close (120,status='keep')	   
       end do
    end do
    fname='clsm/soil_param.first'       
    open (10,file=fname,action='read',       &
         form='formatted',status='old')              

    fname='clsm/cti_stats.dat'           
    open (11,file=fname,action='read',        &
         form='formatted',status='old')              

    fout='clsm/ar.new'               
    open (20,file=fout,action='write',        &
         form='formatted',status='unknown')          

    fout='clsm//bf.dat'               
    open (30,file=fout,action='write',        &
         form='formatted',status='unknown')          

    fout='clsm//ts.dat'               
    open (40,file=fout,action='write',        &
         form='formatted',status='unknown')        

    if (error_file) then 
       fout='clsm/ar_rmse.dat'           
       open (21,file=fout,action='write',        &
            form='formatted',status='unknown')

       fout='clsm/bf_rmse.dat'           
       open (31,file=fout,action='write',        &
            form='formatted',status='unknown')

       fout='clsm/bad_sat_param.tiles'
       open (41,file=fout,action='write',        &
            form='formatted',status='unknown')    
    endif
    fout='clsm/soil_param.dat'
    open (42,file=fout,action='write',        &
         form='formatted',status='unknown')    
    read (11,*) n ! read off nbcatch

    allocate (TOPMEAN (1:nbcatch))
    allocate (TOPVAR  (1:nbcatch))
    allocate (TOPSKEW (1:nbcatch))
    allocate (ARS1 (1:nbcatch))
    allocate (ARS2 (1:nbcatch))
    allocate (ARS3 (1:nbcatch))
    allocate (ARA1 (1:nbcatch))
    allocate (ARA2 (1:nbcatch))
    allocate (ARA3 (1:nbcatch))
    allocate (ARA4 (1:nbcatch))
    allocate (ARW1 (1:nbcatch))
    allocate (ARW2 (1:nbcatch))
    allocate (ARW3 (1:nbcatch))
    allocate (ARW4 (1:nbcatch))
    allocate (BF1 (1:nbcatch))
    allocate (BF2 (1:nbcatch))
    allocate (BF3 (1:nbcatch))
    allocate (TSA1 (1:nbcatch))
    allocate (TSA2 (1:nbcatch))
    allocate (TSB1 (1:nbcatch))
    allocate (TSB2 (1:nbcatch))
    allocate (TABERR1 (1:nbcatch))
    allocate (TABERR2 (1:nbcatch))
    allocate (TABERR3 (1:nbcatch))
    allocate (TABERR4 (1:nbcatch))
    allocate (NORMERR1 (1:nbcatch))
    allocate (NORMERR2 (1:nbcatch))
    allocate (NORMERR3 (1:nbcatch))
    allocate (NORMERR4 (1:nbcatch))
    allocate (BEE        (1:nbcatch))
    allocate (PSIS      (1:nbcatch))
    allocate (POROS     (1:nbcatch))
    allocate (COND      (1:nbcatch))
    allocate (WPWET     (1:nbcatch))
    allocate (soildepth (1:nbcatch))
    allocate (soil_class_top (1:nbcatch))
    allocate (soil_class_com (1:nbcatch))
    allocate (tindex2        (1:nbcatch))
    allocate (pfaf2          (1:nbcatch))

    do n=1,nbcatch

       read(11,'(i10,i8,5(1x,f8.4))') tindex1,pfaf1,meanlu,stdev,minlu,maxlu,coesk                                  
       read(10,*) tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n), &
            BEE(n),PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n)                        

       if(tindex1.ne.tindex2(n))then
          write(*,*)'Warnning 1: tindex mismatched'                        
          stop
       endif

       if(tile_pfs(n).ne.pfaf2(n)) then
          write(*,*)'Warnning 1: pfafstetter mismatched' 
          stop
       endif

       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
          TOPMEAN(n) = meanlu
       else
          TOPMEAN(n) = 0.961*meanlu-1.957                       
       endif

       TOPVAR(n)  = stdev*stdev                                
       TOPSKEW(n) = coesk*stdev*stdev*stdev                   

       if (TOPVAR(n) .eq. 0. .or. coesk .eq. 0. .or. topskew(n) .eq. 0.) then                       
          write(*,*) 'Problem: undefined values:'         
          write(*,*) TOPMEAN(n),TOPVAR(n),coesk,minlu,maxlu
          stop
       endif
    END DO

    inquire(file='clsm/catch_params.nc4', exist=file_exists)

    if(file_exists) then
       status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid) ; VERIFY_(STATUS)
       allocate (parms4file (1:nbcatch, 1:25))
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ),(/1/),(/nbcatch/),BEE  (:))      ; VERIFY_(STATUS)
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'COND' ),(/1/),(/nbcatch/),COND (:))      ; VERIFY_(STATUS)
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'POROS'),(/1/),(/nbcatch/),POROS(:))      ; VERIFY_(STATUS)
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ),(/1/),(/nbcatch/),PSIS (:))      ; VERIFY_(STATUS)
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'WPWET'),(/1/),(/nbcatch/),WPWET(:))      ; VERIFY_(STATUS)
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR'),(/1/),(/nbcatch/),soildepth (:)) ; VERIFY_(STATUS)
       parms4file (:,12) = BEE      (:)
       parms4file (:,16) = COND     (:)
       parms4file (:,18) = POROS    (:)
       parms4file (:,19) = PSIS     (:)
       parms4file (:,24) = WPWET    (:)
       parms4file (:,25) = soildepth(:)
    endif

    rewind(10)

    allocate(low_ind(n_threads))
    allocate(upp_ind(n_threads))
    low_ind(1)         = 1
    upp_ind(n_threads) = nbcatch

    if (running_omp)  then
       do i=1,n_threads-1

          upp_ind(i)   = low_ind(i) + (nbcatch/n_threads) - 1 
          low_ind(i+1) = upp_ind(i) + 1

       end do
    end if


    !$OMP PARALLELDO DEFAULT(NONE)                          &
    !$OMP SHARED( BEE, PSIS,POROS,COND,WPWET,soildepth,     &
    !$OMP        TOPMEAN, TOPVAR, TOPSKEW,                  &
    !$OMP        ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,       &
    !$OMP        ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,         &
    !$OMP        tsa1, tsa2,tsb1, tsb2,                     &
    !$OMP        taberr1,taberr2,normerr1,normerr2,         &
    !$OMP        taberr3,taberr4,normerr3,normerr4,         &
    !$OMP        gwatdep,gwan,grzexcn,gfrc,soil_class_com,  &
    !$OMP        n_threads, low_ind, upp_ind )              &
    !$OMP PRIVATE(k,li,ui,n,i,watdep,wan,rzexcn,frc,ST,AC,  &
    !$OMP COESKEW,profdep)

    do k=1,n_threads

       li = low_ind(k)
       ui = upp_ind(k)

       do n=li,ui
          !      if ((n == 877).or.(n == 880).or.(n == 881)) then
          !         print *,n
          !      endif
          !      print *,n
          !      pause
          !        c Gamma distribution
          CALL TGEN (                              &
               TOPMEAN(n),TOPVAR(n),TOPSKEW(n),     &
               ST,AC,COESKEW)

          !c      write(*,*) 'tgen4 ok'

          !c Areal fractioning parameters
          !      print *,'tileid:' ,n
          CALL SAT_PARAM(                                              &
               BEE(n),PSIS(n),POROS(n),COND(n),              & 
               WPWET(n), ST, AC, COESKEW,n,                  &
               soildepth(n),                                 &
               ars1(n),ars2(n),ars3(n),                      &
               ara1(n),ara2(n),ara3(n),ara4(n),              &
               arw1(n),arw2(n),arw3(n),arw4(n),              &
               taberr1(n),taberr2(n),taberr3(n),taberr4(n),  &
               normerr1(n),normerr2(n),normerr3(n),normerr4(n))


          CALL BASE_PARAM(                                  &
               BEE(n),PSIS(n),POROS(n),COND(n),              &
               ST, AC,                                       &
               bf1(n),bf2(n),bf3(n),                         &
               taberr1(n),taberr2(n),normerr1(n),normerr2(n) &
               )

          profdep=soildepth(n)/1000.
          profdep=amax1(1.,profdep)
          if (grzdep .gt. .75*profdep) then
             i=1
          else
             i=2 
          end if

	  watdep (:,:) =  gwatdep (:,:,soil_class_com(n),i)   
	  wan    (:,:) =  gwan    (:,:,soil_class_com(n),i)
	  rzexcn (:,:) =  grzexcn (:,:,soil_class_com(n),i)
	  frc    (:,:) =  gfrc    (:,:,soil_class_com(n),i)

          CALL TS_PARAM(                       &
               BEE(n),PSIS(n),POROS(n),         &
               ST, AC,                          &
               watdep,wan,rzexcn,frc,           &
               tsa1(n),tsa2(n),tsb1(n),tsb2(n)  &
               )

       END DO
    END DO
    !$OMP ENDPARALLELDO
    tile_pick = 0

    DO n=1,nbcatch
       if((arw1(n).ne.9999.).and.(ars1(n).ne.9999.))then
          if(tile_pick(soil_class_com(n)) == 0)  tile_pick(soil_class_com(n)) = n
       endif
    end do

    DO n=1,nbcatch
       !c Third subroutine for the parameters related to the transfers 
       !c to the water table
       !
       ! Writing the parameters, in the same order as in catchment.def
       !      if((ars1(n).lt.0.).and.(ars2(n).le.0.3).and.(ars3(n).le.0.04).and.(arw1(n).ne.9999.))then
       if((arw1(n).ne.9999.).and.(ars1(n).ne.9999.))then   
          write(20,'(i10,i8,f5.2,11(2x,e14.7))')   &
               tindex2(n),pfaf2(n),gnu,   &
               ars1(n),ars2(n),ars3(n),                   &
               ara1(n),ara2(n),ara3(n),ara4(n),           &
               arw1(n),arw2(n),arw3(n),arw4(n) 
          write(30,'(i10,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(n),bf2(n),bf3(n)
          write(40,'(i10,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
               tsa1(n),tsa2(n),tsb1(n),tsb2(n)

          write(42,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),    &
               BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n)  

          if (allocated (parms4file)) then
             parms4file (n, 1) = ara1(n)
             parms4file (n, 2) = ara2(n)
             parms4file (n, 3) = ara3(n)
             parms4file (n, 4) = ara4(n)
             parms4file (n, 5) = ars1(n)
             parms4file (n, 6) = ars2(n)
             parms4file (n, 7) = ars3(n)
             parms4file (n, 8) = arw1(n)
             parms4file (n, 9) = arw2(n)
             parms4file (n,10) = arw3(n)
             parms4file (n,11) = arw4(n)
             parms4file (n,13) = bf1(n)
             parms4file (n,14) = bf2(n)
             parms4file (n,15) = bf3(n)
             parms4file (n,17) = gnu
             parms4file (n,20) = tsa1(n)
             parms4file (n,21) = tsa2(n)
             parms4file (n,22) = tsb1(n)
             parms4file (n,23) = tsb2(n)
          endif
       else

          if(preserve_soiltype) then    
             picked=.false.
             ! Group3
             !     category 1  : Sand
             !     category 2  : Loamy Sand
             !     category 3  : Sandy Loam
             !     category 8  : Silty Clay Loam
             ! Group2
             !     category 4  : Silt Loam
             !     category 5  : Silt
             !     category 6  : Loam
             !     category 7  : Sandy Clay Loam    
             ! Group1
             !     category 9  : Clay Loam
             !     category 10 : Sandy Clay
             !     category 11 : Silty Clay
             !     category 12 : Clay

             if ((soil_class_com(n)>=9).and.(soil_class_com(n)<=12)) then	
                group=1
             else if ((soil_class_com(n)>=4).and.(soil_class_com(n)<=7)) then
                group=2
             else
                group=3
             endif

             if(tile_pick(soil_class_com(n)) > 0) then
                k = tile_pick(soil_class_com(n))
                picked=.true.
                if (error_file) then
                   write (41,*)n,k
                endif
                write(20,'(i10,i8,f5.2,11(2x,e14.7))')   &
                     tindex2(n),pfaf2(n),gnu,   &
                     ars1(k),ars2(k),ars3(k),                   &
                     ara1(k),ara2(k),ara3(k),ara4(k),           &
                     arw1(k),arw2(k),arw3(k),arw4(k) 
                ars1(n)=ars1(k)
                ars2(n)=ars2(k)
                ars3(n)=ars3(k)
                ara1(n)=ara1(k)
                ara2(n)=ara2(k)
                ara3(n)=ara3(k)
                ara4(n)=ara4(k)
                arw1(n)=arw1(k)
                arw2(n)=arw2(k)
                arw3(n)=arw3(k)
                arw4(n)=arw4(k)

                write(30,'(i10,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
                write(40,'(i10,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
                     tsa1(k),tsa2(k),tsb1(k),tsb2(k)             
                write(42,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),    &
                     BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k)  
                if (allocated (parms4file)) then
                   parms4file (n, 1) = ara1(k)
                   parms4file (n, 2) = ara2(k)
                   parms4file (n, 3) = ara3(k)
                   parms4file (n, 4) = ara4(k)
                   parms4file (n, 5) = ars1(k)
                   parms4file (n, 6) = ars2(k)
                   parms4file (n, 7) = ars3(k)
                   parms4file (n, 8) = arw1(k)
                   parms4file (n, 9) = arw2(k)
                   parms4file (n,10) = arw3(k)
                   parms4file (n,11) = arw4(k)
                   parms4file (n,12) = BEE(k)
                   parms4file (n,13) = bf1(k)
                   parms4file (n,14) = bf2(k)
                   parms4file (n,15) = bf3(k)
                   parms4file (n,16) = COND(k)
                   parms4file (n,17) = gnu
                   parms4file (n,18) = POROS(k)
                   parms4file (n,19) = PSIS(k)
                   parms4file (n,20) = tsa1(k)
                   parms4file (n,21) = tsa2(k)
                   parms4file (n,22) = tsb1(k)
                   parms4file (n,23) = tsb2(k)
                   parms4file (n,24) = wpwet    (k)
                   parms4file (n,25) = soildepth(k)
                endif
             else

                do k =n-1,1,-1 

                   if (group == 1) then
                      if ((soil_class_com(k)>=9).and.(soil_class_com(k)<=12))picked=.true.
                   endif

                   if (group == 2) then
                      if ((soil_class_com(k)>=4).and.(soil_class_com(k)<=7))	picked=.true.
                   endif

                   if (group == 3) then
                      if (((soil_class_com(k)>=1).and.(soil_class_com(k)<=3)).or. &
                           (soil_class_com(k)==8)) picked=.true.
                   endif

                   if (picked) then
                      if (error_file) then
                         write (41,*)n,k
                      endif

                      write(20,'(i10,i8,f5.2,11(2x,e14.7))')   &
                           tindex2(n),pfaf2(n),gnu,   &
                           ars1(k),ars2(k),ars3(k),                   &
                           ara1(k),ara2(k),ara3(k),ara4(k),           &
                           arw1(k),arw2(k),arw3(k),arw4(k) 
                      write(30,'(i10,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
                      write(40,'(i10,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
                           tsa1(k),tsa2(k),tsb1(k),tsb2(k)             
                      write(42,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),    &
                           BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k)  
                      ars1(n)=ars1(k)
                      ars2(n)=ars2(k)
                      ars3(n)=ars3(k)
                      ara1(n)=ara1(k)
                      ara2(n)=ara2(k)
                      ara3(n)=ara3(k)
                      ara4(n)=ara4(k)
                      arw1(n)=arw1(k)
                      arw2(n)=arw2(k)
                      arw3(n)=arw3(k)
                      arw4(n)=arw4(k)

                      if (allocated (parms4file)) then
                         parms4file (n, 1) = ara1(k)
                         parms4file (n, 2) = ara2(k)
                         parms4file (n, 3) = ara3(k)
                         parms4file (n, 4) = ara4(k)
                         parms4file (n, 5) = ars1(k)
                         parms4file (n, 6) = ars2(k)
                         parms4file (n, 7) = ars3(k)
                         parms4file (n, 8) = arw1(k)
                         parms4file (n, 9) = arw2(k)
                         parms4file (n,10) = arw3(k)
                         parms4file (n,11) = arw4(k)
                         parms4file (n,12) = BEE(k)
                         parms4file (n,13) = bf1(k)
                         parms4file (n,14) = bf2(k)
                         parms4file (n,15) = bf3(k)
                         parms4file (n,16) = COND(k)
                         parms4file (n,17) = gnu
                         parms4file (n,18) = POROS(k)
                         parms4file (n,19) = PSIS(k)
                         parms4file (n,20) = tsa1(k)
                         parms4file (n,21) = tsa2(k)
                         parms4file (n,22) = tsb1(k)
                         parms4file (n,23) = tsb2(k)
                         parms4file (n,24) = wpwet    (k)
                         parms4file (n,25) = soildepth(k)
                      endif
                      exit
                   endif

                   if((k==1) .and. (.not. picked)) then
                      print *,'Warning ar.new is bad at n=',n
                      stop
                   endif
                end do
             endif


             !        write(30,'(i8,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(n),bf2(n),bf3(n)
             !        write(40,'(i8,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
             !             tsa1(n),tsa2(n),tsb1(n),tsb2(n)
          else

             dist_save = 1000000.
             k = 0
             do i = 1,nbcatch
                if(i /= n) then
                   if((ars1(i).ne.9999.).and.(arw1(i).ne.9999.)) then

                      tile_distance = (tile_lon(i) - tile_lon(n)) * (tile_lon(i) - tile_lon(n)) + &
                           (tile_lat(i) - tile_lat(n)) * (tile_lat(i) - tile_lat(n))
                      if(tile_distance < dist_save) then
                         k = i
                         dist_save = tile_distance
                      endif
                   endif
                endif
             enddo
             write (41,*)n,k
             write(20,'(i10,i8,f5.2,11(2x,e14.7))')   &
                  tindex2(n),pfaf2(n),gnu,   &
                  ars1(k),ars2(k),ars3(k),                   &
                  ara1(k),ara2(k),ara3(k),ara4(k),           &
                  arw1(k),arw2(k),arw3(k),arw4(k) 
             write(30,'(i10,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)
             write(40,'(i10,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,    &
                  tsa1(k),tsa2(k),tsb1(k),tsb2(k)
             write(42,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.3)') tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),    &
                  BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k)  
             if (allocated (parms4file)) then
                parms4file (n, 1) = ara1(k)
                parms4file (n, 2) = ara2(k)
                parms4file (n, 3) = ara3(k)
                parms4file (n, 4) = ara4(k)
                parms4file (n, 5) = ars1(k)
                parms4file (n, 6) = ars2(k)
                parms4file (n, 7) = ars3(k)
                parms4file (n, 8) = arw1(k)
                parms4file (n, 9) = arw2(k)
                parms4file (n,10) = arw3(k)
                parms4file (n,11) = arw4(k)
                parms4file (n,12) = BEE(k)
                parms4file (n,13) = bf1(k)
                parms4file (n,14) = bf2(k)
                parms4file (n,15) = bf3(k)
                parms4file (n,16) = COND(k)
                parms4file (n,17) = gnu
                parms4file (n,18) = POROS(k)
                parms4file (n,19) = PSIS(k)
                parms4file (n,20) = tsa1(k)
                parms4file (n,21) = tsa2(k)
                parms4file (n,22) = tsb1(k)
                parms4file (n,23) = tsb2(k)
                parms4file (n,24) = wpwet    (k)
                parms4file (n,25) = soildepth(k)
             endif
          endif
       endif


       if (error_file) then
          write(21,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),taberr3(n),taberr4(n), &
               normerr1(n),normerr2(n),normerr3(n),normerr4(n)
          write(31,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),normerr1(n),normerr2(n)
       endif

    END DO

    !      Write(*,*) 'END COMPUTING MODEL PARA'

    close(10,status='keep')
    close(20,status='keep')
    close(30,status='keep')
    close(40,status='keep')
    close(11,status='keep')
    close(12,status='keep')
    close(42,status='keep')
    if (error_file) then
       close(21,status='delete')
       close(31,status='delete')
       close(41,status='keep')
    endif

    if(file_exists) then
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA1' ) ,(/1/),(/nbcatch/), parms4file (:, 1)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA2' ) ,(/1/),(/nbcatch/), parms4file (:, 2)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA3' ) ,(/1/),(/nbcatch/), parms4file (:, 3)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA4' ) ,(/1/),(/nbcatch/), parms4file (:, 4)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARS1' ) ,(/1/),(/nbcatch/), parms4file (:, 5)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARS2' ) ,(/1/),(/nbcatch/), parms4file (:, 6)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARS3' ) ,(/1/),(/nbcatch/), parms4file (:, 7)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW1' ) ,(/1/),(/nbcatch/), parms4file (:, 8)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW2' ) ,(/1/),(/nbcatch/), parms4file (:, 9)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW3' ) ,(/1/),(/nbcatch/), parms4file (:,10)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW4' ) ,(/1/),(/nbcatch/), parms4file (:,11)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ) ,(/1/),(/nbcatch/), parms4file (:,12)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BF1'  ) ,(/1/),(/nbcatch/), parms4file (:,13)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BF2'  ) ,(/1/),(/nbcatch/), parms4file (:,14)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BF3'  ) ,(/1/),(/nbcatch/), parms4file (:,15)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'COND' ) ,(/1/),(/nbcatch/), parms4file (:,16)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'GNU'  ) ,(/1/),(/nbcatch/), parms4file (:,17)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'POROS') ,(/1/),(/nbcatch/), parms4file (:,18)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ) ,(/1/),(/nbcatch/), parms4file (:,19)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSA1' ) ,(/1/),(/nbcatch/), parms4file (:,20)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSA2' ) ,(/1/),(/nbcatch/), parms4file (:,21)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSB1' ) ,(/1/),(/nbcatch/), parms4file (:,22)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSB2' ) ,(/1/),(/nbcatch/), parms4file (:,23)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'WPWET') ,(/1/),(/nbcatch/), parms4file (:,24)) ; VERIFY_(STATUS)
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR') ,(/1/),(/nbcatch/), parms4file (:,25)) ; VERIFY_(STATUS)
       STATUS = NF_CLOSE (NCID) ; VERIFY_(STATUS)
       DEALLOCATE (parms4file)
    endif

  END SUBROUTINE create_model_para

  !--------------------------------------------------------------------

  SUBROUTINE create_model_para_woesten (Maskfile, nbcatch, tile_lon, tile_lat, tile_pfs)

    character(*), intent(in) :: MaskFile
    integer,      intent(in) :: nbcatch
    real,         intent(in) :: tile_lon(:), tile_lat(:)
    integer,      intent(in) :: tile_pfs(:)
    ! -----------------------------------------------

    real, allocatable, dimension (:)  :: a_sand,a_clay,a_silt,a_oc,  &
         atile_sand,atile_clay, grav_vec, soc_vec,&
         poc_vec,a_sand_surf,a_clay_surf,wpwet_surf,poros_surf, pmap

    integer i,j,n,k, tindex1,pfaf1
    integer soil_gswp
    real meanlu,stdev,minlu,maxlu,coesk,rzdep
    real minlat,maxlat,minlon,maxlon
    real,allocatable, dimension (:) ::   &
         BEE, PSIS,POROS,COND,WPWET,soildepth
    REAL, allocatable, dimension(:) :: TOPMEAN, TOPVAR, TOPSKEW
    REAL ST(NAR), AC(NAR),COESKEW
    REAL, allocatable, dimension (:) ::   &
         ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4, &
         ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,   &
         tsa1, tsa2,tsb1, tsb2,               &
         taberr1,taberr2,normerr1,normerr2,   &
         taberr3,taberr4,normerr3,normerr4

    integer, allocatable, dimension (:) :: soil_class_com,tindex2,pfaf2, &
         soil_class_top
    real watdep(nwt,nrz),wan(nwt,nrz),rzexcn(nwt,nrz),frc(nwt,nrz)
    real, allocatable, dimension  (:,:,:) :: &
         gwatdep,gwan,grzexcn,gfrc
    real :: wtdep,wanom,rzaact,fracl,profdep,dist_save,     &
         ncells_top, ncells_top_pro,ncells_sub_pro,tile_distance
    character*512 :: pathout,fname,fout,losfile
    CHARACTER*512 :: version,resoln,continent
    character*6  ::rdep,ext
    integer :: iwt,irz,group
    logical :: picked
    logical :: file_exists
    REAL, ALLOCATABLE, DIMENSION (:,:) :: parms4file
    integer :: ncid, status

    ! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
    !
    ! NOTE: "!$" is for conditional compilation
    !
    logical :: running_omp = .false.
    !
    !$ integer :: omp_get_thread_num, omp_get_num_threads
    !
    integer :: n_threads=1, li, ui
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

    !c-------------------------------------------------------------------------

    ! SoilClasses-SoilHyd-TauParam.dat and SoilClasses-SoilHyd-TauParam.peatmap differ
    ! only in the parameters for the peat class #253.  The file *.peatmap contains
    ! the PEATCLSM parameters from Table 2 of Bechtold et al. 2019 (doi:10.1029/2018MS001574).
    !
    ! Note: K_s = COND*exp(-zks*gnu)   ==> with zks=2 and gnu=1, K_s = 0.135335*COND
    !
    !         K_s      COND     [m/s]
    ! NLv4  7.86e-7   5.81e-6
    ! NLv5  3.79e-6   2.80e-5   <== note *typo* in Table 2 of Bechtold et al. 2019, which erroneously lists K_s=2.8e-5

    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)

    if(use_PEATMAP) then 
       fname = trim(MAKE_BCS_INPUT_DIR)//'/land/soil/SOIL-DATA/SoilClasses-SoilHyd-TauParam.peatmap' 
    else
       fname = trim(MAKE_BCS_INPUT_DIR)//'land/soil/SOIL-DATA/SoilClasses-SoilHyd-TauParam.dat'
    endif
    open (11, file=trim(fname), form='formatted',status='old', &
         action = 'read')
    read (11,'(a)')fout        ! read header line

    losfile =trim(MAKE_BCS_INPUT_DIR)//'/land/soil/soil_water_loss/v2/loss_pd_top/loss_perday_rz1m_'

    allocate (a_sand (1:n_SoilClasses))
    allocate (a_silt (1:n_SoilClasses))
    allocate (a_clay (1:n_SoilClasses))
    allocate (a_oc   (1:n_SoilClasses))
    allocate (gwatdep(1:nwt,1:nrz,1:n_SoilClasses))
    allocate (gwan   (1:nwt,1:nrz,1:n_SoilClasses))
    allocate (grzexcn(1:nwt,1:nrz,1:n_SoilClasses))
    allocate (gfrc   (1:nwt,1:nrz,1:n_SoilClasses))

    do n =1,n_SoilClasses

       ! read sand/clay/orgC for class n defined in SoilClasses-SoilHyd-TauParam.*

       read (11,'(4f7.3)')a_sand(n),a_clay(n),a_silt(n),a_oc(n)
       write (fout,'(i2.2,i2.2,i4.4)')nint(a_sand(n)),nint(a_clay(n)),nint(100*a_oc(n))

       ! open and read loss parameter file for class n (defined through sand/clay/orgC)

       if(n == n_SoilClasses .and. use_PEATMAP) then 
          open (120,file=trim(losfile)//trim(fout)//'.peat',  &
               form='formatted',status='old')
       else
          open (120,file=trim(losfile)//trim(fout),  &
               form='formatted',status='old')
       endif

       do iwt=1,nwt
          do irz=1,nrz
             read(120,2000) wtdep,wanom,rzaact,fracl
2000         format(1x,4e16.8)
             gwatdep(iwt,irz,n)= wtdep
             gwan(iwt,irz,n)   = wanom
             grzexcn(iwt,irz,n)= rzaact
             gfrc(iwt,irz,n)   = amin1(fracl,1.)
          enddo
       enddo
       close (120,status='keep')	   
    end do
    close (11,status='keep')  
    deallocate (a_sand,a_silt,a_clay,a_oc)

    ! open files for *reading*

    fname='clsm/soil_param.first'       
    open (10,file=fname,action='read',       &
         form='formatted',status='old')              

    fname='clsm/cti_stats.dat'           
    open (11,file=fname,action='read',        &
         form='formatted',status='old')  

    ! open files for *writing*

    fout='clsm/ar.new'               
    open (20,file=fout,action='write',        &
         form='formatted',status='unknown')          

    fout='clsm//bf.dat'               
    open (30,file=fout,action='write',        &
         form='formatted',status='unknown')          

    fout='clsm//ts.dat'               
    open (40,file=fout,action='write',        &
         form='formatted',status='unknown')        

    if (error_file) then 
       fout='clsm/ar_rmse.dat'           
       open (21,file=fout,action='write',        &
            form='formatted',status='unknown')

       fout='clsm/bf_rmse.dat'           
       open (31,file=fout,action='write',        &
            form='formatted',status='unknown')

       fout='clsm/bad_sat_param.tiles'
       open (41,file=fout,action='write',        &
            form='formatted',status='unknown')  

    endif

    fout='clsm/soil_param.dat'
    open (42,file=fout,action='write',        &
         form='formatted',status='unknown')       

    read (11,*) n    ! read off header line (number of tiles) -- cti_stats.dat

    allocate (TOPMEAN (1:nbcatch))
    allocate (TOPVAR  (1:nbcatch))
    allocate (TOPSKEW (1:nbcatch))
    allocate (ARS1 (1:nbcatch))
    allocate (ARS2 (1:nbcatch))
    allocate (ARS3 (1:nbcatch))
    allocate (ARA1 (1:nbcatch))
    allocate (ARA2 (1:nbcatch))
    allocate (ARA3 (1:nbcatch))
    allocate (ARA4 (1:nbcatch))
    allocate (ARW1 (1:nbcatch))
    allocate (ARW2 (1:nbcatch))
    allocate (ARW3 (1:nbcatch))
    allocate (ARW4 (1:nbcatch))
    allocate (BF1 (1:nbcatch))
    allocate (BF2 (1:nbcatch))
    allocate (BF3 (1:nbcatch))
    allocate (TSA1 (1:nbcatch))
    allocate (TSA2 (1:nbcatch))
    allocate (TSB1 (1:nbcatch))
    allocate (TSB2 (1:nbcatch))
    allocate (TABERR1 (1:nbcatch))
    allocate (TABERR2 (1:nbcatch))
    allocate (TABERR3 (1:nbcatch))
    allocate (TABERR4 (1:nbcatch))
    allocate (NORMERR1 (1:nbcatch))
    allocate (NORMERR2 (1:nbcatch))
    allocate (NORMERR3 (1:nbcatch))
    allocate (NORMERR4 (1:nbcatch))
    allocate (BEE        (1:nbcatch))
    allocate (PSIS      (1:nbcatch))
    allocate (POROS     (1:nbcatch))
    allocate (COND      (1:nbcatch))
    allocate (WPWET     (1:nbcatch))
    allocate (soildepth (1:nbcatch))
    allocate (soil_class_top (1:nbcatch))
    allocate (soil_class_com (1:nbcatch))
    allocate (tindex2        (1:nbcatch))
    allocate (pfaf2          (1:nbcatch))
    allocate (atile_clay     (1:nbcatch))
    allocate (atile_sand     (1:nbcatch))
    allocate (grav_vec       (1:nbcatch))
    allocate (soc_vec        (1:nbcatch))
    allocate (poc_vec        (1:nbcatch))
    allocate (a_sand_surf    (1:nbcatch))
    allocate (a_clay_surf    (1:nbcatch))
    allocate (wpwet_surf     (1:nbcatch))
    allocate (poros_surf     (1:nbcatch))
    allocate (pmap           (1:nbcatch))

    do n=1,nbcatch

       ! read cti_stats.dat 

       read(11,'(i10,i8,5(1x,f8.4))') tindex1,pfaf1,meanlu,stdev  &
            ,minlu,maxlu,coesk                                  

       ! read soil_param.first
       !
       ! WARNING: Immediately after the present do loop, BEE, COND, POROS, PSIS, WPWET, and 
       !          soildepth will be read again (and thus overwritten) with the values from 
       !          the catch_params.nc4 file.  It is unclear if the values in soil_param.first
       !          and catch_params.nc4 differ. See comments below.

       read(10,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4)') &
            tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),      &
            BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n),       &
            grav_vec(n),soc_vec(n),poc_vec(n),                            &
            a_sand_surf(n),a_clay_surf(n),atile_sand(n),atile_clay(n)                     
       if(tindex1.ne.tindex2(n))then
          write(*,*)'Warnning 1: tindex mismatched'                        
          stop
       endif

       if(tile_pfs(n).ne.pfaf2(n)) then
          write(*,*)'Warnning 1: pfafstetter mismatched' 
          stop
       endif
       if((use_PEATMAP).and.(soil_class_top(n) == 253)) then
          meanlu = 9.3
          stdev  = 0.12
          minlu  = 8.5
          maxlu  = 11.5
          coesk  = 0.25
       endif

       if (index(MaskFile,'GEOS5_10arcsec_mask') /= 0) then
          TOPMEAN(n) = meanlu
       else
          TOPMEAN(n) = 0.961*meanlu-1.957                       
       endif

       TOPVAR(n)  = stdev*stdev                                
       TOPSKEW(n) = coesk*stdev*stdev*stdev                   

       if ( TOPVAR(n) .eq. 0. .or. coesk .eq. 0.            &
            .or. topskew(n) .eq. 0.) then                       
          write(*,*) 'Problem: undefined values:'         
          write(*,*) TOPMEAN(n),TOPVAR(n),coesk,            &
               minlu,maxlu
          stop
       endif
    END DO    ! n=1,nbcatch

    inquire(file='clsm/catch_params.nc4', exist=file_exists)

    if(file_exists) then

       ! Read BEE, COND, POROS, PSIS, WPWET, and soildepth from nc4 file.
       ! It is unclear if parameters in nc4 file differ from those in soil_param.first, which were read 
       ! in the do loop just above.
       ! Probably, the parameters differ by roundoff because soil_param.first is an ASCII file and 
       ! catch_params.nc4 is a netcdf file.  Consequently, the parameters from the nc4 file are used
       ! in the calculation of the ar.new, bf.dat, and ts.dat parameters, which comes next.
       ! To maintain consistency between the parameters in soil_param.first and soil_param.dat where 
       ! no changes are needed, soil_param.first needs to be read again below (so as to overwrite
       ! the values from the nc4 file).  
       ! Why the parameters from the nc4 file are read here in the first place remains a mystery.
       ! Removing this read, however, will (almost certainly) result in non-zero-diff changes
       ! for existing bcs datasets.
       ! - reichle, 28 April 2022        

       status = NF_OPEN ('clsm/catch_params.nc4', NF_WRITE, ncid) ; VERIFY_(STATUS)
       allocate (parms4file (1:nbcatch, 1:25))
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ) ,(/1/),(/nbcatch/), BEE  (:)) ; VERIFY_(STATUS) 
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'COND' ) ,(/1/),(/nbcatch/), COND (:)) ; VERIFY_(STATUS) 
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'POROS') ,(/1/),(/nbcatch/), POROS(:)) ; VERIFY_(STATUS) 
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ) ,(/1/),(/nbcatch/), PSIS (:)) ; VERIFY_(STATUS)    
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'WPWET') ,(/1/),(/nbcatch/), WPWET(:)) ; VERIFY_(STATUS) 
       status = NF_GET_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR') ,(/1/),(/nbcatch/), soildepth (:)) ; VERIFY_(STATUS) 
       parms4file (:,12) = BEE      (:)
       parms4file (:,16) = COND     (:)
       parms4file (:,18) = POROS    (:)
       parms4file (:,19) = PSIS     (:)
       parms4file (:,24) = wpwet    (:)
       parms4file (:,25) = soildepth(:)
    endif

    rewind(10)                ! soil_param.first (so soil_param.first can be read again below...)

    allocate(low_ind(n_threads))
    allocate(upp_ind(n_threads))
    low_ind(1)         = 1
    upp_ind(n_threads) = nbcatch

    if (running_omp)  then
       do i=1,n_threads-1

          upp_ind(i)   = low_ind(i) + (nbcatch/n_threads) - 1 
          low_ind(i+1) = upp_ind(i) + 1

       end do
    end if


    !$OMP PARALLELDO DEFAULT(NONE)                          &
    !$OMP SHARED( BEE, PSIS,POROS,COND,WPWET,soildepth,     &
    !$OMP        TOPMEAN, TOPVAR, TOPSKEW,                  &
    !$OMP        ARS1,ARS2,ARS3, ARA1,ARA2,ARA3,ARA4,       &
    !$OMP        ARW1,ARW2,ARW3,ARW4,bf1, bf2, bf3,         &
    !$OMP        tsa1, tsa2,tsb1, tsb2,                     &
    !$OMP        taberr1,taberr2,normerr1,normerr2,         &
    !$OMP        taberr3,taberr4,normerr3,normerr4,         &
    !$OMP        gwatdep,gwan,grzexcn,gfrc,soil_class_com,  &
    !$OMP        n_threads, low_ind, upp_ind, use_PEATMAP ) &
    !$OMP PRIVATE(k,li,ui,n,i,watdep,wan,rzexcn,frc,ST,AC,  &
    !$OMP COESKEW,profdep)

    do k=1,n_threads

       li = low_ind(k)
       ui = upp_ind(k)

       do n=li,ui

          CALL TGEN (                              &
               TOPMEAN(n),TOPVAR(n),TOPSKEW(n),     &
               ST,AC,COESKEW)

          ! compute areal fractioning parameters (ar.new)

          CALL SAT_PARAM(                                              &
               BEE(n),PSIS(n),POROS(n),COND(n),              & 
               WPWET(n), ST, AC, COESKEW,n,                  &
               soildepth(n),                                 &
               ars1(n),ars2(n),ars3(n),                      &
               ara1(n),ara2(n),ara3(n),ara4(n),              &
               arw1(n),arw2(n),arw3(n),arw4(n),              &
               taberr1(n),taberr2(n),taberr3(n),taberr4(n),  &
               normerr1(n),normerr2(n),normerr3(n),normerr4(n))

          ! compute base flow parameters (bf.dat)

          CALL BASE_PARAM(                                  &
               BEE(n),PSIS(n),POROS(n),COND(n),              &
               ST, AC,                                       &
               bf1(n),bf2(n),bf3(n),                         &
               taberr1(n),taberr2(n),normerr1(n),normerr2(n) &
               )


          watdep (:,:) =  gwatdep (:,:,soil_class_com(n))   
          wan    (:,:) =  gwan    (:,:,soil_class_com(n))
          rzexcn (:,:) =  grzexcn (:,:,soil_class_com(n))
          frc    (:,:) =  gfrc    (:,:,soil_class_com(n))

          ! compute time scale parameters (rzexc-catdef) (ts.dat)

          CALL TS_PARAM(                       &
               BEE(n),PSIS(n),POROS(n),         &
               ST, AC,                          &
               watdep,wan,rzexcn,frc,           &
               tsa1(n),tsa2(n),tsb1(n),tsb2(n)  &
               )

          if(soil_class_com(n) == 253 .and. use_PEATMAP) then

             ! Michel Bechtold paper - PEATCLSM_fitting_CLSM_params.R produced these data values.

             ars1(n) = -7.9514018e-03
             ars2(n) = 6.2297356e-02 
             ars3(n) = 1.9187240e-03                   
             ara1(n) = 8.9551220e+00 
             ara2(n) = 9.8149664e+02 
             ara3(n) = 8.9551220e+00 
             ara4(n) = 9.8149664e+02 
             arw1(n) = 9.9466055e-03 
             arw2(n) = 1.0881960e-02 
             arw3(n) = 1.5309287e-05 
             arw4(n) = 1.0000000e-04 

             bf1(n) = 4.6088086e+02  
             bf2(n) = 1.4237401e-01  
             bf3(n) = 6.9803000e+00

             tsa1(n) = -2.417581e+00  
             tsa2(n) = -4.784762e+00  
             tsb1(n) = -3.700285e-03  
             tsb2(n) = -2.392484e-03

          endif
       END DO
    END DO
    !$OMP ENDPARALLELDO


    ! ----------------------------------------------------------------------------------------
    !
    ! write ar.new, bf.dat, ts.dat, and soil_param.dat

    DO n=1,nbcatch

       ! Read soil_param.first again...; this is (almost certainly) needed to maintain consistency
       ! between soil_param.first and soil_param.dat, see comments above.

       read(10,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4, f8.4)') &
            tindex2(n),pfaf2(n),soil_class_top(n),soil_class_com(n),         &
            BEE(n), PSIS(n),POROS(n),COND(n),WPWET(n),soildepth(n),       &
            grav_vec(n),soc_vec(n),poc_vec(n),                            &
            a_sand_surf(n),a_clay_surf(n),atile_sand(n),atile_clay(n) ,   &
            wpwet_surf(n),poros_surf(n), pmap(n)

       ! This revised if block replaces the complex, nested if block commented out above

       if ( (ars1(n)==9999.) .or. (arw1(n)==9999.) ) then 

          ! some parameter values are no-data --> find nearest tile k with good parameters

          dist_save = 1000000.
          k = 0
          do i = 1,nbcatch
             if(i /= n) then
                if((ars1(i).ne.9999.).and.(arw1(i).ne.9999.)) then

                   tile_distance = (tile_lon(i) - tile_lon(n)) * (tile_lon(i) - tile_lon(n)) + &
                        (tile_lat(i) - tile_lat(n)) * (tile_lat(i) - tile_lat(n))
                   if(tile_distance < dist_save) then
                      k = i
                      dist_save = tile_distance
                   endif
                endif
             endif
          enddo
          ! record in file clsm/bad_sat_param.tiles
          write (41,*)n,k        ! n="bad" tile, k=tile from which parameters are taken

          ! Overwrite parms4file when filling in parameters from neighboring tile k.
          ! For "good" tiles, keep parms4file as read earlier from catch_params.nc4,
          ! which is why this must be done within the "then" block of the "if" statement.
          ! This is necessary for backward 0-diff compatibility of catch_params.nc4.

          parms4file (n,12) = BEE(k)
          parms4file (n,16) = COND(k)
          parms4file (n,18) = POROS(k)
          parms4file (n,19) = PSIS(k)
          parms4file (n,24) = wpwet(k)
          parms4file (n,25) = soildepth(k)

       else

          ! nominal case, all parameters are good

          k = n

       end if

       ! for current tile n, write parameters of tile k into ar.new (20), bf.dat (30), ts.dat (40), 
       !   and soil_param.dat (42)

       write(20,'(i10,i8,f5.2,11(2x,e14.7))')          &
            tindex2(n),pfaf2(n),gnu,                   &
            ars1(k),ars2(k),ars3(k),                   &
            ara1(k),ara2(k),ara3(k),ara4(k),           &
            arw1(k),arw2(k),arw3(k),arw4(k) 

       write(30,'(i10,i8,f5.2,3(2x,e13.7))')tindex2(n),pfaf2(n),gnu,bf1(k),bf2(k),bf3(k)

       write(40,'(i10,i8,f5.2,4(2x,e13.7))')tindex2(n),pfaf2(n),gnu,                      &
            tsa1(k),tsa2(k),tsb1(k),tsb2(k)

       write(42,'(i10,i8,i4,i4,3f8.4,f12.8,f7.4,f10.4,3f7.3,4f7.3,2f10.4, f8.4)')         &
            tindex2(n),pfaf2(n),soil_class_top(k),soil_class_com(k),                      &
            BEE(k), PSIS(k),POROS(k),COND(k),WPWET(k),soildepth(k),                       &
            grav_vec(k),soc_vec(k),poc_vec(k),                                            &
            a_sand_surf(k),a_clay_surf(k),atile_sand(k),atile_clay(k) ,                   &
            wpwet_surf(k),poros_surf(k), pmap(k)

       ! record ar.new, bf.dat, and ts.dat parameters for later writing into catch_params.nc4

       if (allocated (parms4file)) then
          parms4file (n, 1) = ara1(k)
          parms4file (n, 2) = ara2(k)
          parms4file (n, 3) = ara3(k)
          parms4file (n, 4) = ara4(k)
          parms4file (n, 5) = ars1(k)
          parms4file (n, 6) = ars2(k)
          parms4file (n, 7) = ars3(k)
          parms4file (n, 8) = arw1(k)
          parms4file (n, 9) = arw2(k)
          parms4file (n,10) = arw3(k)
          parms4file (n,11) = arw4(k)  
          parms4file (n,13) = bf1(k)
          parms4file (n,14) = bf2(k)
          parms4file (n,15) = bf3(k)
          parms4file (n,17) = gnu
          parms4file (n,20) = tsa1(k)
          parms4file (n,21) = tsa2(k)
          parms4file (n,22) = tsb1(k)
          parms4file (n,23) = tsb2(k)
       endif

       if (error_file) then
          write(21,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),taberr3(n),taberr4(n), &
               normerr1(n),normerr2(n),normerr3(n),normerr4(n)
          write(31,*)tindex2(n),pfaf2(n),taberr1(n),taberr2(n),normerr1(n),normerr2(n)
       endif

    END DO              ! n=1,nbcatch

    !      Write(*,*) 'END COMPUTING MODEL PARA'

    close(10,status='keep')
    close(11,status='keep')
    close(20,status='keep')
    close(30,status='keep')
    close(40,status='keep')
    close(42,status='keep')


    if (error_file) then
       close(21,status='delete')
       close(31,status='delete')
       close(41,status='keep')
    endif

    if(file_exists) then
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA1' ) ,(/1/),(/nbcatch/), parms4file (:, 1)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA2' ) ,(/1/),(/nbcatch/), parms4file (:, 2)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA3' ) ,(/1/),(/nbcatch/), parms4file (:, 3)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARA4' ) ,(/1/),(/nbcatch/), parms4file (:, 4)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARS1' ) ,(/1/),(/nbcatch/), parms4file (:, 5)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARS2' ) ,(/1/),(/nbcatch/), parms4file (:, 6)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARS3' ) ,(/1/),(/nbcatch/), parms4file (:, 7)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW1' ) ,(/1/),(/nbcatch/), parms4file (:, 8)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW2' ) ,(/1/),(/nbcatch/), parms4file (:, 9)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW3' ) ,(/1/),(/nbcatch/), parms4file (:,10)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'ARW4' ) ,(/1/),(/nbcatch/), parms4file (:,11)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BEE'  ) ,(/1/),(/nbcatch/), parms4file (:,12)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BF1'  ) ,(/1/),(/nbcatch/), parms4file (:,13)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BF2'  ) ,(/1/),(/nbcatch/), parms4file (:,14)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'BF3'  ) ,(/1/),(/nbcatch/), parms4file (:,15)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'COND' ) ,(/1/),(/nbcatch/), parms4file (:,16)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'GNU'  ) ,(/1/),(/nbcatch/), parms4file (:,17)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'POROS') ,(/1/),(/nbcatch/), parms4file (:,18)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'PSIS' ) ,(/1/),(/nbcatch/), parms4file (:,19)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSA1' ) ,(/1/),(/nbcatch/), parms4file (:,20)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSA2' ) ,(/1/),(/nbcatch/), parms4file (:,21)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSB1' ) ,(/1/),(/nbcatch/), parms4file (:,22)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'TSB2' ) ,(/1/),(/nbcatch/), parms4file (:,23)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'WPWET') ,(/1/),(/nbcatch/), parms4file (:,24)) ; VERIFY_(STATUS) 
       status = NF_PUT_VARA_REAL(NCID,NC_VarID(NCID,'DP2BR') ,(/1/),(/nbcatch/), parms4file (:,25)) ; VERIFY_(STATUS) 
       STATUS   = NF_CLOSE (NCID) ; VERIFY_(STATUS)
       DEALLOCATE (parms4file)
    endif

  END SUBROUTINE create_model_para_woesten


  !---------------------------------------------------------------------
  
  SUBROUTINE TS_PARAM(                                &
       BEE,PSIS,POROS,                &
       VALX, PX,                      &
       watdep,wan,rzexcn,frc,         &
       tsa1,tsa2,tsb1,tsb2            &
       )

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c                                                                         c
    !c Given pre-computed 1-D relationships between a "local" root zone excess c 
    !c and a "local" catchment deficit, the timescale of the bulk vertical     c
    !c transfer between the two bulk prognostic variables is computed using    c
    !c the distribution of the local deficit established from the distribution c
    !c of the topographic index, then an approximated function of catdef and   c
    !c rzex is derived.                                                        c
    !c                                                                         c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    INTEGER NAR0
    REAL, intent (in) :: BEE, PSIS, POROS
    REAL, intent (in) :: VALX(NAR), PX(NAR)
    real, intent (inout) :: watdep(nwt,nrz),wan(nwt,nrz),  &
         rzexcn(nwt,nrz),frc(nwt,nrz)
    real, intent (out) ::  tsa1, tsa2 ,tsb1, tsb2

    integer :: tex,iwt,irz,n,idep,k, index1,i0
    REAL VALX0(NAR), PX0(NAR),sumta,sumta2,timean,zbar, rzw
    REAL :: term1, term2, sumdef, suma, frcsat,rzexc, rzact
    real zdep(nar),def(nar),wrz(nar),wbin(500),rze(nar)
    real catd(2,2),tsc(2,2), satfrc,sumfrac,sumz,frac
    real, parameter ::  frcmax = .041
    real  wtdep,wanom,rzaact,fracl,profdep,rzdep

    !      logical bug

    !c----------------------------------------------------------------
    !c Is loss.dat compatible with rzdep = 0.49 ???

    rzdep = grzdep

    !c Convert fractions to "per-hour" values
    do iwt=1,nwt
       do irz=1,nrz
          frc(iwt,irz)=1.-((1.-frc(iwt,irz))**(1./24.))
       enddo
    enddo

    nar0=0
    do n=1,nar
       if (px(n) .ne. 0.) then
          nar0=nar0+1
          valx0(nar0)=valx(n)
          px0(nar0)=px(n)
       endif
    enddo

    sumta=0.
    sumta2=0.
    suma=0.
    do n=1,nar0
       sumta=sumta+px0(n)*valx0(n)
       sumta2=sumta2+px0(n)*valx0(n)*valx0(n)
       suma=suma+px0(n)
    enddo

    timean=sumta/suma

    !c**** Loop over two water table depths
    do idep=1,2
       if(idep.eq.1) zbar=1.5 ! zbar in meters
       if(idep.eq.2) zbar=2.0

       !c**** Compute array of water table depths:
       do k=1,nar0
          term1=(1/gnu)*(valx0(k)-timean)
          zdep(k)=zbar-term1
          if(zdep(k) .lt. 0.) zdep(k)=0.
       enddo
       !c            write(*,*)"  End water table depth"
       !c**** Compute array of moisture deficits:
       do k=1,nar0
          term1=(psis-zdep(k))/psis
          term1=term1**(1.-1./bee)
          term2=-psis*(bee/(bee-1.))*(term1-1.)
          def(k)=poros*(zdep(k)-term2)
       enddo

       !c**** Add deficits to produce catdef:
       sumdef=0.
       do k=1,nar0
          sumdef=sumdef+def(k)*px0(k)*1000.
       enddo
       !c            write(*,*)"  End catchment deficit"
       !c**** Compute array of root zone moisture (degree of wetness in root zone):
       do k=1,nar0

          if(zdep(k).eq.0.) then
             wrz(k)=1.
          elseif(zdep(k)-rzdep.lt.0.) then
             term1=((psis-zdep(k))/psis)**(1.-1./bee)
             wrz(k)=(-psis/zdep(k))*(bee/(bee-1.))   &
                  *(term1-1.)
             frcsat=1.-zdep(k)/rzdep
             wrz(k)=(1.-frcsat)*wrz(k)+frcsat*1.
          else
             term1=((psis-zdep(k))/psis)**(1.-1./bee)
             term2=((psis-zdep(k)+rzdep)/psis)    &
                  **(1.-1./bee)
             wrz(k)=(-psis/rzdep)*(bee/(bee-1.))  &
                  *(term1-term2)
          endif
       enddo

       !c       Loop over two root zone excess values:
       do irz=1,2
          if(irz.eq.1) rzexc=-0.1*poros
          if(irz.eq.2) rzexc=0.1*poros

          !c       Determine actual root zone excess
          rzact=0.
          do k=1,nar0
             rze(k)=rzexc
             rzw=wrz(k)*poros
             if(rzw+rze(k) .gt. poros) rze(k)=poros-rzw
             if(rzw+rze(k) .lt. 0.) rze(k)=rzw
             rzact=rzact+rze(k)*px0(k)
          enddo
          !c            write(*,*)"  End root zone excess"
          !c       Compute the average timescale

          satfrc=0.
          do k=1,nar0
             if(zdep(k).lt.0.) satfrc=satfrc+px0(k)
          enddo

          sumfrac=0.
          sumz=0.
          do k=1,nar0
             sumz=sumz+zdep(k)*px0(k)
             if(zdep(k) .lt. 1.) frac=frcmax
             if(zdep(k) .ge. 1.) then
                index1=1+int(((zdep(k)*100.)-99)/5.)
                if(index1.gt.nwt) index1 = nwt
                frac=amin1(frc(index1,1),frcmax)
                do i0=2,nrz
                   if(rze(k) .ge. rzexcn(index1,i0))  &
                        frac=amin1(frc(index1,i0),frcmax)
                enddo
             endif
             sumfrac=sumfrac+frac*px0(k)
          enddo
          !c            write(*,*)"  End average time scale"
          catd(idep,irz)=sumdef
          tsc(idep,irz)=sumfrac

       enddo
    enddo

    tsb1=(alog(tsc(2,2))-alog(tsc(1,2)))/(catd(2,2)-catd(1,2))
    tsb2=(alog(tsc(2,1))-alog(tsc(1,1)))/(catd(2,1)-catd(1,1))
    tsa1=alog(tsc(2,2))-tsb1*catd(2,2)
    tsa2=alog(tsc(2,1))-tsb2*catd(2,1)

  END SUBROUTINE TS_PARAM

  !*********************************************************************
  
  SUBROUTINE BASE_PARAM(                                  &
       BEE,PSIS,POROS,COND,               &
       VALX, PX,                          &
       bf1,bf2,bf3,                       &
       taberr1,taberr2,normerr1,normerr2  &
       )

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c                                                                   c
    !c New way to get baseflow: we parametrize the relationship between  c 
    !c catdef and zbar (two parameters bf1 and bf2).                     c
    !c Then, in the LSM/catchment.f/base.f, we use the original relation c 
    !c from TOPMODEL to infer baseflow from catdef and the mean of the   c
    !c topographic index (topmean=bf3, a third parameter).               c
    !c                                                                   c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    INTEGER IDMAX,i1,i2,i,icount

    REAL, intent (in) :: BEE, PSIS,POROS,COND,VALX(NAR),PX(NAR)
    real zbar(nbdep),catdef(nbdep),bflow(nbdep)
    real, intent (out) :: bf1,bf2,bf3,taberr1,taberr2,normerr1,normerr2
    integer :: n,idep
    real suma,sumta,timean

    real catfit(nbdep),bfit(nbdep),dfit(nbdep),catmean,bfmean
    real catref(nbdep),bref(nbdep)
    real err1, err2
    !      logical, intent (in) :: bug

    sumta=0.
    suma=0.
    do n=1,nar
       sumta=sumta+px(n)*valx(n)
       suma=suma+px(n)
    enddo
    timean=sumta/suma
    bf3 = timean

    !c**** Loop over water table depths

    do idep=1,nbdep

       !c           write(*,*) 'idep=',idep

       CALL BASIDEP(                  &
            IDEP,                      &
            BEE,PSIS,POROS,COND,       &
            VALX,PX,TIMEAN,SUMA,       &
            ZBAR,CATDEF,BFLOW)

    enddo


    i1=10   ! zbar= 0 m
    i2=35   ! zbar= 2.5 m

    bf2=zbar(i2)*SQRT(catdef(i1))               &
         /(SQRT(catdef(i2))-SQRT(catdef(i1)))
    bf1=catdef(i1)/(bf2*bf2)

    if (bf1 .le. 0) write(*,*) 'bf1 le 0 for i=',i
    if (bf2 .le. 0) write(*,*) 'bf2 le 0 for i=',i

    !c Errors: Root mean square errors: only for points where catdef GT 0.5mm

    do idep=1,nbdep
       catref(idep)=0.
       bref(idep)=0.
    enddo
    catmean=0.
    bfmean=0.
    icount=0
    do idep=1,nbdep
       if (catdef(idep) .gt. lim) then
          icount=icount+1
          catref(icount)=catdef(idep)
          bref(icount)=bflow(idep)
          catfit(icount)=bf1*(zbar(idep)+bf2)             &
               *(zbar(idep)+bf2)
          dfit(icount)=SQRT(catdef(idep)/bf1)-bf2
          bfit(icount)=cond*exp(-timean-gnu*dfit(icount)) &
               /gnu
          catmean=catmean+catdef(idep)
          bfmean=bfmean+bflow(idep)
       endif
    enddo
    catmean=catmean/icount
    bfmean=bfmean/icount
    if (icount.gt.1) then
       call RMSE(catref,catfit,icount,err1)
       call RMSE(bref,bfit,icount,err2)

       taberr1=err1
       taberr2=err2
       normerr1=err1/catmean
       normerr2=err2/bfmean
    endif
    !c---------------------------------------------------------------------

  END SUBROUTINE BASE_PARAM

  ! ************************************************************************
  
  SUBROUTINE BASIDEP(                          &
       IDEP,                     &
       BEE,PSIS,POROS,COND,      &
       VALX,PX,TIMEAN,SUMA,      &
       ZBAR,CATDEF,BFLOW)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c                                                                      c
    !c This program returns the eight parameters for the areal fractioning  c
    !c                                                                      c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    INTEGER, intent (in) :: idep
    integer nref, nind,nmax,indmin,locmax,shift,ord,locmin,ordref,width,k
    REAL, intent (in) :: BEE, PSIS, POROS, COND,VALX(NAR), PX(NAR), &
         suma,timean
    real :: dx,sumdef,dz
    real, intent (out) ::  catdef(nbdep),bflow(nbdep),zbar(idep)
    real term1,term2,sum
    real zdep(nar),locdef(nar)
    !      logical bug

    !c-------------------------------------------------------------------------
    !c integral(f(x)dx)=1. for a pdf
    !c here px=f(x)dx

    dx=valx(1)-valx(2)

    if (bug) write(*,*) 'IDEP=',IDEP,' dx=',dx, 'gnu=',gnu

    !c the loops over idmax and nbdep are initiated in sta_params4.f

    zbar(idep)=float(idep-10)*slice ! zdep in meters

    !c**** Compute array of water table depths:
    do k=1,nar
       term1=(1/gnu)*(valx(k)-timean)
       zdep(k)=AMAX1(0.,zbar(idep)-term1)
    enddo

    !c variable change must be reflected in dx
    dz=dx/gnu

    if (bug) write(*,*) 'basidep: ok1'

    !c**** Compute array of moisture deficits:
    do k=1,nar
       term1=(psis-zdep(k))/psis
       term1=term1**(1.-1./bee)
       term2=-psis*(bee/(bee-1.))*(term1-1.)
       locdef(k)=zdep(k)-term2
    enddo

    !c**** Add deficits to produce catdef:
    sumdef=0.
    do k=1,nar
       sumdef=sumdef+locdef(k)*px(k)
    enddo
    catdef(idep)=poros*1000.*sumdef/suma

    if (bug) write(*,*) 'basidep: ok2'

    bflow(idep)=cond*exp(-timean-gnu*zbar(idep))/gnu

    if (bug) write(*,*) 'basidep: ok3'

  END SUBROUTINE BASIDEP

  !*****************************************************************************
  
  SUBROUTINE SAT_PARAM(                                                   &
       BEE,PSIS,POROS,COND,                               &
       WPWET,VALX, PX, COESKEW,PFC,                       &
       soildepth,                                         &
       ARS1,ARS2,ARS3,                                    &
       ARA1,ARA2,ARA3,ARA4,                               &
       ARW1,ARW2,ARW3,ARW4,                               &
       taberr1,taberr2,taberr3,taberr4,                   &
       normerr1,normerr2,normerr3,normerr4,               &
       DBG_UNIT)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c                                                                      c
    !c This program returns the eleven parameters for the areal fractioning c
    !c                                                                      c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    INTEGER, intent (in) :: pfc
    REAL, intent (in) :: BEE,PSIS,POROS,COND,WPWET, &
         VALX(NAR), PX(NAR)
    REAL, intent (in) :: soildepth, COESKEW
    REAL, intent (inout) :: ARS1,ARS2,ARS3,                                 &
         ARA1,ARA2,ARA3,ARA4,                               &
         ARW1,ARW2,ARW3,ARW4,                               &
         taberr1,taberr2,taberr3,taberr4,                   &
         normerr1,normerr2,normerr3,normerr4
    INTEGER idep,n,k,i,icount,iref
    integer nar0
    integer nref, nind,nmax,indmin,locmax,shift,ord,locmin
    integer loc1,loc2,loc3,loc0,flag
    REAL VALX0(NAR), PX0(NAR)
    integer :: adjust,loc2save,inc,dec
    real sumta,suma,timean,upval,loval,profdep
    real rjunk,rjunk2
    integer, intent (in), optional :: DBG_UNIT
    real catdef(nbdep),wmin(nbdep),ar1(nbdep),aa(nbdep),aabis(nbdep)
    real ar2(nbdep),ar3(nbdep),swsrf2(nbdep),swsrf3(nbdep),rzeq(nbdep)
    real zbar0,catdef0,wmin0,RZDEP,wminsave(nbdep)

    real x1,x2,x3,x4,w1,w1_0,w2,w3,w4,ref1
    real y0,f1,f2,f3,g1,g2,g3,df,dg,dx,bf,bg,delta,z1,z2

    real nar1(nbdep),nar2(nbdep),nmean2(nbdep),neq(nbdep)
    real shape, nwm, area1,cdi,nar3(nbdep),nmean3
    real err1,err2,err3,err4,sum
    real tabact(nbdep),tabfit(nbdep)

    integer :: mp,isvd,j,first_loop
    !      REAL*8, allocatable :: A(:,:),AP(:,:)
    !      REAL*8, allocatable :: B(:)
    REAL*8, allocatable, target :: A(:,:)
    REAL*8, allocatable, target :: B(:)
    REAL*8, pointer             :: AP(:,:)
    REAL*8, pointer             :: BP(:)
    REAL*8 V(3,3),W(3),ANS(3),sdmax,sdmin,wbrac

    real :: cdcr1,cdcr2,term1,term2,zmet
    logical :: smooth,ars_svd_loop
    logical, parameter ::  bug=.false.
    logical, parameter :: SingValDecomp = .true.
    integer, parameter :: nl=4, nr=4, m=4, NP=50
    real :: savgol_coeff(NP)  
    integer :: savgol_ind(NP)
    integer :: nbdepl,istart

    ref1 = 100.
    !      print *,'PFC', pfc   
    if (bug) write(*,*) 'starting sat_param'

    if(SingValDecomp) then
       savgol_ind(1)=0 
       j=3
       do i=2, nl+1
          savgol_ind(i)=i-j
          j=j+2
       end do

       j=2
       do i=nl+2, nl+nr+1
          savgol_ind(i)=i-j
          j=j+2
       end do
       call savgol(savgol_coeff,nl+nr+1,nl,nr,0,m)
    endif

    profdep = soildepth
    rzdep =grzdep
    profdep=profdep/1000.
    profdep=amax1(1.,profdep)
    if (rzdep .gt. .75*profdep) then
       rzdep=0.75*profdep
    end if

    zmet=profdep
    term1=-1.+((psis-zmet)/psis)**  &
         ((bee-1.)/bee)
    term2=psis*bee/(bee-1)
    cdcr1=1000.*poros*(zmet-(-term2*term1))
    cdcr2=(1-wpwet)*poros*1000.*zmet
    !c mean of the topographic index distribution

    nar0=0
    do n=1,nar
       if (px(n) .ne. 0.) then
          nar0=nar0+1
          valx0(nar0)=valx(n)
          px0(nar0)=px(n)
       endif
    enddo

    sumta=0.
    suma=0.
    do n=1,nar0
       sumta=sumta+px0(n)*valx0(n)
       suma=suma+px0(n)
    enddo
    timean=sumta/suma

    if (bug) write(*,*) 'ok 0: sumta,suma,nar0=',sumta,suma,nar0

    !c**** Loop over water table depths

    do idep=1,nbdep

       CALL FUNCIDEP(                                       &
            NAR0,IDEP,                              &
            BEE,PSIS,POROS,COND,RZDEP,WPWET,        &
            VALX0,PX0,COESKEW,TIMEAN,SUMA,          &
            CATDEF,AR1,WMIN,AA,AABIS,               &
            AR2,AR3,SWSRF2,SWSRF3,RZEQ)             
    enddo

    nbdepl = 100
    if(catdef(50) > cdcr1 + 20.) nbdepl = 50
    if(soildepth > 6500.)  nbdepl = nbdep

    if (bug) write(*,*) 'funcidep loop ok'

    !c**** for wmin's adjustment, we need an estimate of its limit toward INF
    adjust =0
    ZBAR0=10.
    CALL FUNCZBAR(                                        &
         NAR0,ZBAR0,                           &
         BEE,PSIS,POROS,COND,RZDEP,WPWET,      &
         VALX0,PX0,COESKEW,TIMEAN,SUMA,        &
         CATDEF0,WMIN0)

    if (bug) write(*,*) 'funczbar ok'

    if (wmin0 == 0.9999900) then
       do idep=1,nbdep-1
          if(catdef(idep).le.cdcr1+10.) then
             if((wmin(idep) - wmin(idep +1)) > -0.01) then 
                wmin0=wmin(idep)
             endif
          endif
       enddo
       wmin0 = 0.1*(nint(wmin0*100000.)/10000) -0.02		    
    endif

    if(present(dbg_unit)) then
       write (dbg_unit,*) nbdep,nbdepl,wmin0,cdcr1,cdcr2
       write (dbg_unit,*) catdef
       write (dbg_unit,*) ar1
       write (dbg_unit,*) wmin
    endif

    !c**** AR1 adjustment: 3 points + limit in INF = 0.

    if (bug) write(*,*) 'STARTING AR1'

    ! Singular value decomposition               
    loc1=1
    loc3=nbdepl
    loc2=loc3	 

    do idep = 1,loc2
       if(ar1(idep) < 1.e-10) then	
          loc3 = idep - 1
          exit
       endif
    end do

    first_loop = 0      
    ars_svd_loop = .TRUE.
    DO while (ars_svd_loop)

       first_loop = first_loop + 1
       mp = loc3-loc1+1

       allocate(A(mp,3))
       allocate(AP(mp,3))
       allocate(B(mp))

       a=0.
       ap=0.
       b=0.
       v=0.
       w=0.
       ans=0.

       do isvd=loc1,loc3
          A(isvd-loc1+1,1)=catdef(isvd)
          A(isvd-loc1+1,2)=-catdef(isvd)*ar1(isvd)
          A(isvd-loc1+1,3)=-ar1(isvd)*((catdef(isvd))**2.)
          B(isvd-loc1+1)=ar1(isvd)-1.
       end do

       ap = a
       call svdcmp(ap,mp,3,w,v)
       sdmax=0.
       do j=1,3
          if(w(j).gt.sdmax)sdmax=w(j)
       end do
       sdmin=sdmax*1.0e-6
       do j=1,3
          if(w(j).lt.sdmin)w(j)=0.
       end do

       call svbksb(ap,w,v,mp,3,b,ans)

       ars1 = real(ans(1))
       ars2 = real(ans(2))
       ars3 = real(ans(3))  

       flag=0
       call curve1(ars1,ars2,ars3,cdcr2,flag)
       deallocate (A, AP, B)

       IF(FLAG == 1) THEN
          LOC3 = NBDEP
          LOC1 =1
          IF(first_loop > 1) ars_svd_loop=.FALSE.
       ELSE
          ars_svd_loop=.FALSE.  
       ENDIF
    END DO

    IF (FLAG.EQ.1) then

       flag=0
       loc1=1
       do idep=1,nbdepl
          if (catdef(idep) .le. 20.) loc1=idep
       enddo

       loc3=1
       do idep=1,nbdepl -1
          if ((ar1(idep) >= 0.0001).and.(catdef(idep) <= cdcr1)) loc3=idep + 1
       enddo

       if (loc3.le.loc1+1) then
          loc1=MIN(loc3-4,loc1-4)
          loc1=MAX(1,loc1)
       endif

       !c below is what was used for no regression, but it's not equivalent to the 
       !c IDL program
       loc2=AINT(loc1-1+(loc3-loc1)*3./5.)+1

       w1=ar1(loc1)
       w2=ar1(loc2)
       w3=ar1(loc3)

       if(w3.eq.0.)then
95        loc3=loc3-1
          if(loc3.eq.loc2)loc2=loc2-1
          w3=ar1(loc3)
          w2=ar1(loc2)
          if(w3.eq.0.)goto 95
       endif
       w4=0.

       if((loc1.ge.loc2).or.(loc2.ge.loc3))then
          loc1=10
          loc2=14
          loc3=18
       endif

115    x1=catdef(loc1)
       x2=catdef(loc2)
       x3=catdef(loc3)
       w1=ar1(loc1)
       w2=ar1(loc2)
       w3=ar1(loc3)

       if (bug) then
          write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
          write(*,*) 'x1,x2,x3=',x1,x2,x3
          write(*,*) 'w1,w2,w3=',w1,w2,w3
       endif

       y0=w4
       f1=(1.-w1)/(w1-y0)/x1
       f2=(1.-w2)/(w2-y0)/x2
       f3=(1.-w3)/(w3-y0)/x3
       g1=(1.-y0)/(w1-y0)
       g2=(1.-y0)/(w2-y0)
       g3=(1.-y0)/(w3-y0)
       df=f2-f1
       dg=g2-g1
       dx=x2-x1
       bf=f1-x1*df/dx
       bg=g1-x1*dg/dx

       ars1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
       ars2 = bf+ars1*bg
       ars3 = (df+ars1*dg)/dx

       delta=ars2*ars2-4*ars3
       upval=1.+200.*ars1
       loval=1.+200.*ars2+40000.*ars3
       z1=0.
       z2=0.

       if (delta .ge. 0.) then !if 8
          z1=(-ars2-SQRT(delta))/2./ars3
          z2=(-ars2+SQRT(delta))/2./ars3 
       endif

       if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
            (z2 .gt. 0. .and. z1 .lt. cdcr1) .or.  &
            ((upval/loval).lt.-.01)) then   !if 7
          z1=0.
          z2=0.
          if (loc1 .eq. 10) then 
             loc1=1
1         else  
             loc1=1
             do idep=1,nbdepl
                if (catdef(idep) .gt. 60.) then
                   loc1=idep
                   if(loc1.ge.loc3-1)then
                      !                           write(*,*)'Loc1 exceeded loc3 in 2nd attempt'
                      loc1=loc3-5
                   endif
                   goto 46
                endif
             enddo
          endif
46        loc2=loc1+AINT(float(loc3-loc1)*3./5.)+1
          if(loc2.ge.loc3)loc2=loc3-1
          loc2save=loc2
          INC=1
          DEC=0

47        w1=ar1(loc1)
          w2=ar1(loc2)
          x1=catdef(loc1)
          x2=catdef(loc2)

          if (bug) then
             write(*,*) 'z1,z2=',z1,z2,' -> ar1, 2nd try'
             write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
             write(*,*) 'x1,x2,x3=',x1,x2,x3
             write(*,*) 'w1,w2,w3=',w1,w2,w3
          endif

          f1=(1.-w1)/(w1-y0)/(x1 + 1.e-20)
          f2=(1.-w2)/(w2-y0)/(x2 + 1.e-20)
          g1=(1.-y0)/(w1-y0 + 1.e-20 )
          g2=(1.-y0)/(w2-y0 + 1.e-20)
          df=f2-f1
          dg=g2-g1
          dx=x2-x1
          bf=f1-x1*df/dx
          bg=g1-x1*dg/dx

          ars1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx  + 1.e-10)
          ars2 = bf+ars1*bg
          ars3 = (df+ars1*dg)/dx
          delta=ars2*ars2-4*ars3
          upval=1.+200.*ars1
          loval=1.+200.*ars2+40000.*ars3

          if (delta .ge. 0.) then   !if 6
             z1=(-ars2-SQRT(delta))/2./ars3
             z2=(-ars2+SQRT(delta))/2./ars3
          end if

          if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
               (z2 .gt. 0. .and. z1 .lt. cdcr1) .or.   &
               ((upval/loval).lt.-.01)) then  !if 5
             !c Sarith ---
             z1=0.
             z2=0.
             IF(INC.EQ.1)loc2=loc2+1
             IF(DEC.EQ.1)LOC2=LOC2-1
             if(inc.eq.1)then   !if 4
                if(loc2.ge.loc3)then   !if 3
                   !                     WRITE(*,*)'INCREASING LOC2 FAILED'
                   INC=0
                   DEC=1
                   loc2=loc2save
                else
                   adjust=ADJUST+1
                   goto 47
                end if    !if 3
             endif    !if 4

             if(dec.eq.1)then   !if 2  
                if(loc2.eq.loc1)then  !if 1
                   !                     WRITE(*,*)'Decreasing too failed'
                   INC=1
                   DEC=0
                   ars1=9999. !ars1old
                   ars2=9999. !ars2old
                   ars3=9999. !ars3old
                   !                     write(*,*) 'AR1: PROBLEM for pfc=',pfc
                else
                   adjust=ADJUST+1
                   !c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                   goto 47
                end if   !if 1
             endif  !if 2 
          endif     !if 5
          !c               endif    !if 6
       endif            !if 7

       !c         endif  !if 8
       flag=0
       call curve1(ars1,ars2,ars3,cdcr2,flag)

       IF (FLAG.EQ.1)then
          !            WRITE(*,*)'Curve problem in the catchment pfc=',pfc
          ars1=9999. 
          ars2=9999. 
          ars3=9999. 
          !                     write(*,*) 'Pick values from icatch-1'
          flag=0
       end if
    endif

    adjust=0

    if (bug) write(*,*) 'ar1 adjustment ok'

    !c**** WMIN adjustment: 3 points + limit in INF = wmin0

    if (bug) write(*,*) 'STARTING WMIN'

    w4=wmin0
    y0=w4

    !         write(*,*) 'wmin=',(wmin(idep),idep=1,50)

    loc1=1
    do idep=1,nbdepl
       if (catdef(idep) <= 10.) loc1=idep
    enddo

    loc3=1
    do idep=1,nbdepl - 2
       if ((wmin(idep) >= wmin0).and.(catdef(idep) <= cdcr1)) loc3=idep + 2
    enddo

    loc2=loc1 + 2
    do idep=1,nbdepl -1 
       if ((wmin(idep) >= wmin0).and.(catdef(idep) <= cdcr1/2.))loc2=idep + 1
    enddo

    !c For global catch         
    INC=1
    DEC=0

    if(loc3.eq.loc2)loc2=loc2-2
    if(loc2 <= loc1) loc1= loc1-2
44  loc2save=loc2
    if(loc1 < 1) then
       loc1 =1
       loc2 =2 
       loc3 =3
    endif

    w1=wmin(loc1)
    w2=wmin(loc2)
    w3=wmin(loc3)
    x1=catdef(loc1)
    x2=catdef(loc2)
    x3=catdef(loc3)

    if (bug) then
       write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
       write(*,*) 'x1,x2,x3=',x1,x2,x3
       write(*,*) 'w1,w2,w3,w4=',w1,w2,w3,w4
    endif

    f1=(1.-w1)/(w1-y0)/x1
    f2=(1.-w2)/(w2-y0)/x2
    f3=(1.-w3)/(w3-y0)/x3
    g1=(1.-y0)/(w1-y0)
    g2=(1.-y0)/(w2-y0)
    g3=(1.-y0)/(w3-y0)
    df=f2-f1
    dg=g2-g1
    dx=x2-x1
    bf=f1-x1*df/dx
    bg=g1-x1*dg/dx

    arw1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
    arw2 = bf+arw1*bg
    arw3 = (df+arw1*dg)/dx
    arw4 = y0

    !c wmin=arw4+(1.-arw4)*(1.+arw1*catdef(idep))
    !c     /(1.+arw2*catdef(idep)+arw3*catdef(idep)*catdef(idep))
    !c we want to check the roots of the denominator

    delta=arw2*arw2-4*arw3

    if (delta .ge. 0.) then !if 8

       z1=(-arw2-SQRT(delta))/2./arw3
       z2=(-arw2+SQRT(delta))/2./arw3

       if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.           &
            (z2 .gt. 0. .and. z1 .lt. cdcr1)) then !if 7

          w1_0=w1
          w1=(1.+w1_0)/2.
          x1=x1/4.

          !                  if (gnu .eq. 3.26/1.5) then 
          !                        w1=(1.+w1_0)/3.               ! already difficult
          !                        w3=wmin(nint(cdcr1))                ! with gnu=3.26
          !                        x3=catdef(nint(cdcr1))
          !                        f3=(1.-w3)/(w3-y0)/x3
          !                        g3=(1.-y0)/(w3-y0)
          !                  endif

          f1=(1.-w1)/(w1-y0)/x1
          g1=(1.-y0)/(w1-y0)
          df=f2-f1
          dg=g2-g1
          dx=x2-x1
          bf=f1-x1*df/dx
          bg=g1-x1*dg/dx

          if (bug) then
             write(*,*) 'z1,z2=',z1,z2,' -> wmin, 2nd try'
             write(*,*) 'loc1,loc2,loc3=',loc1,loc2,loc3
             write(*,*) 'x1,x2,x3=',x1,x2,x3
             write(*,*) 'w1,w2,w3=',w1,w2,w3
             write(*,*) 'wmin0=',wmin0
          endif

          arw1 = -(f3-bf-x3*df/dx)/(g3-bg-x3*dg/dx + 1.e-10)
          arw2 = bf+arw1*bg
          arw3 = (df+arw1*dg)/dx
          arw4 = y0

          delta=arw2*arw2-4*arw3

          if (delta .ge. 0.) then  !if 6 
             z1=(-arw2-SQRT(delta))/2./arw3
             z2=(-arw2+SQRT(delta))/2./arw3

             if ((z1 .gt. 0. .and. z1 .lt. cdcr1) .or.         &
                  (z2 .gt. 0. .and. z1 .lt. cdcr1)) then    !if 5  
                !c Sarith ---
                IF(INC.EQ.1)loc2=loc2+1
                IF(DEC.EQ.1)LOC2=LOC2-1
                if(inc.eq.1)then   !if 4
                   if(loc2.eq.loc3)then   !if 3
                      !                     WRITE(*,*)'INCREASING LOC2 FAILED: WMIN'
                      INC=0
                      DEC=1
                      loc2=loc2save
                   else
                      adjust=ADJUST+1
                      !c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                      goto 44
                   end if    !if 3
                endif    !if 4
                if(dec.eq.1)then   !if 2  
                   if(loc2.eq.loc1)then  !if 1
                      !                     WRITE(*,*)'Decreasing too failed: WMIN'
                      INC=1
                      DEC=0

                      arw1=9999. 
                      arw2=9999. 
                      arw3=9999. 
                      arw4=9999. 

                   else
                      adjust=ADJUST+1
                      !c                        write(*,*)'ADJUSTING AR1 CYCLE =',ADJUST
                      goto 44
                   end if   !if 1
                endif  !if 2                         
             endif !if 5     
          endif   !if 6

       endif !if 7
    endif !if 8
    adjust=0
    !         endif ! pfc=12821
    flag=0

    call curve2(arw1,arw2,arw3,arw4,cdcr1,WPWET,flag)

    IF (FLAG.EQ.1) THEN
       arw1=9999. !arw1old
       arw2=9999. !arw2old
       arw3=9999. !arw3old
       arw4=9999. !arw4old
       flag=0
    endif

    if(arw1==9999.) then 
       ! Singular Value Decomposition

       w4=wmin0
       y0=w4

       loc1=1
       loc3=nbdepl

       mp = loc3-loc1+1   

       if(mp.lt.3)then

          write(*,*)'WMIN Note: not sufficient points MP = ',mp
          print *,w4,cdcr1,catdef(loc3),wmin(loc3)
          arw1 = 9999.
          arw2 = 9999.
          arw3 = 9999.
          arw4 = 9999.
       else

          mp = 1
          istart =1
          w4 = wmin(istart)

          if(w4 <=0) then
             do idep=2,nbdepl
                if(wmin(idep) > 0.) istart = idep
                if(wmin(idep) > 0.) exit
             enddo
          endif

          w4 = wmin(istart)

          do idep=istart+1,nbdepl
             !	       if(wmin(idep).lt.w4) then
             if((wmin(idep) - w4).lt.0.0005) then
	       	w4 = wmin(idep)
		mp = mp +1
             endif
          enddo
          loc3 = mp   
          allocate(A(mp,3))
          allocate(AP(mp,3))
          allocate(B(mp))
          allocate(BP(mp))             
          smooth = .false.
          do idep=istart,nbdepl-1
             if(catdef(idep).le.cdcr1+10.) then
                if((wmin(idep) - wmin(idep +1)) < -0.01) smooth = .true.   
             endif
          enddo
          if(smooth) then
             wminsave = wmin
             ! Apply filter to input data
             do i=istart, nbdepl-nr
                wmin(i)=0.
                do j=1, nl+nr+1
                   if (i+savgol_ind(j).gt.0) then  !skip left points that do not exist
                      wmin(i)=wmin(i)+savgol_coeff(j)*wminsave(i+savgol_ind(j))
                   endif
                end do
             enddo
             wmin (istart:istart+4) = wminsave (istart:istart+4)

          endif

          j = 1
          w4 = wmin(istart)
          do isvd=1,size(wmin)
             if (j <= mp) then 
                if(isvd == 1) then
                   wbrac=(wmin(isvd + istart -1)-y0)/(1.-y0 + 1.e-20)
                   A(j,1)=catdef(isvd + istart -1)
                   A(j,2)=-catdef(isvd + istart -1)*wbrac
                   A(j,3)=-wbrac*((catdef(isvd + istart -1))**2.)
                   B(j)=wbrac-1.
                   j = j + 1
                else
                   if((wmin(isvd + istart -1).lt.w4).and.(wmin(isvd + istart -1).gt.y0)) then
                      wbrac=(wmin(isvd + istart -1)-y0)/(1.-y0 + 1.e-20)
                      A(j,1)=catdef(isvd + istart -1)
                      A(j,2)=-catdef(isvd + istart -1)*wbrac
                      A(j,3)=-wbrac*((catdef(isvd + istart -1))**2.)
                      B(j)=wbrac-1.
                      w4 = wmin(isvd + istart -1)
                      j = j + 1 
                   endif
                endif
             endif
          end do

          j = j -1 
          mp = j
          ap => a (1:j,:)
          bp => b (1:j)
          ap(j,1) = catdef(nbdep)
          ap(j,2) = 0.
          ap(j,3) = 0.
          bp (j) = -1.

          call svdcmp(ap,mp,3,w,v)

          sdmax=0.
          do j=1,3
             if(w(j).gt.sdmax)sdmax=w(j)
          end do

          sdmin=sdmax*1.0e-6
          do j=1,3
             if(w(j).lt.sdmin)w(j)=0.
          end do

          call svbksb(ap,w,v,mp,3,bp,ans)

          arw1 = real(ans(1))
          arw2 = real(ans(2))
          arw3 = real(ans(3))
          arw4 = y0

          !c wmin=arw4+(1.-arw4)*(1.+arw1*catdef(idep))
          !c     /(1.+arw2*catdef(idep)+arw3*catdef(idep)*catdef(idep))
          !c we want to check the roots of the denominator

          adjust=0         
          flag=0

          call curve2(arw1,arw2,arw3,arw4,cdcr1,WPWET,flag)

          IF (FLAG.EQ.1) THEN
             !            WRITE(*,*)'Curve2 problem in the catchment:pfc=',pfc

             arw1 = 9999.
             arw2 = 9999.
             arw3 = 9999.
             arw4 = 9999.

             flag=0
          end if
          deallocate (A,  B )
          NULLIFY    (AP, BP)
       end if
    endif

    if(present(dbg_unit)) then
       write (dbg_unit,*) ars1,ars2,ars3
       write (dbg_unit,*) arw1,arw2,arw3,arw4 
    endif

    if (bug) write(*,*) 'wmin adjustment ok'

    !c**** SHAPE PARAMETER ADJUSTMENT: with a straight if coeskew > 0.25
    !c                                 with 2 segments if not

    if (bug) write(*,*) 'STARTING SHAPE'

    x3=catdef(nbdepl)
    w3=aa(nbdepl)
    x1=0.

    if (coeskew .lt. 0.25) then
       w1=0.1
       loc2=20
       do idep=1,nbdepl
          if (catdef(idep) .gt. ref1) then
             loc2=idep
             goto 45
          endif
       enddo
45     x2=catdef(loc2)
       w2=aabis(loc2)
       ara1 = (w1-w2)/(x1-x2)
       ara2 = w1-ara1*x1
       ara3 = (w2-w3)/(x2-x3)
       ara4 = w2-ara3*x2
    else
       w1=1.
       x2=x1
       w2=w1
       ara3 = (w2-w3)/(x2-x3)
       ara4 = w2-ara3*x2
       ara1 = ara3
       ara2 = ara4       
    endif

    if (bug) write(*,*) 'x1,w1,x2,w2,x3,w3',x1,w1,x2,w2,x3,w3

    !**** RMSE checking: on ar1, ar2, swsrf2 and rzeq

    do idep=1,nbdepl
       if(catdef(idep) <= cdcr1) then
          nar1(idep)=AMIN1(1.,AMAX1(0.,(1.+ars1*catdef(idep)) &
               /(1.+ars2*catdef(idep)                          &
               +ars3*catdef(idep)*catdef(idep))))                    

          nwm=AMIN1(1.,AMAX1(0.,arw4+(1.-arw4)*                &
               (1.+arw1*catdef(idep))                          &
               /(1.+arw2*catdef(idep)                          &
               +arw3*catdef(idep)*catdef(idep))))

          !c we have to first determine if there is one or two segments
          if (ara1 .ne. ara3) then
             cdi=(ara4-ara2)/(ara1-ara3)
          else
             cdi=0.
          endif

          if (catdef(idep) .ge. cdi) then
             shape=ara3*catdef(idep)+ara4
          else
             shape=ara1*catdef(idep)+ara2
          endif
          shape =AMIN1(40.,shape)
          area1=exp(-shape*(1.-nwm))*(shape*(1.-nwm)+1.)

          !c the threshold for truncation problems is higher than the "usual"
          !c E-8 to E-10, because it plays together with the uncertainties coming 
          !c from the approximation of the parameters nwm, nar1 and shape.
          if (area1 .ge. 1.-1.E-8) then
             nar1(idep)=1.
             nar2(idep)=0.
             nar3(idep)=0.
             nmean2(idep)=0.  
             nmean3=0.       
             neq(idep)=1.
          else

             if (nwm .gt. wpwet) then
                nar2(idep)=1.-nar1(idep)
             else
                nar2(idep)=AMAX1(0.,((shape*(wpwet-nwm)+1.)        &
                     *exp(-shape*(wpwet-nwm))                       &
                     - (shape*(1.-nwm)+1.)*exp(-shape*(1.-nwm)))   &
                     * (1.-nar1(idep))/(1.-area1))                  
             endif

             nar3(idep)=1.-nar1(idep)-nar2(idep)                          

             if (nar3(idep) .lt. 1.E-8) then ! for nwm le wpwet           

                nmean2(idep)=AMAX1(0.,AMIN1(1.,(nwm + 2./shape +    &
                     shape*exp(-shape*(1.-nwm))*                    &
                     (nwm+nwm/shape-1.-2./shape-2./(shape*shape))) &
                     /(1.-area1)))
                nmean3=0.

             else

                !c WARNING: I think the two values below are false. 
                !c But it is never used in this context, because nwm > wpwet !!
                nmean2(idep)=AMAX1(0.,AMIN1(1.,-shape*(exp(-shape*&
                     (wpwet-nwm))* (nwm*wpwet                      &
                     +nwm/shape-wpwet*wpwet                        &
                     -2.*wpwet/shape-2./(shape*shape))             &
                     - exp(-shape*(1.-nwm))*                       &
                     (nwm+nwm/shape-1.-2./shape-2./(shape*shape)))& 
                     * (1.-nar1(idep))/(1.-area1) / (nar2(idep)+1.e-20)))        

                nmean3=AMAX1(0.,AMIN1(1.,(nwm+2./shape +          &
                     shape*exp(-shape*(wpwet-nwm))*                &
                     (nwm*wpwet+nwm/shape-wpwet                   &
                     *wpwet-2.*wpwet/shape                         &
                     -2./(shape*shape))) * (1.-nar1(idep))         &
                     /(1.-area1)/(nar3(idep) + 1.e-20)))                       
             endif

             neq(idep)=nar1(idep)+nar2(idep)*nmean2(idep)          &
                  +nar3(idep)*nmean3

             if (area1 .ge. 1.-1.E-5) then
                nmean2(idep)=1.  
                nmean3=0.       
                neq(idep)=1.
             endif

          endif
       endif
    enddo

    if (bug) write(*,*) 'shape adjustment ok'
    !c
    !c RMSE

    !c ERR1
    icount=0
    iref=0
    sum=0.
    do i=1,nbdepl
       if(catdef(i) <= cdcr1) then
          tabact(i)=0.
          tabfit(i)=0.
       endif
    enddo

    do i=1,nbdepl
       if(catdef(i) <= cdcr1) then 
          if (catdef(i) .gt. lim) then
             icount=icount+1
             sum=sum+ar1(i)
             tabfit(icount)=nar1(i)
             tabact(icount)=ar1(i)
          endif
       endif
    enddo

    if(icount.gt.1) then
       sum=sum/icount
       call RMSE(tabact,tabfit,icount,err1)
       taberr1=err1
       normerr1=err1/sum
    endif
    !c ERR2
    icount=0
    iref=0
    sum=0.
    do i=1,nbdepl
       if(catdef(i) <= cdcr1) then
          tabact(i)=0.
          tabfit(i)=0.
       endif
    enddo

    do i=1,nbdepl
       if(catdef(i) <= cdcr1) then
          if (catdef(i) .gt. lim) then
             icount=icount+1
             sum=sum+ar2(i)
             tabfit(icount)=nar2(i)
             tabact(icount)=ar2(i)
          endif
       endif
    enddo

    if(icount.gt.1) then
       sum=sum/icount
       call RMSE(tabact,tabfit,icount,err2)
       taberr2=err2
       normerr2=err2/sum
    endif

    !c ERR3
    icount=0
    iref=0
    sum=0.
    do i=1,nbdep
       if(catdef(i) <= cdcr1) then
          tabact(i)=0.
          tabfit(i)=0.
       endif
    enddo

    do i=1,nbdepl
       if(catdef(i) <= cdcr1) then
          if (catdef(i) .gt. lim) then
             icount=icount+1
             sum=sum+swsrf2(i)
             tabfit(icount)=nmean2(i)
             tabact(icount)=swsrf2(i)
          endif
       endif
    enddo

    if(icount.gt.1) then
       sum=sum/icount
       call RMSE(tabact,tabfit,icount,err3)
       taberr3=err3
       normerr3=err3/sum
    endif
    !c ERR4
    icount=0
    iref=0
    sum=0.
    do i=1,nbdepl
       tabact(i)=0.
       tabfit(i)=0.
    enddo

    do i=1,nbdepl
       if(catdef(i) <= cdcr1) then
          if (catdef(i) .gt. lim) then
             icount=icount+1
             sum=sum+rzeq(i)
             tabfit(icount)=neq(i)
             tabact(icount)=rzeq(i)
          endif
       endif
    enddo

    if(icount.gt.1) then 
       sum=sum/icount
       call RMSE(tabact,tabfit,icount,err4)
       taberr4=err4
       normerr4=err4/sum
    endif
  END SUBROUTINE SAT_PARAM
  !
  
  ! ******************************************************************
  
  !c
  SUBROUTINE CURVE1(ars1,ars2,ars3,cdcr2,flag)
    REAL ars1,ars2,ars3,y,x,yp,cdcr2
    INTEGER i,flag
    !c
    yp=1.
    if (abs(ars1+ars2+ars3).le.1.e25) then
       do i=0,CEILING(cdcr2)
          x=float(i)
          if(x > cdcr2) x = cdcr2
          y=(1.+ars1*x)/(1.+ars2*x+ars3*x*x + 1.e-20)
          if((y.gt.0.0).and.(((yp -y) .lt. -1.e-4).or.(y.gt.1.)))then
             flag=1
             goto 99
          endif
          yp=y
       end do
99     continue
    else
       flag=1
    endif

  end SUBROUTINE CURVE1


  ! ******************************************************************

  SUBROUTINE CURVE2(arw1,arw2,arw3,arw4,cdcr1,WPWET,flag)
    REAL arw1,arw2,arw3,arw4,y,x,yp,cdcr1, wpwet
    INTEGER i,flag
    !c
    yp=1.
    if (abs(arw1+arw2+arw3+arw4).le.1.e25) then
       do i=0,CEILING(cdcr1)
          x=float(i)
          if(x > cdcr1) x = cdcr1
          y=arw4+(1.-arw4)*(1.+arw1*x)/(1.+arw2*x+arw3*x*x + 1.e-20)
          if ((y .lt. wpwet).or.((yp -y) .lt. -1.e-4).or.(y.gt.1.)) then
             flag=1
             goto 99
          endif
          yp=y
       end do
99     continue
    else
       flag=1
    endif
  end SUBROUTINE CURVE2

  ! ******************************************************************
  
  subroutine tgen (                         &
       TOPMEAN,TOPVAR,TOPSKEW,               &
       STO,ACO,COESKEW)

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                         c
    ! The difference between tgen4 and tgen3 is that tgen4 deals with arrays  c
    ! of topmean, topvar and topskew and 2-dim arrays of st and ac.           c
    !                                                                         c
    ! This routine determine the theoretical gamma distribution for the       c
    ! soil-topographic indexes (Sivapalan et al., 1987), knowing the three    c
    ! first moments, the min and the max of the observed topographic indexes  c
    ! in a given catchment.                                                   c
    !                                                                         c
    ! Routine from Dave Wolock.                                               c
    ! Modified by Agnes (11-06-98): we don't use min and max anymore, and     c
    ! this strongly improves the behavior for negative skewnesses. It also    c
    ! improves in general the matching of the moments.                        c
    !                                                                         c
    ! We also add a correction on the skewness to have gamma distributions    c
    ! that start and end from the x-axis. It is based on the fact that if     c
    ! TOPETA=1, the gamma is an exponential distribution, and if TOPETA<1,    c
    ! then the gamma distribution increases towards the infinite when x       c
    ! decreases towards 0.                                                    c
    ! To eliminate some numerical pb due to teh discretization of the gamma   c
    ! distribution, we choose skewness=MAX(MIN(1.9, skewness),-1.6)           c
    !                                                                         c
    ! WE MAY NEED TO COMPUTE IN DOUBLE RESOLUTION !!!! BECAUSE OF THE SMALL   c
    ! BIN WIDTH
    !                                                                         c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

    real, parameter :: VALMAX=50.
    REAL, intent (in) :: TOPMEAN,TOPVAR,TOPSKEW
    REAL, intent (out) :: COESKEW
    REAL, dimension (NAR), intent (out) :: STO,ACO

    INTEGER I
    REAL ST(NAR),AC(NAR)
    REAL TOPETA,TOPLAM,TOPSCAL,GAMLN,SCALE,ACLN
    real cumac, cum2,cum3

    !-------------------------------------------------------------------------

    ! topmean is the mean of the ln(a/tanB) distribution
    ! topvar is the variance (2nd moment centered around the mean) of the ...
    ! topskew is the skew (3rd moment centerd around the mean) of the ...
    ! compute the coefficient of skew or skewness (coeskew)

    COESKEW=TOPSKEW/TOPVAR**1.5
    if (coeskew .ge. 0.) then
       COESKEW=AMAX1(0.005, AMIN1(1.9, COESKEW))
    else
       COESKEW=AMAX1(-1.6, AMIN1(-0.005, COESKEW))
    endif

    ! compute the gamma parameters, eta (topeta) and lambda (toplam), and topscal
    ! which is the translation parameter

    TOPETA=4./COESKEW**2
    TOPLAM=SQRT(TOPETA)/SQRT(TOPVAR)
    TOPSCAL=TOPMEAN-TOPETA/TOPLAM

    ! evaluate the gamma function

    CALL GAMMLN (TOPETA,GAMLN)

    CUMAC=0.0

    ! compute the frequency distribution of ln(a/tanB)
    ! st(i) are the values of ln(a/tanB)
    ! ac(i) are the relative frequency values (they should sum to 1)

    DO I=1,NAR

       ST(I)=(FLOAT(I)-0.95)*(VALMAX-TOPSCAL)/FLOAT(NAR)+TOPSCAL
       SCALE=ST(I)-TOPSCAL

       ! below is the logarithmic form of the gamma distribution; this is required 
       ! because the numerical estimate of the logarithm of the gamma function 
       ! is more stable than the one of the gamma function.

       ACLN=TOPETA*ALOG(TOPLAM)+(TOPETA-1.)*ALOG(SCALE)  &
            -TOPLAM*SCALE-GAMLN

       IF(ACLN.LT.-10.) THEN
          AC(I)=0.
       ELSE
          AC(I)=EXP(ACLN)
       ENDIF

       CUMAC=CUMAC+AC(I)

    ENDDO

    ! we want the relative frequencies to sum 1.

    IF (CUMAC.eq.0.) THEN
       !            write(*,*) 'distrib sum=',CUMAC
       stop
    endif
    CUM2=0.
    DO I=1,NAR
       AC(I) = AC(I) / CUMAC
       CUM2=CUM2+AC(I)
    ENDDO

    ! if the real distribution of the topographic indices is negativeley skewed, 
    ! we symetrize the gamma distribution (depending on coeskew**2 and always 
    ! positively skewed), centering on topmean, which preserves topmean and
    ! topvar, and re-establishes a negative skewness.

    IF (COESKEW.LT.0.) then

       do i=1,nar
          STO(I)=2.*TOPMEAN-ST(I)
          ACO(I)=AC(I)

       enddo
    ELSE
       !            if (n .eq. idmax) then
       !               write(*,*) 'last catchment'
       !            endif
       do i=1,nar
          STO(I)=ST(-I+NAR+1)
          ACO(I)=AC(-I+NAR+1)
       enddo
    ENDIF

    !         sum=0.
    !         do i=1,nar
    !            sum=sum+sto(i)*aco(i)
    !         end do

    !         sum=0.
    !         do i=1,nar
    !            sum=sum+aco(i)
    !         end do


  END subroutine tgen

  ! ********************************************************************
  
  SUBROUTINE GAMMLN (XX,GAMLN)

    DOUBLE PRECISION :: COF(6),STP,HALF,ONE,FPF,X,TMP,SER
    REAL, intent(in) :: XX 
    REAL, intent(out) :: GAMLN
    integer :: j

    DATA COF /76.18009173D0,-86.50532033D0,24.01409822D0,    &
         -1.231739516D0,.120858003D-2,-.536382D-5/
    STP = 2.50662827465D0
    HALF= 0.5D0
    ONE = 1.0D0
    FPF = 5.5D0

    X=XX-ONE
    TMP=X+FPF
    TMP=(X+HALF)*LOG(TMP)-TMP
    SER=ONE

    DO  J=1,6
       X=X+ONE
       SER=SER+COF(J)/X
    END DO

    GAMLN=TMP+LOG(STP*SER)

  END SUBROUTINE GAMMLN

  ! ********************************************************************
  
  SUBROUTINE FUNCIDEP(                                         &
       NAR0,IDEP,                              &!I
       BEE,PSIS,POROS,COND,RZDEP,WPWET,        &!I
       VALX,PX,COESKEW,TIMEAN,SUMA,            &!I
       CATDEF,AR1,WMIN,AA,AABIS,               &!O
       AR2,AR3,SWSRF2,SWSRF3,RZEQ)             

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c                                                                      c
    !c This program returns the eight parameters for the areal fractioning  c
    !c                                                                      c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    integer, intent (in) :: NAR0,idep 
    REAL, intent (in) :: BEE, PSIS, POROS, COND, RZDEP, WPWET, COESKEW
    REAL, intent (inout) ::  VALX(NAR), PX(NAR),TIMEAN,SUMA
    !      logical, intent(in) :: bug
    real, dimension (nbdep), intent (inout) :: CATDEF,AR1,WMIN,AA,   &
         AABIS,AR2,AR3,SWSRF2,SWSRF3,RZEQ 
    INTEGER :: width, nref, nind,nmax,indmin,locmax,shift,ord,locmin,ordref
    integer :: indimax10,indmin0,k,n,n1,n2
    real dx,zbar

    real test,term1,term2,sum
    real zdep(nar),locdef(nar),wrz(nar),frcunsat
    real valtest(nbdep,nar),ptest(nbdep,nar),denstest(nbdep,nar)
    real dtest(nbdep,nar),cump
    real x1,x2,y1,y2,wa,wb
    real densaux(nar),densaux2(nar),densmax,aux10
    real :: dz, sumdef
    !c-------------------------------------------------------------------------

    !c integral(f(x)dx)=1. for a pdf
    !c here px=f(x)dx
    dx=valx(1)-valx(2)

    if (bug) write(*,*) 'IDEP=',IDEP,' dx=',dx

    !c the loops over idmax and nbdep are initiated in sta_params4.f

    zbar=float(idep-10)*slice ! zdep in meters

    !c**** Compute array of water table depths:
    do k=1,nar0
       term1=(1/gnu)*(valx(k)-timean)
       zdep(k)=AMAX1(0.,zbar-term1)
    enddo

    !c variable change must be reflected in dx
    dz=dx/gnu

    if (bug) write(*,*) 'funcidep: ok1'

    !c**** Compute array of moisture deficits:
    do k=1,nar0
       term1=(psis-zdep(k))/psis
       term1=term1**(1.-1./bee)
       term2=-psis*(bee/(bee-1.))*(term1-1.)
       locdef(k)=zdep(k)-term2
    enddo

    !c**** Add deficits to produce catdef:
    sumdef=0.
    do k=1,nar0
       sumdef=sumdef+locdef(k)*px(k)
    enddo
    catdef(idep)=poros*1000.*sumdef/suma

    if (bug) write(*,*) 'funcidep: ok2'

    !c**** Compute array of root zone moisture (degree of wetness in root zone):
    do k=1,nar0              
       term1=((psis-zdep(k))/psis)              &
            **(1.-1./bee)
       if(zdep(k).le.0.) then
          wrz(k)=1.
       elseif(zdep(k)-rzdep.lt.0.) then
          term2=(-psis/zdep(k))*(bee/(bee-1.))  &
               *(term1-1.)
          frcunsat=zdep(k)/rzdep
          wrz(k)=frcunsat*term2+(1.-frcunsat)*1.
       else
          term2=((psis-zdep(k)+rzdep)           &
               /psis)**(1.-1./bee)
          wrz(k)=(-psis/rzdep)*(bee/            &
               (bee-1.))*(term1-term2)
       endif

    enddo

    if (bug) write(*,*) 'funcidep: ok3'

    !c**** compute the densities and dx
    !c**** we use a usefull property that is due to the construction of the 
    !c**** gamma distribution in tgen3.f : this distribution is continuous, 
    !c**** with decreasing values on ln(a/tanb) when n goes from 1 to nar0

    !c first we gather in the same bin all the bins with values ge 1 
    nref=1
    nind=1
    ptest(idep,1)=0.
    do k=1,nar0
       if (wrz(k) .eq. 1.) then
          nref=nref+1
          ptest(idep,1) = ptest(idep,1) + px(k)
       endif
    enddo
    if (nref .gt. 1) then
       nind=2
       valtest(idep,1)=1.
    endif
    nmax=nar0-nref+nind
    if (bug) write(*,*) 'nmax,nind,nar0,nref=',nmax,nind,nar0,nref

    !c definition of the probabilities ptest
    if (nmax .eq. 1) then     ! all the bins have values ge 1 
       dtest(idep,1) = 0.0001
       ptest(idep,1) = 1.
    else                      ! distribution in ar2/ar3
       do n=0,nmax-nind
          valtest(idep,nind+n)=wrz(nref+n)
          ptest(idep,nind+n)=px(nref+n)
       enddo

       !c we have to define dtest, the size of each bin
       if (nmax .eq. 2) then
          dtest(idep,2) = valtest(idep,1)-valtest(idep,2)
          dtest(idep,1) = dtest(idep,2)/2.
       else                   ! nmax .gt. 2
          do n=2,nmax-1                     
             dtest(idep,n)=(valtest(idep,n-1)-valtest(idep,n+1))/2.   
          enddo
          dtest(idep,1) = dtest(idep,2)/2.
          dtest(idep,nmax) = dtest(idep,nmax-1)
       endif
    endif

    if (bug) write(*,*) 'funcidep: ok4'

    !c we can now define the probability density: denstest=ptest/dtest
    !c where ptest is the probability and dtest the size of the bin
    do n=1,nmax
       if (ptest(idep,n) .eq. 0.) then
          denstest(idep,n)=0.
       else
          denstest(idep,n)=ptest(idep,n)/dtest(idep,n)
       endif
    enddo

    if (bug) write(*,*) 'funcidep: ok5'

    !c NOW we can estimate the parameters for the approximated distrib
    !c from the actual distrib

    !c 1. AR1=saturated area and AR2 and AR3 + averages of the RZ wetness 
    !c    in the different fractions

    ar1(idep)=0.
    ar2(idep)=0.
    ar3(idep)=0.
    swsrf3(idep)=0.
    swsrf2(idep)=0.
    rzeq(idep)=0.

    if(valtest(idep,1).eq.1.) ar1(idep)=dtest(idep,1)*denstest(idep,1)

    if (nmax .gt. 1) then 
       do n=nind,nmax           
          if (valtest(idep,n) .lt. wpwet) then
             ar3(idep)=ar3(idep)+denstest(idep,n)*dtest(idep,n)
             swsrf3(idep)=swsrf3(idep)+valtest(idep,n)*         &
                  denstest(idep,n)*dtest(idep,n)
          else
             ar2(idep)=ar2(idep)+denstest(idep,n)*dtest(idep,n)
             swsrf2(idep)=swsrf2(idep)+valtest(idep,n)*         &
                  denstest(idep,n)*dtest(idep,n)
          endif
       enddo
    endif

    test=ar1(idep)+ar2(idep)+ar3(idep)
    if (test .gt. 1.+1.e-5 .or. test .lt. 1.-1.e-5) then
       !         write(*,*) 'PROBLEM at depth ',zbar
       !         write(*,*) '  ar1+ar2+ar3=',test
       !         write(*,*) '  ar1=',ar1(idep),' ar2=',ar2(idep),' ar3=', &
       !             ar3(idep)
    endif

    ar1(idep)=ar1(idep)/test
    ar2(idep)=ar2(idep)/test
    ar3(idep)=ar3(idep)/test
    if (ar2(idep) .ne. 0.) swsrf2(idep)=swsrf2(idep)/ar2(idep)
    if (ar3(idep) .ne. 0.) swsrf3(idep)=swsrf3(idep)/ar3(idep)

    rzeq(idep)=ar1(idep)+ar2(idep)*swsrf2(idep)+ar3(idep)*swsrf3(idep)

    if (bug) write(*,*) 'funcidep: ok6'

    !c 2. Maximum density -> shape parameter 
    !c                    -> wmin 

    locmax=3
    shift=15
    ordref=1
    do n=1,nmax
       densaux2(n)=denstest(idep,n)
    enddo

    if (nmax .ge. shift*2) then

       !c we start with sliding mean to facilitate the search for the maximum

       ord=MIN(ordref,nmax/shift)

       call smtot(densaux2,nmax,ord,densaux) 
       !	 print *,nmax,ord,shift,densaux(shift-14),shift-14,size(densaux)
       do n=nmax,shift,-1
          if (densaux(n) .gt. densaux(n-1) .and.     &
               densaux(n) .gt. densaux(n-2) .and.     &
               densaux(n) .gt. densaux(n-3) .and.     &
               densaux(n) .gt. densaux(n-4) .and.     &
               densaux(n) .gt. densaux(n-5) .and.     &
               densaux(n) .gt. densaux(n-6) .and.     &
               densaux(n) .gt. densaux(n-7) .and.     &
               densaux(n) .gt. densaux(n-8) .and.     &
               densaux(n) .gt. densaux(n-9) .and.     &
               densaux(n) .gt. densaux(n-10) .and.    &
               densaux(n) .gt. densaux(n-11) .and.    &
               densaux(n) .gt. densaux(n-12) .and.    &
               densaux(n) .gt. densaux(n-13) .and.    &
               densaux(n) .gt. densaux(n-14))then ! .and.    &
             !                densaux(n) .gt. densaux(n-15)) then
             locmax=n
             goto 30
          endif
       enddo

    else

       aux10=-9999.
       indimax10=3
       do n=1,nmax
          if (densaux2(n) .gt. aux10) then
             aux10=densaux2(n)
             indimax10=n
          endif
       enddo
       locmax=MAX(3,indimax10)
       ! add protection here in case nmax <3 . why 3 ?
       if (locmax > nmax) locmax = nmax
    endif      ! if (nmax .ge. shift+1)
30  densmax=denstest(idep,locmax)
    aa(idep)=exp(1.)*densmax

    if (bug) write(*,*) 'funcidep: ok7'

    !c WMIN=lowest value where the density is strictly gt densmax/100.

    indmin=1
    indmin0=0
    do n=1,nmax
       if (denstest(idep,n) .gt. 0.) indmin0=n
       if (denstest(idep,n) .gt. densmax/100. .and.    &
            valtest(idep,n) .lt. valtest(idep,locmax)) indmin=n
    enddo
    if (indmin .eq.0) indmin=indmin0

    if (indmin .le. 2) then
       wmin(idep) = 0.99999
    else
       x1=valtest(idep,indmin)
       wmin(idep)=x1
    endif

    if (bug) write(*,*) 'funcidep: ok8; first wmin=',wmin(idep)

    !c for negative or low coeskew the previous wmin doesn't give good results...
    !c wmin is higher !!!

    if (coeskew .lt. 1. ) then

       if (locmax .gt. 3 .and. indmin .ge. locmax+4) then
          n2=MAX(locmax+1,(indmin-locmax)/2+locmax)
          x2=valtest(idep,n2)
          y2=denstest(idep,n2)
          n1=locmax
          x1=valtest(idep,n1)
          y1=denstest(idep,n1)
          wa=(y2-y1)/(x2-x1)
          wb=y1-wa*x1
          wmin(idep)=AMAX1(wmin(idep),-wb/wa)
       endif

       !c wmin is even higher in some cases !!!
       if (coeskew .lt. 0.2 ) wmin(idep)=wmin(idep)+0.01

    endif

    if (bug) write(*,*) 'funcidep: ok9; 2nd wmin=',wmin(idep)

    if (valtest(idep,locmax) .le. wmin(idep)) then ! doesn't make sense
       wmin(idep)=valtest(idep,locmax)-dx
    endif
    aabis(idep)=1./(valtest(idep,locmax)-wmin(idep)+1.e-20)

    if (bug) write(*,*) 'funcidep: ok10'

  END SUBROUTINE FUNCIDEP
  
  ! ********************************************************************
  
  SUBROUTINE FUNCZBAR(                                     &   
       NAR0,ZBAR,                            &
       BEE,PSIS,POROS,COND,RZDEP,WPWET,      &
       VALX,PX,COESKEW,TIMEAN,SUMA,          &
       CATDEF,WMIN)                          

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c                                                                      c
    !c This program returns the eight parameters for the areal fractioning  c
    !c                                                                      c
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    INTEGER , intent (in) :: NAR0
    integer nref,nind,nmax,indmin,locmax,shift,ord,locmin,ordref
    integer indimax10,indmin0
    REAL, intent (in) :: BEE, PSIS, POROS, COND, RZDEP, WPWET, COESKEW
    REAL, intent (inout) ::  VALX(NAR), PX(NAR),TIMEAN,SUMA,zbar
    real, intent (inout) :: catdef,wmin

    REAL  dx,dz,sumdef
    real term1,term2
    real zdep(nar),locdef(nar),wrz(nar),frcunsat
    real valtest(nar),ptest(nar),denstest(nar),dtest(nar)
    real x1,x2,y1,y2,wa,wb
    integer n1,n2,k,n
    real densaux(nar),densaux2(nar),densmax,aux10

    !c-------------------------------------------------------------------------
    !c integral(f(x)dx)=1. for a pdf
    !c here px=f(x)dx
    dx=valx(1)-valx(2)

    !c**** Compute array of water table depths:
    do k=1,nar0
       term1=(1/gnu)*(valx(k)-timean)
       zdep(k)=AMAX1(0.,zbar-term1)
    enddo

    !c variable change must be reflected in dx
    dz=dx/gnu

    !c**** Compute array of moisture deficits:
    do k=1,nar0
       term1=(psis-zdep(k))/psis
       term1=term1**(1.-1./bee)
       term2=-psis*(bee/(bee-1.))*(term1-1.)
       locdef(k)=zdep(k)-term2
    enddo

    !c**** Add deficits to produce catdef:
    sumdef=0.
    do k=1,nar0
       sumdef=sumdef+locdef(k)*px(k)
    enddo
    catdef=poros*1000.*sumdef/suma

    !c**** Compute array of root zone moisture (degree of wetness in root zone):
    do k=1,nar0              
       term1=((psis-zdep(k))/psis)  &
            **(1.-1./bee)
       if(zdep(k).le.0.) then
          wrz(k)=1.
       elseif(zdep(k)-rzdep.lt.0.) then
          term2=(-psis/zdep(k))*(bee/(bee-1.))   &
               *(term1-1.)
          frcunsat=zdep(k)/rzdep
          wrz(k)=frcunsat*term2+(1.-frcunsat)*1.
       else
          term2=((psis-zdep(k)+rzdep)     &
               /psis)**(1.-1./bee)
          wrz(k)=(-psis/rzdep)*(bee/      &
               (bee-1.))*(term1-term2)
       endif
    enddo

    !c**** compute the densities and dx
    !c**** we use a usefull property that is due to the construction of the 
    !c**** gamma distribution in tgen3.f : this distribution is continuous, 
    !c**** with decreasing values on ln(a/tanb) when n goes from 1 to nar0
    !c first we gather in the same bin all the bins with values ge 1 
    nref=1
    nind=1
    ptest(1)=0.
    do k=1,nar0
       if (wrz(k) .eq. 1.) then
          nref=nref+1
          ptest(1) = ptest(1) + px(k)
       endif
    enddo
    if (nref .gt. 1) then
       nind=2
       valtest(1)=1.
    endif
    nmax=nar0-nref+nind

    !c definition of the probabilities ptest
    if (nmax .eq. 1) then     ! all the bins have values ge 1 
       dtest(1) = 0.0001
       ptest(1) = 1.
    else                      ! distribution in ar2/ar3
       do n=0,nmax-nind
          valtest(nind+n)=wrz(nref+n)
          ptest(nind+n)=px(nref+n)
       enddo

       !c we have to define dtest, the size of each bin
       if (nmax .eq. 2) then
          dtest(2) = valtest(1)-valtest(2)
          dtest(1) = dtest(2)/2.
       else                   ! nmax .gt. 2
          do n=2,nmax-1
             dtest(n)=(valtest(n-1)-valtest(n+1))/2.            
          enddo
          dtest(1) = dtest(2)/2.
          dtest(nmax) = dtest(nmax-1)
       endif
    endif

    !c we can now define the probability density: denstest=ptest/dtest
    !c where ptest is the probability and dtest the size of the bin
    do n=1,nmax
       if (ptest(n) .eq. 0.) then
          denstest(n)=0.
       else
          denstest(n)=ptest(n)/dtest(n)
       endif
    enddo

    !c NOW we can estimate the parameters for the approximated distrib
    !c from the actual distrib

    !c 2. Maximum density -> shape parameter 
    !c                    -> wmin 

    locmax=3
    shift=15
    ordref=1
    do n=1,nmax
       densaux2(n)=denstest(n)
    enddo

    if (nmax .ge. shift*2) then

       !c we start with sliding mean to facilitate the search for the maximum

       ord=MIN(ordref,nmax/shift)
       call smtot(densaux2,nmax,ord,densaux)

       do n=nmax,shift,-1
          if (densaux(n) .gt. densaux(n-1) .and.         &
               densaux(n) .gt. densaux(n-2) .and.         &
               densaux(n) .gt. densaux(n-3) .and.         &
               densaux(n) .gt. densaux(n-4) .and.         &
               densaux(n) .gt. densaux(n-5) .and.         &
               densaux(n) .gt. densaux(n-6) .and.         &
               densaux(n) .gt. densaux(n-7) .and.         &
               densaux(n) .gt. densaux(n-8) .and.         &
               densaux(n) .gt. densaux(n-9) .and.         &
               densaux(n) .gt. densaux(n-10) .and.        &
               densaux(n) .gt. densaux(n-11) .and.        &
               densaux(n) .gt. densaux(n-12) .and.        &
               densaux(n) .gt. densaux(n-13) .and.        &
               densaux(n) .gt. densaux(n-14)) then ! .and.        &
             !densaux(n) .gt. densaux(n-15)) then
             locmax=n
             goto 30
          endif
       enddo

    else

       aux10=-9999.
       indimax10=3
       do n=1,nmax
          if (densaux2(n) .gt. aux10) then
             aux10=densaux2(n)
             indimax10=n
          endif
       enddo
       locmax=MAX(3,indimax10)
       ! in case nmax < 3. why hard coded 3?
       if(locmax > nmax) locmax = nmax
    endif      ! if (nmax .ge. shift+1) 

30  densmax=denstest(locmax)

    !c WMIN=lowest value where the density is strictly gt densmax/100.

    indmin=1
    indmin0=0
    do n=1,nmax
       if (denstest(n) .gt. 0.) indmin0=n
       if (denstest(n) .gt. densmax/100. .and.         &
            valtest(n) .lt. valtest(locmax)) indmin=n
    enddo
    if (indmin .eq. 0) indmin=indmin0

    if (indmin .le. 2) then
       wmin = 0.99999
    else
       x1=valtest(indmin)
       wmin=x1
    endif

    !c for negative or low coeskew the previous wmin doesn't give good results...
    !c wmin is higher !!!

    if (coeskew .lt. 1. ) then

       if (locmax .gt. 3 .and. indmin .ge. locmax+4) then

          n2=MAX(locmax+1,(indmin-locmax)/2+locmax)
          x2=valtest(n2)
          y2=denstest(n2)
          n1=locmax
          x1=valtest(n1)
          y1=denstest(n1)
          wa=(y2-y1)/(x2-x1)
          wb=y1-wa*x1
          wmin=AMAX1(wmin,-wb/wa)
       endif

       !c wmin is even higher in some cases !!!
       if (coeskew .lt. 0.2 ) wmin=wmin+0.01

    endif

  END SUBROUTINE FUNCZBAR

  ! ******************************************************************
  
  SUBROUTINE RMSE(XX,YY,LEN,ERROR)

    !c---------------------------------------------------------------------------
    !c Computes the root-mean square error ERROR between two one-dimensional
    !c random variables XX and YY of same length LEN
    !c---------------------------------------------------------------------------

    INTEGER, intent (in) :: LEN
    REAL, intent (in) ::  XX(LEN),YY(LEN)
    REAL, intent (out) :: ERROR
    INTEGER :: I

    !c---------------------------------------------------------------------------     
    error=0.
    do i=1,len
       if(abs(xx(i)-yy(i)) >=1.e-10) then
          error=error+(xx(i)-yy(i))*(xx(i)-yy(i))
       endif
    enddo
    error=SQRT(error/float(len))

  END SUBROUTINE RMSE

  ! ******************************************************************
  
  SUBROUTINE SMTOT(XX,LEN,ORD,YY)
    
    !c---------------------------------------------------------------------------
    !c Runs a sliding average of order ORD through the one-dimensional array XX
    !c of length LEN and returns the smoothed YY
    !!c---------------------------------------------------------------------------

    INTEGER, intent(in)  :: LEN

    INTEGER :: ORD,WIDTH,i,ini,n,fin  ! replaced var name "end" w/ "fin" to fix auto-indent, reichle, 24 Dec 2024

    REAL,    intent(in)  :: XX(NAR)
    REAL,    intent(out) :: YY(NAR)
    
    !c---------------------------------------------------------------------------     
    do i=1,nar
       yy(i)=0.
    enddo
    
    width=ord*2+1
    if (width .gt. len/2) then
       write(*,*) 'the order for the sliding average is too large !!!'
       write(*,*) 'regard with the length of the array to be smoothed'
       stop
    endif

    do i=1,len
       ini=MAX(1,i-ord)
       fin=MIN(len,i+ord)
       yy(i)=0.
       do n=ini,fin
          yy(i)=yy(i)+xx(n)
       enddo
       yy(i)=yy(i)/(fin-ini+1)
    enddo
    
  END SUBROUTINE SMTOT

  ! -----------------------------------------------------------------------------------

  subroutine RegridRaster(Rin,Rout)

    ! primitive regridding of integer values from 2-dim array Rin to 2-dim array Rout
    !
    ! If Rout is higher-resolution than Rin, result should be fine:
    !   An Rout grid cell is assigned the value of the Rin grid cell that 
    !   contains the center of the Rout grid cell (oversampling). 
    ! If Rin is higher-resolution than Rout, result is questionable:
    !   An Rout grid cell is assigned the value of the Rin grid cell that is 
    !    near the *corner* of the Rout grid cell. See notes below.

    integer, intent(IN)  :: Rin( :,:)
    integer, intent(OUT) :: Rout(:,:)

    REAL(REAL64) :: xx, yy
    integer      :: i, j, ii, jj
    integer      :: Nx_in, Ny_in, Nx_out, Ny_out

    Nx_in  = size(Rin ,1)
    Ny_in  = size(Rin ,2) 

    Nx_out = size(Rout,1)
    Ny_out = size(Rout,2) 

    !if ( (Nx_in==Nx_out) .and. (Ny_in==Ny_out) ) then
    if (.false.) then    

       ! avoid loop through output grid cells

       Rout = Rin       ! [??] MAY NOT BE 0-DIFF B/C OF MIXED-MODE ARITHMETIC IN LOOP!?!?!?

    else

       ! NOTE: float() yields real*4 but xx was declared real*8

       xx = Nx_in/float(Nx_out)      ! WARNING: mixed mode arithmentic!!! 
       yy = Ny_in/float(Ny_out)      ! WARNING: mixed mode arithmentic!!! 

       do j=1,Ny_out

          ! NOTE: When Rin is finer resolution than Rout, the below use of  
          !          ii = (i-1)*xx + 1                          (1a)
          !          jj = (j-1)*yy + 1                          (1b)
          !       implies that Rout(i,j) is assigned the Rin(ii,jj) value near a corner of 
          !       the (ii,jj) output grid cell, which effectively results in a shift of the 
          !       data by 1/2 of the width of the output grid cell.  This shift could
          !       presumably minimized by using 
          !          ii = NINT( (i-1)*xx + xx/2 )               (2a)
          !          jj = NINT( (j-1)*yy + yy/2 )               (2b)
          !
          !       HOWEVER, equations (2a) and (2b) are preferable when Rout is finer resolution
          !       than Rin, in which case Rout should just be oversampling of Rin.

          jj = (j-1)*yy + 1          ! WARNING: mixed mode arithmetic!!!  Note implied "floor()" operator.
          do i=1,Nx_out
             ii = (i-1)*xx + 1       ! WARNING: mixed mode arithmetic!!!  Note implied "floor()" operator.
             Rout(i,j) = Rin(ii,jj)  
          end do
       end do

    end if

  end subroutine RegridRaster

  ! -----------------------------------------------------------------------------------
  
  subroutine RegridRaster1(Rin,Rout)

    ! same as RegridRaster() but for gridded integer*1 values

    integer*1, intent(IN)  :: Rin( :,:)
    integer*1, intent(OUT) :: Rout(:,:)

    REAL(REAL64) :: xx, yy
    integer      :: i, j, ii, jj
    integer      :: Nx_in, Ny_in, Nx_out, Ny_out

    Nx_in  = size(Rin ,1)
    Ny_in  = size(Rin ,2) 

    Nx_out = size(Rout,1)
    Ny_out = size(Rout,2) 

    !if ( (Nx_in==Nx_out) .and. (Ny_in==Ny_out) ) then
    if (.false.) then    

       Rout = Rin

    else

       xx = Nx_in/float(Nx_out)
       yy = Ny_in/float(Ny_out)

       do j=1,Ny_out
          jj = (j-1)*yy + 1
          do i=1,Nx_out
             ii = (i-1)*xx + 1
             Rout(i,j) = Rin(ii,jj)
          end do
       end do

    end if

  end subroutine RegridRaster1

  ! -----------------------------------------------------------------------------------

  subroutine RegridRaster2(Rin,Rout)

    ! same as RegridRaster() but for gridded integer*2 values

    integer(kind=2), intent(IN)  :: Rin( :,:)
    integer(kind=2), intent(OUT) :: Rout(:,:)

    REAL(REAL64) :: xx, yy
    integer      :: i, j, ii, jj
    integer      :: Nx_in, Ny_in, Nx_out, Ny_out

    Nx_in  = size(Rin ,1)
    Ny_in  = size(Rin ,2) 

    Nx_out = size(Rout,1)
    Ny_out = size(Rout,2) 

    !if ( (Nx_in==Nx_out) .and. (Ny_in==Ny_out) ) then
    if (.false.) then    

       Rout = Rin

    else

       xx = Nx_in/float(Nx_out)
       yy = Ny_in/float(Ny_out)

       do j=1,Ny_out
          jj = (j-1)*yy + 1
          do i=1,Nx_out
             ii = (i-1)*xx + 1
             Rout(i,j) = Rin(ii,jj)
          end do
       end do

    end if

  end subroutine RegridRaster2

  ! -----------------------------------------------------------------------------------

  subroutine RegridRasterReal(Rin,Rout)

    ! same as RegridRaster() but for gridded real values

    real, intent(IN)  :: Rin( :,:)
    real, intent(OUT) :: Rout(:,:)

    REAL(REAL64) :: xx, yy
    integer      :: i, j, ii, jj
    integer      :: Nx_in, Ny_in, Nx_out, Ny_out

    Nx_in  = size(Rin ,1)
    Ny_in  = size(Rin ,2) 

    Nx_out = size(Rout,1)
    Ny_out = size(Rout,2) 

    !if ( (Nx_in==Nx_out) .and. (Ny_in==Ny_out) ) then
    if (.false.) then    

       Rout = Rin

    else

       xx = Nx_in/float(Nx_out)
       yy = Ny_in/float(Ny_out)

       do j=1,Ny_out
          jj = (j-1)*yy + 1
          do i=1,Nx_out
             ii = (i-1)*xx + 1
             Rout(i,j) = Rin(ii,jj)
          end do
       end do

    end if

  end subroutine RegridRasterReal

  !---------------------------------------------------------------------
  
  SUBROUTINE svbksb(u,w,v,m,n,b,x) 

    INTEGER m,mp,n,np,NMAX 
    REAL*8 b(m),u(m,n),v(n,n),w(n),x(n) 
    PARAMETER (NMAX=500)  !Maximum anticipated value of n
    !------------------------------------------------------------------------------------------- 
    ! Solves A "A^" . X = B for a vector X, where A is specified by the arrays u, w, v as returned by 
    ! svdcmp. m and n are the dimensions of a, and will be equal for square matrices. b(1:m) is 
    ! the input right-hand side. x(1:n) is the output solution vector. No input quantities are 
    ! destroyed, so the routine may be called sequentially with different b's. 
    !-------------------------------------------------------------------------------------------

    INTEGER i,j,jj 
    REAL*8 s,tmp(NMAX) 
    do j=1,n !Calculate UTB. 
       s=0. 
       if(w(j).ne.0.)then !Nonzero result only if wj is nonzero. 
          do i=1,m 
             s=s+u(i,j)*b(i) 
          end do
          s=s/(w(j) + 1.d-20) !This is the divide by wj . 
       endif
       tmp(j)=s 
    end do
    do j=1,n !Matrix multiply by V to get answer. 
       s=0. 
       do jj=1,n 
          s=s+v(j,jj)*tmp(jj) 
       end do
       x(j)=s 
    end do
    return 
  END SUBROUTINE svbksb

  !---------------------------------------------------------------------

  SUBROUTINE svdcmp(a,m,n,w,v) 

    INTEGER m,n,NMAX 
    REAL*8, intent (inout)  :: a(m,n)
    REAL*8, intent (out) :: v(n,n),w(n) 
    PARAMETER (NMAX=500)  !Maximum anticipated value of n. 
    !-------------------------------------------------------------------------------------- 
    ! Given a matrix A(1:m,1:n), this routine computes its singular value decomposition, 
    ! A = U . W . Vt. The matrix U replaces A on output. The diagonal matrix of singular 
    ! values W is output as a vector W(1:n). The matrix V (not the transpose Vt) is output 
    ! as V(1:n,1:n). 
    !--------------------------------------------------------------------------------------

    INTEGER i,its,j,jj,k,l,nm 
    REAL*8 anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX) 
    real*8, parameter :: EPS=epsilon(1.0d0)
    g=0.d0  !Householder reduction to bidiagonal form. 
    scale=0.d0 
    anorm=0.d0 
    c =0.d0 
    f =0.d0 
    g =0.d0 
    h =0.d0 
    s =0.d0 
    x =0.d0 
    y =0.d0 
    z =0.d0 
    rv1=0.d0 
    w = 0.d0 
    v = 0.d0 
    do i=1,n 
       l=i+1 
       rv1(i)=scale*g 
       g=0.d0 
       s=0.d0 
       scale=0.d0 
       if(i.le.m)then 
          do k=i,m 
             scale=scale+abs(a(k,i)) 
          end do
          if(scale.ne.0.d0)then 
             do k=i,m 
                a(k,i)=a(k,i)/scale 
                s=s+a(k,i)*a(k,i) 
             end do
             f=a(i,i) 
             g=-dsign(dsqrt(s),f) 
             h=f*g-s 
             a(i,i)=f-g 
             do j=l,n 
                s=0.d0 
                do k=i,m 
                   s=s+a(k,i)*a(k,j) 
                end do
                f=s/h 
                do k=i,m 
                   a(k,j)=a(k,j)+f*a(k,i) 
                end do
             end do
             do k=i,m 
                a(k,i)=scale*a(k,i) 
             end do
          endif
       endif
       w(i)=scale *g 
       g=0.d0 
       s=0.d0 
       scale=0.d0 
       if((i.le.m).and.(i.ne.n))then 
          do k=l,n 
             scale=scale+abs(a(i,k)) 
          end do
          if(scale.ne.0.d0)then 
             do k=l,n 
                a(i,k)=a(i,k)/scale 
                s=s+a(i,k)*a(i,k) 
             end do
             f=a(i,l) 
             g=-sign(sqrt(s),f) 
             h=f*g-s 
             a(i,l)=f-g 
             do k=l,n 
                rv1(k)=a(i,k)/h 
             end do
             do j=l,m 
                s=0.d0 
                do k=l,n 
                   s=s+a(j,k)*a(i,k) 
                end do
                do k=l,n 
                   a(j,k)=a(j,k)+s*rv1(k) 
                end do
             end do
             do k=l,n 
                a(i,k)=scale*a(i,k) 
             end do
          endif
       endif
       anorm=max(anorm,(abs(w(i))+abs(rv1(i)))) 
    end do !do i=1,n

    do i=n,1,-1 !Accumulation of right-hand transformations. 
       if(i.lt.n)then 
          if(g.ne.0.d0)then 
             do j=l,n       !Double division to avoid possible underflow. 
                v(j,i)=(a(i,j)/a(i,l))/g 
             end do
             do j=l,n 
                s=0.d0 
                do k=l,n 
                   s=s+a(i,k)*v(k,j) 
                end do
                do k=l,n 
                   v(k,j)=v(k,j)+s*v(k,i) 
                end do
             end do
          endif
          do j=l,n 
             v(i,j)=0.d0 
             v(j,i)=0.d0 
          end do
       endif
       v(i,i)=1.d0
       g=rv1(i) 
       l=i 
    end do

    do i=min(m,n),1,-1 !Accumulation of left-hand transformations. 
       l=i+1 
       g=w(i) 
       do j=l,n 
          a(i,j)=0.d0 
       end do
       if(g.ne.0.d0)then 
          g=1.d0/g 
          do j=l,n 
             s=0.d0 
             do k=l,m 
                s=s+a(k,i)*a(k,j) 
             end do
             f=(s/a(i,i))*g 
             do k=i,m 
                a(k,j)=a(k,j)+f*a(k,i) 
             end do
          end do
          do j=i,m 
             a(j,i)=a(j,i)*g 
          end do
       else
          do j= i,m 
             a(j,i)=0.d0 
          end do
       endif
       a(i,i)=a(i,i)+1.d0 
    end do

    do k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over 
       !singular values, and over allowed iterations. 
       do its=1,30 
          do l=k,1,-1 !Test for splitting. 
             nm=l-1 !Note that rv1(1) is always zero.
             if( abs(rv1(l)) <= EPS*anorm ) goto 2 
             if( abs(w(nm) ) <= EPS*anorm ) goto 1  
          end do
1         c=0.d0 !Cancellation of rv1(l), if l > 1. 
          s=1.d0 
          do i=l,k 
             f=s*rv1(i) 
             rv1(i)=c*rv1(i) 
             if( abs(f) <= EPS*anorm ) goto 2 
             g=w(i) 
             h=pythag(f,g) 
             w(i)=h 
             h=1.d0/h 
             c= (g*h) 
             s=-(f*h) 
             do j=1,m 
                y=a(j,nm) 
                z=a(j,i) 
                a(j,nm)=(y*c)+(z*s) 
                a(j,i)=-(y*s)+(z*c) 
             end do
          end do
2         z=w(k) 
          if(l.eq.k)then   !Convergence. 
             if(z.lt.0.d0)then !Singular value is made nonnegative. 
                w(k)=-z 
                do j=1,n 
                   v(j,k)=-v(j,k) 
                end do
             endif
             goto 3 
          endif
          if(its.eq.30) print *, 'no convergence in svdcmp' 
          !             if(its.ge.4)  print *, 'its = ',its
          x=w(l) !Shift from bottom 2-by-2 minor. 
          nm=k-1 
          y=w(nm) 
          g=rv1(nm) 
          h=rv1(k) 
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y) 
          g=pythag(f,1.d0) 
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x 
          c=1.d0 !Next QR transformation: 
          s=1.d0 
          do j=l,nm 
             i=j+1 
             g=rv1(i) 
             y=w(i) 
             h=s*g 
             g=c*g 
             z=pythag(f,h) 
             rv1(j)=z 
             c=f/z 
             s=h/z 
             f= (x*c)+(g*s) 
             g=-(x*s)+(g*c) 
             h=y*s 
             y=y*c 
             do jj=1,n 
                x=v(jj,j) 
                z=v(jj,i) 
                v(jj,j)= (x*c)+(z*s) 
                v(jj,i)=-(x*s)+(z*c) 
             end do
             z=pythag(f,h) 
             w(j)=z !Rotation can be arbitrary if z = 0. 
             if(z.ne.0.d0)then 
                z=1.d0/z 
                c=f*z 
                s=h*z 
             endif
             f= (c*g)+(s*y) 
             x=-(s*g)+(c*y) 
             do jj=1,m 
                y=a(jj,j) 
                z=a(jj,i) 
                a(jj,j)= (y*c)+(z*s) 
                a(jj,i)=-(y*s)+(z*c) 
             end do
          end do !j=l;nm 
          rv1(l)=0.d0 
          rv1(k)=f 
          w(k)=x 
       end do !its=1,30
3      continue 
    end do !k=n,1,-1 
    return 
  END SUBROUTINE svdcmp
  !
  ! ________________________________________________________________________________
  !     
  REAL*8 FUNCTION pythag(a,b) 
    REAL*8 a,b 
    !Computes sqrt(a**2 + b**2) without destructive underflow or overflow.
    REAL*8 absa,absb 
    absa=abs(a) 
    absb=abs(b) 
    if(absa.gt.absb)then 
       pythag=absa*sqrt(1.+(absb/absa)**2) 
    else 
       if(absb.eq.0.)then 
          pythag=0. 
       else
          pythag=absb*sqrt(1.+(absa/absb)**2) 
       endif
    endif
    return 
  END FUNCTION pythag
  !
  ! ________________________________________________________________________________
  !     

  SUBROUTINE savgol(c,np,nl,nr,ld,m)

    INTEGER ld,m,nl,np,nr,MMAX 
    real c(np) 
    PARAMETER (MMAX=6)
    !-------------------------------------------------------------------------------------------- 
    !USES lubksb,ludcmp given below. 
    !Returns in c(1:np), in wrap-around order (see reference) consistent with the argument respns 
    !in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward 
    !(past) data points used, while nr is the number of rightward (future) data points, making 
    !the total number of data points used nl +nr+1. ld is the order of the derivative desired 
    !(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also 
    !equal to the highest conserved moment; usual values are m = 2 or m = 4. 
    !--------------------------------------------------------------------------------------------
    INTEGER d,icode,imj,ipj,j,k,kk,mm,indx(MMAX+1) 
    real fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
    if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX  & 
         .or.nl+nr.lt.m) pause ' Bad args in savgol.' 
    do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit. 
       sum=0. 
       if(ipj.eq.0) sum=1. 
       do k=1,nr 
          sum=sum+dfloat(k)**ipj 
       end do
       do k=1,nl 
          sum=sum+dfloat(-k)**ipj 
       end do
       mm=min(ipj,2*m-ipj) 
       do imj=-mm,mm,2 
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum 
       end do
    end do

    call ludcmp(a,m+1,MMAX+1,indx,d,icode)    !Solve them: LU decomposition. 

    do j=1,m+1 
       b(j)=0. 
    end do
    b(ld+1)=1.      !Right-hand side vector is unit vector, depending on which derivative we want. 

    call lubksb(a,m+1,MMAX+1,indx,b)   !Backsubstitute, giving one row of the inverse matrix. 

    do kk=1,np                         !Zero the output array (it may be bigger than the number 
       c(kk)=0.                         !of coefficients).  
    end do
    do k=-nl,nr                        !Each Savitzky-Golay coefficient is the dot product 
       sum=b(1)                         !of powers of an integer with the inverse matrix row. 
       fac=1. 
       do mm=1,m 
          fac=fac*k 
          sum=sum+b(mm+1)*fac 
       end do
       kk=mod(np-k,np)+1                !Store in wrap-around order. 
       c(kk)=sum 
    end do
    return 
  END SUBROUTINE savgol

  !***************************************************************
  !* Given an N x N matrix A, this routine replaces it by the LU *
  !* decomposition of a rowwise permutation of itself. A and N   *
  !* are input. INDX is an output vector which records the row   *
  !* permutation effected by the partial pivoting; D is output   *
  !* as -1 or 1, depending on whether the number of row inter-   *
  !* changes was even or odd, respectively. This routine is used *
  !* in combination with LUBKSB to solve linear equations or to  *
  !* invert a matrix. Return code is 1, if matrix is singular.   *
  !***************************************************************
  Subroutine LUDCMP(A,N,NP,INDX,D,CODE)
    INTEGER, PARAMETER :: NMAX=100
    REAL, PARAMETER :: TINY=1E-12
    real  AMAX,DUM, SUM, A(NP,NP),VV(NMAX)
    INTEGER CODE, D, INDX(N),NP,N,I,J,K,IMAX

    D=1; CODE=0

    DO I=1,N
       AMAX=0.
       DO J=1,N
          IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
          CODE = 1
          RETURN
       END IF
       VV(I) = 1. / AMAX
    END DO ! i loop

    DO J=1,N
       DO I=1,J-1
          SUM = A(I,J)
          DO K=1,I-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
       END DO ! i loop
       AMAX = 0.
       DO I=J,N
          SUM = A(I,J)
          DO K=1,J-1
             SUM = SUM - A(I,K)*A(K,J) 
          END DO ! k loop
          A(I,J) = SUM
          DUM = VV(I)*ABS(SUM)
          IF(DUM.GE.AMAX) THEN
             IMAX = I
             AMAX = DUM
          END IF
       END DO ! i loop  

       IF(J.NE.IMAX) THEN
          DO K=1,N
             DUM = A(IMAX,K)
             A(IMAX,K) = A(J,K)
             A(J,K) = DUM
          END DO ! k loop
          D = -D
          VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(ABS(A(J,J)) < TINY) A(J,J) = TINY

       IF(J.NE.N) THEN
          DUM = 1. / A(J,J)
          DO I=J+1,N
             A(I,J) = A(I,J)*DUM
          END DO ! i loop
       END IF
    END DO ! j loop

    RETURN
  END Subroutine LUDCMP


  !******************************************************************
  !* Solves the set of N linear equations A . X = B.  Here A is     *
  !* input, not as the matrix A but rather as its LU decomposition, *
  !* determined by the routine LUDCMP. INDX is input as the permuta-*
  !* tion vector returned by LUDCMP. B is input as the right-hand   *
  !* side vector B, and returns with the solution vector X. A, N and*
  !* INDX are not modified by this routine and can be used for suc- *
  !* cessive calls with different right-hand sides. This routine is *
  !* also efficient for plain matrix inversion.                     *
  !******************************************************************
  Subroutine LUBKSB(A,N,NP,INDX,B)
    INTEGER :: II,I,J,LL,N,NP
    real  SUM, A(NP,NP),B(N)
    INTEGER INDX(N)

    II = 0

    DO I=1,N
       LL = INDX(I)
       SUM = B(LL)
       B(LL) = B(I)
       IF(II.NE.0) THEN
          DO J=II,I-1
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       ELSE IF(SUM.NE.0.) THEN
          II = I
       END IF
       B(I) = SUM
    END DO ! i loop

    DO I=N,1,-1
       SUM = B(I)
       IF(I < N) THEN
          DO J=I+1,N
             SUM = SUM - A(I,J)*B(J)
          END DO ! j loop
       END IF
       B(I) = SUM / A(I,I)
    END DO ! i loop

    RETURN
  END Subroutine LUBKSB
  
  !
  ! ====================================================================
  !
  
  INTEGER FUNCTION center_pix (x,y,x0,y0,z0,ext_point)

    real, dimension (:), intent(in   ) :: x,y

    real,                intent(inout) :: x0,y0,z0
    logical,             intent(in   ) :: ext_point

    ! ------------------------------------------------------

    real, allocatable, dimension (:,:) :: length_m
    real, allocatable, dimension (:)   :: length
    
    integer :: i,j,npix,ii

    real :: zi, zj

    npix = size (x)
    allocate (length_m (1:npix,1:npix))
    allocate (length   (1:npix))
    length_m =0.
    length   =0.

    do i = 1,npix
       zi = 100. - x(i) - y(i)
       if (.not. ext_point) then
          x0 = x(i)
          y0 = y(i)
          z0 = zi
       endif

       do j = i,npix
          zj = 100. - x(j) - y(j)
          !      length_m (i,j) = abs (x(j) - x0) + &
          !            abs (y(j) - y0) +  abs (zj - z0)
          !
          length_m (i,j) = ((x(j) - x0)*(x(j) - x0) &
               +  (y(j) - y0)*(y(j) - y0) &
               +  (zj - z0)*(zj - z0))**0.5
          length_m (j,i) = length_m (i,j)
       end do
       length (i) = sum(length_m (i,:))
    end do

    center_pix = minloc(length,dim=1)

  END FUNCTION center_pix

  !
  !----------------------------------------------------------
  !

  INTEGER FUNCTION soil_class (min_perc)

    ! Function returns a unique soil class [1-100], 

    type(mineral_perc), intent (in)  :: min_perc

    ! ------------------------------------------------
    
    integer :: clay_row, sand_row, silt_row

    clay_row = ceiling((100.- min_perc%clay_perc)/10.)
    if(clay_row == 0 ) clay_row = 1
    if(clay_row == 11) clay_row = 10

    sand_row = ceiling((min_perc%sand_perc)/10.)
    if(sand_row == 0 ) sand_row = 1
    if(sand_row == 11) sand_row = 10

    silt_row = ceiling((min_perc%silt_perc)/10.)
    if(silt_row == 0 ) silt_row = 1
    if(silt_row == 11) silt_row = 10

    if(clay_row == 1) soil_class=1

    if(clay_row > 1) soil_class=   &
         (clay_row - 1)*(clay_row - 1) + (clay_row - sand_row) + silt_row

  end FUNCTION soil_class

  ! -----------------------------------------------------------------------------------

  SUBROUTINE REFORMAT_VEGFILES

    character*512 :: tmp_string
    integer :: n_tiles
    real, dimension (:), allocatable :: var_array
    character*512 :: header
    real :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
    integer :: month

    tmp_string = 'mkdir -p '//'clsm/g5fmt'
    call execute_command_line(tmp_string)
    tmp_string = '/bin/mv '//'clsm/lai.dat ' //'clsm/g5fmt/.'
    call execute_command_line(tmp_string) 
    tmp_string = '/bin/mv '//'clsm/green.dat ' //'clsm/g5fmt/.'
    call execute_command_line(tmp_string) 

    open (10,file='clsm/g5fmt/lai.dat'  , form = 'unformatted',   &
         convert='little_endian',status='old',action='read' )
    open (11,file='clsm/g5fmt/green.dat', form = 'unformatted',   &
         convert='little_endian',status='old',action='read' )

    open (20,file='clsm/lai.dat', form = 'unformatted',   &
         convert='big_endian',status='unknown',action='write' )
    open (21,file='clsm/green.dat', form = 'unformatted', &
         convert='big_endian',status='unknown',action='write' )

    open (30,file='clsm/catchment.def', form = 'formatted',status='old',action='read' )
    read (30,*) n_tiles
    close(30,status='keep')

    allocate (var_array (1:n_tiles))

    read (10) header
    read (10) var_array
    read (11) header
    read (11) var_array

    do month  =1,12

       read (10) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
       read (10) var_array
       print '(12f3.0,f4.00,f2.0,a6,2f6.2)',a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,'LAI ',minval(var_array),maxval(var_array)
       write (20)var_array(:) 
       read (11) a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
       read (11) var_array
       print '(12f3.0,f4.00,f2.0,a6,2f6.2)',a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,'GREEN ',minval(var_array),maxval(var_array) 
       write (21)var_array(:) 
    end do

  END SUBROUTINE REFORMAT_VEGFILES

  !
  ! --------------------------------------------------------
  !
  ! SUBROUTINE compute_stats (ndata,cti_val,mu,sig,sk)
  !
  ! Subroutine not used as of 24 Dec 2024; removed by reichle, 24 Dec 2024
  !
  ! -----------------------------------------------------------------------------------
  
  SUBROUTINE ascat_r0 (nc,nr, ntiles,tile_id, z0)

    ! 1) ASCAT roughness 
    ! /discover/nobackup/adarmeno/projects/k14/arlems-roughness.x3600_y1800_t1.nc4

    integer,                      intent(in)  :: nc, nr
    integer,                      intent(in)  :: ntiles  
    INTEGER,                      intent(in)  :: tile_id(:,:)

    real, pointer, dimension (:), intent(inout) :: z0

    integer  , parameter               :: N_lon_ascat = 3600, N_lat_ascat = 1800
    integer                            :: i,j, status, varid, ncid
    REAL, ALLOCATABLE, dimension (:)   :: count_pix
    REAL, ALLOCATABLE, dimension (:,:) :: z0_grid, data_grid
    character*512                      :: fout

    ! READ ASCAT source data and regrid
    ! ---------------------------------

    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    status  = NF_OPEN (trim(MAKE_BCS_INPUT_DIR)//'/land/misc/roughness_length/v1/arlems-roughness.x3600_y1800_t1.nc4', NF_NOWRITE, ncid)

    allocate (z0_grid   (1 : NC         , 1 : NR))
    allocate (data_grid (1 : N_lon_ascat, 1 : N_lat_ascat)) 

    status  = NF_INQ_VARID (ncid,'roughness',VarID) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_REAL (ncid,VarID, (/1,1,1/),(/N_lon_ascat, N_lat_ascat,1/), data_grid) ; VERIFY_(STATUS)

    call RegridRasterReal(data_grid, z0_grid)

    status = NF_CLOSE(ncid)

    ! Grid to tile
    ! ------------

    ! Reading tile-id raster file

    allocate (z0        (1:NTILES))
    allocate (count_pix (1:NTILES))

    z0        = 0.
    count_pix = 0.

    do j = 1,nr
       do i = 1, nc
          if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

             ! z0 0. < 0.1
             if((z0_grid(i,j) >= 2.0e-6).and.(z0_grid(i,j) <= 0.1)) then 
                z0 (tile_id(i,j)) = z0 (tile_id(i,j)) + z0_grid(i,j)
                count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
             endif

          endif
       end do
    end do

    where (count_pix > 0.) z0 = z0/count_pix
    where (z0 == 0.)       z0 = 2.0e-6

    deallocate (count_pix)
    deallocate (z0_grid)

  END SUBROUTINE ascat_r0

  ! ----------------------------------------------------------------------------------------------------------------------------
  
  SUBROUTINE jpl_canoph (nc,nr, ntiles, tile_id, z2)

    ! 1) JPL Canopy Height 
    ! /discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/veg/veg_height/v1/Simard_Pinto_3DGlobalVeg_JGR.nc4

    integer,                     intent(in)    :: nc, nr, ntiles
    integer,                     intent(in)    :: tile_id(:,:)
    real, pointer, dimension(:), intent(inout) :: z2

    ! ----------------------------------------------------------

    integer  , parameter                  :: N_lon_jpl = 43200, N_lat_jpl = 21600

    ! ----------------------------------------------------------

    integer                               :: i,j, status, varid, ncid
    REAL,    ALLOCATABLE, dimension (:)   :: count_pix
    INTEGER, ALLOCATABLE, dimension (:,:) :: data_grid, z2_grid
    character*512                         :: fout

    ! READ JPL source data files and regrid
    ! -------------------------------------

    call get_environment_variable ("MAKE_BCS_INPUT_DIR",MAKE_BCS_INPUT_DIR)
    status  = NF_OPEN (trim(MAKE_BCS_INPUT_DIR)//'/land/veg/veg_height/v1/Simard_Pinto_3DGlobalVeg_JGR.nc4', NF_NOWRITE, ncid)

    allocate (z2_grid   (1 : NC         , 1 : NR))
    allocate (data_grid (1 : N_lon_jpl, 1 : N_lat_jpl)) 

    status  = NF_INQ_VARID (ncid,'CanopyHeight',VarID) ; VERIFY_(STATUS)
    status  = NF_GET_VARA_INT (ncid,VarID, (/1,1/),(/N_lon_jpl, N_lat_jpl/), data_grid) ; VERIFY_(STATUS)

    call RegridRaster(data_grid, z2_grid)

    status = NF_CLOSE(ncid)

    ! Grid to tile
    ! ------------

    ! Reading tile-id raster file

    allocate (z2        (1:NTILES))
    allocate (count_pix (1:NTILES))

    z2        = 0.
    count_pix = 0.

    do j = 1,nr
       do i = 1, nc
          if((tile_id(i,j).gt.0).and.(tile_id(i,j).le.NTILES)) then

             if(z2_grid(i,j) >= 0.) then 
                z2 (tile_id(i,j)) = z2 (tile_id(i,j)) + real (z2_grid(i,j))
                count_pix (tile_id(i,j)) = count_pix (tile_id(i,j)) + 1. 
             endif

          endif
       end do
    end do

    where (count_pix >   0.) z2 = z2/count_pix
    where (z2        < 0.01) z2 = 0.01            ! to ensure Z2 >= MIN_VEG_HEIGHT

    deallocate (count_pix)
    deallocate (z2_grid)

  END SUBROUTINE jpl_canoph

  ! ----------------------------------------------------------------------
  
  integer function NC_VarID (NCFID, VNAME) 

    integer, intent (in)      :: NCFID
    character(*), intent (in) :: VNAME
    integer                   :: status

    STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,NC_VarID)
    IF (STATUS .NE. NF_NOERR) &
         CALL HANDLE_ERR(STATUS, trim(VNAME))  

  end function NC_VarID

  ! -----------------------------------------------------------------------

  SUBROUTINE HANDLE_ERR(STATUS, Line)

    INTEGER,      INTENT (IN) :: STATUS
    CHARACTER(*), INTENT (IN) :: Line

    IF (STATUS .NE. NF_NOERR) THEN
       PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
       STOP 'Stopped'
    ENDIF

  END SUBROUTINE HANDLE_ERR

  ! -----------------------------------------------------------------------------------

END module rmTinyCatchParaMod
