#define VERIFY_(A)   IF(A/=0)THEN;PRINT *,'ERROR AT LINE ', __LINE__;STOP;ENDIF
#define ASSERT_(A)   if(.not.A)then;print *,'Error:',__FILE__,__LINE__;stop;endif

module module_irrig_params
! Version 1 - Sarith Mahanama
! Version 2 - Stefano Casirati - LAI min-max obtained from LAI climatology boundary conditions   

  use rmTinyCatchParaMod,  ONLY : RegridRaster,regridrasterreal
  use process_hres_data,   ONLY : get_country_codes
  use MAPL
  use irrigation_module,   ONLY : NCROPS => IRRG_NCROPS
  
  implicit none

  INCLUDE 'netcdf.inc'

  private

  public :: create_irrig_params

contains

  subroutine create_irrig_params (nc, nr, gfile)

    implicit none

    integer      , intent (in) :: nc, nr
    character(*) , intent (in) :: gfile
    REAL,          PARAMETER :: UNDEF = -9999., UNDEFG = 1.e15
    
    ! GRIPC data
    ! ----------
    
    integer,       parameter :: NX_gripc = 86400
    integer,       parameter :: NY_gripc = 43200,  NY_GripcData = 36000
    character*300, parameter :: GRIPC_file = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/crop_fraction_data/v1/irrigtype_salmon2013.flt'
    real, allocatable, dimension (:) :: IGRIPC, RGRIPC, PGRIPC, NGRIPC
    
    ! MIRCA2000 data
    ! --------------
    
    integer, parameter       :: NX_mirca = 4320
    integer, parameter       :: NY_mirca = 2160
    integer, parameter       :: NMON = 12, STRLEN = 20
    real,    parameter       :: DXY_mirca= 360./REAL(NX_mirca)
    real,    parameter       :: lat1_mirca =   90.0 - DXY_mirca / 2.0    !1st grid center lat
    real,    parameter       :: lon1_mirca = -180.0 + DXY_mirca / 2.0   !1st grid center lon 
    character*300, parameter :: MIRCA_pathIrr = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/crop_fraction_data/v1/irrigated/crop_' 
    character*300, parameter :: MIRCA_pathRain = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/crop_fraction_data/v1/rainfed/crop_' 
    real, allocatable, dimension(:,:,:) :: MIFRAC, MRFRAC

    ! Global Irrigated Area data (GIA)
    ! --------------------------------
    
    integer,       parameter :: NX_GIA = 43200
    integer,       parameter :: NY_GIA = 21600,  NY_GIAData = 18000  
    character*300, parameter :: GIA_file = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/irrigation_class/v1/global_irrigated_areas.nc4'
    real, allocatable, dimension (:) :: GIAFRAC

    ! LAI data
    ! --------
    
    integer,       parameter :: NX_LAI = 86400
    integer,       parameter :: NY_LAI = 43200
    character*300, parameter :: LAI_file = '/discover/nobackup/projects/lis/LS_PARAMETERS/MODIS/MCD15A2H.006/MCD15A2H.006_LAI_YYYY'

    ! Irrigation Method
    ! -----------------
    character*300, parameter :: IM_path = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/country_code_IMethod/v1/'
      
    ! Global/Local variables
    ! ----------------------
    
    integer, allocatable, dimension (:,:) :: tile_id
    integer                               :: i,j, NTILES, STATUS, tindex1,pfaf1, NCOutID
    real,    allocatable, dimension (:)   :: tile_lon, tile_lat
    real                                  :: minlon,maxlon,minlat,maxlat
    
    !The codes for the 26 crop classes are as follows:
    !1       Wheat                
    !2       Maize                
    !3       Rice                 
    !4       Barley               
    !5       Rye                  
    !6       Millet               
    !7       Sorghum              
    !8       Soybeans             
    !9       Sunflower            
    !10      Potatoes             
    !11      Cassava              
    !12      Sugar cane           
    !13      Sugar beet           
    !14      Oil palm             
    !15      Rape seed / Canola   
    !16      Groundnuts / Peanuts 
    !17      Pulses               
    !18      Citrus               
    !19      Date palm            
    !20      Grapes / Vine        
    !21      Cotton               
    !22      Cocoa                
    !23      Coffee               
    !24      Others perennial     
    !25      Fodder grasses       
    !26      Others annual            
    
    character(len=STRLEN),dimension(ncrops) :: cname = (/"Wheat               ", &
         "Maize               ", "Rice                ", "Barley              ", &
         "Rye                 ", "Millet              ", "Sorghum             ", &
         "Soybeans            ", "Sunflower           ", "Potatoes            ", &
         "Cassava             ", "Sugar cane          ", "Sugar beet          ", &
         "Oil palm            ", "Rape seed /Canola   ", "Groundnuts/ Peanuts ", &
         "Pulses              ", "Citrus              ", "Date palm           ", &
         "Grapes / Vine       ", "Cotton              ", "Cocoa               ", &
         "Coffee              ", "Other sperennial    ", "Fodder grasses      ", &
         "Others annual       "/)

    ! (1) Reading rst file and NTILES
    ! -------------------------------
    
    open (10,file= trim(gfile)//'.rst',status='old',action='read',  &
         form='unformatted',convert='little_endian')
    allocate (tile_id    (1:nc,1:nr))         
    
    do j=1,nr
       read(10)tile_id(:,j)
    end do
    close (10,status='keep')
    
    open (10,file='clsm/catchment.def',status='old',action='read', form='formatted')
    read (10, *) ntiles
    allocate (tile_lon (1:NTILES))
    allocate (tile_lat (1:NTILES))
    do i = 1, NTILES
       read (10,*) tindex1,pfaf1,minlon,maxlon,minlat,maxlat
       tile_lon(i) = (minlon + maxlon)/2.
       tile_lat(i) = (minlat + maxlat)/2.
    end do
    close (10, status = 'keep')    

    call OpenFile

    ! (3) Process GIA
    ! -----------------

    allocate (GIAFRAC  (NTILES))
    call ReadProcess_GIA (NC, NR, NTILES, tile_id, GIAFRAC)

    ! (4) Process GRIPC and LAI (Min-Max)
    ! ----------------------------------

    allocate (IGRIPC  (NTILES))
    allocate (RGRIPC  (NTILES))
    allocate (PGRIPC  (NTILES))
    allocate (NGRIPC  (NTILES))
    call ReadProcess_GRIPC (NC, NR, NTILES, tile_id, IGRIPC, RGRIPC, PGRIPC, NGRIPC)

    ! (1) Process MIRCA2000
    ! ---------------------
    
    allocate (MIFRAC (NTILES, 12, NCROPS))
    allocate (MRFRAC (NTILES, 12, NCROPS))
    call ReadProcess_MIRCA (NC, NR, NTILES, tile_id, MIFRAC, MRFRAC)

    ! (5) Create parameter file for the model
    ! ---------------------------------------

    call MergeData (NTILES)

    return

    contains

      ! ===========================================================================================
      
      SUBROUTINE OpenFile 

        implicit none
        integer           :: i, m, n,l, lid, mid, cid, vid, sid
        integer, dimension(8)               :: date_time_values
        character (22)                      :: time_stamp
        character (len=STRLEN)              :: ThisCrop
        real                                :: abm_int, peatf_r, gdp_r, hdm_r
        real, dimension (:), allocatable    :: field_cap
        
        status = NF_CREATE ('clsm/irrig.dat'  , NF_NETCDF4, NCOutID) ; VERIFY_(STATUS)
        status = NF_DEF_DIM(NCOutID, 'tile'        , NTILES, lid)    ; VERIFY_(STATUS)
        status = NF_DEF_DIM(NCOutID, 'unknown_dim1', NCROPS, cid)    ; VERIFY_(STATUS)
        status = NF_DEF_DIM(NCOutID, 'unknown_dim2', 2,      mid)    ; VERIFY_(STATUS)
        status = NF_DEF_DIM(NCOutID, 'strlen'      , STRLEN, sid)    ; VERIFY_(STATUS)
        
        ! GIA -> GRIPC MERGE
        ! ------------------
        
        status = NF_DEF_VAR(NCOutID, 'RAINFEDFRAC' , NF_FLOAT, 1 ,(/lid/), vid)      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('fraction of rainfed cropland'),    &
             'fraction of rainfed cropland')                                       ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                       ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)       ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_IRRIGFRAC'   , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('fraction of irrigated cropland'),  &
             'fraction of irrigated cropland')                                     ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                       ; VERIFY_(STATUS)   
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)       ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_PADDYFRAC'   , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('fraction of paddy cropland'),      &
        'fraction of paddy cropland')                                         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                       ; VERIFY_(STATUS)    
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)       ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'FIELDCAP'   , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('soil field capacity'),      &
        'soil field capacity')                                         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 5,'m3/m3')                       ; VERIFY_(STATUS)    
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)       ; VERIFY_(STATUS)
        
        ! GIA- GRIPC -> MIRCA
        ! -------------------
        
        status = NF_DEF_VAR(NCOutID, 'CROPCLASSNAME'   , NF_CHAR, 2 ,(/sid, cid/), vid) ; VERIFY_(STATUS)   
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('Crop Class Name'),      &
             'Crop Class Name')                                                    ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_CROPIRRIGFRAC'   , NF_FLOAT, 2 ,(/lid, cid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('Crop irrigated fraction'),      &
             'Crop irrigated fraction')                      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                       ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)       ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'CROPRAINFEDFRAC'   , NF_FLOAT, 2 ,(/lid, cid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('Crop rainfed fraction'),      &
             'Crop rainfed fraction')                        ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                       ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)       ; VERIFY_(STATUS)
        
        ! Crop calendar
        ! -------------
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_DOY_PLANT'   , NF_FLOAT, 3 ,(/lid, mid, cid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('DOY start planting'),      &
             'DOY start planting')                      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 4,'days')                     ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_DOY_HARVEST'   , NF_FLOAT, 3 ,(/lid, mid, cid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('DOY end harvesting'),      &
             'DOY end harvesting')                      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 4,'days')                     ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'RAINFEDPLANT'   , NF_FLOAT, 3 ,(/lid, mid, cid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('DOY start planting'),      &
             'DOY start planting')                      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 4,'days')                     ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'RAINFEDHARVEST'   , NF_FLOAT, 3,(/lid, mid, cid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('DOY end harvesting'),      &
             'DOY end harvesting')                      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 4,'days')                     ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        
        ! IRRIG TYPE
        ! ----------
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_TYPE'   , NF_FLOAT, 2 ,(/lid, cid/), vid)  ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name',  &
             LEN_TRIM('Preferred Irrig Type : Concurrent (0) SPRINKLER(1) DRIP(2) FLOOD(3) AVOID (negative)'),      &
             'Preferred Irrig Type : Concurrent (0) SPRINKLER(1) DRIP(2) FLOOD(3) AVOID (negative)')            ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                        ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_IRRIGFRAC_SPR' , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('fraction of sprinkler irrigation'),  &
             'fraction of sprinkler irrigation')                                    ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                        ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)          ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)        ; VERIFY_(STATUS)       
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_IRRIGFRAC_DRP' , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('fraction of drip irrigation'),  &
             'fraction of drip irrigation')                                         ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                        ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)          ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)        ; VERIFY_(STATUS)       
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_IRRIGFRAC_FRW' , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('fraction of flood irrigation'),  &
             'fraction of flood irrigation')                                        ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                        ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)          ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)        ; VERIFY_(STATUS)    
        
        ! LAI
        ! ---
        status = NF_DEF_VAR(NCOutID, 'IRRG_LAIMIN'   , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('Minimum LAI irrigated crops'),      &
             'Minimum LAI irrigated crops')                      ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                        ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)          ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)        ; VERIFY_(STATUS)
        
        status = NF_DEF_VAR(NCOutID, 'IRRG_LAIMAX'   , NF_FLOAT, 1 ,(/lid/), vid) ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'long_name', LEN_TRIM('Maximum LAI irrigated crops'),      &
             'Maximum LAI irrigated crops')                        ; VERIFY_(STATUS)
        status = NF_PUT_ATT_TEXT(NCOutID, vid, 'units', 1,'-')                        ; VERIFY_(STATUS)       
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'missing_value', NF_REAL,1,  UNDEFG)   ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'add_offset', NF_REAL,1,  0.)          ; VERIFY_(STATUS)
        status = NF_PUT_ATT_REAL(NCOutID, vid, 'scale_factor', NF_REAL,1,  1.)        ; VERIFY_(STATUS)    
        
        call date_and_time(VALUES=date_time_values)
        
        write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
             date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
             date_time_values(5),':',date_time_values(6),':',date_time_values(7)
        
        status = NF_PUT_ATT_TEXT(NCOutID, NF_GLOBAL, 'CreatedBy', LEN_TRIM('Sarith Mahanama, Stefano Casirati'),   &
             'Sarith Mahanama, Stefano Casirati') ; VERIFY_(STATUS) 
        status = NF_PUT_ATT_TEXT(NCOutID, NF_GLOBAL, 'Contact'   , LEN_TRIM('sarith.p.mahanama@nasa.gov'),                  &
             'sarith.p.mahanama@nasa.gov')
        status = NF_PUT_ATT_TEXT(NCOutID, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))
        
        status = NF_ENDDEF(NCOutID)  
        
        DO n = 1, NCROPS
           i = LEN_TRIM(cname(n))
           ThisCrop = trim(cname(n))
           status = NF_PUT_VARA_text(NCOutID,VarID(NCOutID,'CROPCLASSNAME') ,(/1,n/),(/i,1/),ThisCrop (1:i)) ; VERIFY_(STATUS) 
        END DO

        ! Put field capacity
        
        open (10,file='clsm/CLM4.5_abm_peatf_gdp_hdm_fc',  &
             form='formatted',status='unknown',action = 'read')
        allocate (field_cap(1:NTILES))
 
        do n = 1, NTILES
           read (10,'(2I10, i3, f8.4, f8.2, f10.2, f8.4)' ) i, vid, abm_int, peatf_r, gdp_r, hdm_r, field_cap(n)
        end do

        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'FIELDCAP'  ) ,(/1/),(/NTILES/),field_cap ) ; VERIFY_(STATUS)
        close (10, status = 'keep')
        deallocate (field_cap)

        
      END SUBROUTINE OpenFile

      ! -----------------------------------------------------------------------------------------
    
      SUBROUTINE MergeData (NTILES)
    
        implicit none   
        integer, intent (in)                :: NTILES
        real, allocatable, dimension (:,:)  :: MI, MR
        real                                :: MICROP, MRCROP, MIRICEA, MRRICEA, MICROPA, MRCROPA, DF, SF, FF, ITYPE(3)
        integer                             :: i, j, m, n,l, t                
        real, dimension (:), ALLOCATABLE    :: sprinkler, drip, flood    
        integer :: nc, day1, dayL, day1_2, dayL_2
        integer, dimension (12) :: fmonth, fmonth2, fmonth3
        integer, dimension (12)  :: DOY_MidMonth, DOY_BegMonth, DOY_EndMonth
        logical, dimension (4)  :: found = .false.
        integer, allocatable , dimension (:) :: crop_mons
        real, allocatable, dimension (:,:,:) :: IRRG_DOY_PLANT, IRRG_DOY_HARVEST, RAINFEDPLANT, RAINFEDHARVEST
        real, allocatable, dimension (:,:)   :: IRRG_TYPE
        
        data DOY_BegMonth / 1, 32, 60,  91, 121, 152, 182, 213, 244, 274, 305, 335/          
        data DOY_MidMonth /15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349/       
        data DOY_EndMonth /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366/   
        
        ALLOCATE (FLOOD    (1:NTILES))
        ALLOCATE (SPRINKLER(1:NTILES))
        ALLOCATE (DRIP     (1:NTILES))
        CALL ReadProcess_IMethod (NTILES, sprinkler, drip, flood)        
        
        ! MERGING PROCEDURE  
        ! =================
        !                                                                                                                                  CROP CALENDAR
        ! GIA                 GIA-GRIPC                           GIA-GRIPC-MIRCA                                                          -------------
        ! ---                 ---------                           ---------------                                            
        !                                                                                                                   CAL_STEP 1 : USE MIRCA monthly crop fractions :   
        !               -> NO  -> 100% iCrop              -> YES  scale to match GIA-GRIPC                                1) Year around : 366
        !              |                                 |       1) MIRCA type 3 is paddy                                     2) if not year around check whether 1 or 2 seasons
        !              |                                 |       2) Scale fractions to match                                  
        !              |                                 |          the sum of individual crops                             CAL_STEP 2 : Look for gaps in calendar for GIA-GRIPC-MIRCA and
        ! GIA -> GRIPC_|                       -> MIRCA -|          to GIA-GRIPC                                               use nearest neighbor with similar crop type
        !        HAVE  |                      |   iCrop  |      
        !              |                      |          |                -> YES Scale to match GIA-GRIPC
        !              |      |-> iCrop       |          |               |       1) MIRCA type 3 is paddy          
        !              |      |               |          |               |       2) Scale fractions to match    
        !              |      |               |          |               |          the sum of individual crops 
        !              |-> YES|-> Paddy       |          |-> NO ->MIRCA -|          to GIA-GRIPC          
        !                     |               |                   rCrop  |       3) set MIRCA rCrop to zero.
        !                     |               |                          |
        !                     |-> rCrop       |                          |-> NO Plant Wheat GIA-GRIPC frac    
        !                                     |                          
        !                                     |                          
        !                                     |           -> YES Max (GIA-GRIPC, MIRCA rCrop)                           
        !                                     |          |       1) MIRCA type 3 is paddy       
        !                                     |          |       2) Scale fractions to match    
        !                                     |          |          the sum of individual crops 
        !                                     | ->MIRCA -|          to GIA-GRIPC-MIRCA          
        !                                         rCrop  |       
        !                                                |-> NO  Plant Wheat 
        
        allocate (MI                 (1 : NTILES, 1 : NCROPS))
        allocate (MR                 (1 : NTILES, 1 : NCROPS))    
        allocate (IRRG_DOY_PLANT     (1 : NTILES, 1 : 2, 1 : NCROPS))
        allocate (IRRG_DOY_HARVEST   (1 : NTILES, 1 : 2, 1 : NCROPS))
        allocate (RAINFEDPLANT       (1 : NTILES, 1 : 2, 1 : NCROPS))
        allocate (RAINFEDHARVEST     (1 : NTILES, 1 : 2, 1 : NCROPS))
        allocate (IRRG_TYPE          (1 : NTILES, 1 : NCROPS))
        
        MI     = 0.
        MR     = 0.
        IRRG_TYPE = 0
    
        IRRG_DOY_PLANT     = 998
        IRRG_DOY_HARVEST   = 998
        RAINFEDPLANT       = 998
        RAINFEDHARVEST     = 998
        
        ! Compute annual maximum fractions from MIRCA monthly fractions : 
        !      (1) crop specific and (2) paddy and irrigated crops, seperately for rainfed and irrigated 
        
        DO m = 1,12              
           DO I = 1, NTILES
              IF((maxval(MIFRAC (I,m,:)) > 0.).OR.(maxval(MRFRAC (I,m,:)) > 0.)) then                   
                 DO N = 1, NCROPS                         
                    ! Maximum Crop Fraction over 12 months                         
                    IF(MI (I,n) < MIFRAC (I,m,n)) MI (I,N) = MIFRAC (I,m,n)
                    IF(MR (I,n) < MRFRAC (I,m,n)) MR (I,N) = MRFRAC (I,m,n)                         
                 END DO
              ENDIF
           END DO
        END DO
        
        ! CROP CALENDAR CAL_STEP1 : Compute plant/harvest dates and create HYBRID of GRIPC and MIRCA fractions
        ! ----------------------------------------------------------------------------------------------------
        
        DO I = 1, NTILES             
           IF((maxval(MI (I,:)) > 0.).OR.(maxval(MR (I,:)) > 0.).OR.(IGRIPC (I) > 0.).OR.(PGRIPC (I) > 0.).OR.(RGRIPC (I) > 0.)) THEN
              
              ! Crop planting/Harvesting days
              ! -----------------------------
              ! OUTPPUTS IRRG_DOY_PLANT, IRRG_DOY_HARVEST, RAINFEDPLANT, RAINFEDHARVEST
              
              forall (m=1:12) fmonth3(m) = m
              
              DO N = 1, NCROPS
                 
                 fmonth  = 0.
                 fmonth2 = 0.
                 
                 DO t = 1,2
                    
                    if (t == 1) nc = count (MIFRAC (i,:,n) > 0.)
                    if (t == 2) nc = count (MRFRAC (i,:,n) > 0.)
                    
                    if(nc > 0) then
                       
                       if(nc == 12) then 
                          ! year around
                          
                          day1 = 1
                          dayL = 366
                          day1_2 =998
                          dayL_2 =998
                          
                       else
                          
                          fmonth  = 0.
                          fmonth2 = 0.
                          day1 = 998
                          dayL = 998
                          day1_2 = 998
                          dayL_2 = 998
                          
                          if (t == 1) forall (m=1:12) fmonth(m) = ceiling (MIFRAC (i,m,n)) 
                          if (t == 2) forall (m=1:12) fmonth(m) = ceiling (MRFRAC (i,m,n)) 
                          
                          fmonth2(1) = 1
                          do m = 2,12
                             if(fmonth(m) == fmonth(m-1)) then 
                                fmonth2 (m) = fmonth2(m-1)
                             else 
                                fmonth2 (m) = fmonth2(m-1) + 1
                             endif
                          end do
                          
                          if(maxval (fmonth2) > 3) then
                             
                             ! This crop grows in 2 seasons
                             ! ............................
                             
                             allocate (crop_mons (1:NC))
                             crop_mons = pack(fmonth3, mask = (fmonth > 0.))
                             found = .false.
                             
                             if(fmonth(1) == 1) then
                                if(fmonth(12) == 0) then
                                   ! Season begins on Jan 1
                                   day1   = DOY_BegMonth(crop_mons(1))
                                   found(1) = .true.
                                   do m = 1, nc-1
                                      if((crop_mons(m+1) - crop_mons(m)) > 1) then
                                         dayL   = DOY_EndMonth(crop_mons(m))
                                         day1_2   = DOY_BegMonth(crop_mons(m+1))
                                         dayL_2   = DOY_EndMonth(crop_mons(nc))
                                         found(2) = .true.
                                         found(3) = .true.
                                         found(4) = .true.
                                         exit
                                      endif
                                   enddo                                   
                                else                                
                                   ! season one begins in the fall
                                   do m = 1, nc-1
                                      if((crop_mons(m+1) - crop_mons(m)) > 1) then
                                         if(.not.found(2)) then
                                            dayL   = DOY_EndMonth(crop_mons(m))
                                            day1_2 = DOY_BegMonth(crop_mons(m+1))
                                            found(2) = .true.
                                            found(3) = .true.
                                         elseif (.not.found(4)) then
                                            found(4) = .true.
                                            found(1) = .true.
                                            dayL_2 = DOY_EndMonth(crop_mons(m))
                                            day1   = DOY_BegMonth(crop_mons(m+1))
                                         endif
                                      endif
                                   end do
                                endif
                             else
                                
                                ! season 1 brings in the spring                                
                                day1   = DOY_BegMonth(crop_mons(1))
                                found(1) = .true.
                                do m = 1, nc-1
                                   if((crop_mons(m+1) - crop_mons(m)) > 1) then
                                      dayL   = DOY_EndMonth(crop_mons(m))
                                      day1_2   = DOY_BegMonth(crop_mons(m+1))
                                      dayL_2   = DOY_EndMonth(crop_mons(nc))
                                      found(2) = .true.
                                      found(3) = .true.
                                      found(4) = .true.
                                      exit
                                   endif
                                enddo
                             endif
                             deallocate (crop_mons)                             

                          else
                             
                             ! Single crop season
                             ! ..................
                             if((fmonth(1) == 0).and.(fmonth(12) == 0)) then
                                day1 = DOY_BegMonth (maxloc(fmonth, 1))
                                dayL = DOY_EndMonth (maxloc(fmonth2, 1)-1)
                             else
                                if((fmonth(1) == 1).and.(fmonth(12) == 1)) then
                                   day1 = DOY_BegMonth (maxloc(fmonth2, 1,mask=(fmonth2 > 2)))
                                   dayL = DOY_EndMonth (maxloc(fmonth2, 1,mask=(fmonth2 == 2))-1)
                                endif
                                if((fmonth(1) == 0).and.(fmonth(12) == 1)) then
                                   day1 = DOY_BegMonth (maxloc(fmonth2, 1))
                                   dayL = DOY_EndMonth (12)
                                endif
                                if((fmonth(1) == 1).and.(fmonth(12) == 0)) then
                                   day1 = DOY_BegMonth (1)
                                   dayL = DOY_EndMonth (maxloc(fmonth2, 1)-1)
                                endif
                             endif
                          endif
             
                          if (t == 1) then
                             IRRG_DOY_PLANT   (I,1,N) = day1
                             IRRG_DOY_PLANT   (I,2,N) = day1_2
                             IRRG_DOY_HARVEST (I,1,N) = dayL
                             IRRG_DOY_HARVEST (I,2,N) = dayL_2
                          else
                             RAINFEDPLANT   (I,1,N) = day1
                             RAINFEDPLANT   (I,2,N) = day1_2
                             RAINFEDHARVEST (I,1,N) = dayL
                             RAINFEDHARVEST (I,2,N) = dayL_2
                          endif
                       endif
                    endif
                 END DO
              END DO
              
              ! 1. Main Fractions (OUTPUT) :
              !    1.1 IRRG_IRRIGFRAC   : The maximum value between (1) GRIPC irrigfrac, and (2) sum of MIRCA monthly crop frations without rice
              !    1.2 IRRG_PADDYFRAC   : The maximum value between (1) GRIPC paddyfrac, and (2) monthly rice fractions from MIRCA
              !    1.3 RAINFEDFRAC : The maximum value between (1) GRIPC rainfedfrac, and (2) sum of MIRCA monthly crop frations
              !    1.4 MI (I,CROPS) : Irrigated crop fractions with rice is the 3rd slice crops = 3
              !    1.5 MR (I,CROPS) : rainfed crop fractions with rice is the 3rd slice crops = 3
              
              MICROP = 0.
              MRCROP = 0.                                    
              MIRICEA= 0.
              MRRICEA= 0.
              MICROPA= 0.
              MRCROPA= 0.
              
              DO N = 1, NCROPS
                 IF (n == 3) THEN
                    IF(MIRICEA < MI (I,n)) MIRICEA = MI (I,n)
                    IF(MRRICEA < MR (I,n)) MRRICEA = MR (I,n)
                 ELSE
                    MICROP = MICROP + MI (I,n)
                    MRCROP = MRCROP + MR (I,n)
                 ENDIF
              END DO
              
              IF(MICROPA < MICROP) MICROPA = MICROP
              IF(MRCROPA < MRCROP) MRCROPA = MRCROP
              
              ! GIA-GRIPC
              ! .........
              
              IF(GIAFRAC (I) <= 0 ) THEN
                 ! MASK OUT non-irrigated per GIA
                 RGRIPC (I) = RGRIPC (I) + IGRIPC (I) + PGRIPC (I)
                 IGRIPC (I) = 0.                 
                 PGRIPC (I) = 0.
              ELSE
                 IF ((IGRIPC (I)  + PGRIPC (I)) < 0.) THEN
                    ! GRIPC does not have data
                    PGRIPC (I)  = 0.
                    IGRIPC (I)  = GIAFRAC (I)
                 ELSE
                    ! GRIPC HAVE DATA
                    MICROP = PGRIPC (I) + IGRIPC (I) + 1.e-15 !  RGRIPC (I) 
                    IF (GIAFRAC (I) > MICROP) THEN
                       PGRIPC (I)  = PGRIPC (I) * GIAFRAC (I) / MICROP
                       IGRIPC (I)  = IGRIPC (I) * GIAFRAC (I) / MICROP
                    ENDIF
                 ENDIF
              ENDIF

              if ((RGRIPC(I) + IGRIPC(I) + PGRIPC(I) + NGRIPC(I)) > 0.) then
                 RGRIPC (I) = RGRIPC (I) /(RGRIPC(I) + IGRIPC(I) + PGRIPC(I) + NGRIPC(I))
                 NGRIPC(I)  = NGRIPC(I)  /(RGRIPC(I) + IGRIPC(I) + PGRIPC(I) + NGRIPC(I))
              endif
              
              ! GIA-GRIPC-MIRCA IRRIGATED CROPS
              ! ...............................
              
              IF (IGRIPC(I)  == 0) THEN
                 DO N = 1, NCROPS
                    IF(N /= 3) MI (I,N) = 0. 
                 END DO
              ELSE
                 !  IF(MICROPA > IGRIPC (I)) THEN
                 !  
                 !     ! MIRCA is the larger fraction
                 !     IGRIPC (I) = MICROPA
                 !     !                   CALL STOPIT (1, IGRIPC (I), MICROPA, MI (I,:))
                 !  ELSE
                 
                 IF(MICROPA > 0.) THEN
                    
                    ! MIRCA has data too, thus scale crop fractions to match GIA-GRIPC (i.e. GIA)
                    DO N = 1, NCROPS
                       IF(N /= 3) MI (I,N) = MI (I,N) * IGRIPC (I) /  MICROPA
                    END DO
                    !                      CALL STOPIT (2, IGRIPC (I), MICROPA, MI (I,:))
                    
                 ELSE
                    
                    ! MIRCA does not have data but GRIPC has
                    IF(MRCROPA > 0.) THEN
                       
                       ! Looks like MIRCA rainfed frac has data
                       DO N = 1, NCROPS
                          IF(N /= 3) THEN
                             MI (I,N) = MR (I,N) * IGRIPC (I) /  MRCROPA
                             MR (I,N) = 0.
                             IRRG_DOY_PLANT   (I,1,N) = RAINFEDPLANT   (I,1,N)
                             IRRG_DOY_PLANT   (I,2,N) = RAINFEDPLANT   (I,2,N)
                             IRRG_DOY_HARVEST (I,1,N) = RAINFEDHARVEST (I,1,N)
                             IRRG_DOY_HARVEST (I,2,N) = RAINFEDHARVEST (I,2,N)
                          ENDIF
                       END DO
                       !                         CALL STOPIT (3, IGRIPC (I), MICROPA, MI (I,:))
                       MRCROPA = 0.
                   
                    ELSE
                       
                       ! MIRCA irrigated and rainfed do not have data plant some wheat
                       MI (I,1) =  IGRIPC (I)
                       IRRG_DOY_PLANT   (I,1,1) = 999 
                       IRRG_DOY_PLANT   (I,2,1) = 0
                       IRRG_DOY_HARVEST (I,1,1) = 999
                       IRRG_DOY_HARVEST (I,2,1) = 0
                       !                         CALL STOPIT (4, IGRIPC (I), MICROPA, MI (I,:))
                    ENDIF
                    !   ENDIF
                 ENDIF
              ENDIF
              
              !                CALL STOPIT (5, IGRIPC (I), MICROPA, MI (I,:))
              
              ! GIA-GRIPC-MIRCA PADDY
              ! .....................
              IF(PGRIPC (I) == 0.) THEN
                 MI (I,3) = 0.
              ELSE
                 !   IF(MIRICEA > PGRIPC (I)) THEN
                 !   
                 !      ! MIRCA is the larger fraction
                 !      PGRIPC (I) = MIRICEA 
                 !      
                 !   ELSE
                 
                 IF(MIRICEA > 0.) THEN
                    
                    ! MIRCA has data too, thus scale crop fractions to match GRIPC
                    MI (I,3) = MI (I,3) * PGRIPC (I) /  MIRICEA
                    
                 ELSE
                    
                    ! MIRCA does not have data but GRIPC has
                    IF(MRRICEA > 0.) THEN
                       
                       ! Looks like MIRCA rainfed frac has rice
                       MI (I,3) = MR (I,3) * PGRIPC (I) /  MRRICEA
                       MR (I,3) = 0.
                       MRRICEA = 0.
                       IRRG_DOY_PLANT   (I,1,3) = RAINFEDPLANT   (I,1,3)
                       IRRG_DOY_PLANT   (I,2,3) = RAINFEDPLANT   (I,2,3)
                       IRRG_DOY_HARVEST (I,1,3) = RAINFEDHARVEST (I,1,3)
                       IRRG_DOY_HARVEST (I,2,3) = RAINFEDHARVEST (I,2,3)
                    ELSE
                       
                       ! MIRCA irrigated and rainfed do not have data plant rice to PGRIPC
                       MI (I,3) =  PGRIPC (I)
                       
                       ! Get crop planting days for the nearest neighbor later
                       
                       IRRG_DOY_PLANT   (I,1,3) = 999
                       IRRG_DOY_PLANT   (I,2,3) = 0
                       IRRG_DOY_HARVEST (I,1,3) = 999
                       IRRG_DOY_HARVEST (I,2,3) = 0
                       
                    ENDIF
                 ENDIF
                 !   ENDIF
              ENDIF
              
              ! GIA-GRIPC-MIRCA Rainfed CROPS
              ! .............................
              
              IF( RGRIPC (I) == 0.) THEN
                 MR (I,:) = 0.
              ELSE
                 
                 !   IF(MRCROPA + MRRICEA > RGRIPC (I)) THEN
                 !   
                 !      ! MIRCA is the larger fraction
                 !      RGRIPC (I) = MRCROPA + MRRICEA
                 !      
                 !   ELSE
                 IF(MRCROPA + MRRICEA > 0.) THEN
                    
                    ! MIRCA has data too, thus scale crop fractions to match GRIPC
                    DO N = 1, NCROPS
                       IF(N /= 3) MR (I,N) = MR (I,N) * RGRIPC (I) / (MRCROPA + MRRICEA)
                    END DO
                    
                    IF (MRRICEA > 0.) MR (I,3) = MR (I,3)  * RGRIPC (I) / (MRCROPA + MRRICEA)
                    
                 ELSE
                    
                    ! MIRCA does not have data but GRIPC has
                    ! MIRCA rainfed do not have data plant some wheat
                    MR (I,1) =  RGRIPC (I)   
                    RAINFEDPLANT   (I,1,1) = 999
                    RAINFEDPLANT   (I,2,1) = 0
                    RAINFEDHARVEST (I,1,1) = 999
                    RAINFEDHARVEST (I,2,1) = 0                    
                 ENDIF
              ENDIF
              ! ENDIF
              
              ! IRRG_TYPE
              
              DO N = 1, NCROPS
                 FF = FLOOD     (I)
                 SF = SPRINKLER (I)
                 DF = DRIP      (I)
                 IF(N ==  3) THEN  ! Rice
                    SF = 0.
                    DF = 0.
                    FF = 1.
                    IRRG_TYPE (I, N) = 3 ! Always flood
                 ENDIF
                 IF(N == 10) IRRG_TYPE (I, N) = -1 ! never sprinkler
                 IF(N == 22) IRRG_TYPE (I, N) = -1 ! never sprinkler
                 !IF(N == 10) SF = 0.                                   ! Date palm
                 !IF(N == 22) SF = 0.                                   ! Cocoa
                 !ITYPE = (/SF, DF, FF/)
                 !IRRG_TYPE (I, N) = maxloc(ITYPE, 1)
              END DO
           ENDIF
        END DO
        
        ! CROP CALENDAR CAL_STEP2 Update missing plant harvest dates
        ! ----------------------------------------------------------
        
        DO I = 1, NTILES
           DO N = 1,3,2
              
              ! fill missing crop plant/harvest DOYs in irrigated crops
              ! .......................................................
              
              IF((IRRG_DOY_PLANT(I,1,N) == 999).AND.(MI (I,N) > 0.)) THEN
                 l = getNeighbor (I,day_in = IRRG_DOY_PLANT  (:,1,N))
                 IRRG_DOY_PLANT  (I,1,N) = IRRG_DOY_PLANT  (l,1,N)
                 IRRG_DOY_HARVEST(I,1,N) = IRRG_DOY_HARVEST(l,1,N)
                 if(N == 1) then
                    IF(RAINFEDPLANT(I,1,N) == 999) THEN
                       RAINFEDPLANT  (I,1,N) = IRRG_DOY_PLANT  (l,1,N)
                       RAINFEDHARVEST(I,1,N) = IRRG_DOY_HARVEST(l,1,N)         
                    ENDIF
                 endif
              endif
              
              IF(N == 1) THEN
                 ! fill missing crop plant/harvest DOYs in rainfed crops
                 ! .....................................................
                 ! temperorily commented out to save time, because we don't irrigate here anyway
                 ! IF((RAINFEDPLANT(I,1,N) == 999).AND.(MR (I,N) > 0.)) THEN
                 !    print *,'RAINFEDPLANT(I,1,N)',I,N, MR (I,N)
                 !    l = getNeighbor (I,day_in = IRRG_DOY_PLANT  (:,1,N))
                 !    RAINFEDPLANT  (I,1,N) = IRRG_DOY_PLANT  (l,1,N)
                 !    RAINFEDHARVEST(I,1,N) = IRRG_DOY_HARVEST(l,1,N)   
                 ! endif
              ENDIF
           END DO
           if(  ((IRRG_DOY_PLANT  (I,1,1) == 999).and.(MI (i,1) > 0.))   .OR.            &
                ((IRRG_DOY_PLANT  (I,1,3) == 999).and.(MI (i,3) > 0.))          ) then
              print *, i, IRRG_DOY_PLANT  (I,1,1:3), MI (i,1:3)
              stop
           endif
        END DO
        
        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_IRRIGFRAC'  ) ,(/1/),(/NTILES/),IGRIPC ) ; VERIFY_(STATUS)
        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_PADDYFRAC'  ) ,(/1/),(/NTILES/),PGRIPC ) ; VERIFY_(STATUS)
        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'RAINFEDFRAC') ,(/1/),(/NTILES/),RGRIPC ) ; VERIFY_(STATUS)      
         
        do n = 1,NCROPS
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_CROPIRRIGFRAC'    ) ,(/1,n  /),(/NTILES,1  /), MI              (:,n)  ) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'CROPRAINFEDFRAC'  ) ,(/1,n  /),(/NTILES,1  /), MR              (:,n)  ) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_TYPE'        ) ,(/1,n  /),(/NTILES,1  /), IRRG_TYPE       (:,n)  ) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_DOY_PLANT'   ) ,(/1,1,n/),(/NTILES,1,1/), IRRG_DOY_PLANT  (:,1,n)) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_DOY_PLANT'   ) ,(/1,2,n/),(/NTILES,1,1/), IRRG_DOY_PLANT  (:,2,n)) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_DOY_HARVEST' ) ,(/1,1,n/),(/NTILES,1,1/), IRRG_DOY_HARVEST(:,1,n)) ; VERIFY_(STATUS)         
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_DOY_HARVEST' ) ,(/1,2,n/),(/NTILES,1,1/), IRRG_DOY_HARVEST(:,2,n)) ; VERIFY_(STATUS)
           
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'RAINFEDPLANT'     ) ,(/1,1,n/),(/NTILES,1,1/), RAINFEDPLANT    (:,1,n)) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'RAINFEDPLANT'     ) ,(/1,2,n/),(/NTILES,1,1/), RAINFEDPLANT    (:,2,n)) ; VERIFY_(STATUS)
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'RAINFEDHARVEST'   ) ,(/1,1,n/),(/NTILES,1,1/), RAINFEDHARVEST  (:,1,n)) ; VERIFY_(STATUS)          
           status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'RAINFEDHARVEST'   ) ,(/1,2,n/),(/NTILES,1,1/), RAINFEDHARVEST  (:,2,n)) ; VERIFY_(STATUS)          
        end do
        
        status = NF_CLOSE(NCOutID) 
    
      END SUBROUTINE MergeData

      !----------------------------------------------------------------------------------------

      SUBROUTINE ReadProcess_IMethod (NTILES, f_sprink, f_drip, f_flood)
        
        implicit none
        
        integer, INTENT (IN)                :: NTILES
        real, dimension (:), INTENT(INOUT)  :: f_sprink, f_drip, f_flood
        integer                             :: i,j,n, k, N_METHOD, cnt_code, st_code
        character*2                         :: ST_NAME
        character*3                         :: CNT_ABR    
        integer, parameter                  :: N_STATES = 50, N_COUNTRY = 256
        real                                :: s_dum, d_dum, f_dum
        real, dimension (:), allocatable    :: us_sprink, us_drip, us_flood, us_tarea
        real, dimension (:), allocatable    :: sprink, drip, flood, tarea
        integer,      dimension (:),pointer :: index_range
        character*2,  dimension (:),pointer :: ST_NAME_ABR 
        character*3,  dimension (:),pointer :: CNT_NAME_ABR
        
        ! Read state fractions
        print *, '   '
        print *, '.........................................................................'   
        print *, 'PROCESSING IRRIGATION METHOD DATA '
        
        open (10, file = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/country_code_IMethod/v1/US_IMethod.2015'    , form = 'formatted', status ='old', action = 'read')
        open (11, file = '/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/country_code_IMethod/v1/Global_IMethod.data', form = 'formatted', status ='old', action = 'read')    
        
        READ (11, *) N_METHOD
        
        call get_country_codes (index_range=index_range, ST_NAME_ABR = ST_NAME_ABR, CNT_NAME_ABR = CNT_NAME_ABR)

        allocate (sprink   (0:N_COUNTRY))
        allocate (drip     (0:N_COUNTRY))
        allocate (flood    (0:N_COUNTRY))
        allocate (tarea    (0:N_COUNTRY))
        allocate (us_sprink(1:N_STATES ))
        allocate (us_drip  (1:N_STATES ))
        allocate (us_flood (1:N_STATES ))
        allocate (us_tarea (1:N_STATES ))
        
        sprink = 0.
        drip   = 0.
        flood  = 0.
        tarea  = 0.
        us_sprink = 0.
        us_drip   = 0.
        us_flood  = 0.
        us_tarea  = 0.
        
        do i = 1, N_METHOD
           read (11, *) CNT_ABR,s_dum, d_dum, f_dum
           do k = 1, N_COUNTRY
              if(CNT_ABR == CNT_NAME_ABR(k)) then
                 sprink(index_range(k)) = s_dum
                 drip  (index_range(k)) = d_dum
                 flood (index_range(k)) = f_dum
              endif
           end do
        end do

        tarea  = sprink + drip + flood
        where (tarea > 0.) 
           sprink = sprink / tarea
           drip   = drip   / tarea
           flood  = flood  / tarea
        endwhere

        do i = 1, N_STATES
           read (10, *) ST_NAME,s_dum, d_dum, f_dum
           do k = 1, N_STATES
              if(ST_NAME == ST_NAME_ABR(k)) then
                 us_sprink(k) = s_dum
                 us_drip  (k) = d_dum
                 us_flood (k) = f_dum
              endif
           end do
        end do
                
        us_tarea  = us_sprink + us_drip + us_flood
        where (us_tarea > 0.)
           us_sprink = us_sprink / us_tarea
           us_drip   = us_drip   / us_tarea
           us_flood  = us_flood  / us_tarea
        endwhere
        
        close (10, status = 'keep')
        close (11, status = 'keep')

        ! map irrig method fractions
        
        open (10,file='clsm/country_and_state_code.data',  &
             form='formatted',status='old', action = 'read')    
        
        ! allocate and initialize
                
        f_flood  = 1.
        f_sprink = 0.
        f_drip   = 0.
        
        tile_loop : do n = 1, NTILES
           read (10, '(i10, 2I4)') j, cnt_code, st_code 
           if (cnt_code < 257) then
              if(tarea (cnt_code) > 0.) then
                 f_flood (n) = flood (cnt_code)
                 f_sprink(n) = sprink(cnt_code)
                 f_drip  (n) = drip  (cnt_code)
              endif
           endif
           
           ! overwrite with US state methods
           if(st_code /= 999) then
              f_flood (n) = us_flood (st_code)
              f_sprink(n) = us_sprink(st_code)
              f_drip  (n) = us_drip  (st_code)                     
           endif
        end do tile_loop

        close (10, status = 'keep')
        
        CALL update_IMethod_bycounty (NTILES, f_sprink, f_drip, f_flood)

        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_IRRIGFRAC_SPR') ,(/1/),(/NTILES/), f_sprink) ; VERIFY_(STATUS)
        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_IRRIGFRAC_DRP'     ) ,(/1/),(/NTILES/), f_drip  ) ; VERIFY_(STATUS)
        status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_IRRIGFRAC_FRW'    ) ,(/1/),(/NTILES/), f_flood ) ; VERIFY_(STATUS)
        
        print *, 'DONE PROCESSING IRRIGATION METHOD DATA '

        deallocate (us_sprink, us_drip, us_flood, us_tarea, sprink, drip, flood, tarea)

      END SUBROUTINE ReadProcess_IMethod
 
      !----------------------------------------------------------------------------------------

      SUBROUTINE update_IMethod_bycounty (NTILES, f_sprink, f_drip, f_flood)

        implicit none
        integer, INTENT (IN)                :: NTILES
        real, dimension (:), INTENT(INOUT)  :: f_sprink, f_drip, f_flood
        integer,       parameter            :: NX_cb = 43200, NY_cb = 21600,  NY_cbData = 10800
        integer,       parameter            :: cb_states = 72, cb_county = 900, cb_countyUS = 3220
        integer                             :: i,j, n, status, ncid, I0(1), j0(1),SS, CCC
        real,    dimension(:,:),allocatable :: SFR, DFR, FFR
        integer, dimension  (:),allocatable :: GEOID
        integer, dimension(:,:),allocatable :: POLYID
        real (kind =8)                      :: XG(NX_cb),YG(NY_cb), y0, x0, dxy
        integer                             :: ii(NX_cb),jj(NY_cb)
        
        allocate (SFR (1:cb_county,1:cb_states))
        allocate (DFR (1:cb_county,1:cb_states))
        allocate (FFR (1:cb_county,1:cb_states))
        allocate (GEOID   (1:cb_countyUS))
        allocate (POLYID(1:NX_cb,1:NY_cb))

        POLYID = -9999

        status = NF_OPEN ('/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/irrigation/fraction_drip_flood_sprinkler/v1/cb_2015_us_county_30arcsec.nc4',NF_NOWRITE, ncid) ; VERIFY_(STATUS)    
        do j = 1, NY_cbData
           status = NF_GET_VARA_INT(NCID,VarID(NCID,'POLYID') ,(/1,j/),(/NX_cb, 1/), POLYID (:,NY_cb - j + 1)) ; VERIFY_(STATUS) ! reading north to south
        end do
        do j = 1, cb_states
           status = NF_GET_VARA_REAL(NCID,VarID(NCID,'IRRG_IRRIGFRAC_SPR') ,(/1,j/),(/cb_county, 1/), SFR (:,j)) ; VERIFY_(STATUS)
           status = NF_GET_VARA_REAL(NCID,VarID(NCID,'IRRG_IRRIGFRAC_DRP'     ) ,(/1,j/),(/cb_county, 1/), DFR (:,j)) ; VERIFY_(STATUS)
           status = NF_GET_VARA_REAL(NCID,VarID(NCID,'IRRG_IRRIGFRAC_FRW'    ) ,(/1,j/),(/cb_county, 1/), FFR (:,j)) ; VERIFY_(STATUS) 
        end do
        status = NF_GET_VARA_INT(NCID,VarID(NCID,'GEOID'    ) ,(/1/),(/cb_countyUS/), GEOID) ; VERIFY_(STATUS) 
        status = NF_CLOSE(NCID) ; VERIFY_(STATUS)

        dxy = 360.d0/NX_cb
        do i = 1, NX_cb 
           xg(i) = (i-1)*dxy -180.d0 + dxy/2.d0
        end do
        do i = 1, NY_cb
           yg(i) = (i-1)*dxy -90.d0 + dxy/2.d0
        end do

        do n = 1, NTILES
           
           x0 = dble (tile_lon(n))
           y0 = dble (tile_lat(n))
           II = 0
           JJ = 0
           WHERE ((xg >= x0).and.(xg < x0 + dxy)) II = 1
           WHERE ((yg >= y0).and.(yg < y0 + dxy)) JJ = 1
           
           I0 = FINDLOC(II,1)
           J0 = FINDLOC(JJ,1)

           if((POLYID(I0(1), J0(1)) >= 1).AND.(POLYID(I0(1), J0(1)) <= cb_countyUS)) then
              SS = GEOID(POLYID(I0(1),J0(1))) / 1000
              CCC= GEOID(POLYID(I0(1),J0(1))) - SS*1000
              f_sprink (n) = SFR (CCC,SS)
              f_drip   (n) = DFR (CCC,SS)
              f_flood  (n) = FFR (CCC,SS)
           endif
           
        END DO
        
        deallocate (SFR, DFR, FFR, GEOID, POLYID)
        
      END SUBROUTINE update_IMethod_bycounty
      
      !----------------------------------------------------------------------------------------

      SUBROUTINE ReadProcess_MIRCA (NC, NR, NTILES, tile_id, MIFRAC, MRFRAC)

        implicit none
    
        INTEGER, INTENT (IN)    :: NTILES, NC, NR
        INTEGER, INTENT (IN), DIMENSION(:,:):: tile_id 
        REAL,DIMENSION(:,:,:),INTENT(INOUT) :: MIFRAC, MRFRAC
        character*2                         :: TT
        integer                             :: i, j, n, m 
        real,dimension (:), allocatable     :: read_ir , read_rn, cnt_pix1, cnt_pix2
        real,dimension (:,:,:), allocatable :: mon_ir  , mon_rn
        real, pointer, dimension (:,:)      :: var_raster1, var_raster2
        real,parameter                      :: radius = MAPL_radius, pi = MAPL_PI
        real                                :: D2R, latc, lonc, area
        
        ! --------- VARIABLES FOR *OPENMP* PARALLEL ENVIRONMENT ------------
        !
        ! NOTE: "!$" is for conditional compilation
        !
        logical :: running_omp = .false.
        !
        !$ integer :: omp_get_thread_num, omp_get_num_threads
        !
        integer :: n_threads=1, li, ui, t_count, n_used
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
        
        MIFRAC = 0.
        MRFRAC = 0.

        n_used = MIN(n_threads, NCROPS/2)
        if(NTILES > 6000000) n_used = 1 ! otherwise multiple threads will run out of virtual memory
        allocate(low_ind(n_used))
        allocate(upp_ind(n_used))
        low_ind(1)         = 1
        upp_ind(n_used) = NCROPS
        
        if (running_omp)  then
           do i=1,n_used-1  
              upp_ind(i)   = low_ind(i) + (NCROPS/n_used) - 1 
              low_ind(i+1) = upp_ind(i) + 1
           end do
        end if

        D2R= PI/180.
                
        !$OMP PARALLEL DEFAULT (NONE) &
        !$OMP COPYIN (read_ir, read_rn,var_raster1, var_raster2, mon_ir, mon_rn, cnt_pix1, cnt_pix2) &
        !$OMP SHARED( n_used, low_ind, upp_ind, tile_id, NTILES, MIFRAC,MRFRAC, NC, NR, D2R) &
        !$OMP PRIVATE(n,i,j,m,tt,t_count, read_ir, read_rn,var_raster1, var_raster2, mon_ir, mon_rn, cnt_pix1, cnt_pix2, latc, area) 
        
        allocate (read_ir     (1:NX_mirca*12))
        allocate (read_rn     (1:NX_mirca*12))
        allocate (mon_ir      (1:NX_mirca,1:NY_mirca,1:12)) 
        allocate (mon_rn      (1:NX_mirca,1:NY_mirca,1:12))       
        allocate (var_raster1 (1:nc,1:nr))
        allocate (var_raster2 (1:nc,1:nr))
        allocate (cnt_pix1    (1:NTILES))
        allocate (cnt_pix2    (1:NTILES))
        
        !$OMP DO
        
        DO t_count = 1,n_used
           
           CROP_TYPE : DO n = low_ind(t_count),upp_ind(t_count)
              
              write (tt, '(i2.2)') n
              
              print *, '   '
              print *, '.........................................................................'
              print *, 'PROCESSING MIRCA : crop_', tt,'_irrigated_12.flt' 
              print *, 'PROCESSING MIRCA : crop_', tt,'_rainfed_12.flt' 
              
              mon_ir = UNDEF
              mon_rn = UNDEF
              
              open (50 + t_count, file = trim(MIRCA_pathIrr)//tt//'_irrigated_12.flt', action = 'read',        &
                   form = 'unformatted', access='direct', recl=NX_mirca*12)
              open (100 + t_count, file = trim(MIRCA_pathRain)//tt//'_rainfed_12.flt',   action = 'read',        &
                   form = 'unformatted', access='direct', recl=NX_mirca*12)
              
              mirca_rows : do j = 1, NY_mirca
                 read(50  + t_count,rec= NY_mirca - J + 1) read_ir
                 read(100 + t_count,rec= NY_mirca - J + 1) read_rn
                 latc = lat1_mirca - (j-1) * DXY_MIRCA
                 area = (sin(d2r*(latc+0.5*dxy_mirca)) - sin(d2r*(latc-0.5*dxy_mirca)))*(dxy_mirca*d2r)
                 area = area * radius * radius / 10000. ! in ha
                 
                 mirca_COLS : do i = 1, NX_MIRCA
                    mon_ir (i,j,:) = read_ir ((i-1)*12 + 1:  (i-1)*12 + 12)/area
                    mon_rn (i,j,:) = read_rn ((i-1)*12 + 1:  (i-1)*12 + 12)/area
                 end do mirca_cols
                 
              end do mirca_rows
              
              close (50  + t_count, status ='keep')
              close (100 + t_count, status ='keep')
              
              !  Grid 2 tile
              !  -----------
              do m = 1,12
                 
                 var_raster1 = 0.
                 var_raster2 = 0.
               
                 call RegridRasterReal(mon_ir(:,:,m),var_raster1)
                 call RegridRasterReal(mon_rn(:,:,m),var_raster2)
                 
                 cnt_pix1 = 0.
                 cnt_pix2 = 0.
                 
                 do j = 1,nr
                    do i = 1,nc                  
                       if((var_raster1 (i,j) > 0.).and.(tile_id (i,j) >= 1).AND.(tile_id (i,j) <= NTILES)) then
                          MIFRAC (tile_id(i,j),m,n) = MIFRAC (tile_id(i,j),m,n) + var_raster1 (i,j)
                          cnt_pix1 (tile_id(i,j))   = cnt_pix1 (tile_id(i,j)) + 1. 
                       endif
                       if((var_raster2 (i,j) > 0.).and.(tile_id (i,j) >= 1).AND.(tile_id (i,j) <= NTILES)) then
                          MRFRAC (tile_id(i,j),m,n) = MRFRAC (tile_id(i,j),m,n) + var_raster2 (i,j)
                          cnt_pix2 (tile_id(i,j))   = cnt_pix2 (tile_id(i,j)) + 1. 
                       endif
                    end do
                 end do
                 
                 do i = 1, NTILES
                    if(cnt_pix1(i) > 0.) MIFRAC (i,m,n) = MIFRAC (i,m,n)/cnt_pix1(i)
                    if(cnt_pix2(i) > 0.) MRFRAC (i,m,n) = MRFRAC (i,m,n)/cnt_pix2(i)
                 end do
               
              end do
           END DO CROP_TYPE
        END DO
        !$OMP END DO
        !$OMP END PARALLEL   
        
      print *,'DONE MIRCA PROCESSING'
    END SUBROUTINE ReadProcess_MIRCA
    
    !----------------------------------------------------------------------------------------
    
    SUBROUTINE ReadProcess_GRIPC (NC, NR, NTILES, tile_id,IGRIPC, RGRIPC, PGRIPC, NGRIPC)
      
      INTEGER, INTENT (IN)    :: NTILES, NC, NR
      INTEGER, INTENT (IN), DIMENSION(:,:):: tile_id 
      REAL,DIMENSION(:),INTENT(INOUT)     :: IGRIPC, RGRIPC, PGRIPC, NGRIPC
      real, allocatable     :: var_in (:,:), tot_cnt (:), min_cnt(:), max_cnt(:)           
      integer               :: i,j, n, r, status, DOY, NCID, NCIDW,xid, yid,vid
      integer, pointer      :: iraster  (:,:)
      character*3           :: DDD
      integer*2, allocatable, dimension (:,:) :: Lai_clim
      real,      allocatable, dimension (:,:) :: clim_min, clim_max, clim_lai
      real, allocatable, dimension (:)        :: LAI_MIN, LAI_MAX, lai
      real,allocatable :: dum,yr,mn,dy,nt
      logical                                 :: write_lai = .false.
    
    !V2: Min Max LAI from LAI Climatology for consistence    
     
      allocate( var_in(NX_gripc,NY_gripc)) 
      var_in = UNDEF
      
      open ( 10, file = trim(GRIPC_file), form = 'unformatted', access='direct', recl=(NX_gripc))
      
      !- Read input file::
      
      do j = 1, NY_gripcdata
         r = NY_gripc -j + 1
         read(10,rec=j) var_in(:, r)
         do i = 1, NX_gripc
            if( var_in(i, r) == 0. ) var_in(i, r) = -9999.
            if( var_in(i, r) == 4. ) var_in(i, r) = -9999.
         end do
      end do
      close( 10 )
      
      allocate(iraster(NX_gripc,NY_gripc),stat=STATUS); VERIFY_(STATUS)
      call RegridRaster(tile_id,iraster)
      
      allocate (tot_cnt  (1:ntiles))
      allocate (min_cnt  (1:ntiles))
      allocate (max_cnt  (1:ntiles))
      
      RGRIPC  = 0.
      IGRIPC  = 0.
      PGRIPC  = 0.
      tot_cnt = 0.
      min_cnt = 0.
      max_cnt = 0.
      allocate (LAI_MIN (NTILES))
      allocate (LAI_MAX (NTILES))
      allocate (lai (NTILES))
      LAI_MIN = -9999.
      LAI_MAX = -9999.
     
      do j = 1,NY_gripc
         do i =  1,NX_gripc
            if((iraster (i,j) >=1).and.(iraster (i,j) <=ntiles)) then
               tot_cnt (iraster (i,j)) = tot_cnt (iraster (i,j)) + 1.
               if (var_in(i,j) == 1) RGRIPC(iraster (i,j)) = RGRIPC(iraster (i,j)) + 1.
               if (var_in(i,j) == 2) IGRIPC(iraster (i,j)) = IGRIPC(iraster (i,j)) + 1.
               if (var_in(i,j) == 3) PGRIPC(iraster (i,j)) = PGRIPC(iraster (i,j)) + 1.
               if (var_in(i,j) == 4) NGRIPC(iraster (i,j)) = NGRIPC(iraster (i,j)) + 1.

            endif
         end do
      end do
      
      RGRIPC = RGRIPC / tot_cnt
      IGRIPC = IGRIPC / tot_cnt
      PGRIPC = PGRIPC / tot_cnt
      NGRIPC = NGRIPC / tot_cnt

   allocate(dum)
   allocate(yr)
   allocate(mn)
   allocate(dy)
   allocate(nt)

        open (43,file='clsm/lai.dat',      &
        form='unformatted',status='unknown',convert='little_endian',action='read')
        LAI_MAX=-9999.
        LAI_MIN=9999.
        do j = 1, 48
                read(43) yr,mn,dy,dum,dum,dum,yr,mn,dy,dum,dum,dum,nt,dum
	        read(43) lai
                do i = 1, NTILES
                      	if (lai(i)>LAI_MAX(i)) LAI_MAX(i)=lai(i)
		        if (lai(i)<LAI_MIN(i)) LAI_MIN(i)=lai(i)
                end do
        end do
        close (43)
        where (LAI_MIN == 9999.) LAI_MIN=-9999.

      status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_LAIMAX'  ) ,(/1/),(/NTILES/),LAI_MAX ) ; VERIFY_(STATUS)
      status = NF_PUT_VARA_REAL(NCOutID,VarID(NCOutID,'IRRG_LAIMIN'  ) ,(/1/),(/NTILES/),LAI_MIN  ) ; VERIFY_(STATUS)
      
      deallocate (var_in, iraster, min_cnt, max_cnt, tot_cnt, LAI_MIN, LAI_MAX, lai, yr, mn, dy, dum, nt)

      print *,'DONE PROCESSING GRIPC and LAI'
      
    END SUBROUTINE ReadProcess_GRIPC

    ! ----------------------------------------------------------------------

    SUBROUTINE ReadProcess_GIA (NC, NR, NTILES, tile_id, GIAFRAC)
    
      implicit none

      INTEGER, INTENT (IN)    :: NTILES, NC, NR
      INTEGER, INTENT (IN), DIMENSION(:,:):: tile_id 
      REAL,DIMENSION(NTILES),INTENT(INOUT)  :: GIAFRAC

      integer, allocatable, dimension (:,:) :: var_in 
      integer, allocatable, dimension (:,:) :: irrig
      real,dimension (:), allocatable       :: cnt_pix1
      integer    :: i,j, status, NCID

      ! Read GIA data
    ! -------------

      allocate( var_in(NX_GIA,NY_GIA)) 
      var_in = UNDEF

      print *, 'PROCESSING GIA : ', trim (GIA_file)
      status = NF_OPEN (trim (GIA_file),NF_NOWRITE, NCID) ; VERIFY_(STATUS)    
      
      do j = NY_GIAData, 1, -1
         status = NF_GET_VARA_INT(NCID,VarID(NCID,'IrrigClass') ,(/1,j/),(/NX_GIA, 1/), var_in (:,j + 3600 )) ; VERIFY_(STATUS)    

      end do
      
      status = NF_CLOSE(NCID) ; VERIFY_(STATUS)
      allocate (irrig (1:NC, 1:NR))
      irrig = undef
      
      if (NC /= NX_GIA) then
         call RegridRaster (var_in, irrig)
      else
         irrig = var_in
      endif
      
      ! Compute Fractions on tiles
      ! --------------------------
            
      allocate (cnt_pix1 (NTILES))
      cnt_pix1 = 0.
      GIAFRAC  = 0.
      
      geos_rows : do j = 1, NR
         geos_cols : do i = 1, NC
            if ((tile_id(i,j) > 0).and.(tile_id (i,j) <= NTILES)) then
               cnt_pix1(tile_id(i,j)) = cnt_pix1(tile_id(i,j)) + 1.
               if((irrig (i,j) > 0) .AND. (irrig (i,j) < 4)) then
                  GIAFRAC(tile_id (i,j)) = GIAFRAC(tile_id (i,j)) + 1.
               endif
            endif
         enddo geos_cols
      end do geos_rows
      where (cnt_pix1 > 0) GIAFRAC = GIAFRAC / cnt_pix1
      deallocate (irrig, var_in)
      deallocate (cnt_pix1)
      print *,'DONE PROCESSING GIA'
      RETURN
    END SUBROUTINE ReadProcess_GIA
      
    ! ----------------------------------------------------------------------

    integer function VarID (NCFID, VNAME) 
      
      integer, intent (in)      :: NCFID
      character(*), intent (in) :: VNAME
      integer                   :: status
      
      STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
      IF (STATUS .NE. NF_NOERR) &
           CALL HANDLE_ERR(STATUS, trim(VNAME))  
      
    end function VarID
    
    ! -----------------------------------------------------------------------
    
    SUBROUTINE HANDLE_ERR(STATUS, Line)
      
      INTEGER,      INTENT (IN) :: STATUS
      CHARACTER(*), INTENT (IN) :: Line
      
      IF (STATUS .NE. NF_NOERR) THEN
         PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
         STOP 'Stopped'
      ENDIF
      
    END SUBROUTINE HANDLE_ERR
    
    ! -----------------------------------------------------------------------------
   
     integer function getNeighbor (tid_in, lai_in, day_in)
     
       implicit none
       integer, intent (in)                :: tid_in
       real, optional, dimension (NTILES)  :: lai_in, day_in
       integer                             :: i, nplus    
       logical                             :: tile_found
       logical, allocatable, dimension (:) :: mask
       integer, allocatable, dimension (:) :: sub_tid
       real   , allocatable, dimension (:) :: sub_lon, sub_lat, rev_dist
       real                                :: dw, min_lon, max_lon, min_lat, max_lat
       integer, allocatable, dimension (:) :: TILEID
              
       allocate (mask   (1:  NTILES))
       allocate (TILEID (1:  NTILES))
       forall (i=1:NTILES) TILEID (i) = i
              
       dw = 0.5
       getNeighbor = -9999
       
       ZOOMOUT : do  
          
          tile_found = .false. 
          
          ! Min/Max lon/lat of the working window
          ! -------------------------------------
          
          min_lon = MAX(tile_lon (tid_in) - dw, -180.)
          max_lon = MIN(tile_lon (tid_in) + dw,  180.)
          min_lat = MAX(tile_lat (tid_in) - dw,  -90.)
          max_lat = MIN(tile_lat (tid_in) + dw,   90.) 
          
          mask = .false.
          if(present (lai_in)) then
             mask =  ((tile_lat >= min_lat .and. tile_lat <= max_lat).and.(tile_lon >= min_lon .and. tile_lon <= max_lon).and.(lai_in >= 0.))
          endif
          if(present (day_in)) then
             mask =  ((tile_lat >= min_lat .and. tile_lat <= max_lat).and.(tile_lon >= min_lon .and. tile_lon <= max_lon).and.(day_in < 998))
          endif
          nplus =  count(mask = mask)
          
          if(nplus < 0) then
             dw = dw + 0.5
             CYCLE
          endif
          
          allocate (sub_tid (1:nplus))
          allocate (sub_lon (1:nplus))
          allocate (sub_lat (1:nplus))
          allocate (rev_dist  (1:nplus))
          
          sub_tid = PACK (TILEID  , mask= mask) 
          sub_lon = PACK (tile_lon, mask= mask)
          sub_lat = PACK (tile_lat, mask= mask)
          
          ! compute distance from the tile
          
          sub_lat = sub_lat * MAPL_PI/180.
          sub_lon = sub_lon * MAPL_PI/180.
          
          SEEK : if(getNeighbor < 0) then
             
             rev_dist  = 1.e20
             
             do i = 1,nplus
                
                rev_dist(i) = haversine(to_radian(tile_lat(tid_in)), to_radian(tile_lon(tid_in)), &
                     sub_lat(i), sub_lon(i))
                
             end do
             
             FOUND : if(minval (rev_dist) < 1.e19) then
                if(present (lai_in)) then
                   if(lai_in(sub_tid(minloc(rev_dist,1))) >= 0.) then
                      getNeighbor = sub_tid(minloc(rev_dist,1)) 
                      tile_found = .true.
                   endif
                endif
                if(present (day_in)) then
                   if(day_in(sub_tid(minloc(rev_dist,1))) < 998) then
                      getNeighbor = sub_tid(minloc(rev_dist,1)) 
                      tile_found = .true.
                   endif
                endif
             endif FOUND
             
          endif SEEK
          
          deallocate (sub_tid, sub_lon, sub_lat, rev_dist)
          
          if(tile_found) GO TO 100
          
          ! if not increase the window size
          dw = dw + 0.5
          
       end do ZOOMOUT
       
100    continue
          
       deallocate (mask)
   
     end function getNeighbor

     ! *****************************************************************************

     ! IRRGRR - duplicates same in Utils/mk_restarts/getids.F90
     
      function to_radian(degree) result(rad)
   
        real,intent(in) :: degree
        real :: rad
   
        rad = degree*MAPL_PI/180.
   
      end function to_radian
   
      ! *****************************************************************************

      ! IRRGRR - duplicates same in Utils/mk_restarts/getids.F90
      
      real function haversine(deglat1,deglon1,deglat2,deglon2)
        ! great circle distance -- adapted from Matlab 
        real,intent(in) :: deglat1,deglon1,deglat2,deglon2
        real :: a,c, dlat,dlon,lat1,lat2
        real,parameter :: radius = MAPL_radius
        
   !     dlat = to_radian(deglat2-deglat1)
   !     dlon = to_radian(deglon2-deglon1)
        !     lat1 = to_radian(deglat1)
   !     lat2 = to_radian(deglat2)
        dlat = deglat2-deglat1
        dlon = deglon2-deglon1
        lat1 = deglat1
        lat2 = deglat2     
        a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
        if(a>=0. .and. a<=1.) then
           c = 2*atan2(sqrt(a),sqrt(1-a))
           haversine = radius*c / 1000.
        else
           haversine = 1.e20
        endif
      end function
      
      ! *****************************************************************************
    end subroutine create_irrig_params
    
  END module module_irrig_params

! -----------------------------------------------------------------------------

!PROGRAM irrig_model
!
!  use module_irrig_params
!  
!  call create_irrig_params (43200, 21600, 'rst/SMAP_EASEv2_M36_964x406')
!  
!END PROGRAM irrig_model
