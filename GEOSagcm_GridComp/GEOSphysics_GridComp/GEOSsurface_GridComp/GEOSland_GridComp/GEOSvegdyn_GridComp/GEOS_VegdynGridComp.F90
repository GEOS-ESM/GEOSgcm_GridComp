
!  $Id$

#include "MAPL_Generic.h"


!=============================================================================
module GEOS_VegdynGridCompMod

!BOP

! !MODULE: GEOS_Vegdyn -- child to the "Land" gridded component.  

!DESCRIPTION:
!   {\tt GEOS\_Vegdyn} is a gridded component that performs the
!   necessary interpolation to provide refreshed values of the 
!   dynamic vegetation values prescribed by external data/observations.\\
!
! Added MODIS_LAI as an IMPORT.
! Exports from this routine are the instaneous values of the
! vegetation parameters on tilespace to be used in other components
! of the land subroutine.  All exports and imports are stored on the
! tile grid inherited from the parent routine.\\
! 
! I. Parameter Class 1: Time AND spatially dependent parameters 
! from a binary data file\\
! 
! Current list: LAI, GRN, NDVI \\
! 
! The gridded component stores the surrounding observations of 
! each parameter in the internal state.  If the run method 
! discovers that the current internal state does not contain the 
! observed values required to interpolate the values at the current 
! time, it performs the required i/o to refresh the values of 
! the internal state.  The first iteration of the run method 
! always has to fill the values.  No restart is required by this 
! gridded component for these parameters.  (A restart *is* now
! required for Vegetation Class 3 \\
!
! INTERNALS: ITY, Z2CH, ASCATZ0\\
!
! EXPORTS:  LAI, GRN, NDVI\\
!
! !USES: 

  use ESMF
  use MAPL
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

  integer :: MODIS_DVG
  integer, parameter		     :: NTYPS = MAPL_NumVegTypes
  real,    dimension(   NTYPS)       :: VGRT
  ! real,    dimension(   NTYPS)       :: VGZ2   
  character(len=ESMF_MAXSTR)         :: GRIDNAME, MODIS_PATH
  
  data VGRT  / 19700., 7000., 9400., 7000., 7000., 14000./
  ! commented out legacy look-up table for veg heights, which are now always from bcs via restarts, - reichle, 17 March 2020
  ! data VGZ2 / 35.0, 20.0, 17.0, 0.6, 0.5, 0.6/ ! Dorman and Sellers (1989)   
  
contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS, OFFLINE_MODE
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    type(ESMF_Config)                       :: SCF
    character(len=ESMF_MAXSTR)              :: SURFRC

! Local derived type aliases

    type(MAPL_MetaComp),pointer             :: MAPL=>null()

!=============================================================================

! Begin...

!------------------------------------------------------------
! Get my name and set-up traceback handle
!------------------------------------------------------------

    call ESMF_GridCompGet(GC                             ,&
                          NAME=COMP_NAME                 ,&
                          RC=STATUS )

    VERIFY_(STATUS)

    Iam = trim(COMP_NAME) // 'SetServices'

! -----------------------------------------------------------
! Set the Run entry point
! -----------------------------------------------------------

    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get the configuration
! -----------------------------------------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get the intervals
! -----------------------------------------------------------

    !call MAPL_GetResource ( MAPL,DT, Label="RUN_DT:", RC=STATUS)
    !VERIFY_(STATUS)

    !RUN_DT = nint(DT)

! -----------------------------------------------------------
! Get experiment configuration parameters
! -----------------------------------------------------------
    
    call MAPL_GetResource (MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='MODIS_LAI:', value=MODIS_DVG, DEFAULT=0, __RC__ ) ; VERIFY_(STATUS)

    if (MODIS_DVG == 1) then
       call ESMF_ConfigGetAttribute (SCF, value= MODIS_PATH, label='MODIS_PATH:',  &
            DEFAULT='/discover/nobackup/projects/gmao/ssd/land/l_data/LandBCs_files_for_mkCatchParam/V001/MCD15A2H.006/IAV_smoothed/', __RC__ )
       call MAPL_GetResource (MAPL, OFFLINE_MODE, Label="CATCHMENT_OFFLINE:", DEFAULT=0, RC=STATUS)
       if (OFFLINE_MODE == 0) then
          call MAPL_GetResource( MAPL, GRIDNAME, label='AGCM_GRIDNAME:', RC=STATUS ) ; VERIFY_(STATUS)
       else
          call MAPL_GetResource( MAPL, GRIDNAME, label='GRIDNAME:',      RC=STATUS ) ; VERIFY_(STATUS)
       endif    
    endif
    call ESMF_ConfigDestroy      (SCF, __RC__)
    
! -----------------------------------------------------------
! At the moment, this will refresh when the land parent 
! needs to refresh.
!
!    call ESMF_ConfigGetFloat ( CF, DT, Label=trim(COMP_NAME)//&
!    "_DT:", default=DT, RC=STATUS)
!     VERIFY_(STATUS)
!
!    MY_STEP = nint(DT)
!
! -----------------------------------------------------------

! -----------------------------------------------------------
! Set the state variable specs.
! -----------------------------------------------------------

!BOS

! -----------------------------------------------------------
!   Import States
! -----------------------------------------------------------

    IF (MODIS_DVG == 2) THEN

       call MAPL_AddImportSpec(gc, &
            short_name = "MODIS_LAI", &
            LONG_NAME  = 'MODIS Leaf Area Index',                     &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsTileOnly,                           &
            VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
       VERIFY_(STATUS)
      
    ENDIF

! -----------------------------------------------------------
! Internal State 
! -----------------------------------------------------------

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'ITY'                                     ,&
         LONG_NAME  = 'vegetation_type'			        ,&
         UNITS      = '1'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'Z2CH'                                    ,&
         LONG_NAME  = 'vegetation_height'			,&
         UNITS      = 'm'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddInternalSpec(GC                                ,&
         SHORT_NAME = 'ASCATZ0'                                 ,&
         LONG_NAME  = 'ASCAT_roughness_length'	                ,&
         UNITS      = 'm'                                       ,&
         DIMS       = MAPL_DimsTileOnly                         ,&
         VLOCATION  = MAPL_VLocationNone                        ,&
         FRIENDLYTO = trim(COMP_NAME)                           ,&
         RC=STATUS  )
    VERIFY_(STATUS)  

! -----------------------------------------------------------
! These are variables that are considered time-independent
! and are stored and retrieved as-is
! -----------------------------------------------------------

! -----------------------------------------------------------
! Export Variables
! -----------------------------------------------------------

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'LAI'                                     ,&
       LONG_NAME  = 'leaf_area_index'                         ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'GRN'                                     ,&
       LONG_NAME  = 'greeness_fraction'			      ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  	 

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'ROOTL'                                   ,&
       LONG_NAME  = 'root_length_density'                     ,&
       UNITS      = 'm+2'                                     ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  	 

    call MAPL_AddExportSpec(GC                                ,&
       SHORT_NAME = 'NDVI'                                    ,&
       LONG_NAME  = 'normalized_difference_vegetation_index'  ,&
       UNITS      = '1'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )

    VERIFY_(STATUS)  

!EOS

!------------------------------------------------------------
! Set generic init and final methods
!------------------------------------------------------------

    call MAPL_GenericSetServices(GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

! -----------------------------------------------------------
! RUN -- Run method for the vegdyn component
! -----------------------------------------------------------

  subroutine RUN (GC,IMPORT, EXPORT, CLOCK, RC )

! -----------------------------------------------------------
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC    
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    type(ESMF_Clock),    intent(inout) :: CLOCK
    integer, optional,   intent(  out) :: RC

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

    type (MAPL_MetaComp),     pointer   :: MAPL=>null()
    type (ESMF_State       )            :: INTERNAL

! INTERNAL pointers
 
    real, dimension(:), pointer :: ITY
    real, dimension(:), pointer :: Z2CH
    real, dimension(:), pointer :: ASCATZ0

! IMPORT Pointers

    real, save, dimension(:), pointer :: MODIS_LAI

! EXPORT pointers 

    real, dimension(:), pointer :: LAI
    real, dimension(:), pointer :: GRN
    real, dimension(:), pointer :: ROOTL
    real, dimension(:), pointer :: NDVI
  
! Time attributes and placeholders

    type(ESMF_Time)         :: CURRENT_TIME, MODIS_TIME
    type(ESMF_Time), save   :: MODIS_RING
    type(ESMF_TimeInterval) :: M8, TIME_DIFF

! Others

    character(len=ESMF_MAXSTR)         :: LAIFile
    character(len=ESMF_MAXSTR)         :: GRNFile
    character(len=ESMF_MAXSTR)         :: NDVIFile
    character(len=ESMF_MAXSTR)         :: LAItpl
    character(len=ESMF_MAXSTR)         :: GRNtpl
    character(len=ESMF_MAXSTR)         :: NDVItpl
    integer                            :: NUM_LDAS_ENSEMBLE, ens_id_width
    integer                            :: MOD_DOY,DOY,MFDOY,CUR_YY,CUR_MM,CUR_DD, &
                                          MOD_YY, MOD_MM, MOD_DD
    logical, save                      :: first = .true.
    logical                            :: b4_modis_date = .false.
    integer                            :: MODIS_FIRSTDATE = 20020704
    
! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet(GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
  
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
! -----------------------------------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS)
    VERIFY_(STATUS) 

    call MAPL_GetResource ( MAPL, NUM_LDAS_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, ens_id_width, Label="ENS_ID_WIDTH:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

! -----------------------------------------------------------
! Get file names from configuration
! -----------------------------------------------------------

    if(NUM_LDAS_ENSEMBLE > 1) then
       !comp_name should be vegdynxxxx....
       call MAPL_GetResource(MAPL, LAItpl, label = 'LAI'//comp_name(7:7+ens_id_width-1)//'_FILE:', &
            RC=STATUS )
       call MAPL_GetResource(MAPL, GRNtpl, label = 'GREEN'//comp_name(7:7+ens_id_width-1)//'_FILE:', &
            RC=STATUS )
       call MAPL_GetResource(MAPL, NDVItpl, label = 'NDVI'//comp_name(7:7+ens_id_width-1)//'_FILE:', &
            RC=STATUS )

       if (STATUS/=ESMF_SUCCESS) then
          call MAPL_GetResource(MAPL, LAItpl, label = 'LAI_FILE:', &
               default = '../input/lai%s.data', RC=STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, GRNtpl, label = 'GREEN_FILE:', &
               default = '../input/green%s.data', RC=STATUS )
          VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, NDVItpl, label = 'NDVI_FILE:', &
                default = '../input/ndvi%s.data', RC=STATUS )
       endif

       call ESMF_CFIOStrTemplate(LAIFILE, LAItpl,'GRADS', xid=comp_name(7:7+ens_id_width-1), stat=status)
       VERIFY_(STATUS)
       call ESMF_CFIOStrTemplate(GRNFILE, GRNtpl,'GRADS', xid=comp_name(7:7+ens_id_width-1), stat=status)
       VERIFY_(STATUS)
       call ESMF_CFIOStrTemplate(NDVIFILE,NDVItpl,'GRADS',xid=comp_name(7:7+ens_id_width-1), stat=status)
       VERIFY_(STATUS)
    else
       call MAPL_GetResource(MAPL, LAIFILE, label = 'LAI_FILE:', &
         default = 'lai.dat', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, GRNFILE, label = 'GREEN_FILE:', &
         default = 'green.dat', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, NDVIFILE, label = 'NDVI_FILE:', &
          default = 'ndvi.dat', RC=STATUS )
       VERIFY_(STATUS)
   endif

! get pointers to internal variables
! ----------------------------------
  
    call MAPL_GetPointer(INTERNAL,      ITY,      'ITY' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL,      Z2CH,     'Z2CH', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL,   ASCATZ0,  'ASCATZ0', RC=STATUS)
    VERIFY_(STATUS)

! get pointers to EXPORTS
! -----------------------

    call MAPL_GetPointer(EXPORT, LAI,   'LAI',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, GRN,   'GRN',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ROOTL, 'ROOTL',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, NDVI,   'NDVI',  RC=STATUS)
    VERIFY_(STATUS)

! Do the lai greeness and ndvi interpolation
! ------------------------------------------

    call ESMF_ClockGet  ( CLOCK, currTime=CURRENT_TIME, RC=STATUS )
    VERIFY_(STATUS)

    IF (MODIS_DVG == 0) THEN
       
       ! read lai_clim_IMxJM.data file from BCSDIR
       call MAPL_ReadForcing(MAPL,'LAI',LAIFILE,CURRENT_TIME,LAI,ON_TILES=.true.,RC=STATUS)
       VERIFY_(STATUS)
       
    ELSE IF (MODIS_DVG == 1) THEN

       ! read interannually varying LAI data on tile space
       call ESMF_TimeIntervalSet(M8, h=24*8, rc=status ) ; VERIFY_(STATUS)
       MOD_YY = MODIS_FIRSTDATE / 10000
       MOD_MM = (MODIS_FIRSTDATE - MOD_YY*10000) / 100
       MOD_DD = MODIS_FIRSTDATE - (MOD_YY*10000 + MOD_MM*100)
       call ESMF_TimeSet (MODIS_TIME, yy=MOD_YY, mm=MOD_MM, dd=MOD_DD, rc=status) ; VERIFY_(STATUS)
       call ESMF_TimeGet (MODIS_TIME, DayOfYear=MFDOY, RC=STATUS)                 ; VERIFY_(STATUS)

       if (first) then          
          call ESMF_TimeGet (CURRENT_TIME, YY = CUR_YY, MM = CUR_MM, DD = CUR_DD, DayOfYear=DOY, RC=STATUS); VERIFY_(STATUS)
          MOD_DOY = modis_date (DOY)
          call ESMF_TimeSet(MODIS_TIME, yy=CUR_YY, mm=CUR_MM, dd=CUR_DD, rc=status) ; VERIFY_(STATUS)
          call ESMF_TimeIntervalSet(TIME_DIFF, h=24*(DOY -MOD_DOY), rc=status )     ; VERIFY_(STATUS)
          MODIS_RING = MODIS_TIME - TIME_DIFF
          MODIS_RING = MODIS_RING + M8
          
          ALLOCATE (MODIS_LAI (1:SIZE(ITY)))
          if((MOD_YY*1000 + MFDOY) > (CUR_YY*1000 + MOD_DOY)) b4_modis_date = .true.
          call read_modis_lai (MAPL,CUR_YY, MOD_DOY, MODIS_LAI, b4_modis_date)
          first = .false.
       endif

       if (CURRENT_TIME ==  MODIS_RING) then
          call ESMF_TimeGet (CURRENT_TIME, YY = CUR_YY, DayOfYear=DOY, RC=STATUS); VERIFY_(STATUS)
          MOD_DOY = modis_date (DOY)
          if((MOD_YY*1000 + MFDOY) > (CUR_YY*1000 + MOD_DOY)) b4_modis_date = .true.
          call read_modis_lai (MAPL,CUR_YY, DOY, MODIS_LAI, b4_modis_date)
          if (DOY < 361) then
             MODIS_RING = CURRENT_TIME + M8
          else
             call ESMF_TimeSet(MODIS_TIME, yy=CUR_YY+1, mm=1, dd=1, rc=status) ; VERIFY_(STATUS)
             MODIS_RING = MODIS_TIME
          endif
       endif
       
       LAI = MODIS_LAI
       LAI = min(7., max(0.0001, LAI))
       
       
    ELSE IF(MODIS_DVG == 2) THEN

       ! import via ExtData
       call MAPL_GetPointer(IMPORT,MODIS_LAI, 'MODIS_LAI', RC=STATUS) ; VERIFY_(STATUS)       
       LAI = MODIS_LAI
       LAI = min(7., max(0.0001, LAI))
              
    ENDIF

    call MAPL_ReadForcing(MAPL,'GRN',GRNFILE,CURRENT_TIME,GRN,ON_TILES=.true.,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_ReadForcing(MAPL,'NDVI',NDVIFILE,CURRENT_TIME,NDVI,ON_TILES=.true.,RC=STATUS)
    VERIFY_(STATUS)

! Vegetation types used to index into tables
! Root length density no longer depends on time of year
! -----------------------------------------------------

    ROOTL = VGRT(nint(ITY))
    
!  All done
! ---------

    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  contains

    ! ---------------------------------------------------------------------------
    
    integer function modis_date (DOY) result (MOD_DOY)
      
      implicit none
      integer, intent(in) :: DOY
      integer, parameter  :: N_MODIS_DATES = 46
      integer, dimension (N_MODIS_DATES), parameter ::  &
           MODIS_DOYS = (/                              &
           1  ,  9, 17, 25, 33, 41, 49, 57, 65,         &
           73 , 81, 89, 97,105,113,121,129,137,         &
           145,153,161,169,177,185,193,201,209,         &
           217,225,233,241,249,257,265,273,281,         &
           289,297,305,313,321,329,337,345,353,361/)
      integer :: i
      
      if (DOY < MODIS_DOYS(N_MODIS_DATES)) then
         do i = 1, N_MODIS_DATES
            if (MODIS_DOYS(i) > DOY) exit
         end do         
         MOD_DOY = MODIS_DOYS(i-1)
      else
         MOD_DOY = MODIS_DOYS(N_MODIS_DATES)
      endif
            
    end function modis_date
    
    ! ---------------------------------------------------------------------------

    subroutine read_modis_lai (MAPL,CUR_YY, MOD_DOY,MODIS_LAI, b4_modis_date) 

      implicit none
      integer, intent (in)                     :: CUR_YY, MOD_DOY
      logical, intent (in)                     :: b4_modis_date
      real, dimension (:), intent (inout)      :: MODIS_LAI
      type(MAPL_MetaComp),pointer, intent (in) :: MAPL
      type(ESMF_Grid)                          :: TILEGRID
      type(MAPL_LocStream)                     :: LOCSTREAM
      integer, pointer                         :: mask(:)
      integer                                  :: status, unit
      character*300                            :: filename
      CHARACTER(len=7)                         :: YYYYDoY
      logical                                  :: file_exists
      
      if(b4_modis_date) then
         WRITE (YYYYDoY,'(a4,i3.3)') 'YYYY',MOD_DOY
      else
         WRITE (YYYYDoY,'(i4.4,i3.3)') CUR_YY,MOD_DOY
      endif
      
      filename = trim(MODIS_PATH)//'/'//trim(GRIDNAME)//'/lai_data.'//YYYYDoY
      inquire(file=filename, exist=file_exists)
      if (.not. file_exists) then
          _ASSERT(.FALSE.,'Missing : '//trim(filename))
      endif
      call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)            ; VERIFY_(STATUS)
      call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS); VERIFY_(STATUS)
      call MAPL_TileMaskGet(tilegrid,  mask, rc=status)              ; VERIFY_(STATUS)
      unit = GETFILE(trim(filename), form="unformatted", RC=STATUS)  ; VERIFY_(STATUS)
      call MAPL_VarRead(unit,tilegrid,MODIS_LAI,mask=mask,RC=STATUS) ; VERIFY_(STATUS)
      call FREE_FILE(unit, RC=STATUS)                                ; VERIFY_(STATUS)
      
    end subroutine read_modis_lai
    
    ! ---------------------------------------------------------------------------
  end subroutine RUN

end module GEOS_VegdynGridCompMod
