#include "MAPL_Generic.h"
module catch_wrap_stateMod
  use ESMF
  use MAPL
  
  private
  public :: T_CATCH_STATE
  public :: CATCH_WRAP
  public :: T_CATCHCN_STATE
  public :: CATCHCN_WRAP
  public :: surface_params_to_wrap_state
  
  type T_CATCH_STATE
     type (ESMF_FieldBundle)  :: Bundle
     !
     logical :: LDAS_CORRECTOR
     !
     ! CATCH_OFFLINE:
     !    0: DEFAULT for GCM,      (WW,CH,CM,CQ,FR) are required in Catchment restart file 
     !    1: DEFAULT for GEOSldas, (WW,CH,CM,CQ,FR) are optional
     !    2: Option  for GEOSldas, (WW,CH,CM,CQ,FR) are optional for input restart but will be in output
     !                               restart; select when using GEOSldas to create restarts for the GCM.
     ! see also GEOSldas repo: src/Applications/LDAS_App/GEOSldas_LDAS.rc
     integer :: CATCH_OFFLINE
     !
     integer :: CATCH_SPINUP
     !
     ! some (but not all) resource parameters from GEOS_SurfaceGridComp.rc:
     integer :: USE_ASCATZ0, Z0_FORMULATION, AEROSOL_DEPOSITION, N_CONST_LAND4SNWALB
     integer :: CHOOSEMOSFC, MOSFC_EXTRA_DERIVS_OFFL_LAND, SNOW_ALBEDO_INFO
     real    :: SURFLAY  
     real    :: FWETC, FWETL
     logical :: USE_FWET_FOR_RUNOFF
     integer :: RUN_IRRIG, IRRIG_METHOD
  end type T_CATCH_STATE
  
  type CATCH_WRAP
     type (T_CATCH_STATE), pointer :: PTR
  end type CATCH_WRAP
  
  type, extends(T_CATCH_STATE) :: T_CATCHCN_STATE
     ! resource parameters from GEOS_SurfaceGridComp.rc:     
     integer :: ATM_CO2, PRESCRIBE_DVG
     real    :: CO2
     integer :: CO2_YEAR_IN
     real    :: DTCN            
  end type T_CATCHCN_STATE
  
  type CATCHCN_WRAP
    type (T_CATCHCN_STATE), pointer :: PTR
  end type CATCHCN_WRAP

contains 

  subroutine surface_params_to_wrap_state(statePtr, scf, rc)

    ! obtain resource variables from "SURFRC" file; they are ultimately stored in CATCH_INTERNAL_STATE or CATCHCN_INTERNAL_STATE
    
    class(T_CATCH_STATE), pointer,  intent(inout) :: statePtr
    type(ESMF_Config),              intent(inout) :: SCF    
    integer,              optional, intent(  out) :: rc

    real :: FWETC_default, FWETL_default
    integer:: status, ii
    
    ! ************************************************* !
    !  For documentation, see GEOS_SurfaceGridComp.rc.  !
    ! ************************************************* !
    
    call MAPL_GetResource(    SCF, statePtr%SURFLAY,                  label='SURFLAY:',                  DEFAULT=50.,           __RC__ )
    call MAPL_GetResource(    SCF, statePtr%USE_ASCATZ0,              label='USE_ASCATZ0:',              DEFAULT=0,             __RC__ )
    call MAPL_GetResource(    SCF, statePtr%CHOOSEMOSFC,              label='CHOOSEMOSFC:',              DEFAULT=1,             __RC__ )

    ! MOSFC_EXTRA_DERIVS_OFFL_LAND: Resource parameter for *offline* (LDAS) mode.
    !
    !  Over *land*, use derivatives of exchange coeffs w.r.t. temp. & humidity.
    !    
    !    0 : None,               default for Helfand.  
    !    1 : Analytical derivs,  default for Louis,    *not* available for Helfand.
    !    2 : Numerical  derivs. 
    !    3 : Numerical  derivs,  via virtual temp.,    *not* available for Helfand,  same as 2 but faster than 2.
    !
    !  Runtimes: Helfand takes ~10 times longer than Louis.  In offline mode, Helfand consumes
    !            about as much CPU as Catchment.  Option 2 triples the runtime of the MOSFC scheme.
    !            Option 3 doubles the runtime of the Louis scheme.  

    if (statePtr%CATCH_OFFLINE==0) then

       statePtr%MOSFC_EXTRA_DERIVS_OFFL_LAND = 0        ! must be 0 for GCM

    else

       ! offline (LDAS) mode;  default for MOSFC_EXTRA_DERIVS_OFFL_LAND depends on CHOOSEMOSFC (Louis or Helfand)
       
       if     (statePtr%CHOOSEMOSFC==0) then
          
          ! Louis
          call MAPL_GetResource( SCF, statePtr%MOSFC_EXTRA_DERIVS_OFFL_LAND,  label='MOSFC_EXTRA_DERIVS_OFFL_LAND:',  DEFAULT=1,             __RC__ )
          ! make sure parameter value is allowed
          ii = statePtr%MOSFC_EXTRA_DERIVS_OFFL_LAND ; _ASSERT(ii>=0 .and. ii<=3, 'unknown MOSFC_EXTRA_DERIVS_OFFL_LAND for Louis  ')
          
       elseif (statePtr%CHOOSEMOSFC==1) then
          
          ! Helfand
          call MAPL_GetResource( SCF, statePtr%MOSFC_EXTRA_DERIVS_OFFL_LAND,  label='MOSFC_EXTRA_DERIVS_OFFL_LAND:',  DEFAULT=0,             __RC__ )
          ! make sure parameter value is allowed (analytical derivs not implemented for Helfand)
          ii = statePtr%MOSFC_EXTRA_DERIVS_OFFL_LAND ; _ASSERT(ii==0 .or. ii==2, 'unknown MOSFC_EXTRA_DERIVS_OFFL_LAND for Helfand')   
          
       else
          
          _ASSERT(.FALSE.,'unknown CHOOSEMOSFC')
          
       end if

    end if

    ! for CatchCN, must have MOSFC_EXTRA_DERIVS_OFFL_LAND<=1     (numerical derivatives not yet implemented for CatchCN)
    
    select type (statePtr)
    type is (T_CATCHCN_STATE) ! CATCHCN
       
       _ASSERT( statePtr%MOSFC_EXTRA_DERIVS_OFFL_LAND<=1, 'selected choice for MOSFC_EXTRA_DERIVS_OFFL_LAND not yet implemented for CatchCN')
       
    end select

    ! -------------------------
    
    call MAPL_GetResource(    SCF, statePtr%USE_FWET_FOR_RUNOFF,      label='USE_FWET_FOR_RUNOFF:',      DEFAULT=.FALSE.,       __RC__ )
    call MAPL_GetResource(    SCF, statePtr%Z0_FORMULATION,           label='Z0_FORMULATION:',           DEFAULT=4,             __RC__ )
    
    if (.NOT. statePtr%USE_FWET_FOR_RUNOFF) then
       FWETC_default = 0.02
       FWETL_default = 0.02
    else
       FWETC_default = 0.005   ! NOT ready for science!
       FWETL_default = 0.025   ! NOT ready for science!
    endif
    
    call MAPL_GetResource(    SCF, statePtr%FWETC,                    label='FWETC:',                    DEFAULT=FWETC_default, __RC__ )
    call MAPL_GetResource(    SCF, statePtr%FWETL,                    label='FWETL:',                    DEFAULT=FWETL_default, __RC__ )
    call MAPL_GetResource(    SCF, statePtr%SNOW_ALBEDO_INFO,         label='SNOW_ALBEDO_INFO:',         DEFAULT=0,             __RC__ )
    call MAPL_GetResource(    SCF, statePtr%N_CONST_LAND4SNWALB,      label='N_CONST_LAND4SNWALB:',      DEFAULT=0,             __RC__ )
    call MAPL_GetResource(    SCF, statePtr%AEROSOL_DEPOSITION,       label='AEROSOL_DEPOSITION:',       DEFAULT=0,             __RC__ )
    call MAPL_GetResource(    SCF, statePtr%RUN_IRRIG,                label='RUN_IRRIG:',                DEFAULT=0,             __RC__ )
    call MAPL_GetResource(    SCF, statePtr%IRRIG_METHOD,             label='IRRIG_METHOD:',             DEFAULT=0,             __RC__ )
    
    select type (statePtr)
    type is (T_CATCHCN_STATE) ! CATCHCN
       
       call MAPL_GetResource( SCF, statePtr%DTCN,                     label='DTCN:',                     DEFAULT=5400.,         __RC__ )
       call MAPL_GetResource( SCF, statePtr%ATM_CO2,                  label='ATM_CO2:',                  DEFAULT=2,             __RC__ )
       call MAPL_GetResource( SCF, statePtr%PRESCRIBE_DVG,            label='PRESCRIBE_DVG:',            DEFAULT=0,             __RC__ )
       call MAPL_GetResource( SCF, statePtr%CO2,                      label='CO2:',                      DEFAULT=350.e-6,       __RC__ )
       call MAPL_GetResource( SCF, statePtr%CO2_YEAR_IN,              label='CO2_YEAR:',                 DEFAULT=-9999,         __RC__ )
       
    end select

    _RETURN(_SUCCESS)

  end subroutine surface_params_to_wrap_state 
  
end module catch_wrap_stateMod

! ========================= EOF ======================================================
