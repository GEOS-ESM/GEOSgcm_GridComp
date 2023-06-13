#include "MAPL_Generic.h"
module CATCHCN_INTERNAL_TYPE_Mod
  use ESMF
  use MAPL

  private
  public :: CATCHCN_WRAP
  public :: T_CATCHCN_STATE
  public :: surface_params_to_internal

  type T_CATCHCN_STATE
    ! CATCHCN_OFFLINE:
    !    0: DEFAULT for GCM,      (WW,CH,CM,CQ,FR) are required in CATCHCNment restart file 
    !    1: DEFAULT for GEOSldas, (WW,CH,CM,CQ,FR) are optional
    !    2: Option  for GEOSldas, (WW,CH,CM,CQ,FR) are optional for input restart but will be in output
    !                               restart; select when using GEOSldas to create restarts for the GCM.
    ! see also GEOSldas repo: src/Applications/LDAS_App/GEOSldas_LDAS.rc
    integer :: CATCHCN_OFFLINE
    integer :: RUN_IRRIG, USE_ASCATZ0, IRRIG_METHOD, AEROSOL_DEPOSITION, N_CONST_LAND4SNWALB
    integer :: ATM_CO2, PRESCRIBE_DVG,CHOOSEMOSFC
    real    :: SURFLAY              ! Default (Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params
                                !         (Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params
    real    :: CO2
    integer :: CO2_YEAR_IN          ! years when atmospheric carbon dioxide concentration increases, starting from 1850
    real    :: DTCN                 ! Time step for carbon/nitrogen routines in CatchmentCN model (default 5400)
    real    :: FWETC, FWETL
    logical :: USE_FWET_FOR_RUNOFF 
  end type T_CATCHCN_STATE

  type CATCHCN_WRAP
    type (T_CATCHCN_STATE), pointer :: PTR
  end type CATCHCN_WRAP

contains 

  subroutine surface_params_to_internal(wrap, scf, rc)
    type(CATCHCN_WRAP), intent(inout) :: wrap
    type(ESMF_Config), intent(inout) :: SCF    
    integer, optional, intent(out) :: rc
    type(T_CATCHCN_STATE), pointer :: statePtr
    real :: FWETC_default, FWETL_default
    integer:: status, SNOW_ALBEDO_INFO

    statePtr => wrap%ptr
    call MAPL_GetResource (SCF, statePtr%SURFLAY,             label='SURFLAY:',             DEFAULT=50.,     __RC__ )
    call MAPL_GetResource (SCF, statePtr%USE_ASCATZ0,         label='USE_ASCATZ0:',         DEFAULT=0,       __RC__ )
    call MAPL_GetResource (SCF, statePtr%RUN_IRRIG,           label='RUN_IRRIG:',           DEFAULT=0,       __RC__ )
    call MAPL_GetResource (SCF, statePtr%IRRIG_METHOD,        label='IRRIG_METHOD:',        DEFAULT=0,       __RC__ )
    call MAPL_GetResource (SCF, statePtr%CHOOSEMOSFC,         label='CHOOSEMOSFC:',         DEFAULT=1,       __RC__ )
    call MAPL_GetResource (SCF, statePtr%USE_FWET_FOR_RUNOFF, label='USE_FWET_FOR_RUNOFF:', DEFAULT=.FALSE., __RC__ )

    if (.NOT. statePtr%USE_FWET_FOR_RUNOFF) then
       FWETC_default = 0.02
       FWETL_default = 0.02
    else
       FWETC_default = 0.005
       FWETL_default = 0.025
    endif
    call MAPL_GetResource (SCF, statePtr%FWETC, label='FWETC:', DEFAULT=FWETC_default, __RC__ )
    call MAPL_GetResource (SCF, statePtr%FWETL, label='FWETL:', DEFAULT=FWETL_default, __RC__ )

    ! SNOW ALBEDO -- so far, only parameterization based on look-up table is implemented for CatchCN
    ! 0 : parameterization based on look-up table 
    ! 1 : MODIS-derived snow albedo (backfilled with global land average snow albedo)
    call MAPL_GetResource (SCF, SNOW_ALBEDO_INFO, label='SNOW_ALBEDO_INFO:', DEFAULT=0, __RC__)
    _ASSERT( SNOW_ALBEDO_INFO==0, "SNOW_ALBEDO_INFO must be 0 for CatchCN")


    ! GOSWIM ANOW_ALBEDO 
    ! 0 : GOSWIM snow albedo scheme is turned off
    ! 9 : i.e. N_CONSTIT in Stieglitz to turn on GOSWIM snow albedo scheme 
    call MAPL_GetResource (SCF, statePtr%N_CONST_LAND4SNWALB, label='N_CONST_LAND4SNWALB:', DEFAULT=0, __RC__ )

    ! Get parameters to zero the deposition rate 
    ! 1: Use all GOCART aerosol values, 0: turn OFF everythying, 
    ! 2: turn off dust ONLY,3: turn off Black Carbon ONLY,4: turn off Organic Carbon ONLY
    ! __________________________________________
    call MAPL_GetResource (SCF, statePtr%AEROSOL_DEPOSITION, label='AEROSOL_DEPOSITION:', DEFAULT=0, __RC__ )

    ! CATCHCN
    call MAPL_GetResource (SCF, statePtr%DTCN, label='DTCN:', DEFAULT=5400. , __RC__ )
    ! ATM_CO2
    ! 0: uses a fix value defined by CO2
    ! 1: CT tracker monthly mean diurnal cycle
    ! 2: CT tracker monthly mean diurnal cycle scaled to match EEA global average CO2
    ! 3: spatially fixed interannually varyiing CMIP from getco2.F90 look up table (AGCM only)
    ! 4: import AGCM model CO2 (AGCM only)
    call MAPL_GetResource (SCF, statePtr%ATM_CO2, label='ATM_CO2:', DEFAULT=2  , __RC__ )

    ! PRESCRIBE_DVG: Prescribe daily LAI and SAI data from an archived CATCHCN simulation 
    ! 0--NO Run CN Model interactively
    ! 1--YES Prescribe interannually varying LAI and SAI
    ! 2--YES Prescribe climatological LAI and SAI
    ! 3--Estimated LAI/SAI using anomalies at the beginning of the foeecast and climatological LAI/SAI
    call MAPL_GetResource (SCF, statePtr%PRESCRIBE_DVG, label='PRESCRIBE_DVG:', DEFAULT=0  , __RC__ )

    ! Global mean CO2 
    call MAPL_GetResource (SCF, statePtr%CO2,         label='CO2:',      DEFAULT=350.e-6, __RC__ )
    call MAPL_GetResource (SCF, statePtr%CO2_YEAR_IN, label='CO2_YEAR:', DEFAULT=  -9999, __RC__ )

    _RETURN(_SUCCESS)
  end subroutine

end module
