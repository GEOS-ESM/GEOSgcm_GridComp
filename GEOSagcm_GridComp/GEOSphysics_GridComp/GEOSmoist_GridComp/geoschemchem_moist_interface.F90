#include "MAPL_Generic.h"
!BOP
!
! !MODULE: geoschemchem_moist_interface 
!
! !DESCRIPTION: This module contains routines and variables to control 
!  the convective washout of GEOS-Chem species by GEOS moist. 
!\\
!\\
! !INTERFACE:
!
MODULE geoschemchem_moist_interface
!
! !USES:
!
USE ESMF
USE MAPL

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: compute_ki_gcc_aerosol
  PUBLIC  :: compute_ki_gcc_gas
  PUBLIC  :: GCCparams
  PUBLIC  :: GCCparams_init
  PUBLIC  :: GCCparams_set
  PUBLIC  :: GCCparams_cleanup
  PUBLIC  :: is_gcc_species
  PUBLIC  :: get_gcc_diagID
  PUBLIC  :: GCCspecies
  PUBLIC  :: GCCmax   
  PUBLIC  :: GCCfsol  
  PUBLIC  :: GCCscav  
  PUBLIC  :: GCCdiag_init 
  PUBLIC  :: GCCdiag_cleanup
  PUBLIC  :: GCCdiag_AddExports
  PUBLIC  :: GCCdiag_FillExports
  PUBLIC  :: GCCdiag_Count
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GCC_E_ICE 
!
! !REMARKS:
!  References:
!  D.J. Jacob, H. Liu, C. Mari, and R. M. Yantosca, "Harvard wet
!  deposition scheme for GMI", Harvard Atmospheric Chemistry Modeling 
!  Group, March 2000.
!
! !REVISION HISTORY:
!  20210312 - christoph.a.keller@nasa.gov - initial version
!  See https://github.com/GEOS-ESM/GEOSgcm_GridComp for full history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  REAL, PARAMETER       :: KC_DEFAULT = 5.e-3  ! s-1  ! Default autoconversion parameter for GEOS-Chem species
!
! MODULE VARIABLES:
!
 ! Type to keep all species washout parameters, passed down from GEOS-Chem.
 TYPE GCCparams_vars
   ! for aerosols
   REAL    :: KcScal1,KcScal2,KcScal3
   ! for gases
   REAL    :: convfaci2g
   REAL    :: retfactor
   REAL    :: liq_and_gas
   REAL    :: online_cldliq 
 END TYPE GCCparams_vars
 TYPE (GCCparams_vars), ALLOCATABLE :: GCCparams(:)

 ! GCCfsol and GCCscav hold the GEOS-Chem moist diagnostics. nGCCfsol and nGCCscav store the number of 
 ! requested diagnostics for a given run. If zero, no of the diagnostics calculations will be performed.
 REAL, ALLOCATABLE              :: GCCfsol(:,:,:,:)
 REAL, ALLOCATABLE              :: GCCscav(:,:,:)
 INTEGER, SAVE                  :: nGCCfsol = -1
 INTEGER, SAVE                  :: nGCCscav = -1


! Variables to compute GEOS-Chem moist exports 
 INTEGER, PARAMETER             :: GCCmax = 28
 CHARACTER(LEN=5), PARAMETER    :: GCCspecies(GCCmax) = &
                                            (/ 'SO2',   &
                                               'SO4',   &
                                               'HNO3',  &
                                               'DST1',  &
                                               'DST2',  &
                                               'DST3',  &
                                               'DST4',  &
                                               'NIT',   &
                                               'SALA',  &
                                               'SALC',  &
                                               'BCPI',  &
                                               'BCPO',  &
                                               'OCPI',  &
                                               'OCPO',  &
                                               'ALD2',  &
                                               'ALK4',  &
                                               'CH2O',  &
                                               'Br2',   &
                                               'HBr',   &
                                               'MAP',   &
                                               'MEK',   &
                                               'MTPA',  &
                                               'MTPO',  &
                                               'MVK',   &
                                               'NH3',   &
                                               'NH4',   &
                                               'PAN',   &
                                               'CFC12'  /)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1          !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GCCdiag_AddExports 
!
! !DESCRIPTION: Create exports related to convective transport/washout of
!  GEOS-Chem species. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GCCdiag_AddExports( GC, RC )
!
! !PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)  :: GC       ! Ref to GridComp
    INTEGER,             INTENT(OUT)    :: RC       ! Success or failure
!
! !REVISION HISTORY:
!  20210312 - christoph.a.keller@nasa.gov - initial version
!  See https://github.com/GEOS-ESM/GEOSgcm_GridComp for full history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER                         :: I
    CHARACTER(LEN=ESMF_MAXSTR)      :: spcname
 
    __Iam__('GCCdiag_AddExports')

    !=======================================================================
    ! Begins here 
    !=======================================================================

    DO I = 1,GCCmax
       spcname = TRIM(GCCspecies(I))

       call MAPL_AddExportSpec(GC,                                             &
         SHORT_NAME='GF_ConvScav_'//TRIM(spcname),                             &
         LONG_NAME ='GEOS-Chem_'//TRIM(spcname)//'_dry_air_tendency_due_to_GF_conv_scav', &
         UNITS     ='kg m-2 s-1',                                              &
         DIMS      = MAPL_DimsHorzOnly,                                        &
          __RC__ )

       call MAPL_AddExportSpec(GC,                                             &
         SHORT_NAME='GF_WetLossConvFrac_'//TRIM(spcname),                      &
         LONG_NAME ='GEOS-Chem_'//TRIM(spcname)//'_fraction_lost_in_GF_convection', &
         UNITS     ='1',                                                       &
         DIMS      = MAPL_DimsHorzVert,                                        &
         VLOCATION = MAPL_VLocationCenter,                                     &
          __RC__ )
    ENDDO

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCCdiag_AddExports 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1          !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GCCdiag_FillExports 
!
! !DESCRIPTION: Fill exports related to convective transport/washout of
!  GEOS-Chem species. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GCCdiag_FillExports( EXPORT, RC )
!
! !PARAMETERS:
!
    TYPE(ESMF_State),    INTENT(INOUT)  :: EXPORT    ! Export state
    INTEGER,             INTENT(OUT)    :: RC        ! Success or failure
!
! !REVISION HISTORY:
!  20210312 - christoph.a.keller@nasa.gov - initial version
!  See https://github.com/GEOS-ESM/GEOSgcm_GridComp for full history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER                            :: I
    CHARACTER(LEN=ESMF_MAXSTR)         :: spcname
    REAL, POINTER, DIMENSION(:,:)      :: GCptr2d
    REAL, POINTER, DIMENSION(:,:,:)    :: GCptr3d
 
    __Iam__('GCCdiag_FillExports')

    !=======================================================================
    ! Begins here 
    !=======================================================================

    ! Loop over all species
    do I=1,GCCmax
       if ( allocated(GCCscav) ) then
          spcname = 'GF_ConvScav_'//TRIM(GCCspecies(I))
          call MAPL_GetPointer(EXPORT, GCptr2d, TRIM(spcname), NotFoundOk=.TRUE., __RC__ )
          if ( associated(GCptr2d) ) then
             IF ( ALLOCATED(GCCscav) ) GCptr2d  = GCCscav(:,:,I)
          endif
       endif
       if ( allocated(GCCfsol) ) then
          spcname = 'GF_WetLossConvFrac_'//TRIM(GCCspecies(I))
          call MAPL_GetPointer(EXPORT, GCptr3d, TRIM(spcname), NotFoundOk=.TRUE., __RC__ )
          if ( associated(GCptr3d) ) then
             IF ( ALLOCATED(GCCfsol) ) GCptr3d  = GCCfsol(:,:,:,I)
          endif
       endif
    enddo

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCCdiag_FillExports 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1          !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GCCdiag_Count
!
! !DESCRIPTION: Count GEOS-Chem moist diagnostics to be filled. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GCCdiag_Count( EXPORT, RC )
!
! !PARAMETERS:
!
    TYPE(ESMF_State),    INTENT(INOUT)  :: EXPORT    ! Export state
    INTEGER,             INTENT(OUT)    :: RC        ! Success or failure
!
! !REVISION HISTORY:
!  20210317 - christoph.a.keller@nasa.gov - initial version
!  See https://github.com/GEOS-ESM/GEOSgcm_GridComp for full history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER                            :: I
    CHARACTER(LEN=ESMF_MAXSTR)         :: spcname
    REAL, POINTER, DIMENSION(:,:)      :: GCptr2d
    REAL, POINTER, DIMENSION(:,:,:)    :: GCptr3d
 
    __Iam__('GCCdiag_FillExports')

    !=======================================================================
    ! Begins here 
    !=======================================================================

    ! Reset counters
    nGCCscav = 0
    nGCCfsol = 0

    ! Loop over all species
    do I=1,GCCmax
       spcname = 'GF_ConvScav_'//TRIM(GCCspecies(I))
       call MAPL_GetPointer(EXPORT, GCptr2d, TRIM(spcname), NotFoundOk=.TRUE., __RC__ )
       if ( associated(GCptr2d) ) then
          nGCCscav = nGCCscav + 1
       endif
       spcname = 'GF_WetLossConvFrac_'//TRIM(GCCspecies(I))
       call MAPL_GetPointer(EXPORT, GCptr3d, TRIM(spcname), NotFoundOk=.TRUE., __RC__ )
       if ( associated(GCptr3d) ) then
          nGCCfsol = nGCCfsol + 1
       endif
    enddo

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCCdiag_Count
!EOC


!------------------------------------------------------------------------------

  SUBROUTINE GCCdiag_init( IM, JM, LM, RC )
    INTEGER,      INTENT(IN )      :: IM, JM, LM 
    INTEGER,      INTENT(OUT)      :: RC   ! Return code
    __Iam__('GCCdiag_init')
    IF ( nGCCscav > 0 ) THEN
       IF ( .NOT. ALLOCATED(GCCscav) ) THEN
          ALLOCATE(GCCscav(IM,JM,GCCmax), STAT=STATUS )
          ASSERT_(STATUS==0)
       ENDIF
       GCCscav(:,:,:) = 0.0
    ENDIF
    IF ( nGCCfsol > 0 ) THEN
       IF ( .NOT. ALLOCATED(GCCfsol) ) THEN
          ALLOCATE(GCCfsol(IM,JM,LM,GCCmax), STAT=STATUS )
          ASSERT_(STATUS==0)
       ENDIF
       GCCfsol(:,:,:,:) = 0.0
    ENDIF

    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GCCdiag_init
!-----------------------------------------------------------------------------------------


  SUBROUTINE GCCdiag_cleanup ( RC )
    INTEGER,      INTENT(OUT)      :: RC   ! Return code
    __Iam__('GCCdiag_cleanup')
    IF ( ALLOCATED(GCCfsol) ) DEALLOCATE(GCCfsol)
    IF ( ALLOCATED(GCCscav) ) DEALLOCATE(GCCscav)
    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE GCCdiag_cleanup
!-----------------------------------------------------------------------------------------


  SUBROUTINE GCCparams_init( KM, RC )
    INTEGER,      INTENT(IN )      :: KM   ! # of species
    INTEGER,      INTENT(OUT)      :: RC   ! Return code
    __Iam__('GCCparams_init')
    ALLOCATE(GCCparams(KM), STAT=STATUS ) 
    ASSERT_(STATUS==0)
    GCCparams(1:KM)%KcScal1        = 1.
    GCCparams(1:KM)%KcScal2        = 1.
    GCCparams(1:KM)%KcScal3        = 1.
    GCCparams(1:KM)%retfactor      = 1.
    GCCparams(1:KM)%online_cldliq  = 1.
    GCCparams(1:KM)%liq_and_gas    = 0.
    GCCparams(1:KM)%convfaci2g     = 0.
    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE GCCparams_init
!-----------------------------------------------------------------------------------------

  SUBROUTINE GCCparams_set( Field, k, RC )

    TYPE(ESMF_Field), INTENT(INOUT) :: Field   ! species field  
    INTEGER,          INTENT(IN )   :: k       ! species index 
    INTEGER,          INTENT(OUT)   :: RC      ! Return code

    LOGICAL            :: isPresent
    REAL, DIMENSION(3) :: Vect_KcScal
    REAL               :: ival

    __Iam__('GCCparams_set')

    ! KC scale factors for GEOS-Chem
    Vect_KcScal(:) = 1.0
    call ESMF_AttributeGet  (FIELD,"SetofKcScalFactors",isPresent=isPresent, __RC__ )
    if (isPresent) then
       call ESMF_AttributeGet  (FIELD,"SetofKcScalFactors",Vect_KcScal, __RC__ )
       GCCparams(k)%KcScal1 = Vect_KcScal(1)
       GCCparams(k)%KcScal2 = Vect_KcScal(2)
       GCCparams(k)%KcScal3 = Vect_KcScal(3)
    ENDIF
    ! Gas-phase washout parameter for GEOS-Chem
    call ESMF_AttributeGet (FIELD,"RetentionFactor",isPresent=isPresent, __RC__ )
    if (isPresent) then
       call ESMF_AttributeGet (FIELD,"RetentionFactor",ival, __RC__ )
       GCCparams(k)%retfactor = ival
    endif
    call ESMF_AttributeGet (FIELD,"LiqAndGas",isPresent=isPresent, __RC__ )
    if (isPresent) then
       call ESMF_AttributeGet (FIELD,"LiqAndGas",ival, __RC__ )
       GCCparams(k)%liq_and_gas = ival
    endif
    call ESMF_AttributeGet (FIELD,"ConvFacI2G",isPresent=isPresent, __RC__ )
    if (isPresent) then
       call ESMF_AttributeGet (FIELD,"ConvFacI2G",ival, __RC__ )
       GCCparams(k)%convfaci2g = ival
    endif
    call ESMF_AttributeGet (FIELD,"OnlineCLDLIQ",isPresent=isPresent, __RC__ )
    if (isPresent) then
       call ESMF_AttributeGet (FIELD,"OnlineCLDLIQ",ival, __RC__ )
       GCCparams(k)%online_cldliq = ival
    endif

    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GCCparams_set
!-----------------------------------------------------------------------------------------

  SUBROUTINE GCCparams_cleanup( RC )
    INTEGER,      INTENT(OUT)      :: RC   ! Return code
    __Iam__('GCCparams_cleanup')
    IF ( ALLOCATED(GCCparams) ) DEALLOCATE(GCCparams) 
    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE GCCparams_cleanup
!-----------------------------------------------------------------------------------------

  FUNCTION is_gcc_species( specname ) RESULT( TrueOrFalse )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    ! Check if this is a GEOS-Chem species, based on its name. Assumes that all GEOS-Chem
    ! species start with 'SPC_'. 
    !EOP
    !=====================================================================================
    character(len=*), intent(in) :: specname   ! species name to check 
    logical                      :: TrueOrFalse
    TrueOrFalse = .FALSE.
    if ( LEN(TRIM(specname))>4 ) then
       if ( TRIM(specname(1:4)) == 'SPC_' ) then
          TrueOrFalse = .TRUE. 
       endif
    endif
  END FUNCTION is_gcc_species
!-----------------------------------------------------------------------------------------

  FUNCTION get_gcc_diagID( specname ) RESULT( diagID )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    ! Return the diagnostics ID for the given species. 
    !EOP
    !=====================================================================================
    character(len=*), intent(in) :: specname   ! species name to check 
    integer                      :: diagID 
    ! local
    integer    :: i
    ! start here
    diagID = -1
    do i = 1,GCCmax
       if( TRIM(specname) == 'SPC_'//TRIM(GCCspecies(i)) ) then
          diagID=i
          exit
       endif
    enddo
  END FUNCTION get_gcc_diagID 
!-----------------------------------------------------------------------------------------


  SUBROUTINE compute_ki_gcc_aerosol( temp, kcscal1, kcscal2, kcscal3, kc_scaled )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    ! Compute the loss rate of a GEOS-Chem aerosol. This follows the parameterization described 
    ! in Jacob et al. 2000, as implemented in GEOS-Chem.
    !EOP
    !=====================================================================================
    implicit none

    real,    intent(in)      :: temp               ! temperature [K]
    real,    intent(in)      :: kcscal1            ! temperature-dependent scale factor for temperature range 1 
    real,    intent(in)      :: kcscal2            ! temperature-dependent scale factor for temperature range 2 
    real,    intent(in)      :: kcscal3            ! temperature-dependent scale factor for temperature range 3 
    real,    intent(out)     :: kc_scaled          ! loss rate [s-1]
    ! parameter
    real, parameter       :: TEMP1  = 237.0        ! K
    real, parameter       :: TEMP2  = 258.0        ! K

    ! start with default kc, then scale based on temperature and aerosol-specific scale factor
    kc_scaled = KC_DEFAULT    
    if ( temp < TEMP1 ) then
       kc_scaled = kc_scaled * kcscal1
    else if ( (temp>=TEMP1) .and. (temp<TEMP2) ) then
       kc_scaled = kc_scaled * kcscal2
    else
       kc_scaled = kc_scaled * kcscal3
    endif

  END SUBROUTINE compute_ki_gcc_aerosol
!-----------------------------------------------------------------------------------------


   SUBROUTINE compute_ki_gcc_gas( temp, press, q, cldh2o, heff, liq_and_gas, convfaci2g, retfactor, online_cldliq, kc_scaled )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Compute the loss rate of a GEOS-Chem gas. This follows the parameterization described 
    !  in Jacob et al. 2000, as implemented in GEOS-Chem.
    !EOP
    !=====================================================================================
     implicit none

     real,    intent(in)      :: temp               ! temperature [K]
     real,    intent(in)      :: press              ! pressure [Pa] 
     real,    intent(in)      :: q                  ! water vapor mixing ratio [kg/kg] 
     real,    intent(in)      :: cldh2o             ! cloud total water [kg/kg] 
     real,    intent(in)      :: heff               ! effective gas/aq Henry constant [-]
     real,    intent(in)      :: liq_and_gas        ! species considers ice and liquid phase?
     real,    intent(in)      :: convfaci2g         ! conversion factor for ice/gas ratio
     real,    intent(in)      :: retfactor          ! retention factor [unitless] 
     real,    intent(in)      :: online_cldliq      ! calculate cloud liquid/ice online or use default GEOS-Chem parameterization
     real,    intent(out)     :: kc_scaled          ! loss rate [s-1]

     ! parameter
     real, parameter       :: T_zero = 273.16  ! K, as in ConvPar_GF_GEOS5 
     real, parameter       :: T_ice  = 250.16  ! K, as in ConvPar_GF_GEOS5
     real, parameter       :: TEMP3  = 248.0   ! K
     real, parameter       :: TEMP4  = 268.0   ! K

     ! local variables
     real            :: fract_liq_f
     real            :: cldliq, cldice, c_h2o
     real            :: i2g, l2g
     real            :: c_tot, f_l, f_i
     real            :: airdens

     ! Compute environmental variables: cloud water and H2O mixing ratio.
     ! This corresponds to the computations done in SETUP_WETSCAV in module
     ! wetscav_mod.F90 in GEOS-Chem.

     ! compute cloud liquid water content and cloud ice water. 
     ! Compute either based on environmental variables or use original GEOS-Chem formulation
     if ( online_cldliq == 1.0 ) then
        ! compute from cloud total water, using formulation as suggested by Saulo Freitas
        fract_liq_f = min(1., (max(0.,(temp-T_ice))/(T_zero-T_ice))**2)
        ! liquid and ice water in kg/kg 
        cldliq = cldh2o * fract_liq_f       ! kg/kg
        cldice = cldh2o * (1.-fract_liq_f)  ! kg/kg
        ! to convert to cm3/cm3, need air density
        airdens = 100.*press/(287.04*temp*(1.+0.608*q))
        cldliq  = cldliq*airdens*1.e-3      ! cm3/cm3
        cldice  = cldice*airdens*1.e-3      ! cm3/cm33

     else
        ! original GEOS-Chem formulation
        IF ( temp >= TEMP4 ) THEN
           cldliq = 1e-6
        ELSE IF ( temp > TEMP3 .and. temp < TEMP4 ) THEN
           cldliq = 1e-6 * ((temp - 248.0) / 20.0 )
        ELSE
           cldliq = 0.0
        ENDIF
        cldliq = max(cldliq,0.0)      ! cm3 H2O/cm3 air
        cldice = max(1e-6-cldliq,0.0) ! cm3 ice/cm3 air
     endif

     ! mixing ratio of H2O [v/v]: compute using Dalton's law
     c_h2o = gcc_e_ice(temp) / press

     ! ice to gas ratio 
     i2g = 0.0
     if ( (liq_and_gas==1.) .and. (c_h2o>0.0) ) then
        i2g = ( cldice / c_h2o ) * convfaci2g
     endif

     ! liquid to gas ratio
     l2g = Heff * cldliq

     ! fraction of species in liquid & ice phases (Eqs. 4, 5, 6, Jacob et al, 2000)
     c_tot = 1.0 + l2g + i2g
     f_l   = l2g / c_tot
     f_i   = i2g / c_tot

     ! compute the rate constant Ki for loss of species from
     ! convective updraft scavenging (Eq. 1, Jacob et al, 2000)
     if ( temp >= TEMP4 ) then
        kc_scaled = KC_DEFAULT * ( f_l + f_l )
     else if ( temp > TEMP3 .and. temp < TEMP4 ) THEN
        kc_scaled = KC_DEFAULT * ( ( retfactor * f_L ) + F_I )
     else
        kc_scaled = KC_DEFAULT * f_i
     endif

   END SUBROUTINE compute_ki_gcc_gas
!-----------------------------------------------------------------------------------------


  FUNCTION gcc_e_ice( temp ) RESULT( vpress )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Calculate saturation vapor pressure over ice at given temperature - adapted from GEOS-Chem
    !  Marti & Mauersberber (GRL '93) formulation of saturation
    !  vapor pressure of ice [Pa] is: log P = A/TK + B
    !EOP
    !=====================================================================================
    real, intent(in) :: temp   ! Temperature [K]
    real             :: vpress ! Saturation vapor pressure [hPa]
    ! parameter
    real, parameter  :: A = -2663.5
    real, parameter  :: B =  12.537
    ! Saturation vap press of Ice [Pa]
    if ( temp <= 1.e-5 ) then
       vpress = 0.0
    else
       vpress = ( 10.**( A/temp + B ) )
    endif

  END FUNCTION gcc_e_ice
!-----------------------------------------------------------------------------------------

END MODULE geoschemchem_moist_interface
