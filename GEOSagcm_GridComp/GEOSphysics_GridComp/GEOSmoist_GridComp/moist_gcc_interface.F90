#include "MAPL_Generic.h"
!BOP
!
! !MODULE: moist_gcc_interface 
!
! !DESCRIPTION: This module contains routines and variables to control 
!  the convective washout of GEOS-Chem species by GEOS moist. This includes 
!  routines to compute the washout fraction of gases and aerosols, using the
!  same paraterizations as used within the GEOS-Chem CTM. These parameterizations
!  require a few additional species attributes (other than the Henry coefficients),
!  which are also obtained from the species bundles.\\
!  In addition to the washout parameterization routines, this module also contains 
!  routines to diagnose the column integrated loss of a species due to convective 
!  washout as well as the washout fraction per grid cell. These diagnostics are
!  made available for all species listed below. More species can be added to this
!  list as needed.\\
!  In order to avoid excessive memory usage / computations, diagnostics are only
!  calculated for the fields requested in HISTORY.rc.\\ 
!\\
!\\
! !INTERFACE:
!
MODULE moist_gcc_interface
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
  PUBLIC  :: get_w_upd_gcc
  PUBLIC  :: henry_gcc
  PUBLIC  :: GCC_AddExports
  PUBLIC  :: GCC_FillExportConvScav
  PUBLIC  :: GCC_FillExportConvFrac
  PUBLIC  :: GCC_check_params
  PUBLIC  :: GCC_get_ndiag
  PUBLIC  :: GCC_get_diagID
  PUBLIC  :: GCC_ConvFrac
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: gcc_e_ice 
!
! !REMARKS:
!  References:
!  D.J. Jacob, H. Liu, C. Mari, and R. M. Yantosca, "Harvard wet
!  deposition scheme for GMI", Harvard Atmospheric Chemistry Modeling 
!  Group, March 2000.
!
! !REVISION HISTORY:
!  20210315 - christoph.a.keller@nasa.gov - initial version
!  See https://github.com/GEOS-ESM/GEOSgcm_GridComp for full history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Default autoconversion parameter for GEOS-Chem species [s-1]
  REAL, PARAMETER       :: KC_DEFAULT_GCC = 5.e-3

! Species for which to provide diagnostics. Need to hardcode this since we create the export during SetServices,
! when we cannot yet import the species information from GEOS-Chem... 
 INTEGER, PARAMETER             :: GCCmax = 141
 CHARACTER(LEN=8), PARAMETER    :: GCCspecies(GCCmax) = &
                                            (/ 'ACTA',  &
                                               'AERI',  &
                                               'ALD2',  &
                                               'AONITA', &
                                               'AROMP4', &
                                               'AROMP5', &
                                               'ATOOH',  &
                                               'BALD',   &
                                               'BCPI',  &
                                               'BCPO',  &
                                               'BENZP', &
                                               'Br2',   &
                                               'BrCl',  &
                                               'BrSALA', &
                                               'BrSALC', &
                                               'BZCO3H', &
                                               'BZPAN',  &
                                               'CH2O',  &
                                               'CSL',   &
                                               'DST1',  &
                                               'DST2',  &
                                               'DST3',  &
                                               'DST4',  &
                                               'EOH',   &
                                               'ETHLN', &
                                               'ETHN',  &
                                               'ETHP',  &
                                               'ETP',   &
                                               'GLYC',  &
                                               'GLYX',  &
                                               'H2O2',  &
                                               'HBr',   &
                                               'HC5A',  &
                                               'HCl',   &
                                               'HCOOH', &
                                               'HI',    &
                                               'HMHP',  &
                                               'HMML',  &
                                               'HMS',   &
                                               'HNO3',  &
                                               'HOBr',  &
                                               'HOCl',  &
                                               'HOI',   &
                                               'HONIT', &
                                               'HPETHNL', &
                                               'I2',    &
                                               'I2O2',  &
                                               'I2O3',  &
                                               'I2O4',  &
                                               'IBr',   &
                                               'ICHE',  &
                                               'ICl',   &
                                               'ICN',   &
                                               'ICPDH', &
                                               'IDCHP', &
                                               'IDHDP', &
                                               'IDHPE', &
                                               'IEPOXA',&
                                               'IEPOXB',&
                                               'IEPOXD',&
                                               'IHN1',  &
                                               'IHN2',  &
                                               'IHN3',  &
                                               'IHN4',  &
                                               'INDIOL',&
                                               'INPB',  &
                                               'INPD',  &
                                               'IONITA',&
                                               'IONO',  &
                                               'IONO2', &
                                               'ISALA', &
                                               'ISALC', &
                                               'ITCN',  &
                                               'ITHN',  &
                                               'LIMO',  &
                                               'LVOC',  &
                                               'LVOCOA',&
                                               'MACR1OOH',&
                                               'MAP',   &
                                               'MCRDH', &
                                               'MCRENOL', &
                                               'MCRHN', &
                                               'MCRHNB', &
                                               'MCRHP', &
                                               'MCT',   &
                                               'MEK',   &
                                               'MGLY',  &
                                               'MOH',   &
                                               'MONITA', &
                                               'MONITS', &
                                               'MONITU', &
                                               'MP',    &
                                               'MPAN',  &
                                               'MPN',   &
                                               'MSA',   &
                                               'MTPA',  &
                                               'MTPO',  &
                                               'MVK',   &
                                               'MVKDH',  &
                                               'MVKHC',  &
                                               'MVKHCB', &
                                               'MVKHP',  &
                                               'MVKN',   &
                                               'MVKPC',  &
                                               'NH3',   &
                                               'NH4',   &
                                               'NIT',   &
                                               'NITs', &
                                               'NPHEN',  &
                                               'OCPI',  &
                                               'OCPO',  &
                                               'PAN',   &
                                               'pFe',   &
                                               'PHEN',  &
                                               'PP',    &
                                               'PPN',   &
                                               'PROPNN', &
                                               'PRPE',  &
                                               'PRPN',  &
                                               'PYAC',  &
                                               'R4N2',  &
                                               'R4P',   &
                                               'RA3P',  &
                                               'RB3P',  &
                                               'RIPA',  &
                                               'RIPB',  &
                                               'RIPC',  &
                                               'RIPD',  &
                                               'RP',    &
                                               'SALA',  &
                                               'SALAAL', &
                                               'SALACL', &
                                               'SALC',  &
                                               'SALCAL', &
                                               'SALCCL', &
                                               'SO2',   &
                                               'SO4',   & 
                                               'SO4s',  & 
                                               'SOAGX', & 
                                               'SOAIE', & 
                                               'SOAS' /)

   ! Name templates for diagnostics
   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER     :: Prefix_ConvScav = 'GCC_ConvScav_GF_'
   CHARACTER(LEN=ESMF_MAXSTR), PARAMETER     :: Prefix_ConvFrac = 'GCC_ConvFrac_GF_'

   ! Helper arrays to keep track of requested diagnostics. If the diagnostics are to be filled for one
   ! of the species, the corresponding slot in the array below will be filled with the species index
   ! as used by MOIST.
   INTEGER, SAVE     :: ConvScavDiag(GCCmax) = -1
   INTEGER, SAVE     :: ConvFracDiag(GCCmax) = -1

   ! Public array to hold the convective fractions of all requested species. These are 3D fields.
   REAL, ALLOCATABLE :: GCC_ConvFrac(:,:,:,:)

CONTAINS
!EOC

!---------------------------------------------------------------------------------------------------
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Create the GF diagnostics exports for selected GEOS-Chem species 
    !EOP
    !=====================================================================================
  SUBROUTINE GCC_AddExports( GC, RC )
    TYPE(ESMF_GridComp), INTENT(INOUT)  :: GC       ! Ref to GridComp
    INTEGER,             INTENT(OUT)    :: RC       ! Success or failure
    !LOCAL VARIABLES:
    INTEGER                         :: I
    CHARACTER(LEN=ESMF_MAXSTR)      :: spcname

    __Iam__('GCCdiag_AddExports')

    DO I = 1,GCCmax
       spcname = TRIM(GCCspecies(I))

       call MAPL_AddExportSpec(GC,                                             &
         SHORT_NAME=TRIM(Prefix_ConvScav)//TRIM(spcname),                         &
         LONG_NAME ='GEOS-Chem_'//TRIM(spcname)//'_dry_air_tendency_due_to_GF_conv_scav', &
         UNITS     ='kg m-2 s-1',                                              &
         DIMS      = MAPL_DimsHorzOnly,                                        &
          __RC__ )

       call MAPL_AddExportSpec(GC,                                             &
         SHORT_NAME=TRIM(Prefix_ConvFrac)//TRIM(spcname),                  &
         LONG_NAME ='GEOS-Chem_'//TRIM(spcname)//'_fraction_lost_in_GF_convection', &
         UNITS     ='1',                                                       &
         DIMS      = MAPL_DimsHorzVert,                                        &
         VLOCATION = MAPL_VLocationCenter,                                     &
          __RC__ )
    ENDDO

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCC_AddExports
!EOC

!---------------------------------------------------------------------------------------------------
  SUBROUTINE GCC_FillExportConvScav( EXPORT, IM, JM, Arr2d, DiagID, RC ) 
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Fill the GEOS-Chem exports fro the convective scavenging 
    !EOP
    !=====================================================================================
    TYPE(ESMF_State), INTENT(INOUT)   :: EXPORT    
    INTEGER,          INTENT(IN)      :: IM, JM    
    REAL,             INTENT(IN)      :: Arr2d(IM,JM) 
    INTEGER,          INTENT(IN)      :: DiagID
    INTEGER,          INTENT(OUT)     :: RC
    ! Local variables
    INTEGER                           :: I, II, cnt
    CHARACTER(LEN=ESMF_MAXSTR)        :: diagname
    REAL, POINTER, DIMENSION(:,:)     :: GCptr2d

    __Iam__('GCC_FillExportConvScav')

    ! Get species index for the given diagnostics index. The passed diagnostics index
    ! denotes the nth active diagnostics in the list of GEOS-Chem species available for
    ! diagnostics. Need to map it back to the actual index in the full list of GEOS-Chem
    ! species as specified in array GCCspecies
    cnt = 0
    II = -1
    do I = 1,GCCmax
        if ( ConvScavDiag(I) > 0 ) cnt = cnt + 1 
        if ( cnt == DiagID ) then
           II = I
           exit
        endif
    enddo
    ASSERT_(II>0)
    ! Now that species index is known, can construct diagnostics name and fill it
    diagname = TRIM(Prefix_ConvScav)//TRIM(GCCspecies(II))
    call MAPL_GetPointer(EXPORT, GCptr2d, TRIM(diagname), NotFoundOk=.TRUE., __RC__ )
    if ( associated(GCptr2d) ) GCptr2d = Arr2d
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCC_FillExportConvScav
!EOC

!---------------------------------------------------------------------------------------------------
  SUBROUTINE GCC_FillExportConvFrac( EXPORT, DiagID, RC ) 
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Fill the GEOS-Chem exports for the convective fractions.
    !EOP
    !=====================================================================================
    TYPE(ESMF_State), INTENT(INOUT)   :: EXPORT    
    INTEGER,          INTENT(IN)      :: DiagID
    INTEGER,          INTENT(OUT)     :: RC
    ! Local variables
    INTEGER                           :: I, II, cnt
    CHARACTER(LEN=ESMF_MAXSTR)        :: diagname
    REAL, POINTER, DIMENSION(:,:,:)   :: GCptr3d

    __Iam__('GCC_FillExportConvFrac')

    if ( allocated(GCC_ConvFrac) ) then
       ! Get species index for the given diagnostics index. The passed diagnostics index
       ! denotes the nth active diagnostics in the list of GEOS-Chem species available for
       ! diagnostics. Need to map it back to the actual index in the full list of GEOS-Chem
       ! species as specified in array GCCspecies
       cnt = 0
       II = -1
       do I = 1,GCCmax
           if ( ConvFracDiag(I) > 0 ) cnt = cnt + 1 
           if ( cnt == DiagID ) then
              II = I
              exit
           endif
       enddo
       ASSERT_(II>0)
       ! Now that species index is known, can construct diagnostics name and fill it
       diagname = TRIM(Prefix_ConvFrac)//TRIM(GCCspecies(II))
       call MAPL_GetPointer(EXPORT, GCptr3d, TRIM(diagname), NotFoundOk=.TRUE., __RC__ )
       if ( associated(GCptr3d) ) GCptr3d = GCC_ConvFrac(:,:,:,DiagID)
    endif

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCC_FillExportConvFrac
!EOC

!---------------------------------------------------------------------------------------------------
  SUBROUTINE GCC_check_params( EXPORT, k, SpcName, FIELD, is_gcc, Vect_KcScal, retfactor, &
                               liq_and_gas, convfaci2g, online_cldliq, online_vud, use_gocart, &
                               ftemp_threshold, RC )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Check for GEOS-Chem parameters 
    !EOP
    !=====================================================================================
    TYPE(ESMF_State), INTENT(INOUT)   :: EXPORT    ! Export state
    INTEGER,          INTENT(IN)      :: k             ! Species ID in MOIST
    CHARACTER(LEN=*), INTENT(IN)      :: spcname
    TYPE(ESMF_FIELD), INTENT(INOUT)   :: FIELD
    LOGICAL, INTENT(OUT)              :: is_gcc
    REAL, DIMENSION(3), INTENT(OUT)   :: Vect_KcScal
    REAL, INTENT(OUT)                 :: retfactor
    REAL, INTENT(OUT)                 :: liq_and_gas
    REAL, INTENT(OUT)                 :: convfaci2g
    REAL, INTENT(OUT)                 :: online_cldliq
    REAL, INTENT(OUT)                 :: online_vud
    LOGICAL, INTENT(OUT)              :: use_gocart 
    REAL, INTENT(OUT)                 :: ftemp_threshold
    INTEGER, INTENT(OUT)              :: RC       ! Success or failure
    ! Local variables
    CHARACTER(LEN=ESMF_MAXSTR)        :: shortname
    CHARACTER(LEN=ESMF_MAXSTR)        :: diagname
    INTEGER                           :: I
    REAL                              :: doit
    LOGICAL                           :: isPresent
    REAL, POINTER, DIMENSION(:,:)     :: GCptr2d
    REAL, POINTER, DIMENSION(:,:,:)   :: GCptr3d

    __Iam__('GCC_check_params')

    ! Default initial values
    is_gcc          = .FALSE.
    Vect_KcScal(:)  = 1.0
    retfactor       = 1.0
    liq_and_gas     = 0.0
    online_cldliq   = 0.0
    online_vud      = 1.0
    use_gocart      = .FALSE.
    ftemp_threshold = -999.0
    ! check if this is a GEOS-Chem species
    if ( LEN(TRIM(SpcName)) > 4 ) then
    if ( TRIM(SpcName(1:4)) == 'SPC_' ) then
       is_gcc = .TRUE.
       ! KC scale factors for GEOS-Chem
       call ESMF_AttributeGet  (FIELD,"SetofKcScalFactors",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet  (FIELD,"SetofKcScalFactors",Vect_KcScal, __RC__ )
       endif
       ! Gas-phase washout parameter for GEOS-Chem
       call ESMF_AttributeGet (FIELD,"RetentionFactor",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"RetentionFactor",retfactor, __RC__ )
       endif
       call ESMF_AttributeGet (FIELD,"LiqAndGas",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"LiqAndGas",liq_and_gas, __RC__ )
       endif
       call ESMF_AttributeGet (FIELD,"ConvFacI2G",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"ConvFacI2G",convfaci2g, __RC__ )
       endif
       call ESMF_AttributeGet (FIELD,"OnlineCLDLIQ",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"OnlineCLDLIQ",online_cldliq, __RC__ )
       endif
       call ESMF_AttributeGet (FIELD,"OnlineVUD",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"OnlineVUD",online_vud, __RC__ )
       endif
       call ESMF_AttributeGet (FIELD,"UseGOCART",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"UseGOCART",doit, __RC__ )
          if ( doit==1.0 ) use_gocart = .TRUE.
       endif
       call ESMF_AttributeGet (FIELD,"GOCARTfTempThreshold",isPresent=isPresent, __RC__ )
       if (isPresent) then
          call ESMF_AttributeGet (FIELD,"GOCARTfTempThreshold",ftemp_threshold, __RC__ )
       endif
       ! check if exports are requested for this species. If so, store the species 
       ! index as used by MOIST in the corresponding slot in the local diagnostics
       ! counter array
       shortname = SpcName(5:LEN(TRIM(SpcName)))
       diagname = TRIM(Prefix_ConvScav)//TRIM(shortname)
       call MAPL_GetPointer(EXPORT, GCptr2d, TRIM(diagname), NotFoundOk=.TRUE., __RC__ )
       if ( associated(GCptr2d) ) then
          do I=1,GCCmax
             if ( TRIM(GCCspecies(I)) == TRIM(shortname) ) then
                ConvScavDiag(I) = k
                exit
             endif
          enddo   
       endif
       diagname = TRIM(Prefix_ConvFrac)//TRIM(shortname)
       call MAPL_GetPointer(EXPORT, GCptr3d, TRIM(diagname), NotFoundOk=.TRUE., __RC__ )
       if ( associated(GCptr3d) ) then
          do I=1,GCCmax
             if ( TRIM(GCCspecies(I)) == TRIM(shortname) ) then
                ConvFracDiag(I) = k
                exit
             endif
          enddo   
       endif
    end if
    end if

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE GCC_check_params
!EOC

!---------------------------------------------------------------------------------------------------
  FUNCTION GCC_get_ndiag( diagtype ) RESULT ( ndiag )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Returns the number of 'active' diagnostics for type 1 or 2 
    !EOP
    !=====================================================================================
    integer, intent(in) :: diagtype
    integer             :: ndiag
    ! local variables
    integer             :: I
    ! starts here
    ndiag = 0
    if ( diagtype == 1 ) then
       do I=1,GCCmax
           if ( ConvScavDiag(I) > 0 ) ndiag = ndiag + 1
       enddo
    end if
    if ( diagtype == 2 ) then
       do I=1,GCCmax
           if ( ConvFracDiag(I) > 0 ) ndiag = ndiag + 1
       enddo
    end if
  END FUNCTION GCC_get_ndiag 

!---------------------------------------------------------------------------------------------------
  FUNCTION GCC_get_diagID( diagtype, specID ) RESULT ( diagID )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Maps the MOIST species index to the corresponding GEOS-Chem diagnstics list 
    !EOP
    !=====================================================================================
    integer, intent(in) :: diagtype   ! diagnostics type
    integer, intent(in) :: specID     ! species index as counted by MOIST
    integer             :: diagID     ! GEOS-chem diagnostics ID
    ! local variables
    integer             :: I, diagCount
    ! starts here
    diagID = 0
    diagCount = 0
    if ( diagtype == 1 ) then
       do I=1,GCCmax
           if ( ConvScavDiag(I) > 0 ) diagCount = diagCount + 1
           if ( ConvScavDiag(I) == specID ) then
              diagID = diagCount
              exit
           endif
       enddo
    end if
    if ( diagtype == 2 ) then
       do I=1,GCCmax
           if ( ConvFracDiag(I) > 0 ) diagCount = diagCount + 1
           if ( ConvFracDiag(I) == specID ) then
              diagID = diagCount
              exit
           endif
       enddo
    end if
  END FUNCTION GCC_get_diagID

!---------------------------------------------------------------------------------------------------
  FUNCTION henry_gcc( hstar, dhr, ak0, dak, temp ) RESULT( henry_coeff )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Return henry coefficient (liquid to gas) as defined by GEOS-Chem
    !EOP
    !=====================================================================================
    real, intent(in) :: hstar               ! Henry coefficient [M/atm]]
    real, intent(in) :: dhr                 ! temperature dependency of hstar [-d ln(kH) / d(1/T)]
    real, intent(in) :: ak0                 ! pKa value [-]
    real, intent(in) :: dak                 ! temperature dependency of ak0, currently not used
    real, intent(in) :: temp                ! ambient temperature [K]
    real             :: henry_coeff         ! effective gas/aq constant [-] 
    ! parameter
    real*8, parameter :: pH   = 4.5d0
    REAL*8, PARAMETER :: TREF = 298.15d0        ! [K          ]
    REAL*8, PARAMETER :: R    = 8.3144598d0     ! [J K-1 mol-1]
    REAL*8, PARAMETER :: ATM  = 101.325d0       ! [mPa (!)    ]
    ! local variables
    real*8          :: hstar8, dhr8, ak08, temp8, h8
    ! cast all variables to r*8 locally to prevent overflows
    hstar8 = hstar
    dhr8   = dhr
    ak08   = ak0
    temp8  = temp
    ! calculate henry coefficient
    h8 = hstar8 * exp ( dhr8 * (1./temp8 - 1./TREF) ) * R * temp8 / ATM
    if ( ak08 > 0.0d0 ) then
       h8 = h8 * ( 1.0 + 10.0**(pH-ak08) )
    endif
    ! limit henry coefficient to 1.0e30
    henry_coeff = real( min(h8,1.0d+30) )

  END FUNCTION henry_gcc

!---------------------------------------------------------------------------------------------------
   SUBROUTINE compute_ki_gcc_gas( temp, press, q, cldh2o, Heff, liq_and_gas, convfaci2g, retfactor, online_cldliq, kc_scaled, l2g )
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
     real,    intent(in)      :: Heff               ! effective gas/aq Henry constant [-]
     real,    intent(in)      :: liq_and_gas        ! species considers ice and liquid phase?
     real,    intent(in)      :: convfaci2g         ! conversion factor for ice/gas ratio
     real,    intent(in)      :: retfactor          ! retention factor [unitless] 
     real,    intent(in)      :: online_cldliq      ! calculate cloud liquid/ice online or use default GEOS-Chem parameterization
     real,    intent(out)     :: kc_scaled          ! loss rate [s-1]
     real,    intent(out)     :: l2g                ! liquid to gas ratio 

     ! parameter
     real, parameter       :: T_zero = 273.16  ! K, as in ConvPar_GF_GEOS5 
     real, parameter       :: T_ice  = 250.16  ! K, as in ConvPar_GF_GEOS5
     real, parameter       :: TEMP3  = 248.0   ! K
     real, parameter       :: TEMP4  = 268.0   ! K

     ! local variables
     real            :: fract_liq_f
     real            :: cldliq, cldice, c_h2o
     real            :: i2g
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
        cldice  = cldice*airdens*1.e-3      ! cm3/cm3

     else
        ! original GEOS-Chem formulation
        IF ( temp >= TEMP4 ) THEN
           cldliq = 1e-6
        ELSE IF ( temp > TEMP3 .and. temp < TEMP4 ) THEN
           cldliq = 1e-6 * ((temp-TEMP3)/(TEMP4-TEMP3))
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
        kc_scaled = KC_DEFAULT_GCC * ( f_l + f_l )
     else if ( temp > TEMP3 .and. temp < TEMP4 ) THEN
        kc_scaled = KC_DEFAULT_GCC * ( ( retfactor * f_l ) + f_i )
     else
        kc_scaled = KC_DEFAULT_GCC * f_i
     endif

   END SUBROUTINE compute_ki_gcc_gas

!---------------------------------------------------------------------------------------------------
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
    kc_scaled = KC_DEFAULT_GCC
    if ( temp < TEMP1 ) then
       kc_scaled = kc_scaled * kcscal1
    else if ( (temp>=TEMP1) .and. (temp<TEMP2) ) then
       kc_scaled = kc_scaled * kcscal2
    else
       kc_scaled = kc_scaled * kcscal3
    endif

  END SUBROUTINE compute_ki_gcc_aerosol

!-----------------------------------------------------------------------------------------
  FUNCTION get_w_upd_gcc( vud, xland, online_vud ) RESULT( w_upd )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Return updraft vertical velocity to be used for GEOS-Chem convective washout. 
    !EOP
    !=====================================================================================
    real, intent(in) :: vud        ! online updraft velocity [m/s]
    real, intent(in) :: xland      ! land flag (1.-FRLAND): greater value means more water 
    real, intent(in) :: online_vud ! use online vud (1.0) or set vud based on land/water (0.0)
    real             :: w_upd      ! updraft velocity to use 
    ! use environment vud if specified so 
    if ( online_vud == 1.0 ) then
       w_upd = vud
    ! use parameterization otherwise: 10m/s over land, 5m/s over water.
    else
       ! over water
       if ( xland > 0.9 ) then
          w_upd = 5.0
       ! over land
       else
          w_upd = 10.0
       endif
    endif

  END FUNCTION get_w_upd_gcc

!---------------------------------------------------------------------------------------------------
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

!EOC
END MODULE moist_gcc_interface 
