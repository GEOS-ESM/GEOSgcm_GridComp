#include "MAPL_Generic.h"
!BOP
!
! !MODULE: moist_gcc_interface 
!
! !DESCRIPTION: This module contains routines and variables to control 
!  the convective washout of GEOS-Chem species by GEOS moist. 
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
!  20210312 - christoph.a.keller@nasa.gov - initial version
!  See https://github.com/GEOS-ESM/GEOSgcm_GridComp for full history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  REAL, PARAMETER       :: KC_DEFAULT_GCC = 5.e-3  ! s-1  ! Default autoconversion parameter for GEOS-Chem species

CONTAINS
!EOC

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
    !write(*,*) 'henry_gcc (B)',hstar,dhr,ak0,dak,temp,henry_coeff

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
