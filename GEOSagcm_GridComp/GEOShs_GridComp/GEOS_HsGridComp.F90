!  $Id: GEOS_HsGridComp.F90,v 1.27 2007/05/17 13:09:31 f4mjs Exp $

#include "MAPL_Generic.h"

  module GEOS_HSGridCompMod

!BOP
! \renewcommand{\comp}{{\tt GEOS\_HSGridCompMod}}
!
! !MODULE: GEOS_HSGridCompMod
!
! !DESCRIPTION:
!
!  \comp  implements the Held and Suarez (1994) benchmark
!  forcing for global atmospheric dynamical cores.
!  The Williamson, et al. (1998) stratospheric modifications
!  to the benchmark are also implemented, and a mechanism is included to add
!  a localized, gaussian heat source.
!
!
!  {\it Scientific Basis:} 
!
! {\bf The Held-Suarez forcing:}
!The wind tendencies are given by
!$$
!\frac{\partial {\bf v}}{\partial t} = \cdots - k_v(\sigma) \, {\bf v}
!$$
!where
!$$
! k_v(\sigma) = k_f f_1(\sigma)
!$$
!and
!$$
! f_1(\sigma) = \max \left[\,0, \frac{\sigma - \sigma_b}{1 - \sigma_b} \right]
!$$
!
! The heating tendencies are given by:
!$$
!\frac{\partial T}{\partial t} = \cdots - k_T(\phi,\sigma) \left[\, T - T_{eq}(\phi,p) \,\right]
!$$
!$$
!T_{eq} = \max \left[   T_{strat}, 
!   \left(\frac{p}{p_\circ}\right)^\kappa 
!   \left( T_0  - (\Delta T)_y\,\sin^2\!\phi 
!          - (\Delta\theta)_z  \log\left(\frac{p}{p_\circ}\right) \cos^2\!\phi\, 
!   \right)     \right] 
!$$
!where
!$$
!k_T = k_a + (k_s - k_a) \cos^4\!\phi \, f_1(\sigma).
!$$
!
! The parameter values suggested by Held and Suarez (1994) are:
! $$\sigma_b = .7,$$
! $$k_f = 1 \hbox{ day}^{-1},$$
! $$k_a = \frac{1}{40 \hbox{ days}},$$
! $$k_s = \frac{1}{4 \hbox{ days}},$$
! $$(\Delta T)_y = 60  \hbox{K },$$
! $$(\Delta \theta)_z = 10 \hbox{K },$$
! $$T_{strat} = 200 \hbox{K },$$
! $$T_0 = 315 \hbox{K }.$$
!
! {\bf The Williamson Modification:}
!  Williamson's modification in the stratosphere simply alters $T_{eq}$ at upper levels.
!  For $p<p_D$, the constant HS value is replaced with the following distribution:
!$$
!T_{eq} = T_{strat} \,\left[
! \min \left(1,\frac{p}{p_D}\right)^\frac{R_d \gamma_D}{g} +  
! \min \left(1,\frac{p}{p_I}\right)^\frac{R_d \gamma_I}{g} - 1 \right], 
!$$
! where $T_{strat}$ is as in HS,
!$$
!  p_I = p_D - \frac{p_D-p(\sigma=0)}{2}( 1 + \tanh(\alpha\frac{|\phi|-\phi_0}{(\delta\phi)_0})),
!$$
!and defaults for the other parameters are
!  $\gamma_I=-3.345 \times 10^{-3}$ K m$^{-1}$, 
!                $\gamma_D=2.0 \times 10^{-3}$ K m$^{-1}$,
!$p_D = 100 \hbox{ hPa} $,
!$\alpha = 2.65 $,
!$(\delta\phi)_0 = 15^\circ $,
!$\phi_0 = 60^\circ $.
!
! {\bf The Local Heat Source:}
! A local heat form can also be included. It has the form
!$$
!  Q(\lambda,\phi) = Q_{max} e^{ -\frac{1}{2} 
!    \left( \frac{(\phi_h-\phi)^2}{(\delta\phi)_h^2}
!   +       \frac{(\lambda_h-\lambda)^2}{(\delta\lambda)_h^2} \right) }
!     h(p),
!$$
!where $\lambda$ is the longitude, $p_s$ is the surface pressure,
!\begin{equation*}
!     h(p) = \left\{ \begin{array}{ll}
!                     0   &  \mbox{ $p < p_1$} \\

!                     \sin(\pi  \frac{p_s-p}{p_s-p_1}) & \mbox{$ p_1 \le p < p_s$}, \\
!                   \end{array}
!            \right.
!\end{equation*}
!  $Q_{max}$ is in K day$^{-1}$, with a default of zero, and defaults for the
!  other parameters are $\phi_h=\lambda_h=0$, 
!  $(\delta\phi)_h=5^\circ$, $(\delta\lambda)_h=30^\circ$,
!  $p_1=$0 hPa.
!  
!
!
!  {\it Code Implementation:} 
!  
!  
!   All HS parameters are
!  optionally adjustable from the configuration, with published values
!  being the defaults.  Its parameters are also adjustable
!  from the configuration and also default to published values. The Williamson
!  modification can be disabled by choosing $p_D=0$. Finally, we also include
!  a mechanism for adding a localized, gaussian heat source whose strength and shape can
!  be controlled from the configuration. This is zero by default.
!
! {\it References:}
!
!  Held, I. M., and M. J. Suarez, 1994: A proposal for the intercomparison
!  of the dynamical cores of atmospheric general circulation models. {\em Bulletin
!  of the American Meteorological Society}, {\bf 75(10)}, 1825-1830. 
!
!  Williamson, D.L., J. G. Olson, and B. A. Boville, 1998: A comparison of 
!  semi-Lagrangian and Eulerian tropical climate simulations. {\em Mon. Wea. Rev.}, 
!  {\bf 126}, 1001-1012
!

! !USES:

  use ESMF
  use MAPL_Mod
  use hs_oacc_mod
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all

  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:

  public  SetServices

!EOP

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !DESCRIPTION:  This version uses the MAPL\_GenericSetServices, which sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of the MAPL\_MetaComp and putting it in the 
!   gridded component (GC). The MAPL\_MetaComp contains an ESMF state to use as
!   an internal state. This routine describes the contents of this Internal state,
!   as well as of the conventional Import and Export states by making call 
!   to MAPL.
!
!

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional,   intent(  OUT) :: RC  ! return code

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm
    integer                           :: STATUS
    character(len=ESMF_MAXSTR)        :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Set Data Services for the GC
! ----------------------------

!BOP

! !IROUTINE: State Descriptions

! !DESCRIPTION: The component uses all three states (Import, Export
!  and Internal). There is no Private (non-ESMF) Internal state. All
!  three are managed by MAPL. The Internal state contains only invariant,
!  horizontally dependent quantities that are set in Initialize, therefore
!  it never needs checkpointing.


! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DEFAULT    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DEFAULT    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'TEMP',                                      &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DEFAULT    = 300.0,                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'PLE',                                       &
         LONG_NAME  = 'edge_pressure',                             &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
                                                        RC=STATUS  )
     VERIFY_(STATUS)


! !INTERNAL STATE:

    call MAPL_AddInternalSpec(gc,                                  &
         SHORT_NAME = 'SPHI2',                                     &
         LONG_NAME  = 'sine_latitude_squared',                     &
         UNITS      = '1',                                         &
         DEFAULT    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RESTART    = MAPL_RestartSkip,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec(gc,                                  &
         SHORT_NAME = 'CPHI2',                                     &
         LONG_NAME  = 'cosine_latitude_squared',                   &
         UNITS      = '1',                                         &
         DEFAULT    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RESTART    = MAPL_RestartSkip,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec(gc,                                  &
         SHORT_NAME = 'HFCN',                                      &
         LONG_NAME  = 'horizontal_structure_of_localized_heating', &
         UNITS      = '1',                                         &
         DEFAULT    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RESTART    = MAPL_RestartSkip,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec(gc,                                  &
         SHORT_NAME = 'P_I',                                       &
         LONG_NAME  = 'Williamson_interface_pressure',             &
         UNITS      = '1',                                         &
         DEFAULT    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RESTART    = MAPL_RestartSkip,                        &
                                                        RC=STATUS  )
     VERIFY_(STATUS)


! !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'pressure_weighted_air_temperature_tendency',&
         UNITS      = 'Pa K s-1',                                  &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_tendency',                    &
         UNITS      = 'm s-2',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    & 
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_tendency',                   &
         UNITS      = 'm s-2',                                     &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    & 
         SHORT_NAME = 'T_EQ',                                      &
         LONG_NAME  = 'equilibrium_temperature',                   &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    & 
         SHORT_NAME = 'THEQ',                                      &
         LONG_NAME  = 'equilibrium_potential_temperature',         &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'TAUX',                                      &
         LONG_NAME  = 'eastward_surface_stress',                   &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'TAUY',                                      &
         LONG_NAME  = 'northward_surface_stress',                  &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
                                                        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    & 
         SHORT_NAME = 'DISS',                                      &
         LONG_NAME  = 'frictional_dissipation',                    & 
         UNITS      = 'W m-2',                                     &
         DIMS       =  MAPL_DimsHorzOnly,                          &
         VLOCATION  =  MAPL_VLocationNone,                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
!EOP

! Set the Initialize and Run entry points
! --------------------------------------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE, Initialize, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN"   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="INIT"  ,RC=STATUS); VERIFY_(STATUS)

! SetServices clean-up on the way back up through the hierarchy
!--------------------------------------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP
! !IROUTINE: Initialize
!
! !DESCRIPTION: Here we initialize the internal state, which in HS contains
!   two-dimensional invariant arrays used in the forcing calculations. They
!   are done here simply for economy. Initialize calls MAPL\_GenericInitialize.  
!   If the import state needs to be restarted, MAPL will do it 
!   if a restart file is provided. 
!
! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

    use openacc
    
! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

 
!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)    :: IAm 
    integer                       :: STATUS
    character(len=ESMF_MAXSTR)    :: COMP_NAME

! Pointers to geography info in the MAPL MetaComp

    real,                 pointer :: LATS (:,:)
    real,                 pointer :: LONS (:,:)

! Pointers to internals

    real,                 pointer :: SPHI2(:,:)
    real,                 pointer :: CPHI2(:,:)
    real,                 pointer :: HFCN (:,:)
    real,                 pointer :: P_I  (:,:)

! Local variables

    real                          :: DX, DY, X0, Y0
    real                          :: AFAC, PHI0, QMAX
    logical                       :: halting_mode(5)

! Local derived type aliases

    type (MAPL_MetaComp), pointer :: MAPL
    type (ESMF_State   )          :: INTERNAL 

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet       ( GC, name=COMP_NAME,        RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my instance of the MAPL MetaComp
!-----------------------------------

    call MAPL_GetObjectFromGC   ( GC, MAPL,                   RC=STATUS)
    VERIFY_(STATUS)

! Call Generic Initialize 
! -----------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

! Start total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"INIT" )

! Get coordinate information from MAPL_MetaComp
!--------------------------------------------

    call MAPL_Get(MAPL,                &
       LATS          = LATS,           & ! These are in radians
       LONS          = LONS,           & ! These are in radians
       INTERNAL_ESMF_STATE=INTERNAL,   &
                             RC=STATUS )
    VERIFY_(STATUS)

! Get pointers to internal variables
!-----------------------------------

    call MAPL_GetPointer(INTERNAL, SPHI2, 'SPHI2', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CPHI2, 'CPHI2', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, HFCN , 'HFCN' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, P_I  , 'P_I'  , RC=STATUS); VERIFY_(STATUS)

! Initialize geometric factors
!-----------------------------

    SPHI2 = sin(LATS)**2
    CPHI2 = cos(LATS)**2

! Precompute Local heating distribution
!--------------------------------------


!BOR

! !RESOURCE_ITEM: K day-1 :: Local heating 
    call MAPL_GetResource(MAPL,QMAX,'QMAX:',default=.000    , RC=STATUS ) 
    VERIFY_(STATUS)

    if(QMAX /= 0.0) then

! Get local heating parameters from the configuration
!----------------------------------------------------

! !RESOURCE_ITEM: degrees :: Central longitude of local heating, $\lambda_h$
       call MAPL_GetResource(MAPL,X0  ,'X0:'  ,default=0.0     , RC=STATUS)
       VERIFY_(STATUS)                                               
! !RESOURCE_ITEM: degrees :: Central latitude of local heating, $\phi_h$
       call MAPL_GetResource(MAPL,Y0  ,'Y0:'  ,default=0.0     , RC=STATUS)
       VERIFY_(STATUS)                                               
! !RESOURCE_ITEM: degrees :: Longitudinal width of local heating, $(\delta\lambda)_h$
       call MAPL_GetResource(MAPL,DX  ,'DX:'  ,default=30.     , RC=STATUS)
       VERIFY_(STATUS)                                               
! !RESOURCE_ITEM: degrees :: Latitudinal width of local heating, $(\delta\phi)_h$
       call MAPL_GetResource(MAPL,DY  ,'DY:'  ,default=5.0     , RC=STATUS)
       VERIFY_(STATUS)                             

       ! Horizontal structure of optional stationary heat source
       !--------------------------------------------------------

       X0 = X0 * (MAPL_PI/180.0)
       Y0 = Y0 * (MAPL_PI/180.0)
       DX = DX * (MAPL_PI/180.0)
       DY = DY * (MAPL_PI/180.0)

       HFCN =          ((LONS     -X0)/DX)**2
       HFCN = min(HFCN,((LONS+360.-X0)/DX)**2)
       HFCN = min(HFCN,((LONS-360.-X0)/DX)**2)

       HFCN = HFCN +   ((LATS     -Y0)/DY)**2

       where(HFCN < 10.)
          HFCN = exp(-0.5*HFCN)
       elsewhere
          HFCN = 0.0
       end where
    else
       HFCN = 0.0
    endif
                        
! Get Williamson parameters from the configuration
!-------------------------------------------------
                                                
! !RESOURCE_ITEM: 1 :: Williamson's $\alpha$ parameter
    call MAPL_GetResource(MAPL,AFAC,'AFAC:',default=2.65/15., RC=STATUS)
    VERIFY_(STATUS)
! !RESOURCE_ITEM: degrees :: Williamson's $\phi_o$ parameter
    call MAPL_GetResource(MAPL,PHI0,'PHI0:',default=60.     , RC=STATUS)
    VERIFY_(STATUS)                                      
!EOR
                                           
! Horizontal structure of Willianson parameterization
!----------------------------------------------------

    P_I = ( 1.0 + tanh( AFAC*( abs(LATS) - PHI0 ) ) )

! All Done
!---------

    call MAPL_TimerOff(MAPL,"INIT" )
    call MAPL_TimerOff(MAPL,"TOTAL")

    !OACC GCC10 Implementation can trigger floating-point exception
    !per https://gcc.gnu.org/wiki/OpenACC/Implementation%20Status
    !we saw those on PTX load so we deactivate the ieee mode around
    !the init
    call ieee_get_halting_mode(ieee_all, halting_mode)
    call ieee_set_halting_mode(ieee_all, .false.)
    call acc_set_device_type(acc_device_default)
    call ieee_set_halting_mode(ieee_all, halting_mode)

    RETURN_(ESMF_SUCCESS)
  end subroutine INITIALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!BOP

! !IROUTINE: Run

! !DESCRIPTION: The Run method of the HS Gridded Component. It computes
!  tendencies for temperature and winds as described above,
!  as well as a number of
!  diagnostics. These are all optional Exports. The input winds 
!  and/or temperatures may be Friendly. All Friendlies are updated
!  using a time step
!  obtained from 'RUN\_DT:' in the configuration. Both U and V must be
!  Friendly for the wind to be updated. It is an error for one to be Friendly
!  and not the other.  
!
! !INTERFACE:

  subroutine Run( GC, IMPORT, EXPORT, CLOCK, RC )
    include 'mpif.h'

    INTERFACE
        subroutine c_call_hs_oacc(CPHI2, DISS, DTDT, DUDT, DVDT, &
                          	HFCN, P_I, PLE, SPHI2, TAUX, TAUY, T, &
                          	THEQ, T_EQ, U, V, &
                          	DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
                          	FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
                          	SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank) bind(C, name='c_call_hs_oacc')
        
            use, intrinsic :: iso_c_binding

            integer(C_INT),  value :: FRICQ, IM, JM, LM, compType, rank
            ! The logical variables are declared with (C_INT) since (C_BOOL) seems
            ! to be 1 byte, while Fortran (with gfortran) uses 4-byte logicals
            logical(C_INT), value :: FriendlyTemp, FriendlyWind
            real(C_FLOAT),   value :: DAYLEN, DELH, DELV1, DT, GAM_D, GAM_I, P_D, P_1
            real(C_FLOAT),   value :: QMAX, SIG1, T0, TAUA, TAUF, TAUS, TSTRT

            type(C_PTR), value :: CPHI2, HFCN, P_I, SPHI2
            type(C_PTR), value :: DISS, TAUX, TAUY
            type(C_PTR), value :: PLE
            type(C_PTR), value :: DTDT, DUDT, DVDT, T, T_EQ, THEQ, U, V
        end subroutine
    end interface

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

 
!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)        :: IAm 
    integer                           :: STATUS
    character(len=ESMF_MAXSTR)        :: COMP_NAME
    
! Local derived type aliases

    type (MAPL_MetaComp    ), pointer :: MAPL
    type (ESMF_State       )          :: INTERNAL
    type (ESMF_Grid        )          :: GRID
    type (ESMF_Logical     )          :: Friendly

    type(ESMF_VM) :: vm

    integer :: comm, rank, ierr, length
    character * (MPI_MAX_PROCESSOR_NAME) nameNode

! Pointers to inputs

    real, pointer, dimension(:,:,:)   :: PLE, T, U, V

! Pointers to internals

    real, pointer, dimension(:,:  )   :: HFCN, SPHI2, CPHI2, P_I

! Pointers to outputs

    real, pointer, dimension(:,:,:)   :: DTDT, DUDT, DVDT, DTDT_p, DUDT_p, DVDT_p
    real, pointer, dimension(:,:,:)   :: THEQ, T_EQ, THEQ_p, T_EQ_p
    real, pointer, dimension(:,:  )   :: TAUX, TAUY, DISS, TAUX_p, TAUY_p, DISS_p

! Scratch arrays and working pointers

    real, allocatable, dimension(:,:) :: PII
    real, allocatable, dimension(:,:) :: DP
    real, allocatable, dimension(:,:) :: PL
    real, allocatable, dimension(:,:) :: UU
    real, allocatable, dimension(:,:) :: VV
    real, allocatable, dimension(:,:) :: VR
    real, allocatable, dimension(:,:) :: TE
    real, allocatable, dimension(:,:) :: F1
    real, allocatable, dimension(:,:) :: RR
    real, allocatable, dimension(:,:) :: DS
    real, allocatable, dimension(:,:) :: DM
    real, allocatable, dimension(:,:) :: PK

    real,     pointer, dimension(:,:) :: PS
    real,     pointer, dimension(:,:) :: PT

    integer :: L
    integer :: IM, JM, LM
    integer :: FRICQ

    logical :: FriendlyTemp
    logical :: FriendlyWind

    real    :: DT
    real    :: KA  
    real    :: KS  
    real    :: KF

! 8 Held-Suarez parameters

    real    :: SIG1
    real    :: TAUA
    real    :: TAUS
    real    :: TAUF
    real    :: DELV1
    real    :: DELH
    real    :: T0
    real    :: TSTRT

! 3 Localized heating parameters

    real    :: P_1
    real    :: QMAX

! 3 Williamson parameters

    real    :: P_D
    real    :: GAM_I
    real    :: GAM_D
    
    real, parameter :: DAYLEN=86400.

    integer :: run_gtfv3 = 0
!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my MAPL_MetaComp
!---------------------

    call MAPL_GetObjectFromGC( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start the TOTAL timer
!----------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN")

! Get esmf internal state from generic state.
!-------------------------------------------

    call MAPL_Get(MAPL,                  &
         IM = IM, JM=JM, LM=LM,          &
         INTERNAL_ESMF_STATE=INTERNAL,   &
                               RC=STATUS )
    VERIFY_(STATUS)

! Pointers TO INTERNALS
!----------------------

    call MAPL_GetPointer(INTERNAL, SPHI2, 'SPHI2', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CPHI2, 'CPHI2', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, HFCN , 'HFCN' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, P_I  , 'P_I'  , RC=STATUS); VERIFY_(STATUS)

! Pointers TO IMPORTS
!-------------------- 

    call MAPL_GetPointer(IMPORT  , PLE  , 'PLE'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , T    , 'TEMP' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , U    , 'U'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT  , V    , 'V'    , RC=STATUS); VERIFY_(STATUS)

! Pointers TO EXPORTS
!-------------------- 

    call MAPL_GetPointer(EXPORT  , DUDT , 'DUDT' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DVDT , 'DVDT' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DTDT , 'DTDT' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , DISS , 'DISS' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , THEQ , 'THEQ' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , T_EQ , 'T_EQ' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TAUX , 'TAUX' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT  , TAUY , 'TAUY' , RC=STATUS); VERIFY_(STATUS)

! Get parameters from the configuration
!--------------------------------------

! 8 Held-Suarez parameters

!BOR

! !RESOURCE_ITEM: days :: H-S $k_a^{-1}$ parameter
    call MAPL_GetResource(MAPL,TAUA ,'TAUA:' ,default=40. ,     RC=STATUS )
    VERIFY_(STATUS)                  
! !RESOURCE_ITEM: days :: H-S $k_s^{-1}$ parameter
    call MAPL_GetResource(MAPL,TAUS ,'TAUS:' ,default=4.0 ,     RC=STATUS )
    VERIFY_(STATUS)                  
! !RESOURCE_ITEM: days :: H-S $k_f^{-1}$ parameter
    call MAPL_GetResource(MAPL,TAUF ,'TAUF:' ,default=1.0 ,     RC=STATUS )
    VERIFY_(STATUS)
! !RESOURCE_ITEM: days :: H-S $\sigma_b$ parameter
    call MAPL_GetResource(MAPL,SIG1 ,'SIG1:' ,default=0.7 ,     RC=STATUS )
    VERIFY_(STATUS)                  
! !RESOURCE_ITEM: K :: H-S $T_o$ parameter
    call MAPL_GetResource(MAPL,T0   ,'T0:'   ,default=315.,     RC=STATUS )
    VERIFY_(STATUS)                  
! !RESOURCE_ITEM: K :: H-S $(\Delta \theta)_z$ parameter
    call MAPL_GetResource(MAPL,DELV1,'DELV1:',default=10. ,     RC=STATUS )
    VERIFY_(STATUS)                  
! !RESOURCE_ITEM: K :: H-S $(\Delta T)_y$ parameter
    call MAPL_GetResource(MAPL,DELH ,'DELH:' ,default=60. ,     RC=STATUS )
    VERIFY_(STATUS)                  
! !RESOURCE_ITEM: K :: H-S $T_{strat}$ parameter
    call MAPL_GetResource(MAPL,TSTRT,'TSTRT:',default=200.,     RC=STATUS )
    VERIFY_(STATUS)
! !RESOURCE_ITEM: 1 OR 0 :: Controls addition of dissipation heat.
    call MAPL_GetResource(MAPL,FRICQ,'DISSQ:',default=1,        RC=STATUS )
    VERIFY_(STATUS)

! 3 Localized heating parameters

! !RESOURCE_ITEM: Pa :: Local heating $p_1$ parameter
    call MAPL_GetResource(mapl,P_1  ,'P1:'   ,default=0.0 ,     RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(mapl,QMAX ,'QMAX:' ,default=0.0 ,     RC=STATUS )
    VERIFY_(STATUS)                                        

! 3 Williamson parameters

! !RESOURCE_ITEM: K m-1 :: Williamson's $\gamma_I$ parameter
    call MAPL_GetResource(mapl,GAM_I,'GAM_I:',default=-3.345e-3,RC=STATUS)
    VERIFY_(STATUS)                                        
! !RESOURCE_ITEM: K m-1 :: Williamson's $\gamma_D$ parameter
    call MAPL_GetResource(mapl,GAM_D,'GAM_D:',default=.002     ,RC=STATUS)
    VERIFY_(STATUS)                                        
! !RESOURCE_ITEM: Pa :: Williamson's $p_D$ parameter
    call MAPL_GetResource(mapl,P_D  ,'P_D:'  ,default=1.E4     ,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, run_gtfv3, 'RUN_GTFV3:', default=0, RC=STATUS)
    VERIFY_(STATUS)
  
!EOR

! Check for Friendliness
!-----------------------

    Friendly = MAPL_VerifyFriendly(IMPORT,'V',trim(COMP_NAME),RC=STATUS)
    VERIFY_(STATUS)

    FriendlyWind = Friendly == ESMF_TRUE

    Friendly = MAPL_VerifyFriendly(IMPORT,'U',trim(COMP_NAME),RC=STATUS)
    VERIFY_(STATUS)

    ASSERT_((Friendly == ESMF_TRUE) .eqv. FriendlyWind)

    Friendly = MAPL_VerifyFriendly(IMPORT,'TEMP',trim(COMP_NAME),RC=STATUS)
    VERIFY_(STATUS)

    FriendlyTemp = Friendly == ESMF_TRUE

! If we have a friendly, we need a time step
!-------------------------------------------
    
    if(FriendlyTemp .or. FriendlyWind) then
       call MAPL_GetResource(MAPL, DT,  'RUN_DT', RC=STATUS)
       VERIFY_(STATUS)
    else
       DT = huge(0.0)
    end if

! Allocate 10 2D scratch arrays
!------------------------------

    allocate(DP (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(PL (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(UU (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(VV (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(VR (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(TE (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(DS (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(PII(IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(F1 (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(RR (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(DM (IM,JM),STAT=STATUS); VERIFY_(STATUS)
    allocate(PK (IM,JM),STAT=STATUS); VERIFY_(STATUS)

! Get MPI Rank
!-------------

  call ESMF_VMGetCurrent(vm, RC=status)
  call ESMF_VMGet(vm, mpiCommunicator=comm, RC=status)
  call MPI_Comm_rank(comm, rank, ierr)

  if (run_gtfv3 == 0) then
    ! Begin calculations.
    !-------------------

        PS  => PLE(:,:,LM)
        PT  => PLE(:,:, 0)

        PII =  P_D - (P_D - PT)*0.5*P_I

    ! Initialize vertically integrated diagnostics
    !---------------------------------------------

        if(associated(DISS)) DISS = 0.0
        if(associated(TAUX)) TAUX = 0.0
        if(associated(TAUY)) TAUY = 0.0

    ! Loop invariants
    !----------------

        KA   = 1.0/(DAYLEN*TAUA)
        KS   = 1.0/(DAYLEN*TAUS)
        KF   = 1.0/(DAYLEN*TAUF)

        LEVELS: do L = 1,LM

        DP  = (PLE(:,:,L)-PLE(:,:,L-1))
        PL  = (PLE(:,:,L)+PLE(:,:,L-1))*0.5
        DM  = DP / MAPL_GRAV
        PK  = (PL/MAPL_P00)**MAPL_KAPPA 

    ! H&S equilibrium temperature
    !----------------------------

        TE  = PK*( T0 - DELH*SPHI2 - DELV1*CPHI2*log( PL/MAPL_P00 ) )
        TE  = max( TE, TSTRT )

    ! Williamson Stratospheric modifications to equilibrium temperature
    ! -----------------------------------------------------------------

        where( PL < P_D )
            TE  = TSTRT*( min(1.0,PL/P_D)**(MAPL_RGAS*GAM_D/MAPL_GRAV)     &
                        + min(1.0,PL/PII)**(MAPL_RGAS*GAM_I/MAPL_GRAV) - 1 )
        end where

    !  Exports of equilibrium T and Theta
    !------------------------------------

        if(associated(T_EQ)) T_EQ(:,:,L) = TE
        if(associated(THEQ)) THEQ(:,:,L) = TE/PK

    ! Vertical structure of timescales in H&S.
    !---------------------------------------------

        F1  = max(0.0, ( (PL/PS)-SIG1 )/( 1.0-SIG1 ) )

    ! Atmospheric heating from H&S
    !-----------------------------

        RR = (KA + (KS-KA)*F1*CPHI2**2) * (TE-T(:,:,L))

        if(associated(DTDT)) DTDT(:,:,L) =            DP*RR
        if(FriendlyTemp    ) T   (:,:,L) = T(:,:,L) + DT*RR

    ! Wind tendencies
    !----------------

        UU  = -U(:,:,L)*(F1*KF)
        VV  = -V(:,:,L)*(F1*KF)

        if(associated(DUDT)) DUDT(:,:,L) = UU
        if(associated(DVDT)) DVDT(:,:,L) = VV

        if(FriendlyWind) then
            U(:,:,L) = U(:,:,L) + DT*UU
            V(:,:,L) = V(:,:,L) + DT*VV
        end if

    !  Frictional heating from H&S drag
    !----------------------------------

        DS = U(:,:,L)*UU + V(:,:,L)*VV

        if(associated(DISS)) DISS           = DISS        - DS*DM
        if(FRICQ /= 0) then
            if(associated(DTDT)) DTDT(:,:,L) = DTDT(:,:,L) - DS*(DP/MAPL_CP  )
            if(FriendlYTemp    ) T   (:,:,L) = T   (:,:,L) - DS*(DT/MAPL_CP  )
        end if

    !  Surface stresses from vertically integrated H&S surface drag
    !--------------------------------------------------------------

        if(associated(TAUX)) TAUX = TAUX - UU*DM
        if(associated(TAUY)) TAUY = TAUY - VV*DM

    ! Localized heat source, if any
    !------------------------------

        if((associated(DTDT).or.FriendlyTemp) .and. QMAX/=0.0) then
            where(PL > P_1)
                VR = HFCN*(QMAX/DAYLEN)*sin( MAPL_PI*(PS-PL)/(PS-P_1) )
            elsewhere
                VR = 0.
            end where

            if(associated(DTDT)) DTDT(:,:,L) = DTDT(:,:,L) + DP*VR 
            if(FriendlyTemp    ) T   (:,:,L) = T   (:,:,L) + DT*VR
        end if

        enddo LEVELS

        ! if(associated(DISS) .eqv. .FALSE.) then
        !     write(*,*) 'DISS is not associated'
        ! else
        !     DISS_p => DISS
        ! endif

        ! if(associated(TAUX) .eqv. .FALSE.) then
        !     write(*,*) 'TAUX is not associated'
        ! else
        !     TAUX_p => TAUX
        ! endif

        ! if(associated(TAUY) .eqv. .FALSE.) then
        !     write(*,*) 'TAUY is not associated'
        ! else
        !     TAUY_p => TAUY
        ! endif

        ! if(associated(DTDT) .eqv. .FALSE.) then
        !     write(*,*) 'DTDT is not associated'
        ! else
        !     DTDT_p => DTDT
        ! endif

        ! if(associated(DUDT) .eqv. .FALSE.) then
        !     write(*,*) 'DUDT is not associated'
        ! else
        !     DUDT_p => DUDT
        ! endif

        ! if(associated(DVDT) .eqv. .FALSE.) then
        !     write(*,*) 'DVDT is not associated'
        ! else
        !     DVDT_p => DVDT
        ! endif

        ! if(associated(THEQ) .eqv. .FALSE.) then
        !     write(*,*) 'THEQ is not associated'
        ! else
        !     THEQ_p => THEQ
        ! endif

        ! if(associated(T_EQ) .eqv. .FALSE.) then
        !     write(*,*) 'T_EQ is not associated'
        ! else
        !     T_EQ_p => T_EQ
        ! endif

        ! if(associated(CPHI2) .eqv. .FALSE.) write(*,*) 'CPHI2 is not associated'
        ! if(associated(HFCN) .eqv. .FALSE.) write(*,*) 'HFCN is not associated'
        ! if(associated(P_I) .eqv. .FALSE.) write(*,*) 'P_I is not associated'
        ! if(associated(PLE) .eqv. .FALSE.) write(*,*) 'PLE is not associated'
        ! if(associated(SPHI2) .eqv. .FALSE.) write(*,*) 'SPHI2 is not associated'
        ! if(associated(U) .eqv. .FALSE.) write(*,*) 'U is not associated'
        ! if(associated(V) .eqv. .FALSE.) write(*,*) 'V is not associated'
        ! if(associated(T) .eqv. .FALSE.) write(*,*) 'T is not associated'

        ! write(*,*) 'Calling c_call_hs_oacc'
        ! call c_call_hs_oacc(C_LOC(CPHI2), C_LOC(DISS_p), C_LOC(DTDT_p), C_LOC(DUDT_p), C_LOC(DVDT_p), &
        !                     C_LOC(HFCN), C_LOC(P_I), C_LOC(PLE), C_LOC(SPHI2), C_LOC(TAUX_p), C_LOC(TAUY_p), C_LOC(T), &
        !                     C_LOC(THEQ_p), C_LOC(T_EQ), C_LOC(U), C_LOC(V), &
        !                     DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
        !                     FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
        !                     SIG1, TAUA, TAUF, TAUS, TSTRT, T0, 0)
    else
        call MPI_GET_PROCESSOR_NAME(nameNode, length, ierr)
        !write(*,*) 'Rank = ', rank, 'on node', nameNode
        DISS_P => null()
        TAUX_P => null()
        TAUY_P => null()
        THEQ_P => null()
        call held_suarez_oacc(CPHI2, DISS_P, DTDT, DUDT, DVDT, &
                            HFCN, P_I, PLE, SPHI2, TAUX_P, TAUY_P, T, &
                            THEQ_P, T_EQ, U, V, &
                            DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, &
                            FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, &
                            SIG1, TAUA, TAUF, TAUS, TSTRT, T0, 2, rank)
    end if
! Free 10 scratch arrays
!-----------------------

    deallocate(PK ) 
    deallocate(DM ) 
    deallocate(RR ) 
    deallocate(F1 ) 
    deallocate(PII) 
    deallocate(DS ) 
    deallocate(TE ) 
    deallocate(VR ) 
    deallocate(VV ) 
    deallocate(UU )
    deallocate(PL ) 
    deallocate(DP ) 

! Close timers
!-------------

    call MAPL_TimerOff(MAPL,"RUN")
    call MAPL_TimerOff(MAPL,"TOTAL")

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine RUN

end module GEOS_HSGridCompMod



