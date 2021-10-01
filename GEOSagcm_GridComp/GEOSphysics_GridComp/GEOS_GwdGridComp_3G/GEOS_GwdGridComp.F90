#include "MAPL_Generic.h"

module GEOS_GwdGridComp

!BOP

! !MODULE: GEOS_Gwd -- A Module to compute the forcing due to parameterized gravity wave drag

! !DESCRIPTION:
! 
!   {\tt GWD} is a light-weight gridded component to compute the forcing
! due to gravity wave drags. It operates on the ESMF grid that appears in the
! gridded component passed to its {\tt Initialize} method. Unlike
! heavier gridded components, it does not enforce its own grid.
! The only restrictions are that it be a 3-dimensional grid
! in which one dimension is aligned with the vertical coordinate and
! only the horizontal dimensions are decomposed.
!
! The gravity wave drag scheme is based on NCAR WACCM1b gw\_drag routine.
! The scheme includes parameterizations for orographic (stationary) gravity
! waves (Kiehl et al. 1996), and for a spectrum of traveling gravity waves 
!(Sassi et al. 2003; http://acd.ucar.edu/models/WACCM). Both parameteriz-
! ations are based on Lindzen's [1981] formulation. The interested reader 
! is referred to those publications for details of the mathematical
! derivations.
!

! !USES:

  use ESMF
  use MAPL

  use gw_drag, only: gw_intr
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

contains

!BOP
! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:
  subroutine SetServices ( gc, rc )

! !ARGUMENTS:
    type(ESMF_GridComp), intent(inout) :: gc  ! gridded component
    integer, intent(out)               :: rc  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (gc). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: status
    character(len=ESMF_MAXSTR)              :: comp_name

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( gc, NAME=comp_name, _RC)
    Iam = trim(comp_name) // Iam

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_INITIALIZE,  Initialize, _RC)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,  Run, _RC)

! Set the state variable specs.
! -----------------------------
!BOS
#include "DU2G_Export___.h"
#include "DU2G_Import___.h"
!EOS

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, _RC)

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOP

  ! !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component

  ! !INTERFACE:

  subroutine Initialize ( gc, import, export, clock, rc )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component 
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION: The Initialize method of the GWD Physics Gridded Component first 
    !   calls the Initialize method of the children.  Then, if using the NCAR GWD
    !   scheme, calls the initialization routines.

    !EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: status
    character(len=ESMF_MAXSTR)              :: comp_name

! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL

!=============================================================================

   ! Begin...

   ! Get my name and set-up traceback handle
   ! ---------------------------------------

      Iam = 'Initialize'
      call ESMF_GridCompGet( gc, NAME=comp_name, _RC)
      Iam = trim(comp_name) // Iam

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC ( gc, MAPL, _RC)

      ! Call Generic Initialize for GWD gc
      !-----------------------------------

      call MAPL_GenericInitialize ( gc, import, export, clock, _RC)

      ! All done
      !---------

      RETURN_(ESMF_SUCCESS)
   end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: RUN -- Run method for the GWD component

! !INTERFACE:
subroutine RUN ( gc, import, export, clock, rc )

! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component 
  type(ESMF_State),    intent(inout) :: import ! Import state
  type(ESMF_State),    intent(inout) :: export ! Export state
  type(ESMF_Clock),    intent(inout) :: clock  ! The clock
  integer, optional,   intent(  out) :: rc     ! Error code:

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: status
  character(len=ESMF_MAXSTR)          :: comp_name

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: mapl
  type (ESMF_Alarm       )            :: alarm

  integer                             :: im, jm, lm
  integer                             :: pgwv
  real                                :: effgworo, effgwbkg
  real                                :: CDMBGWD1, CDMBGWD2
  real                                :: bgstressmax
  real, pointer, dimension(:,:)       :: lats

! Rayleigh friction parameters

  REAL                                :: H0, HH, Z1, TAU1

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run"
   call ESMF_GridCompGet( gc, name=comp_name, _RC)
   Iam = trim(comp_name) // Iam

! Retrieve the pointer to the state
! -------------------------------
   call MAPL_GetObjectFromGC ( gc, MAPL, _RC)

! Retrieve configuration resource
! -------------------------------
   config = MAPL%get_config(_RC)

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL, &
         im=im, jm=jm, lm=lm,        &
         runalarm=alarm, lats=lats, _RC)

! Gravity wave drag
! -----------------

    call config%get(effgworo,    "EFFGWORO",    default=0.250, _RC)
    call config%get(effgwbkg,    "EFFGWBKG",    default=0.125, _RC)
    call config%get(pgwv,        "PGWV",        default=nint(4*lm/72.0), _RC)
    call config%get(bgstressmax, "BGSTRESSMAX", default=0.9,  _RC)

! Rayleigh friction
! -----------------
    
    call config%get(Z1,   "RAYLEIGH_Z1:",   default=75000.,  _RC)
    call config%get(TAU1, "RAYLEIGH_TAU1:", default=172800., _RC)
    call config%get(H0,   "RAYLEIGH_H0:",   default=7000.,   _RC)
    call config%get(HH,   "RAYLEIGH_HH:",   default=7500.,   _RC)

! If its time, recalculate the GWD tendency
! -----------------------------------------

   if ( ESMF_AlarmIsRinging( ALARM ) ) then
      call ESMF_AlarmRingerOff(ALARM, _RC)
      call state%t_profiler%start('driver')
      call Gwd_Driver(_RC)
      call state%t_profiler%stop('driver')
   endif

   _RETURN(ESMF_SUCCESS)

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Gwd_Driver(rc)
      integer, optional, intent(OUT) :: rc

!  Locals

      character(len=ESMF_MAXSTR)      :: IAm
      integer                         :: status

      type (ESMF_TimeInterval)        :: TINT

!  Pointers from Import state

      real, pointer, dimension(:)      :: PREF
      real, pointer, dimension(:,:)    :: SGH
      real, pointer, dimension(:,:,:)  :: PLE, T, Q, U, V
      !++jtb Array for moist/deepconv heating
      real, pointer, dimension(:,:,:)  :: HT_dpc

!  Pointers to Export state

      real, pointer, dimension(:)      :: PREF_EXP
      real, pointer, dimension(:,:)    :: SGH_EXP
      real, pointer, dimension(:,:,:)  :: PLE_EXP, T_EXP, Q_EXP, U_EXP, V_EXP

      real, pointer, dimension(:,:)    :: CLDSTD
      real, pointer, dimension(:,:)    :: UBAR,    VBAR
      real, pointer, dimension(:,:)    :: UBASE,   VBASE
      real, pointer, dimension(:,:)    :: TAUGWX,  TAUGWY
      real, pointer, dimension(:,:)    :: TAUOROX, TAUOROY
      real, pointer, dimension(:,:)    :: TAUBKGX, TAUBKGY
      real, pointer, dimension(:,:,:)  :: TAUOROXZ,TAUOROYZ,FEOROZ,FEPOROZ
      real, pointer, dimension(:,:,:)  :: TAUBKGXZ,TAUBKGYZ,FEBKGZ,FEPBKGZ
      real, pointer, dimension(:,:)    :: TAUOROXT,TAUOROYT,FEOROT,FEPOROT
      real, pointer, dimension(:,:)    :: TAUOROXS,TAUOROYS,FEOROS,FEPOROS
      real, pointer, dimension(:,:)    :: TAUBKGXT,TAUBKGYT,FEBKGT,FEPBKGT
      real, pointer, dimension(:,:)    :: TAUBKGXS,TAUBKGYS,FEBKGS,FEPBKGS
      real, pointer, dimension(:,:)    :: TAUMSTX, TAUMSTY
      real, pointer, dimension(:,:)    :: KEGWD, KEORO,  KERAY,  KEBKG, KERES
      real, pointer, dimension(:,:)    :: PEGWD, PEORO,  PERAY,  PEBKG, BKGERR

      real, pointer, dimension(:,:,:)  :: DTDT, DUDT, DVDT, TTMGW
      real, pointer, dimension(:,:,:)  :: DTDT_ORO, DUDT_ORO, DVDT_ORO
      real, pointer, dimension(:,:,:)  :: DTDT_BKG, DUDT_BKG, DVDT_BKG
      real, pointer, dimension(:,:,:)  :: DTDT_RAY, DUDT_RAY, DVDT_RAY
      real, pointer, dimension(:,:,:)  :: DTGENBKG, DUGENBKG, DVGENBKG
      
! local variables

      real,              dimension(im,jm,lm  ) :: ZM, PMID, PDEL, RPDEL, PMLN
      real,              dimension(im,jm,lm  ) :: DUDT_ORG, DVDT_ORG, DTDT_ORG
      real,              dimension(im,jm,lm  ) :: DUDT_GWD, DVDT_GWD, DTDT_GWD
      real,              dimension(im,jm,lm  ) :: DUDT_RAH, DVDT_RAH, DTDT_RAH
      real,              dimension(im,jm,lm  ) :: DUDT_TOT, DVDT_TOT, DTDT_TOT
      real,              dimension(im,jm,lm+1) :: PILN,   ZI
      real,              dimension(      lm  ) :: ZREF, KRAY
      real,              dimension(im,jm     ) :: TAUXO_TMP, TAUYO_TMP
      real,              dimension(im,jm     ) :: TAUXB_TMP, TAUYB_TMP
      real,              dimension(im,jm,lm+1) :: TAUXO_3D , TAUYO_3D , FEO_3D, FEPO_3D
      real,              dimension(im,jm,lm+1) :: TAUXB_3D , TAUYB_3D , FEB_3D, FEPB_3D
      real,              dimension(im,jm,lm  ) :: DUBKGSrc , DVBKGSrc , DTBKGSrc
      real,              dimension(im,jm)      :: KEGWD_X, KEORO_X,  KERAY_X,  KEBKG_X, KERES_X
      real,              dimension(im,jm)      :: PEGWD_X, PEORO_X,  PERAY_X,  PEBKG_X, BKGERR_X

      integer                                  :: J, K, L
      real(ESMF_KIND_R8)                       :: DT_R8
      real                                     :: DT     ! time interval in sec

!  Begin...
!----------

      IAm = "Gwd_Driver"

! Get time step
!-------------------------------------------------

      call ESMF_AlarmGet( ALARM, ringInterval=TINT,    _RC)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8,      _RC)

      DT = DT_R8

! Pointers to import/export quantities
!-------------------------------------
#include "Gwd_DeclarePointer___.h"


      CALL PREGEO(im*jm,   lm,   &
                    PLE, LATS,   PMID,  PDEL, RPDEL,     PILN,     PMLN)

! Compute ZM
!-------------

      call GEOPOTENTIAL( im*jm, lm,                  &
           PILN,  PMLN,  PLE,   PMID, PDEL, RPDEL,   &
           T,     Q,     ZI,    ZM                   )

! Do gravity wave drag calculations on a list of soundings
!---------------------------------------------------------

      call state%t_profiler%start('INTR')

          ! Use GEOS GWD    
          call gw_intr   (im*jm,      lm,         DT,                  &
               PGWV,                                                   &
               PLE,       T,          U,          V,      SGH,   PREF, &
               PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
               DUDT_GWD,  DVDT_GWD,   DTDT_GWD,                        &
               DUDT_ORG,  DVDT_ORG,   DTDT_ORG,                        &
               TAUXO_TMP, TAUYO_TMP,  TAUXO_3D,   TAUYO_3D,  FEO_3D,   &
               TAUXB_TMP, TAUYB_TMP,  TAUXB_3D,   TAUYB_3D,  FEB_3D,   &
               FEPO_3D,   FEPB_3D,    DUBKGSrc,   DVBKGSrc,  DTBKGSrc, &
               BGSTRESSMAX, effgworo, effgwbkg,   _RC            )

    call state%t_profiler%stop('INTR')

    CALL POSTINTR(im*jm, lm, DT, H0, HH, Z1, TAU1, &
          PREF,     &
          PDEL,     &
          U,        &
          V,        &
          DUDT_GWD, &
          DVDT_GWD, &
          DTDT_GWD, &
          DUDT_ORG, &
          DVDT_ORG, &
          DTDT_ORG, &

          DUDT_TOT, &
          DVDT_TOT, &
          DTDT_TOT, &
          DUDT_RAH, &
          DVDT_RAH, &
          DTDT_RAH, &
          PEGWD_X,  &
          PEORO_X,  &
          PERAY_X,  &
          PEBKG_X,  &
          KEGWD_X,  &
          KEORO_X,  &
          KERAY_X,  &
          KEBKG_X,  &
          KERES_X,  &
          BKGERR_X  )

!! Tendency diagnostics
!!---------------------
    call copy_if_associated(DUDT, from=DUDT_TOT)
    call copy_if_associated(DVDT, from=DVDT_TOT)
    call copy_if_associated(DTDT, from=DTDT_TOT)

    call copy_if_associated(DUDT_RAY, from=DUDT_RAH)
    call copy_if_associated(DVDT_RAY, from=DVDT_RAH)
    call copy_if_associated(DTDT_RAY, from=DTDT_RAH)

!! KE dIagnostics
!!----------------

    call copy_if_associated(PEGWD,  from=PEGWD_X)
    call copy_if_associated(PEORO,  from=PEORO_X)
    call copy_if_associated(PERAY,  from=PERAY_X)
    call copy_if_associated(PEBKG,  from=PEBKG_X)
    call copy_if_associated(KEGWD,  from=KEGWD_X)
    call copy_if_associated(KEORO,  from=KEORO_X)
    call copy_if_associated(KERAY,  from=KERAY_X)
    call copy_if_associated(KEBKG,  from=KEBKG_X)
    call copy_if_associated(KERES,  from=KERES_X)
    call copy_if_associated(BKGERR, from=BKGERR_X)

!! Tendency diagnostics
!!---------------------

    call copy_if_associated(DUDT_ORO, from=DUDT_ORG)
    call copy_if_associated(DVDT_ORO, from=DVDT_ORG)
    call copy_if_associated(DTDT_ORO, from=DTDT_ORG)

! Orographic stress
!------------------

    call copy_if_associated(TAUGWX,  TAUXO_TMP + TAUXB_TMP)
    call copy_if_associated(TAUGWY,  TAUYO_TMP + TAUYB_TMP)
    call copy_if_associated(TAUOROX, TAUXO_TMP)
    call copy_if_associated(TAUOROY, TAUYO_TMP)
    call copy_if_associated(TAUBKGX, TAUXB_TMP)
    call copy_if_associated(TAUBKGY, TAUYB_TMP)

! Export unweighted T Tendency
!-----------------------------

    if(associated(TTMGW)) then
       if(associated(DTDT )) then
          TTMGW = DTDT
       else
          TTMGW = 0.0
       end if
    end if

! AMM modify T_EXP to be the T AFTER GWD, ie., add the tendency*dt
! (need to do this before DTDT is pressure weighted for the dynamics)
    if(associated(T_EXP   )) T_EXP    = T + DTDT*DT

! DTDT has to be pressure weighted and is all due to frictional heating.
!-----------------------------------------------------------------------

    if(associated(DTDT    )) then
       DTDT = DTDT*PDEL 
    end if

    if(associated(PREF_EXP)) PREF_EXP = PREF
    if(associated(SGH_EXP )) SGH_EXP  = SGH
    if(associated(PLE_EXP )) PLE_EXP  = PLE
    if(associated(Q_EXP   )) Q_EXP    = Q
    if(associated(U_EXP   )) U_EXP    = U
    if(associated(V_EXP   )) V_EXP    = V

! All done
!-----------

    _RETURN(ESMF_SUCCESS)
   end subroutine GWD_DRIVER

  end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine geopotential(pcols  , pver   ,                   &
         piln   , pmln   , pint  , pmid   , pdel   , rpdel  , &
         t      , q      , zi     , zm     )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!
!-----------------------------------------------------------------------

    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: pcols                ! Number of longitudes
    integer, intent(in) :: pver                 ! Number of vertical layers

    real,    intent(in) :: piln (pcols,pver+1)  ! Log interface pressures
    real,    intent(in) :: pmln (pcols,pver)    ! Log midpoint pressures
    real,    intent(in) :: pint (pcols,pver+1)  ! Interface pressures
    real,    intent(in) :: pmid (pcols,pver)    ! Midpoint pressures
    real,    intent(in) :: pdel (pcols,pver)    ! layer thickness
    real,    intent(in) :: rpdel(pcols,pver)    ! inverse of layer thickness
    real,    intent(in) :: t    (pcols,pver)    ! temperature
    real,    intent(in) :: q    (pcols,pver)    ! specific humidity

! Output arguments

    real,    intent(out) :: zi(pcols,pver+1)    ! Height above surface at interfaces
    real,    intent(out) :: zm(pcols,pver)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    logical  :: fvdyn              ! finite volume dynamics
    integer  :: i,k                ! Lon, level indices
    real     :: hkk                ! diagonal element of hydrostatic matrix
    real     :: hkl                ! off-diagonal element
    real     :: tv                 ! virtual temperature
    real     :: tvfac              ! Tv/T

    real, parameter :: ROG     = MAPL_RGAS/MAPL_GRAV
!
!-----------------------------------------------------------------------
!

! Set dynamics flag

    fvdyn = .true.

! The surface height is zero by definition.

    I_LOOP: do i = 1, pcols
       zi(i,pver+1) = 0.0

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k)

       do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics

          if (fvdyn) then
             hkl = piln(i,k+1) - piln(i,k)
             hkk = piln(i,k+1) - pmln(i,k)
          else
             hkl = pdel(i,k) / pmid(i,k)
             hkk = 0.5 * hkl
          end if

! Now compute tv, zm, zi

          tvfac   = 1. + MAPL_VIREPS * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + ROG * tv * hkk
          zi(i,k) = zi(i,k+1) + ROG * tv * hkl
       end do
    end do I_LOOP

    return
  end subroutine geopotential

!----------------------------------------------------------------------- 

  subroutine pregeo(pcols,pver,&
    ple,lats,pmid,pdel,rpdel,piln,pmln)

    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!

    integer, intent(in) :: pcols    ! Number of longitudes
    integer, intent(in) :: pver     ! Number of vertical layers

    real,    intent(in) :: ple (pcols,pver+1)    ! Interface pressures
    real,    intent(in) :: lats(pcols)           ! latitude in radian

! Output arguments
    
    real,    intent(out) :: pmid  (pcols,pver)   ! Midpoint pressures
    real,    intent(out) :: pdel  (pcols,pver)   ! layer thickness
    real,    intent(out) :: rpdel (pcols,pver)   ! inverse of layer thickness
    real,    intent(out) :: piln  (pcols,pver+1) ! Log interface pressures
    real,    intent(out) :: pmln  (pcols,pver)   ! Log midpoint pressures

!
!---------------------------Local variables-----------------------------
!
    integer :: i,k

    real    :: hvsd  ! Efficiency factor

    real, parameter :: PI_GWD  = 4.0*atan(1.0)  ! This is *not* MAPL_PI

!
!-----------------------------------------------------------------------
!

! Form pressure factors
!----------------------

    I_LOOP: DO I = 1, PCOLS
       DO K = 1, PVER
           PMID(I,K) = 0.5*(  PLE(I,K  ) + PLE(I,K+1) )
           PDEL(I,K) =        PLE(I,K+1) - PLE(I,K  )
          RPDEL(I,K) = 1.0 / PDEL(I,K)
           PILN(I,K) = log(   PLE(I,K) )
           PMLN(I,K) = log(  PMID(I,K) ) !
       END DO
       PILN(I,PVER+1)  = log( PLE(I,PVER+1)  )
    END DO I_LOOP

  end subroutine pregeo

  subroutine postintr(pcols,pver,dt, h0, hh, z1, tau1, &
        pref, &
        pdel, &
        u, &
        v, &
        dudt_gwd, &
        dvdt_gwd, &
        dtdt_gwd, &
        dudt_org, &
        dvdt_org, &
        dtdt_org, &

        ! Outputs
        dudt_tot, &
        dvdt_tot, &
        dtdt_tot, &
        dudt_rah, &
        dvdt_rah, &
        dtdt_rah, &
        pegwd, &
        peoro, &
        peray, &
        pebkg, &
        kegwd, &
        keoro, &
        keray, &
        kebkg, &
        keres, &
        bkgerr )
    
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!

    integer, intent(in) :: PCOLS ! Number of longitudes
    integer, intent(in) :: PVER  ! Number of vertical layers
    real,    intent(in) :: DT    ! Time step
    real,    intent(in) :: H0, HH, Z1, TAU1 ! Rayleigh friction parameters

    real,    intent(in) :: PREF(PVER+1)
    real,    intent(in) :: PDEL(PCOLS,PVER)
    real,    intent(in) :: U(PCOLS,PVER)
    real,    intent(in) :: V(PCOLS,PVER)

    real,    intent(in) :: DUDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DVDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DTDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DUDT_ORG(PCOLS,PVER)
    real,    intent(in) :: DVDT_ORG(PCOLS,PVER)
    real,    intent(in) :: DTDT_ORG(PCOLS,PVER)

    real,    intent(out) :: DUDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DVDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DTDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DUDT_RAH(PCOLS,PVER)
    real,    intent(out) :: DVDT_RAH(PCOLS,PVER)
    real,    intent(out) :: DTDT_RAH(PCOLS,PVER)
    real,    intent(out) :: PEGWD(PCOLS)
    real,    intent(out) :: PEORO(PCOLS)
    real,    intent(out) :: PERAY(PCOLS)
    real,    intent(out) :: PEBKG(PCOLS)
    real,    intent(out) :: KEGWD(PCOLS)
    real,    intent(out) :: KEORO(PCOLS)
    real,    intent(out) :: KERAY(PCOLS)
    real,    intent(out) :: KEBKG(PCOLS)
    real,    intent(out) :: KERES(PCOLS)
    real,    intent(out) :: BKGERR(PCOLS)

!
!---------------------------Local variables-----------------------------
!
    integer :: i,k
    real :: zref, kray
!
!-----------------------------------------------------------------------
!

    I_LOOP: DO I = 1, PCOLS
       PEGWD(I)  = 0.0
       PEORO(I)  = 0.0
       PERAY(I)  = 0.0
       PEBKG(I)  = 0.0
       KEGWD(I)  = 0.0
       KEORO(I)  = 0.0
       KERAY(I)  = 0.0
       KEBKG(I)  = 0.0
       KERES(I)  = 0.0
       BKGERR(I) = 0.0

       DO K = 1, PVER 

! Rayleigh friction
!------------------

          ZREF     = H0 * LOG(MAPL_P00/(0.5*(PREF(K)+PREF(K+1))))
          KRAY     = (1.0/TAU1)*( 1.0 - TANH( (Z1-ZREF)/HH ) )
          KRAY     = KRAY/(1+DT*KRAY)

          DUDT_RAH(I,K) = -U(I,K)*KRAY
          DVDT_RAH(I,K) = -V(I,K)*KRAY

          DTDT_RAH(I,K) = - ((U(I,K) + (0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) + &
                             (V(I,K) + (0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K)   ) * (1.0/MAPL_CP)

          DUDT_TOT(I,K) = DUDT_RAH(I,K) + DUDT_GWD(I,K)
          DVDT_TOT(I,K) = DVDT_RAH(I,K) + DVDT_GWD(I,K)
          DTDT_TOT(I,K) = DTDT_RAH(I,K) + DTDT_GWD(I,K)

! KE dIagnostics
!----------------

          PEGWD(I) = PEGWD(I) +  DTDT_TOT(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
          PEORO(I) = PEORO(I) +  DTDT_ORG(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
          PERAY(I) = PERAY(I) +  DTDT_RAH(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
          PEBKG(I) = PEBKG(I) + (DTDT_GWD(I,K)-DTDT_ORG(I,K))*PDEL(I,K)*(MAPL_CP/MAPL_GRAV)

          KEGWD(I) = KEGWD(I) + ((U(I,K)+(0.5*DT)*DUDT_TOT(I,K))*DUDT_TOT(I,K) +   &
                                 (V(I,K)+(0.5*DT)*DVDT_TOT(I,K))*DVDT_TOT(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)

          KEORO(I) = KEORO(I) + ((U(I,K)+(0.5*DT)*DUDT_ORG(I,K))*DUDT_ORG(I,K) +   &
                                 (V(I,K)+(0.5*DT)*DVDT_ORG(I,K))*DVDT_ORG(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)

          KERAY(I) = KERAY(I) + ((U(I,K)+(0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) +   &
                                 (V(I,K)+(0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)

          KEBKG(I) = KEBKG(I) + ((U(I,K)+(0.5*DT)*(DUDT_GWD(I,K) - DUDT_ORG(I,K)))*(DUDT_GWD(I,K) - DUDT_ORG(I,K)) +     &
                                 (V(I,K)+(0.5*DT)*(DVDT_GWD(I,K) - DVDT_ORG(I,K)))*(DVDT_GWD(I,K) - DVDT_ORG(I,K))   ) * &
                                  PDEL(I,K)*(1.0/MAPL_GRAV)
       END DO

       BKGERR(I) = -( PEBKG(I) + KEBKG(I) )
       KERES(I)  =    PEGWD(I) + KEGWD(I) + BKGERR(I)

    END DO I_LOOP

  end subroutine postintr

end module GEOS_GwdGridComp
