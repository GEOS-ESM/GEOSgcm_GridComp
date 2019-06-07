
#define _REALLY_UPDATE_
#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: DynVec -- A Module to dynamics states vector

! !INTERFACE:

module DynVec_GridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
  use fv_statemod, only: INTERP_AGRID_TO_DGRID
  use fv_statemod, only: fv_getpkz
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!

!EOP

! Common block (for now)
  integer,  parameter :: r8 = 8
  integer,  parameter :: r4 = 4

contains

!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

  subroutine SetServices ( GC, RC )

! ! ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    type (MAPL_MetaComp),         pointer   :: MAPL

    integer DYNVEC

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Set the Run entry point
! -----------------------
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   Run1, rc=status)
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

! !IMPORT STATE:

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DUDT',                                      &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DVDT',                                      &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DTDT',                                      &
         LONG_NAME  = 'temperature_analysis_increment',            &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DPEDT',                                     &
         LONG_NAME  = 'edge_pressure_analysis_increment',          &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'QTOT',                                      &
         LONG_NAME  = 'total_specific_humidity',                   &
         UNITS      = 'kg kg-1',                                   &
!        default    = 1.0e-6,                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

! !INTERNAL STATE:

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'PKZ',                                       &
         LONG_NAME  = 'pressure_to_kappa',                         &
         UNITS      = 'Pa$^\kappa$',                               &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'DZ',                                        &
         LONG_NAME  = 'height_thickness',                          &
         UNITS      = 'm',                                         &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

    call MAPL_AddInternalSpec ( gc,                                &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                     &
         PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = 'Pa',                                        &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PKZ',                                       &
         LONG_NAME  = 'pressure_to_kappa',                         &
         UNITS      = 'Pa$^\kappa$',                               &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DZ',                                        &
         LONG_NAME  = 'height_thickness',                          &
         UNITS      = 'm',                                         &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                     &
!        PRECISION  = ESMF_KIND_R8,                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( gc, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! ! IROUTINE: RUN -- Run method for DYNVEC component

! !INTERFACE:

subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: All this component does is return the state of dynamics
!                as read from a restart. This should really live in the
!                dynamics.

!EOP


! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

  type (MAPL_MetaComp),     pointer   :: MAPL
  type (ESMF_State)                   :: INTERNAL

  real(r4), parameter :: KAPPA        = MAPL_KAPPA

  integer                             :: IM, JM, LM
  integer, pointer, dimension(:)      :: iheader

  real(r8), pointer, dimension(:)         :: ak
  real(r8), pointer, dimension(:)         :: bk
  real(r8), pointer, dimension(:,:,:)     :: bkg_u
  real(r8), pointer, dimension(:,:,:)     :: bkg_v
  real(r8), pointer, dimension(:,:,:)     :: bkg_pt
  real(r8), pointer, dimension(:,:,:)     :: bkg_pe
  real(r8), pointer, dimension(:,:,:)     :: bkg_pkz
  real(r8), pointer, dimension(:,:,:)     :: bkg_dz
  real(r8), pointer, dimension(:,:,:)     :: bkg_w

  real(r8), pointer, dimension(:,:,:)     :: bkg_dp
  real(r8), pointer, dimension(:,:,:)     :: bkg_pl

  real(r4), pointer, dimension(:)         :: ak4
  real(r4), pointer, dimension(:)         :: bk4
  real(r4), pointer, dimension(:,:,:)     :: ana_u
  real(r4), pointer, dimension(:,:,:)     :: ana_v
  real(r4), pointer, dimension(:,:,:)     :: ana_pt
  real(r4), pointer, dimension(:,:,:)     :: ana_pe
  real(r4), pointer, dimension(:,:,:)     :: ana_pkz
  real(r4), pointer, dimension(:,:,:)     :: ana_dz
  real(r4), pointer, dimension(:,:,:)     :: ana_w

  real(r8), allocatable, dimension(:,:)   :: dps
  real(r8), allocatable, dimension(:,:)   :: sumq
  real(r8), allocatable, dimension(:,:)   :: lnbkgpe
  real(r8), allocatable, dimension(:,:,:) :: ana_td
  real(r8), allocatable, dimension(:,:,:) :: pkxy,dpek,dpkz
  real(r8), allocatable, dimension(:,:,:) :: dua,dva,dud,dvd

  real, pointer, dimension(:,:,:)     :: du
  real, pointer, dimension(:,:,:)     :: dv
  real, pointer, dimension(:,:,:)     :: dt
  real, pointer, dimension(:,:,:)     :: qtot_bkg
  real, pointer, dimension(:,:,:)     :: dpe

  real, allocatable, dimension(:,:)   :: area
  real, allocatable, dimension(:,:)   :: qint

  character(len=ESMF_MAXSTR)          :: areafn
  character(len=ESMF_MAXSTR)          :: incremental
  type (ESMF_Config)                  :: CF
  type(ESMF_Grid)                     :: grid
  real                                :: sum_pt
  real(r8)  bkg_pdry
  integer L, ll,lu, approach

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

   Iam = "Run1"
   call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // Iam

   if ( MAPL_AM_I_ROOT() ) then
       print *, 'Now running ',trim(Iam)
   endif

! Retrieve the pointer to the state
!----------------------------------

   call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
   VERIFY_(STATUS)

! Local aliases to the state, grid, and configuration
! ---------------------------------------------------

!  call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,    &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                                       RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GC, grid=grid, rc=status)
    VERIFY_(STATUS)

! **********************************************************************
! ****               Get Pointers to BKG Import Data                ****
! **********************************************************************
#if 0
    if ( MAPL_AM_I_ROOT() ) then
       call ESMF_StatePrint(IMPORT)
    end if
#endif

!   Get pointers to import variables
!   --------------------------------
    call MAPL_GetPointer(import,       du, 'DUDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,       dv, 'DVDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,       dt, 'DTDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import,      dpe,'DPEDT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(import, qtot_bkg, 'QTOT',  RC=STATUS)
    VERIFY_(STATUS)

!   Get pointers to internal variables
!   ----------------------------------
    call MAPL_GetPointer(internal,        ak, 'AK',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,        bk, 'BK',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,     bkg_u,  'U',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,     bkg_v,  'V',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,    bkg_pt, 'PT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,    bkg_pe, 'PE', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,   bkg_pkz,'PKZ', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,    bkg_dz, 'DZ', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(internal,     bkg_w,  'W', RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, incremental, label ='INCREMENTAL:', default='NO', rc=STATUS )
    VERIFY_(STATUS)

#ifdef _REALLY_UPDATE_

!   Calculate background dry mass
!   -----------------------------
    allocate( bkg_pl(IM,JM,LM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( bkg_dp(IM,JM,LM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( sumq(IM,JM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( qint(IM,JM),STAT=STATUS )
    VERIFY_(STATUS)
    allocate( area(IM,JM),STAT=STATUS )
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, areafn, label ='AREA_FILENAME:', default='NONE', rc=STATUS )
    VERIFY_(STATUS)

    lu = GETFILE( trim(areafn), form="unformatted", RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_VarRead(lu, grid, area, RC=STATUS)
    VERIFY_(STATUS)
    call FREE_FILE(lu, RC=STATUS)
    VERIFY_(STATUS)

!   bkg_pl  = 0.5*(bkg_pe(:,:,1:LM)+bkg_pe(:,:,0:LM-1))
!   bkg_dp  =      bkg_pe(:,:,1:LM)-bkg_pe(:,:,0:LM-1)
    bkg_pl  = 0.5*(bkg_pe(:,:,2:LM+1)+bkg_pe(:,:,1:LM))
    bkg_dp  =      bkg_pe(:,:,2:LM+1)-bkg_pe(:,:,1:LM)

    sumq = 0.0_r8
    do L=1,LM
       sumq = sumq + qtot_bkg(:,:,L)*bkg_dp(:,:,L)
    enddo
    qint = bkg_pe(:,:,LM)-sumq
    bkg_pdry = 0.0_r8
    call MAPL_AreaMean( bkg_pdry, qint, area, grid, rc=STATUS ); VERIFY_(STATUS)
    if ( MAPL_AM_I_ROOT() ) then
       print *, 'Background Total Dry Mass: ', bkg_pdry
    endif

    deallocate( area )
    deallocate( qint )
    deallocate( sumq )
    deallocate( bkg_dp )
    deallocate( bkg_pl )

!   Handle update to p-to-the-kappa
!   -------------------------------
    allocate(pkxy(IM,JM,LM+1))

    if ( trim(incremental) == 'YES' .or. trim(incremental) == 'yes' ) then

       approach = 2
       if ( approach == 1 ) then

!        Compute incremental change to p-to-the-kappa
!        --------------------------------------------
         allocate(dpek(IM,JM,LM+1))
         allocate(dpkz(IM,JM,LM))
         allocate(lnbkgpe(IM,JM))
         pkxy = bkg_pe**kappa
         dpek = kappa*pkxy*dpe/bkg_pe ! kappa*(bkg_pe)**(kappa-1)*dpe
         do L=1,LM ! note: bkg_pe ranges from 0:LM
            lnbkgpe = log(bkg_pe(:,:,L)/bkg_pe(:,:,L-1)) 
            dpkz(:,:,L) = (    ( dpek (:,:,L+1) -   dpek(:,:,L) )*lnbkgpe &
                          -  (   pkxy (:,:,L+1) -  pkxy (:,:,L) )* &
                             (   dpek (:,:,L+1) * bkg_pe(:,:,L-1)  &
                             -   dpek (:,:,L)   * bkg_pe(:,:,L  )) &
                             / (bkg_pe(:,:,L)   * bkg_pe(:,:,L-1)) &
                          )  / (kappa*(lnbkgpe**2))
         enddo
         deallocate(lnbkgpe)

!        Update potential temperature
!        ----------------------------
         bkg_pt = bkg_pt + (bkg_pt*dpkz - dt) / bkg_pkz !   Convert temperature back to potential temperature

!        Update p-to-the-kappa
!        ---------------------
         bkg_pkz = bkg_pkz + dpkz

!        Update pressures
!        ----------------
         bkg_pe  = bkg_pe  + dpe
         if(any(bkg_pe<0.0d0)) then
           if ( MAPL_AM_I_ROOT() ) then
              print *, 'Update leading to negative pressures, aborting'
           end if
           ASSERT_(.FALSE.)
         endif

         deallocate( dpkz )
         deallocate( dpek )

       else ! approach = 2

!        Update p-to-the-kappa
!        ---------------------
         allocate(dpkz(IM,JM,LM))
         allocate(dpek(IM,JM,LM+1))
         dpek = bkg_pe+dpe
         pkxy = (dpek)**kappa
         do L=1,LM
            dpkz(:,:,L) = ( pkxy(:,:,L+1)-pkxy(:,:,L) ) &
                        / ( kappa*( log(dpek(:,:,L+1))-log(dpek(:,:,L))) )
         enddo
         deallocate(dpek)
         print *, 'after pkz ', minval(dpkz), maxval(dpkz)
         pkxy = bkg_pe**kappa
         do L=1,LM
            bkg_pkz(:,:,L) = ( pkxy(:,:,L+1)-pkxy(:,:,L) ) &
                           / ( kappa*( log(bkg_pe(:,:,L))-log(bkg_pe(:,:,L-1))) )
         enddo
         ! dpkz=ana_pkz
         dpkz = dpkz-bkg_pkz ! dpkz = incremental change in pkz

!        Update potential temperature
!        ----------------------------
         bkg_pt = bkg_pt + (bkg_pt*dpkz - dt) / bkg_pkz !   Convert temperature back to potential temperature

!        Update pressures
!        ----------------
         bkg_pe  = bkg_pe  + dpe
         if(any(bkg_pe<0.0d0)) then
           if ( MAPL_AM_I_ROOT() ) then
              print *, 'Update leading to negative pressures, aborting'
           end if
           ASSERT_(.FALSE.)
         endif
         pkxy = bkg_pe**kappa
         do L=1,LM
            bkg_pkz(:,:,L) = ( pkxy(:,:,L+1)-pkxy(:,:,L) ) &
                         / ( kappa*( log(bkg_pe(:,:,L))-log(bkg_pe(:,:,L-1))) )
         enddo

         deallocate(dpkz)
         endif

    else ! <.not.incremental>

!      Update dry temperature
!      ----------------------
       allocate(ana_td(IM,JM,LM))
       ana_td = bkg_pt*bkg_pkz + dt

!      Update pressures
!      ----------------
       bkg_pe  = bkg_pe  + dpe
       if(any(bkg_pe<0.0d0)) then
         if ( MAPL_AM_I_ROOT() ) then
            print *, 'Update leading to negative pressures, aborting'
         end if
         ASSERT_(.FALSE.)
       endif

!      Update p-to-the-kappa
!      ---------------------
       pkxy = bkg_pe**kappa
       do L=1,LM
          bkg_pkz(:,:,L) = ( pkxy(:,:,L+1)-pkxy(:,:,L) ) &
                         / ( kappa*( log(bkg_pe(:,:,L))-log(bkg_pe(:,:,L-1))) )
       enddo
!      call fv_getPKZ(bkg_pkz,bkg_pt,0.0_r8*bkg_pt,bkg_pe,bkg_dz,.false.)

!      Update potential temperature
!      ----------------------------
       bkg_pt = ana_td / bkg_pkz    !   Convert temperature back to potential temperature

       deallocate(ana_td)
    endif ! <incremental> 

    deallocate(pkxy)

!   Convert A-grid increments to D-grid
!   -----------------------------------
    allocate(dua(im  ,jm  ,lm), & ! U-Wind
             dva(im  ,jm  ,lm)  ) ! V-Wind
    allocate(dud(im  ,jm+1,lm), & ! U-Wind
             dvd(im+1,jm  ,lm)  ) ! V-Wind
    dud=0.0
    dvd=0.0
    dua=du
    dva=dv
    call INTERP_AGRID_TO_DGRID(dua, dva, dud, dvd)

!   Update wind fields
!   ------------------
    bkg_u  = bkg_u  + dud(:,:jm,:)
    bkg_v  = bkg_v  + dvd(:im,:,:)

    deallocate(dua,dva)
    deallocate(dud,dvd)

#endif /* _REALLY_UPDATE_ */

!   Get pointers to export variables
!   --------------------------------
    call MAPL_GetPointer(export,       ak4, 'AK',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,       bk4, 'BK',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,     ana_u,  'U',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,     ana_v,  'V',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,    ana_pt, 'PT',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,    ana_pe, 'PE', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,   ana_pkz,'PKZ', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,    ana_dz, 'DZ', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(export,     ana_w,  'W', RC=STATUS)
    VERIFY_(STATUS)

    if (associated(ak4)) then
       ak4 = ak
    endif
    if (associated(bk4)) then
       bk4 = bk
    endif
    if (associated(ana_u)) then
       ana_u = bkg_u
    end if
    if (associated(ana_v)) then
       ana_v = bkg_v
    end if
    if (associated(ana_pt)) then
       ana_pt = bkg_pt
    end if
    if (associated(ana_pe)) then
       ana_pe = bkg_pe
    end if
    if (associated(ana_pkz)) then
       ana_pkz = bkg_pkz
    end if
    if (associated(ana_dz)) then
       ana_dz = bkg_dz
    end if
    if (associated(ana_w)) then
       ana_w = bkg_w
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine RUN1

end module DynVec_GridCompMod
