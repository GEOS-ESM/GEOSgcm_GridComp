!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_OradGridCompMod -- Implements absorption of solar radiation in the ocean.

! !INTERFACE:

module GEOS_OradGridCompMod

! !USES:

  use ESMF
  use MAPL

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  logical          :: ocean_extData

! !DESCRIPTION:
! 
!   {\tt GEOS\_Orad} is a light-weight gridded component that updates
!      the the solar radiation penetrating the ocean
!

!EOP

   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices, which sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp. The import
!   and internal variables are allocated and initialized by generic.  Here
!   generic is used for the ocean grid.

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer       :: MAPL => null()

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource (MAPL, ocean_extData, Label="OCEAN_EXT_DATA:", DEFAULT=.FALSE., __RC__ ) ! .TRUE. or .FALSE.

!BOS

! !INTERNAL STATE:

   if (ocean_extData) then
     call MAPL_AddImportSpec(GC,                               &
          SHORT_NAME = 'DATA_KPAR',                            &
          LONG_NAME  = 'PAR_extinction_coefficient',           &
          UNITS      = 'm-1',                                  &
          DIMS       = MAPL_DimsHorzOnly,                      &
          VLOCATION  = MAPL_VLocationNone,                     &
          RC=STATUS  )
     VERIFY_(STATUS)
   end if

!  !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'SWHEAT',                            &
        LONG_NAME          = 'solar_heating_rate',                &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'KPAR',                              &
        LONG_NAME          = 'PAR_extinction_coefficient',        &
        units              = 'm-1',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                                 &
        LONG_NAME  = 'cosine_of_the_solar_zenith_angle',              &
        UNITS      = '1',                                             &
        SHORT_NAME = 'COSZ',                                          &
        DIMS       = MAPL_DimsHorzOnly,                               &
        VLOCATION  = MAPL_VLocationNone,                              &
                                                           RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENUVR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENPAR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENUVF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENPAF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DRNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DFNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                 &
        SHORT_NAME = 'FROCEAN',                                       &
        LONG_NAME  = 'ocean_fraction_of_grid_cell',                   &
        UNITS      = '1',                                             &
        DIMS       = MAPL_DimsHorzOnly,                               &
        VLOCATION  = MAPL_VLocationNone,                              &
                                                           RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                 &
        SHORT_NAME = 'DH',                                            &
        LONG_NAME  = 'Layer mass',                                    &
        UNITS      = 'dyn-m',                                         &
        DIMS       = MAPL_DimsHorzVert,                               &
        VLOCATION  = MAPL_VLocationCenter,                            &
                                                           RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENUVR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENPAR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENUVF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'PENPAF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DRNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DFNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWFLX'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME = 'FROCEAN',                                       &
        LONG_NAME  = 'ocean_fraction_of_grid_cell',                   &
        UNITS      = '1',                                             &
        DIMS       = MAPL_DimsHorzOnly,                               &
        VLOCATION  = MAPL_VLocationNone,                              &
                                                           RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME = 'DH',                                            &
        LONG_NAME  = 'Layer mass',                                    &
        UNITS      = 'dyn-m',                                         &
        DIMS       = MAPL_DimsHorzVert,                               &
        VLOCATION  = MAPL_VLocationCenter,                            &
                                                           RC=STATUS  )
     VERIFY_(STATUS)

!EOS

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="RUN"   ,RC=STATUS)
    VERIFY_(STATUS)
  
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

!  All done
!-----------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP

! !IROUTINE: RUN -- First Run stage for the Orad component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the penetrating radiation. Eventually
!   this will do a full radiative transfer calculation.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer   :: MAPL
  type (ESMF_Time)                    :: CurrentTime
  character(len=ESMF_MAXSTR)          :: DATAFILE
  real, pointer, dimension(:,:)       :: KPAR  => null()
  real, pointer, dimension(:,:)       :: Z     => null()
  real, pointer, dimension(:,:)       :: UVR   => null()
  real, pointer, dimension(:,:)       :: PAR   => null()
  real, pointer, dimension(:,:)       :: NIR   => null()
  integer                             :: L
  integer                             :: IM,JM,LM

! pointers to export

   real, pointer, dimension(:,:,:)  :: QSW   => null()
   real, pointer, dimension(:,:  )  :: KPARX => null()

! pointers to import

   real, pointer, dimension(:,:)   :: COSZ      => null()
   real, pointer, dimension(:,:)   :: FR        => null()
   real, pointer, dimension(:,:)   :: PRUVR     => null()
   real, pointer, dimension(:,:)   :: PRPAR     => null()
   real, pointer, dimension(:,:)   :: PRUVF     => null()
   real, pointer, dimension(:,:)   :: PRPAF     => null()
   real, pointer, dimension(:,:)   :: DRNIR     => null()
   real, pointer, dimension(:,:)   :: DFNIR     => null()
   real, pointer, dimension(:,:,:) :: H         => null()
   real, pointer, dimension(:,:)   :: data_kpar => null()

! ponters to export

   real, pointer, dimension(:,:)   :: FRx    => null()
   real, pointer, dimension(:,:)   :: PRUVRx => null()
   real, pointer, dimension(:,:)   :: PRPARx => null()
   real, pointer, dimension(:,:)   :: PRUVFx => null()
   real, pointer, dimension(:,:)   :: PRPAFx => null()
   real, pointer, dimension(:,:)   :: DRNIRx => null()
   real, pointer, dimension(:,:)   :: DFNIRx => null()
   real, pointer, dimension(:,:)   :: SWFLX  => null()
   real, pointer, dimension(:,:,:)   :: Hx   => null()

   real, parameter                 :: KUVR = 0.09

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run"

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

   call MAPL_TimerOn(MAPL,"TOTAL")
   call MAPL_TimerOn(MAPL,"RUN" )

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QSW   , 'SWHEAT'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,KPARX , 'KPAR'    ,    RC=STATUS); VERIFY_(STATUS)

! Our only possible exports are the heating rates
!------------------------------------------------

! Pointers to inputs
!-------------------
   
   call MAPL_GetPointer(IMPORT,COSZ  , 'COSZ'    ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,FR    , 'FROCEAN' ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PRUVF , 'PENUVF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PRPAF , 'PENPAF'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PRUVR , 'PENUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PRPAR , 'PENPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DRNIR , 'DRNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DFNIR , 'DFNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,H     , 'DH'      ,    RC=STATUS); VERIFY_(STATUS)

! Pointers to exports
!-------------------
   
   call GET_POINTER(EXPORT,FRx    , 'FROCEAN' ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,PRUVFx , 'PENUVF'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,PRPAFx , 'PENPAF'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,PRUVRx , 'PENUVR'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,PRPARx , 'PENPAR'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,DRNIRx , 'DRNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,DFNIRx , 'DFNIR'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT,SWFLX , 'SWFLX'  ,    RC=STATUS); VERIFY_(STATUS)
   call GET_POINTER(EXPORT, Hx , 'DH'  ,    RC=STATUS); VERIFY_(STATUS)

! Grid sizes
!-----------
   
   IM = size(H,1)
   JM = size(H,2)
   LM = size(H,3)

   if(associated(KPARX)) then
      KPAR => KPARX
   else
      allocate(KPAR(IM,JM), stat=STATUS)
      VERIFY_(STATUS)
   end if

   allocate(   Z(IM,JM), stat=STATUS)
   VERIFY_(STATUS)
   allocate(UVR (IM,JM), stat=STATUS)
   VERIFY_(STATUS)
   allocate(PAR (IM,JM), stat=STATUS)
   VERIFY_(STATUS)
   allocate(NIR (IM,JM), stat=STATUS)
   VERIFY_(STATUS)

! Get current time from clock
!----------------------------

   call ESMF_ClockGet(CLOCK, currTime=CurrentTime, rc=STATUS)
   VERIFY_(STATUS)

! Get KPAR from data file
!------------------------

   if (.not. ocean_extData) then
     call MAPL_GetResource(MAPL,DATAFILE,LABEL="KPAR_FILE:"     , RC=STATUS)
     VERIFY_(STATUS)

     if (datafile == '/dev/null') then
        kpar = 1.0
     else
        call MAPL_ReadForcing(MAPL,'KPAR',DATAFILE,CURRENTTIME,KPAR, RC=STATUS)
        VERIFY_(STATUS)
     end if
   else
     call MAPL_GetPointer(import, data_kpar, 'DATA_KPAR', __RC__)
     KPAR = data_kpar
   end if

! Use Beer'S Law to compute flux divergence
!------------------------------------------
   UVR=0.0
   PAR=0.0
   NIR=0.0
   where(FR > 0.0)
      UVR = PRUVF+PRUVR
      PAR = PRPAF+PRPAR
      NIR = DRNIR+DFNIR
   endwhere

   Z   = 0.0

   if(associated(QSW)) then
      QSW(:,:,1) = UVR + PAR + NIR

      do L=2,LM   
         Z = Z + H(:,:,L-1)
         QSW(:,:,L  ) = UVR*exp(-KUVR*Z) + PAR*exp(-KPAR*Z)
         QSW(:,:,L-1) = QSW(:,:,L-1) - QSW(:,:,L)
      enddo
      Z = Z + H(:,:,LM)
      QSW(:,:,LM) = QSW(:,:,LM) - (UVR*exp(-KUVR*Z) + PAR*exp(-KPAR*Z))
   end if

   if(.not.associated(KPARX)) deallocate(KPAR)

   deallocate(Z   )
   deallocate(PAR )
   deallocate(UVR )
   deallocate(NIR )

   if (associated(FRx))    FRX = FR
   if (associated(PRUVFx)) PRUVFX = PRUVF
   if (associated(PRPAFx)) PRPAFX = PRPAF
   if (associated(PRUVRx)) PRUVRX = PRUVR
   if (associated(PRPARx)) PRPARX = PRPAR
   if (associated(DRNIRx)) DRNIRX = DRNIR
   if (associated(DFNIRx)) DFNIRX = DFNIR
   if (associated(SWFLX)) SWFLX = DFNIR+DRNIR+PRPAR+PRUVR+PRPAF+PRUVF
   if (associated(Hx)) Hx = H

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN" )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)

 end subroutine RUN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_OradGridCompMod

