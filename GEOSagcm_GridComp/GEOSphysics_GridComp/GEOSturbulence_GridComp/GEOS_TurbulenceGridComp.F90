!   $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_TurbulenceGridCompMod

!BOP

!  !MODULE: GEOS_Turbulence --- An GEOS generic atmospheric turbulence component

! !USES:

  use ESMF
  use GEOS_Mod
  use MAPL
  use LockEntrain
  use shoc
  use edmf_mod, only: run_edmf,mfparams
  use scm_surface, only : surface_layer, surface

#ifdef _CUDA
  use cudafor
#endif

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_TurbulenceGridComp} computes atmospheric tendencies due to turbulence.
!   Its physics is a combination of the first-order scheme of Louis---for stable PBLs
!   and free atmospheric turbulence---with a modified version of the non-local-K
!   scheme proposed by Lock for unstable and cloud-topped boundary layers.
!   In addition to diffusive tendencies, it adds the effects orographic form drag
!   for features with horizontal scales of 2 to 20 km following Beljaars et al. (2003,
!   ECMWF Tech. Memo. 427).  
!
!\vspace{12 pt}
!\noindent
!{\bf Grid Considerations}
!
!   Like all GEOS\_Generic-based components, it works on an inherited 
!   3-dimensional ESMF grid. It assumes that the first two (inner) dimensions span the
!   horizontal and the third (outer) dimension is the vertical. In the horizontal,
!   one or both dimensions can be degenerate, effectively supporting 
!   single-columns (1-D), and slices (2-D). No horizontal dimension needs to be
!   aligned with a particular coordinate. In the vertical, the only assumption
!   is that columns are indexed from top to bottom.
!
!\vspace{12 pt}
!\noindent
!{\bf Methods}
!
!   {\tt GEOS\_TurbulenceGridComp} uses the default Initialize and Finalize methods
!   of GEOS\_Generic. It has a 2-stage Run method that can be used in conjunction with
!   two-stage surface calculations to implement semi-implicit time differencing.
!
!\vspace{12 pt}
!\noindent
!{\bf Time Behavior}
!
!   {\tt GEOS\_TurbulenceGridComp} assumes both run stages will be invoked every 
!   RUN\_DT seconds, where RUN\_DT is required in the configuration. On this interval
!   both run stages will perform diffusion updates using diffusivities found in the
!   internal state.  The diffusivities in the internal state may be refreshed intermitently
!   by specifying MY\_STEP and ACCUMINT in the configuration. Accumulated imports used
!   in the intermittent refreshing are valid only on MY\_STEP intervals. Currently the
!   origin of these intervals is the beginning of the run. Accumulation of these imports
!   is done for a period ACCUMINT prior to the valid time. Both ACCUMINT and MY\_STEP are
!   in seconds.
!
!\vspace{12 pt}
!\noindent
!{\bf Working with Bundles and Friendlies}
!
!   {\tt GEOS\_TurbulenceGridComp} works on bundles of quantities to be diffused
!   and with corresponding bundles of their tendencies, surface values, etc.
!   These bundles may contain an arbitrary number of conservative quantities and
!   no requirements or restrictions are placed on what quantities they contain.
!   Quantities required for the calculation, such as pressures, stability, etc
!   are passed separately from the diffused quantities. Little distinction is made
!   of what is in the bundle, except that needed to decide what diffusivity applies
!   to the quantity and in what form its effects are implemented.
!
!   Quantities to be diffused can be marked as "Friendly-for-diffusion". In that case,
!   {\tt GEOS\_TurbulenceGridComp} directly updates the quantity; otherwise it 
!   merely computes its tendency, placing it in the appropriate bundle and treating
!   the quantity itself as read-only.
!
!   In working with bundled quantities, corresponding fields must appear in the 
!   same order in all bundles. Some of these fields, however, 
!   may be ``empty'' in the sense that the data pointer has not been allocated.
!   
!   {\tt GEOS\_TurbulenceGridComp} works with two bundles; TR in the import
!   state and TRI in the export state. The surface values (SHAT, UHAT, VHAT, QHAT)
!   and their changes (SDEL, UDEL, VDEL, QDEL) are explicit imports, and the
!   surface fluxes (SFLX, UFLX, VFLX, QFLX) and their sensitivities to the surface
!   value (SDFL, UDFL, VDFL, QDFL) are explicit exports. Quantities other than
!   S, U, V and Q have zero surface flux.
!
!   All fields in the export bundles are checked for associated pointers before being
!   updated.
!
!   Fields in the TR bundle can have four attributes:
! \begin{itemize}
! \item FriendlyTo[{\it Component Name}]: default=false --- If true, TR field is updated.
! \item WeightedTendency: default=true --- If true, tendencies (TRI) are pressure-weighted
! \item DiffuseLike: ('S','Q','M') default='S' --- Use mixing coefficients for either
!          heat, moisture or momentum.
! \end{itemize}
!  
!   Only fields in the TR bundle are checked for friendly status. Non-friendly
!   fields in TR and all other bundles are treated with the usual Import/Export
!   rules.
!
!\vspace{12 pt}
!\noindent
!{\bf Other imports and exports}
!
!   In addition to the updates of these bundles, {\tt GEOS\_TurbulenceGridComp} produces
!   a number of diagnostic exports, as well as frictional heating contributions. The latter 
!   are NOT added by {\tt GEOS\_TurbulenceGridComp}, but merely exported to be added
!   elsewhere in the GCM.
!
!\vspace{12 pt}
!\noindent
!{\bf Two-Stage Interactions with the Surface}
!
!   The two-stage scheme for interacting with the surface module is as follows:
! \begin{itemize}
!   \item  The first run stage takes the surface values of the diffused quantities
!      and the surface exchange coefficients as input. These are, of course, on the 
!      grid turbulence is working on.
!   \item  It then does the full diffusion calculation assuming the surface values are
!      fixed, i.e., the explicit surface case. In addition, it also computes derivatives of the
!      tendencies wrt surface values. These are to be used in the second stage.
!   \item The second run stage takes the increments of the surface values as inputs
!      and produces the final results, adding the implicit surface contributions. 
!   \item It also computes the frictional heating due to both implicit and explicit
!       surface contributions.
! \end{itemize}
!
!\vspace{12 pt}
!\noindent
!{\bf GEOS-5 Specific Aspects}
!
!   In GEOS-5, {\tt GEOS\_TurbulenceGridComp} works on the atmosphere's lat-lon grid,
!   while surface quantities are computed during the first run stage of the each of
!   the tiled surface components.  The tiled quantities are properly aggregated to
!   the GEOS-5 lat-lon grid by the first stage of {\tt GEOS\_SurfaceGridComp}, which
!   is called immediately before the first run stage of {\tt GEOS\_TurbulenceGridComp}.
!
!EOP

    logical                             :: dflt_false = .false.
    character(len=ESMF_MAXSTR)          :: dflt_q     = 'Q'

    logical :: DEBUG_TRB

! Generic interfaces for module-scope helpers defined below (after
! `contains`) - see their own doc comments at the definition site for
! why each exists (both replace/stand in for deleted or API-incompatible
! MAPL2-era routines). Fortran requires generic `interface ...
! module procedure` blocks declaring dispatch for same-module procedures
! to live in the module's specification part, not interspersed among the
! procedure bodies themselves.
   interface TURB_MaxMin
      module procedure TURB_MaxMin2
      module procedure TURB_MaxMin3
   end interface TURB_MaxMin

   interface TURB_BundleGetPointerToData
      module procedure TURB_BundleGetPointerToData3
   end interface TURB_BundleGetPointerToData

contains

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!               the Initialize and Finalize services to generic versions. It also
!   allocates our instance of a generic state and puts it in the 
!   gridded component (GC). Here we only set the two-stage run method and
!   declare the data services.
! \newline
! !REVISION HISTORY: 
!   ??Jul2006 E.Novak./Todling - Added output defining TLM/ADM trajectory

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp)                :: GC  ! gridded component
    integer, intent(out)               :: RC  ! return code
!EOP

    integer                            :: STATUS
    integer                            :: DO_WAVES
    integer                            :: DO_SEA_SPRAY
    integer                            :: SCM_SL
    integer                            :: NUMUP

!=============================================================================

! Begin...

    call MAPL_GridCompGetResource( GC, "DEBUG_TRB", DEBUG_TRB, default=.FALSE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompGetResource( GC, "USE_WAVES", DO_WAVES, default=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompGetResource( GC, "USE_SEA_SPRAY", DO_SEA_SPRAY, default=0, RC=STATUS)
    VERIFY_(STATUS)

! Set the Run entry points
! ------------------------
! NOTE: MAPL3 identifies Run phases by name (phase_name=), not by
! call-order-based implicit phase numbering like MAPL2. A parent
! component invokes a specific phase via
! MAPL_GridCompRunChild(gc, "TURBULENCE", phase_name="Run1"/"Run2", RC=...)
! (real precedent: GOCART2G_GridCompMod.F90's own SetServices registers
! its two Run phases the same way, and calls its own children's phases
! the same way).

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run1, phase_name="Run1", RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, Run2, phase_name="Run2", RC=STATUS )
    VERIFY_(STATUS)

! Every declared EXPORT must always be allocated regardless of whether a
! downstream coupler/HISTORY requested it: MAPL2's per-fetch
! `MAPL_GetPointer(EXPORT, ptr, 'NAME', ALLOC=.TRUE.)` (and the
! PDFSHAPE-conditional `ALLOC=PDFALLOC` variant) forced this at fetch
! time; `MAPL_StateGetPointer` (the real MAPL3 replacement) has no
! `alloc=` keyword at all (confirmed: not in
! superstructure/state/StateGetPointer.F90's real interface, despite one
! stale internal call site referencing it - see
! mapl3-porting-notes.md). The real MAPL3 mechanism is component-wide:
! `MAPL_GridCompSet(gc, activate_all_exports=.true.)` forces the
! framework to allocate every EXPORT unconditionally (confirmed via
! OuterMetaComponent/advertise_variable.F90's set_default_activation).
! This is a safe over-approximation of the ~40 individually-ALLOC'd
! MAPL2 exports (a small extra memory cost for fields that were
! genuinely meant to be lazy/optional, in exchange for guaranteed no
! null-pointer crashes on the many exports this component's numerics
! write into unconditionally).
    call MAPL_GridCompSet(GC, activate_all_exports=.true., RC=STATUS)
    VERIFY_(STATUS)

! Get number of EDMF updrafts and SCM single-level flag (both feed
! COND/UNGRIDDED expressions in the ACG-generated includes below, so
! must be set before the #include points)
! -------------------------------------------------------------------
    call MAPL_GridCompGetResource( GC, "EDMF_NUMUP", NUMUP, default=10, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompGetResource( GC, "SCM_SL", SCM_SL, default=0, RC=STATUS)
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------
! Static IMPORT/EXPORT/INTERNAL specs generated by mapl_acg() from
! Turbulence_StateSpecs.rc (see CMakeLists.txt). See
! mapl3-porting-notes.md for the full MAPL2->MAPL3 spec-column mapping
! this was derived from.

#include "Turbulence_Import___.h"

! Dynamically-sized ESMF_FieldBundle imports carrying an arbitrary
! number of diffused-quantity tracers - NOT expressible via the ACG
! YAML (ITEMTYPE column only supports 'F'/'V', no bundle option).
! itemtype=MAPL_STATEITEM_FIELDBUNDLE is the real, working MAPL3
! replacement for MAPL2's DATATYPE=MAPL_BundleItem, confirmed via the
! already-ported GOCART2G_GridCompMod.F90 (its AERO_DP bundle specs).

    call MAPL_GridCompAddSpec(GC,                                          &
         state_intent      = ESMF_STATEINTENT_IMPORT,                      &
         short_name        = 'TR',                                         &
         standard_name     = 'diffused_quantities',                        &
         units             = 'X',                                          &
         dims              = 'xyz',                                        &
         vertical_stagger  = MAPL_VERTICAL_STAGGER_CENTER,                  &
         itemtype          = MAPL_STATEITEM_FIELDBUNDLE,                   &
         restart_mode      = MAPL_RESTART_SKIP,                            &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

#include "Turbulence_Export___.h"

! EDMF_PLUMES_W/_THL/_QT: the ACG `UNGRIDDED`/`ungridded_dim_array`
! column emits a bare `[NUMUP]` integer-array literal, but the real
! `gridcomp_add_spec` dummy argument is `type(UngriddedDim),
! ungridded_dim_array(:)` - a genuine ACG bug (verified by inspecting
! the real generated code and cross-checking MAPL_Generic.F90's actual
! signature; see mapl3-porting-notes.md). Declared manually here instead,
! wrapping NUMUP via the real `MAPL_UngriddedDim(extent, name=, units=)`
! constructor.
    call MAPL_GridCompAddSpec(GC,                                          &
         state_intent      = ESMF_STATEINTENT_EXPORT,                      &
         short_name        = 'EDMF_PLUMES_W',                              &
         standard_name     = 'Vertical_velocity_of_individual_EDMF_plumes',&
         units             = 'm s-1',                                      &
         dims              = 'xyz',                                        &
         vertical_stagger  = MAPL_VERTICAL_STAGGER_EDGE,                   &
         ungridded_dim_array = [MAPL_UngriddedDim(NUMUP, name='plume', units='1')], &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_GridCompAddSpec(GC,                                          &
         state_intent      = ESMF_STATEINTENT_EXPORT,                      &
         short_name        = 'EDMF_PLUMES_THL',                            &
         standard_name     = 'Liquid_water_potential_temperature_of_EDMF_plumes', &
         units             = 'K',                                          &
         dims              = 'xyz',                                        &
         vertical_stagger  = MAPL_VERTICAL_STAGGER_EDGE,                   &
         ungridded_dim_array = [MAPL_UngriddedDim(NUMUP, name='plume', units='1')], &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_GridCompAddSpec(GC,                                          &
         state_intent      = ESMF_STATEINTENT_EXPORT,                      &
         short_name        = 'EDMF_PLUMES_QT',                             &
         standard_name     = 'Total_water_of_individual_EDMF_plumes',      &
         units             = 'kg kg-1',                                    &
         dims              = 'xyz',                                        &
         vertical_stagger  = MAPL_VERTICAL_STAGGER_EDGE,                   &
         ungridded_dim_array = [MAPL_UngriddedDim(NUMUP, name='plume', units='1')], &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

! Dynamically-sized ESMF_FieldBundle exports, one entry per diffused
! tracer in the matching TR import bundle - same bundle-itemtype
! treatment as the IMPORT bundles above.

    call MAPL_GridCompAddSpec(GC,                                          &
         state_intent      = ESMF_STATEINTENT_EXPORT,                      &
         short_name        = 'TRI',                                        &
         standard_name     = 'diffusion_tendencies',                       &
         units             = 'X kg m-2 s-1',                               &
         dims              = 'xyz',                                        &
         vertical_stagger  = MAPL_VERTICAL_STAGGER_CENTER,                  &
         itemtype          = MAPL_STATEITEM_FIELDBUNDLE,                   &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

#include "Turbulence_Internal___.h"

! TKESHOC needs a nonzero fill_value (1e-4, its cold-start minimum TKE)
! - the ACG `FILL` column is currently non-functional (a real, verified
! bug in MAPL_GridCompSpecs_ACG.py's is_printable(), see
! mapl3-porting-notes.md), so declared manually here instead, using the
! `fill_value=` keyword directly (unaffected by the ACG bug - only the
! YAML->ACG path to it is broken).
!
! NOTE: MAPL2 also made TKESHOC (and PDF_A/QT2/QT3, now ADD2EXPORT'd via
! the Internal include above) FRIENDLYTO='DYNAMICS:TURBULENCE' instead
! of just 'TURBULENCE' whenever <COMP_NAME>_DO_SHOC was set - i.e.
! friendly to DYNAMICS too, participating in cross-component
! tracer-advection aggregation, not merely self-exported. MAPL3's
! SERVICE-item-type replacement for that specific cross-component case
! is not yet documented/stable (see
! src/Shared/@MAPL/docs/mapl3/diffs-from-mapl2.md's own "[REVIEW
! NEEDED]" note on this exact mechanism). Deferred - only the
! friendly-to-self (ADD2EXPORT) half is implemented here. TODO: revisit
! once MAPL3 SERVICE semantics for advected-tracer aggregation are
! documented and GEOSsuperdyn's DO_SHOC-aware side (if any) is ported.

    call MAPL_GridCompAddSpec(GC,                                          &
         state_intent      = ESMF_STATEINTENT_INTERNAL,                    &
         short_name        = 'TKESHOC',                                    &
         standard_name     = 'turbulent_kinetic_energy_from_SHOC',         &
         units             = 'm+2 s-2',                                    &
         dims              = 'xyz',                                        &
         vertical_stagger  = MAPL_VERTICAL_STAGGER_CENTER,                  &
         fill_value        = 1.0e-4,                                       &
         add_to_export     = .true.,                                       &
                                                              RC=STATUS  )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

   end subroutine SetServices


!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: Initialize -- Creates the "turbulence_run_alarm" used by Run1
!   to gate the intermittent REFRESH (diffusivity recompute) step.

! !DESCRIPTION: MAPL2's default generic Initialize exposed a component-
!   wide "RUNALARM" (via MAPL_Get(MAPL, RUNALARM=alarm, ...)) that fired
!   on this component's own configured run interval; MAPL3 has no such
!   automatic per-component alarm (MAPL_GridCompGet exposes no RunAlarm
!   output). Since this component never configured a run interval
!   distinct from the model heartbeat, the replacement alarm here rings
!   every heartbeat by default - the same effective cadence as MAPL2's
!   RUNALARM had in the absence of a distinct <NAME>_DT. Matches the
!   pattern established in GEOS_IrradGridComp.F90's own Initialize
!   ("irrad_lw_alarm").

! !INTERFACE:

   subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

     type(ESMF_GridComp)                :: GC
     type(ESMF_State)                   :: IMPORT
     type(ESMF_State)                   :: EXPORT
     type(ESMF_Clock)                   :: CLOCK
     integer,              intent(out)  :: RC

     integer                            :: STATUS
     real                                :: run_dt
     type(ESMF_TimeInterval)            :: run_alarm_interval
     type(ESMF_Alarm)                   :: run_alarm

     call MAPL_ClockGet(CLOCK, dt=run_dt, RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_TimeIntervalSet(run_alarm_interval, s=nint(run_dt), RC=STATUS)
     VERIFY_(STATUS)
     run_alarm = ESMF_AlarmCreate( &
          name="turbulence_run_alarm", &
          clock=CLOCK, &
          ringInterval=run_alarm_interval, &
          sticky=.true., RC=STATUS)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
     _UNUSED_DUMMY(IMPORT)
     _UNUSED_DUMMY(EXPORT)

   end subroutine Initialize

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: RUN1 -- First run stage for the {\tt MAPL_TurbulenceGridComp} component

! !INTERFACE:

   subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp)                :: GC
    type(ESMF_State)                   :: IMPORT
    type(ESMF_State)                   :: EXPORT
    type(ESMF_Clock)                   :: CLOCK
    integer,              intent(out)  :: RC

! !DESCRIPTION: The first run stage of {\tt GEOS\_TurbulenceGridComp} computes the diffusivities,
!   sets-up the matrix for a backward-implicit computation of the surface fluxes,
!   and solves this system for a fixed surface value of the diffused quantity. Run1
!   takes as inputs the surface exchange coefficients (i.e., $\rho |U| C_{m,h,q}$) for
!   momentun, heat, and moisture, as well as the pressure, temperature, moisture, 
!   and winds for the sounding. These are used only for computing the diffusivities
!   and, as explained above,  are not the temperatures, moistures, etc. being diffused.
!
!   The computation of turbulence fluxes for fixed surface values is done at every
!   time step in the contained subroutine {\tt DIFFUSE}; but the computation of 
!   diffusivities and orographic drag coefficients, as well as the set-up of the
!   vertical difference matrix and its LU decomposition
!   can be done intermittently for economy in the contained subroutine  {\tt REFRESH}.
!   The results of this calculation are stored in an internal state. 
!   Run1 also computes the sensitivity of the 
!   atmospheric tendencies and the surface flux to changes in the surface value.
!
!   The diffusivities are computed by calls to {\tt LOUIS\_KS} and {\tt ENTRAIN}, which
!   compute the Louis et al. (1983) and Lock (2000) diffusivities. The Louis 
!   diffusivities are computed for all conditions, and {\tt ENTRAIN} overrides them 
!   where appropriate. Lock can be turned off from the resource file.


!

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (ESMF_State       )            :: INTERNAL 
    type (ESMF_Alarm       )            :: ALARM   
    type (ESMF_Grid        )            :: esmfgrid

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: nn

! Local variables

    real, dimension(:,:,:), pointer     :: AKS, BKS, CKS, DKS
    real, dimension(:,:,:), pointer     :: AKQ, BKQ, CKQ, DKQ
    real, dimension(:,:,:), pointer     :: AKV, BKV, CKV, DKV, EKV, FKV
    real, dimension(:,:,:), pointer     :: PLE, ZLE, SINC
    real, dimension(:,:,:), pointer     :: ZLS, ZLES
    real, dimension(:,:  ), pointer     :: CU, CT, CQ, ZPBL, PHIS
    integer                             :: IM, JM, LM
    real                                :: DT
 
! EDMF-related variables
    real, dimension(:,:,:), pointer    :: AKSS, BKSS, CKSS, YS
    real, dimension(:,:,:), pointer    :: AKQQ, BKQQ, CKQQ, YQV,YQL,YQI
    real, dimension(:,:,:), pointer    :: AKUU, BKUU, CKUU, YU,YV

! SHOC-related variables
    integer                             :: DO_SHOC, SCM_SL
    real, dimension(:,:,:), pointer     :: TKESHOC,TKH,QT2,QT3,WTHV2,WQT_DC, &
                                           PDF_A,MFAW,QTFLXMF,SLFLXMF
   
    real, dimension(:,:), pointer   :: EVAP, SH

! Idealized SCM surface layer variables
    real, dimension(:,:), pointer :: cu_scm, ct_scm, ssurf_scm, qsurf_scm

! Sea spray
    integer :: DO_WAVES
    integer :: DO_SEA_SPRAY
    real, dimension(:,:), pointer :: SH_SPR  => null()
    real, dimension(:,:), pointer :: LH_SPR  => null()
    real, dimension(:,:), pointer :: SH_SPRX => null()
    real, dimension(:,:), pointer :: LH_SPRX => null()


! Begin... 
!---------

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'Run1'

    call MAPL_GridCompTimerStart(GC, "TOTAL")
    call MAPL_GridCompTimerStart(GC, "-RUN1")

! Get parameters from the grid/internal state.
! ---------------------------------------------
! NOTE: MAPL2's MAPL_Get(MAPL, RUNALARM=ALARM, ...) has no MAPL3
! equivalent - MAPL_GridCompGet exposes no RunAlarm output. Instead, a
! "turbulence_run_alarm" is created once in Initialize (below) and
! looked up here by name every call, matching the pattern established
! in GEOS_IrradGridComp.F90's Initialize/Run ("irrad_lw_alarm").

    call MAPL_GridCompGet(GC, grid=esmfgrid, num_levels=LM, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridGet(esmfgrid, IM=IM, JM=JM, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompGetInternalState(GC, INTERNAL, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ClockGetAlarm(CLOCK, alarmname="turbulence_run_alarm", alarm=ALARM, RC=STATUS)
    VERIFY_(STATUS)

! Sea spray
    call MAPL_GridCompGetResource(GC, "USE_WAVES", DO_WAVES, DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompGetResource(GC, "USE_SEA_SPRAY", DO_SEA_SPRAY, DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    if (DO_WAVES/=0 .and. DO_SEA_SPRAY/=0) then
        call MAPL_StateGetPointer(IMPORT, SH_SPR, 'SHFX_SPRAY', RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_StateGetPointer(IMPORT, LH_SPR, 'LHFX_SPRAY', RC=STATUS)
        VERIFY_(STATUS)

        call MAPL_StateGetPointer(EXPORT, SH_SPRX, 'SHFX_SPRAY', RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_StateGetPointer(EXPORT, LH_SPRX, 'LHFX_SPRAY', RC=STATUS)
        VERIFY_(STATUS)

        if (associated(SH_SPRX)) SH_SPRX = SH_SPR
        if (associated(LH_SPRX)) LH_SPRX = LH_SPR
    end if    

! Get all pointers that are needed by both REFRESH and DIFFUSE
!-------------------------------------------------------------

! Get pressure & height structure; this is instantaneous.
!-----------------------------------------------

     call MAPL_StateGetPointer(IMPORT,  PLE,   'PLE',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,  ZLE,   'ZLE',     RC=STATUS)
     VERIFY_(STATUS)

! Get surface exchange coefficients
!----------------------------------

     call MAPL_StateGetPointer(IMPORT,  CU,     'CM',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,  CT,     'CT',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,  CQ,     'CQ',     RC=STATUS)
     VERIFY_(STATUS)

!----- variables needed for SHOC and EDMF -----
    call MAPL_StateGetPointer(IMPORT, SH,   'SH',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(IMPORT, EVAP, 'EVAP',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(IMPORT, WTHV2, 'WTHV2',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(IMPORT, WQT_DC, 'WQT_DC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(IMPORT, PHIS,   'PHIS',    RC=STATUS)
    VERIFY_(STATUS)

!----- Variables for idealized SCM surface layer ------
    call MAPL_GridCompGetResource(GC, "SCM_SL", SCM_SL, default=0, RC=STATUS)
    if (SCM_SL /= 0) then
      call MAPL_StateGetPointer(INTERNAL, cu_scm,    'cu_scm', RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, ct_scm,    'ct_scm', RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, ssurf_scm, 'ssurf_scm', RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, qsurf_scm, 'qsurf_scm', RC=STATUS)
      VERIFY_(STATUS)
    end if

! Get pointers from internal state
!---------------------------------
    call MAPL_StateGetPointer(INTERNAL, AKS,   'AKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, BKS,   'BKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, CKS,   'CKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, DKS,   'DKS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, AKQ,   'AKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, BKQ,   'BKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, CKQ,   'CKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, DKQ,   'DKQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, AKV,   'AKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, BKV,   'BKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, CKV,   'CKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, DKV,   'DKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, EKV,   'EKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, FKV,   'FKV',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, SINC,  'SINC',    RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, ZPBL,  'ZPBL',    RC=STATUS)
    VERIFY_(STATUS)

!----- SHOC-related variables -----
    call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_DO_SHOC", DO_SHOC, &
                           default=0, RC=STATUS)
    call MAPL_StateGetPointer(INTERNAL, TKESHOC,'TKESHOC', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, TKH,    'TKH',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, QT3,    'QT3',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, QT2,    'QT2',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, PDF_A,  'PDF_A',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, MFAW,    'MFAW',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, QTFLXMF, 'QTFLXMF',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, SLFLXMF, 'SLFLXMF',     RC=STATUS)
    VERIFY_(STATUS)
!
! edmf variables
!
    
! a,b,c and rhs for s
    call MAPL_StateGetPointer(INTERNAL, AKSS,   'AKSS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, BKSS,   'BKSS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, CKSS,   'CKSS',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, YS,   'YS',     RC=STATUS)
    VERIFY_(STATUS)
! a,b,c for moisture and rhs for qv,ql,qi    
    call MAPL_StateGetPointer(INTERNAL, AKQQ,   'AKQQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, BKQQ,   'BKQQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, CKQQ,   'CKQQ',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, YQV,   'YQV',     RC=STATUS)
    VERIFY_(STATUS)  
    call MAPL_StateGetPointer(INTERNAL, YQL,   'YQL',     RC=STATUS)
    VERIFY_(STATUS)  
    call MAPL_StateGetPointer(INTERNAL, YQI,   'YQI',     RC=STATUS)
    VERIFY_(STATUS)   
! a,b,c and rhs for wind speed    
    call MAPL_StateGetPointer(INTERNAL, AKUU,   'AKUU',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, BKUU,   'BKUU',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, CKUU,   'CKUU',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, YU,   'YU',     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_StateGetPointer(INTERNAL, YV,   'YV',     RC=STATUS)
    VERIFY_(STATUS)


! Get application's timestep from configuration
!----------------------------------------------

    call MAPL_GridCompGetResource(GC, "RUN_DT", DT, RC=STATUS)
    VERIFY_(STATUS)

! If its time, do the refresh
! ---------------------------

    if ( ESMF_AlarmIsRinging(ALARM, rc=status) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_GridCompTimerStart(GC, "--REFRESHKS")
        call REFRESH(IM,JM,LM,RC=STATUS)
        VERIFY_(STATUS)
       call MAPL_GridCompTimerStop(GC, "--REFRESHKS")
    endif

! Solve the free atmosphere problem
! ---------------------------------

    call MAPL_GridCompTimerStart(GC, "--DIFFUSE")
     call DIFFUSE(IM,JM,LM,RC=STATUS)
     VERIFY_(STATUS)
    call MAPL_GridCompTimerStop(GC, "--DIFFUSE")

!  All done with RUN1
!--------------------

    call MAPL_GridCompTimerStop(GC, "-RUN1")
    call MAPL_GridCompTimerStop(GC, "TOTAL")
    RETURN_(ESMF_SUCCESS)

  contains

!=============================================================================
!=============================================================================

!BOP

! !CROUTINE: REFRESH -- Refreshes diffusivities.

! !INTERFACE:

   subroutine REFRESH(IM,JM,LM,RC)

! !ARGUMENTS:

     integer,           intent(IN)       :: IM,JM,LM
     integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: 
!   {\tt REFRESH} can be called intermittently to compute new values of the 
!   diffusivities. In addition it does all possible calculations that depend
!   only on these. In particular, it sets up the semi-implicit tridiagonal
!   solver in the vertical and does the LU decomposition. It also includes the
!   local effects of orographic drag, so that it to is done implicitly.
!
!   Diffusivities are first computed with the Louis scheme ({\tt LOUIS\_KS}),
!   and then, where appropriate,
!   they are overridden by the Lock values ({\tt ENTRAIN}).
!   Once diffusivities are computed, {\tt REFRESH} sets-up the tridiagonal
!   matrices for the semi-implicit vertical diffusion calculation and performs
!   their $LU$ decomposition. 
!
!   {\tt REFRESH} requires surface exchange coefficients for heat, moisture, and
!   momentum,  The calculations in the interior are also
!   done for momentum, heat, and water diffusion. Heat and water mixing
!   coefficients differ only at the surface, but these affect the entire $LU$
!   decomposition, and so all three decompositions are saved in the internal state. 
!
!   For a conservatively diffused quantity $q$, we have
!   $$
!   \frac{\partial q}{\partial t} = -g \frac{\partial }{\partial p} 
!       \left(\rho K_q \frac{\partial q}{\partial z} \right)
!   $$
!   In finite difference form, using backward time differencing, this becomes
!   $$
!   \begin{array}{rcl}
!   {q^{n+1}_l - q^{n}_l} & = & - \frac{g}{\delta_l p}^*
!     \delta_l \left[
!      \left( \frac{\Delta t \rho K_q}{\delta_l z} \right)^* (\delta_l q)^{n+1} \right]   \\
!   &&\\
!                         & = & - \alpha_l ( \beta_{l+\frac{1}{2}}(q_{l+1}-q_l)^{n+1} - 
!                                            \beta_{l-\frac{1}{2}}(q_l-q_{l-1})^{n+1} ) \\
!   &&\\
!   \alpha_l & = & \frac{g \Delta t}{(p_{l+\frac{1}{2}}-p_{l-\frac{1}{2}})^*} \\
!   &&\\
!   \beta_{l+\frac{1}{2}} & = & \left( \frac{ (\rho K_q)^*_{l+\frac{1}{2}}}{(z_{l+1}-z_{l})^*} \right) \\
!   \end{array}
!   $$
!   where the subscripts denote levels, superscripts denote times, and the $*$ superscript
!   denotes evaluation at the refresh time.
!   The following tridiagonal set is then solved for $q^{n+1}_l$:
!   $$
!   a_l q_{l-1} + b_l q_l + c_l q_{l+1} = q_l
!   $$
!   where
!   $$
!   \begin{array}{rcl}
!   a_l & = & \alpha_l \beta_{l-\frac{1}{2}} \\
!   c_l & = & \alpha_l \beta_{l+\frac{1}{2}} \\
!   b_l & = & 1 - a_l - c_l.
!   \end{array}
!   $$
!   At the top boundary, we assume $K_q=0$, so  $ \beta_{\frac{1}{2}}=0$ and $a_1=0$.
!   At the surface, $ \beta_{L+\frac{1}{2}}= \rho_s |U|_s C_{m,h,q}$, the surface exchange coefficient.
!   

!EOP
 
     character(len=ESMF_MAXSTR)          :: IAm='Refresh'
     integer                             :: STATUS

     character(len=ESMF_MAXSTR)          :: TYPE
     character(len=ESMF_MAXSTR)          :: NAME
     type (ESMF_Field)                   :: FIELD
     type (ESMF_Array)                   :: ARRAY
     type (ESMF_FieldBundle)             :: TR


     real, dimension(:,:,:), pointer     :: TH, U, V, OMEGA, Q, T, RI, DU, RADLW, RADLWC, LWCRT
     real, dimension(:,:  ), pointer     :: AREA, VARFLT
     real, dimension(:,:,:), pointer     :: KH, KM, QLTOT, QITOT, QRTOT, QSTOT, QGTOT, FCLD
     real, dimension(:,:,:), pointer     :: ALM,ALH
     real, dimension(:    ), pointer     :: PREF

     real, dimension(IM,JM,1:LM-1)       :: TVE, RDZ
     real, dimension(IM,JM,LM)           :: THV, TV, Z, DMI, PLO, QL, QI, QA, TSM, USM, VSM
     real, dimension(IM,JM,0:LM)         :: ZL0
     real, dimension(IM,JM)              :: drycblh
     integer, dimension(IM,JM)           :: SMTH_LEV

     real, dimension(:,:,:), pointer     :: EKH, EKM, KHLS, KMLS, KHRAD, KHSFC
     real, dimension(:,:  ), pointer     :: Z0, Z0H, EIS
     real, dimension(:,:  ), pointer     :: BSTAR, USTAR, PPBL, WERAD, WESFC,VSCRAD,KERAD,DBUOY,ZSML,ZCLD,ZRADML,FRLAND,TRINVBS,TRINVFRQ,TRINVDELT
     real, dimension(:,:  ), pointer     :: TCZPBL => null()
     real, dimension(:,:  ), pointer     :: ZPBL2 => null()
     real, dimension(:,:  ), pointer     :: ZPBL10P => null()
     real, dimension(:,:  ), pointer     :: ZPBLHTKE => null()
     real, dimension(:,:,:), pointer     :: TKE => null()
     real, dimension(:,:  ), pointer     :: ZPBLRI => null()
     real, dimension(:,:  ), pointer     :: ZPBLTHV => null()
     real, dimension(:,:  ), pointer     :: ZPBLQV => null()
     real, dimension(:,:  ), pointer     :: ZPBLRFRCT => null()
     real, dimension(:,:  ), pointer     :: SBIFRQ => null()
     real, dimension(:,:  ), pointer     :: SBITOP => null()
     real, dimension(:,:  ), pointer     :: KPBL => null()
     real, dimension(:,:  ), pointer     :: KPBL_SC => null()
     real, dimension(:,:  ), pointer     :: ZPBL_SC => null()                
     real, dimension(:,:  ), pointer     :: WEBRV,VSCBRV,DSIEMS,CHIS,ZCLDTOP,DELSINV,SMIXT,ZRADBS,CLDRF,VSCSFC,RADRCODE

     real, dimension(:,:,:), pointer     :: AKSODT, CKSODT
     real, dimension(:,:,:), pointer     :: AKQODT, CKQODT
     real, dimension(:,:,:), pointer     :: AKVODT, CKVODT

     real, dimension(:,:,:), pointer     :: LSHOC,BRUNTSHOC,ISOTROPY, &
                                            LSHOC1,LSHOC2,LSHOC3, & 
                                            SHOCPRNUM,&
                                            TKEBUOY,TKESHEAR,TKEDISS,TKEDISSx, &
                                            SL2, SL3, W2, W3, WSL, SLQT !, W3CANUTO, QT2DIAG,SL2DIAG,SLQTDIAG
     real, dimension(:,:), pointer       :: edmf_depth, lobukhov
     real                               :: obukhov_denom

! EDMF variables
     real, dimension(:,:,:), pointer     :: edmf_dry_a,edmf_moist_a,edmf_frc, edmf_dry_w,edmf_moist_w, &
                                            edmf_dry_qt,edmf_moist_qt, &
                                            edmf_dry_thl,edmf_moist_thl, &
                                            edmf_dry_u,edmf_moist_u,  &
                                            edmf_dry_v,edmf_moist_v,  &
                                            edmf_moist_qc,edmf_buoyf,edmf_mfx, &
                                            edmf_dmfx, edmf_w2, &
                                            edmf_w3, edmf_wqt, edmf_slqt, & 
                                            edmf_wsl, edmf_qt3, edmf_sl3, &
                                            edmf_entx, edmf_tke,          &
                                            edmf_dqrdt, edmf_dqsdt

   real, dimension(IM,JM,0:LM)          ::  ae3,aw3,aws3,awqv3,awql3,awqi3,awu3,awv3

   real, dimension(IM,JM) :: zpbl_test

   real, dimension(:,:,:,:), pointer    :: EDMF_PLUMES_W, EDMF_PLUMES_THL, EDMF_PLUMES_QT

     logical                             :: ALLOC_TCZPBL, CALC_TCZPBL
     logical                             :: ALLOC_ZPBL2, CALC_ZPBL2
     logical                             :: ALLOC_ZPBL10p, CALC_ZPBL10p
     logical                             :: ALLOC_ZPBLRI, CALC_ZPBLRI
     logical                             :: PDFALLOC

     real                                :: LOUIS_B_KH, LOUIS_B_KM
     real                                :: LOUIS_C_KH, LOUIS_C_KM 
     real                                :: LOUIS_D_KH, LOUIS_D_KM 
     real                                :: ALHFAC, ALMFAC
     real                                :: LAMBDAM, LAMBDAM2
     real                                :: LAMBDAH, LAMBDAH2
     real                                :: ZKHMENV
     real                                :: MINTHICK
     real                                :: MINSHEAR
     real                                :: AKHMMAX
     real                                :: C_B, LAMBDA_B
     logical                             :: USE_EIS
     real                                :: PRANDTLSFC,PRANDTLRAD,BETA_RAD,BETA_SURF,KHRADFAC,TPFAC,ENTRATE_SURF
     real                                :: PCEFF_SURF, VSCALE_SURF, KHSFCFAC_LND, KHSFCFAC_OCN

     real                                :: SMTH_HGT
     integer                             :: I,J,L,LOCK_ON,ITER,MO_MAX_ITER
     integer                             :: KPBLMIN,PBLHT_OPTION

     ! SCM idealized surface-layer parameters
     integer :: SCM_SL          ! 0:    use exchange coefficients from surface grid comp
                                ! else: idealized surface layer specified in AGCM.rc
     integer :: SCM_SL_FLUX     ! 0: prescribed roughness length and surface relative humidity,
                                !    all fluxes from surface layer theory
                                ! 1: prescribed thermodynamic fluxes,
                                !    along with roughness length roughness length and surface relative humidity
                                !    momentum fluxes from surface layer theory
                                ! 2: prescribed thermodynamic fluxes,
                                !    based on SHOBS and LHOBS read from SCM forcing file
                                ! 3: prescribed Monin-Obhkov length,
                                !    along with roughness length and surface relative humidity,
                                !    all fluxes from surface layer theory
                                ! else: use prescribed surface exchange coefficients
     real    :: SCM_SH          ! prescribed surface sensible heat flux (Wm-1) (for SCM_SL_FLUX == 1)
     real    :: SCM_LH          ! prescribed surface latent heat flux (Wm-1) (for SCM_SL_FLUX == 1)
     real    :: SCM_Z0          ! surface roughness length (m)
     real    :: SCM_ZETA        ! Monin-Obkhov length scale (m) (for SCM_SL_FLUX == 3)
     real    :: SCM_RH_SURF     ! Surface relative humidity
     real    :: SCM_TSURF       ! Sea surface temperature (K)
     
     ! SCM idealized surface parameters
     integer :: SCM_SURF      ! 0:    native surface from GEOS
                              ! else: idealized surface with prescribed cooling
     real    :: SCM_DTDT_SURF ! Surface heating rate (Ks-1)
     real, dimension(:,:),   pointer     :: SHOBS, LHOBS

     ! mass-flux constants/parameters
     integer :: DOMF, NumUp
     real    :: L0,L0fac

     real, dimension(IM,JM)    :: L02
     real, dimension(IM,JM,LM) :: QT,THL,SL,EXF

     ! Variables for idealized surface layer     
     real, dimension(IM,JM), target :: bstar_scm, ustar_scm, sh_scm, evap_scm, zeta_scm

     real, dimension(im,jm,0:lm) :: edmfdrya, edmfmoista,     &
                                    edmfdryw, edmfmoistw,     &
                                    edmfdryqt, edmfmoistqt,   &
                                    edmfdrythl, edmfmoistthl, &
                                    edmfdryu, edmfmoistu,     &
                                    edmfdryv, edmfmoistv,     &
                                    edmfmoistqc
     real, dimension(im,jm,lm)   :: zlo, pk, rho
     real, dimension(im,jm)      :: edmfZCLD
     real, dimension(im,jm,0:lm) :: RHOE, RHOAW3, edmf_mf, mfwsl, mfwqt, mftke
     real, dimension(im,jm,lm)   :: buoyf, mfw2, mfw3, mfqt3,     &
                                    mfsl3, mfqt2, mfsl2,   &
                                    mfslqt, edmf_ent

     real                                :: a1,a2
     real,               dimension(IM,JM,LM) :: dum3d,tmp3d,WVP
     real,               dimension(LM+1) :: temparray, htke
     real,               dimension(IM,JM,LM  ) :: tcrib !TransCom bulk Ri
     real,               dimension(LM+1) :: thetav
     real,               dimension(IM,JM,LM+1) :: tmp3de

! variables associated with SHOC
     real, dimension( IM, JM, LM )       :: QPL,QPI
     integer                             :: DO_SHOC, DOPROGQT2
     real                                :: SL2TUNE, QT2TUNE, SLQT2TUNE,          &
                                            SKEW_TGEN, SKEW_TDIS, QT2_TDIS, FREE_ATM_SIGS
     real    :: PDFSHAPE

     real    :: lambdadiss

     integer :: locmax
     real    :: maxkh,minlval
     real, dimension(IM,JM) :: thetavs,thetavh,uv2h,kpbltc,kpbl2,kpbl10p,kpblri
     real    :: maxdthvdz,dthvdz
     real    :: u_shear, v_shear, z_diff
     logical :: has_solid_deck

     ! PBL-top diagnostic
     ! -----------------------------------------

     real, parameter :: tcri_crit = 0.25

     logical :: is_marine_stratocumulus, is_shallow_cumulus
     integer :: cloud_base_idx, cloud_top_idx
     real :: cloud_fraction_mean, cloud_thickness

     real(kind=MAPL_R8), dimension(IM,JM,LM) :: AKX, BKX
     real, dimension(IM,JM,LM) :: DZ, DTM, TM

     logical :: JASON_BELJAARS, JASON_LOUIS, JASON_LOCK, JASON_PBL_SC
     real(kind=MAPL_R8), dimension(IM,JM,LM) :: AERTOT
     real, dimension(:,:,:), pointer     :: S
     integer :: NTR, K, LTOP, LMAX
     real :: maxaero
     real :: eis_stable, kh_thresh

#ifdef _CUDA
     type(dim3) :: Grid, Block
     integer :: blocksize_x, blocksize_y
#endif

! Get tracer bundle for aerosol PBL calculation
!-----------------------------------

    call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)

    call ESMF_FieldBundleGet(TR, fieldCOUNT=NTR, RC=STATUS)
    VERIFY_(STATUS)

! Get Sounding from the import state
!-----------------------------------

     call MAPL_StateGetPointer(IMPORT,     T,       'T', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,     Q,      'QV', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,    TH,      'TH', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,     U,       'U', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,     V,       'V', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, OMEGA,   'OMEGA', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,  AREA,   'AREA',  RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,VARFLT,  'VARFLT', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,  PREF,    'PREF', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, RADLW,   'RADLW', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,RADLWC,  'RADLWC', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, QLTOT,   'QLTOT', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, QITOT,   'QITOT', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, QRTOT,   'QRTOT', RC=STATUS); VERIFY_(STATUS) 
     call MAPL_StateGetPointer(IMPORT, QSTOT,   'QSTOT', RC=STATUS); VERIFY_(STATUS) 
     call MAPL_StateGetPointer(IMPORT, QGTOT,   'QGTOT', RC=STATUS); VERIFY_(STATUS) 
     call MAPL_StateGetPointer(IMPORT,  FCLD,    'FCLD', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, BSTAR,   'BSTAR', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT, USTAR,   'USTAR', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,FRLAND,  'FRLAND', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,    Z0,      'Z0', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,   Z0H,     'Z0H', RC=STATUS); VERIFY_(STATUS)
     call MAPL_StateGetPointer(IMPORT,   EIS,     'EIS', RC=STATUS); VERIFY_(STATUS)


! Get turbulence parameters from configuration
!---------------------------------------------
     if (LM .eq. 72) then
       call MAPL_GridCompGetResource(GC, "JASON_PBL_SC", JASON_PBL_SC, default=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, "JASON_BELJAARS", JASON_BELJAARS, default=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, "JASON_LOUIS", JASON_LOUIS, default=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, "JASON_LOCK", JASON_LOCK, default=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PBLHT_OPTION", PBLHT_OPTION, default=4, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SMTH_HGT", SMTH_HGT, default=0.0, RC=STATUS); VERIFY_(STATUS)
     else
       call MAPL_GridCompGetResource(GC, "JASON_PBL_SC", JASON_PBL_SC, default=.FALSE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, "JASON_BELJAARS", JASON_BELJAARS, default=.FALSE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, "JASON_LOUIS", JASON_LOUIS, default=.FALSE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, "JASON_LOCK", JASON_LOCK, default=.FALSE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PBLHT_OPTION", PBLHT_OPTION, default=5, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SMTH_HGT", SMTH_HGT, default=300.0, RC=STATUS); VERIFY_(STATUS)
     endif

     ! Beljaars
     if (JASON_BELJAARS) then
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_C_B", C_B, default=6.0, RC=STATUS); VERIFY_(STATUS)
     else
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_C_B", C_B, default=-3.0, RC=STATUS); VERIFY_(STATUS)
     endif
     ! Louis
     if (JASON_LOUIS) then
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_B_KH", LOUIS_B_KH, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_C_KH", LOUIS_C_KH, default=7.5, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_D_KH", LOUIS_D_KH, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_B_KM", LOUIS_B_KM, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_C_KM", LOUIS_C_KM, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_D_KM", LOUIS_D_KM, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ALHFAC", ALHFAC, default=1.2, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ALMFAC", ALMFAC, default=1.2, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAM", LAMBDAM, default=160.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAM2", LAMBDAM2, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAH", LAMBDAH, default=160.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAH2", LAMBDAH2, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ZKHMENV", ZKHMENV, default=3000., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_MINTHICK", MINTHICK, default=2.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_MINSHEAR", MINSHEAR, default=0.0030, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDA_B", LAMBDA_B, default=1500., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_AKHMMAX", AKHMMAX, default=500., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_MO_MAX_ITER", MO_MAX_ITER, default=0, RC=STATUS); VERIFY_(STATUS)
     else
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_B_KH", LOUIS_B_KH, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_C_KH", LOUIS_C_KH, default=7.5, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_D_KH", LOUIS_D_KH, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_B_KM", LOUIS_B_KM, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_C_KM", LOUIS_C_KM, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOUIS_D_KM", LOUIS_D_KM, default=5.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ALHFAC", ALHFAC, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ALMFAC", ALMFAC, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAM", LAMBDAM, default=150.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAM2", LAMBDAM2, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAH", LAMBDAH, default=200.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDAH2", LAMBDAH2, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ZKHMENV", ZKHMENV, default=3200., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_MINTHICK", MINTHICK, default=2.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_MINSHEAR", MINSHEAR, default=0.0030, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDA_B", LAMBDA_B, default=1500., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_AKHMMAX", AKHMMAX, default=500., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_MO_MAX_ITER", MO_MAX_ITER, default=0, RC=STATUS); VERIFY_(STATUS)
     endif
     ! Lock
     if (JASON_LOCK) then
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDADISS", LAMBDADISS, default=50.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_KHRADFAC", KHRADFAC, default=0.85, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_KHSFCFAC_LND", KHSFCFAC_LND, default=0.60, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_KHSFCFAC_OCN", KHSFCFAC_OCN, default=0.30, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PRANDTLSFC", PRANDTLSFC, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PRANDTLRAD", PRANDTLRAD, default=0.75, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_BETA_RAD", BETA_RAD, default=0.20, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_BETA_SURF", BETA_SURF, default=0.25, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ENTRATE_SURF", ENTRATE_SURF, default=1.5e-3, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_TPFAC", TPFAC, default=20.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PCEFF_SURF", PCEFF_SURF, default=0.5, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOCK_ON", LOCK_ON, default=1, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_VSCALE_SURF", VSCALE_SURF, default=2.5e-3, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_USE_EIS", USE_EIS, default=.false., RC=STATUS); VERIFY_(STATUS)
     else
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LAMBDADISS", LAMBDADISS, default=15., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_KHRADFAC", KHRADFAC, default=0.8, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_KHSFCFAC_LND", KHSFCFAC_LND, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_KHSFCFAC_OCN", KHSFCFAC_OCN, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PRANDTLSFC", PRANDTLSFC, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PRANDTLRAD", PRANDTLRAD, default=0.75, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_BETA_RAD", BETA_RAD, default=0.15, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_BETA_SURF", BETA_SURF, default=0.10, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_ENTRATE_SURF", ENTRATE_SURF, default=1.5e-3, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_TPFAC", TPFAC, default=0.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_PCEFF_SURF", PCEFF_SURF, default=0.375, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_LOCK_ON", LOCK_ON, default=1, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_VSCALE_SURF", VSCALE_SURF, default=2.5e-3, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_USE_EIS", USE_EIS, default=.false., RC=STATUS); VERIFY_(STATUS)
     endif

     call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_DO_SHOC", DO_SHOC, default=0, RC=STATUS); VERIFY_(STATUS)
     if (DO_SHOC /= 0) then
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_PRNUM", SHOCPARAMS%PRNUM, default=-0.9, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_LAMBDA", SHOCPARAMS%LAMBDA, default=0.25, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_TSCALE", SHOCPARAMS%TSCALE, default=400., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_CK", SHOCPARAMS%CKVAL, default=0.1, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_CEFAC", SHOCPARAMS%CEFAC, default=1.0, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_CESFAC", SHOCPARAMS%CESFAC, default=4., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_LENOPT", SHOCPARAMS%LENOPT, default=3, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_LENFAC1", SHOCPARAMS%LENFAC1, default=5., RC=STATUS); VERIFY_(STATUS)       
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_LENFAC2", SHOCPARAMS%LENFAC2, default=2., RC=STATUS); VERIFY_(STATUS)       
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_LENFAC3", SHOCPARAMS%LENFAC3, default=0.5, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_SHC_BUOY_OPTION", SHOCPARAMS%BUOYOPT, default=2, RC=STATUS); VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, 'PDFSHAPE', PDFSHAPE, DEFAULT = 6.0, RC=STATUS); VERIFY_(STATUS)
     else
       call MAPL_GridCompGetResource(GC, 'PDFSHAPE', PDFSHAPE, DEFAULT = 1.0, RC=STATUS); VERIFY_(STATUS)
     end if

     call MAPL_GridCompGetResource(GC, 'DOPROGQT2', DOPROGQT2, DEFAULT = 1, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'SL2TUNE', SL2TUNE, DEFAULT = 4.0, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'QT2TUNE', QT2TUNE, DEFAULT = 1.0, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'SLQT2TUNE', SLQT2TUNE, DEFAULT = 7.0, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'SKEW_TDIS', SKEW_TDIS, DEFAULT = 900.0, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'SKEW_TGEN', SKEW_TGEN, DEFAULT = 900.0, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'QT2_TDIS', QT2_TDIS, DEFAULT = 900.0, RC=STATUS); VERIFY_(STATUS)
     call MAPL_GridCompGetResource(GC, 'FREE_ATM_SIGS', FREE_ATM_SIGS, DEFAULT = 0.05, RC=STATUS); VERIFY_(STATUS)

! Get pointers from export state...
!-----------------------------------

     PDFALLOC = (PDFSHAPE.ge.5)

     call MAPL_StateGetPointer(EXPORT,      KH,      'KH', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,      KM,      'KM', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,      RI,      'RI', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,      DU,      'DU',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,     EKH,     'EKH', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,     EKM,     'EKM', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    KHLS,    'KHLS', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    KMLS,    'KMLS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   KHSFC,   'KHSFC', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   KHRAD,   'KHRAD', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    PPBL,    'PPBL',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    KPBL,    'KPBL',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    KPBL_SC, 'KPBL_SC', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBL_SC, 'ZPBL_SC',            RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    TCZPBL,  'TCZPBL',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBL2,  'ZPBL2',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBL10p,  'ZPBL10p',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBLHTKE,  'ZPBLHTKE',         RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    TKE,  'TKE',         RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBLRI,  'ZPBLRI',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBLTHV,  'ZPBLTHV',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBLQV,  'ZPBLQV',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZPBLRFRCT, 'ZPBLRFRCT',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    SBIFRQ,  'SBIFRQ',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    SBITOP,  'SBITOP',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   LWCRT,   'LWCRT', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   WERAD,   'WERAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   WESFC,   'WESFC',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   DBUOY,   'DBUOY',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, TRINVBS, 'TRINVBS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, TRINVDELT, 'TRINVDELT',           RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, TRINVFRQ, 'TRINVFRQ',             RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  VSCRAD,  'VSCRAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  VSCsfc,  'VSCSFC',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   KERAD,   'KERAD',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  VSCBRV,  'VSCBRV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   WEBRV,   'WEBRV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    CHIS,    'CHIS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  DSIEMS,  'DSIEMS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZCLD,    'ZCLD', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZSML,    'ZSML', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  ZRADML,  'ZRADML', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  ZRADBS,  'ZRADBS', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, ZCLDTOP, 'ZCLDTOP',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, DELSINV, 'DELSINV',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,RADRCODE,'RADRCODE',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   SMIXT,   'SMIXT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,   CLDRF,   'CLDRF',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,     ALH,     'ALH',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,     ALM,     'ALM',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  AKSODT,  'AKSODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  CKSODT,  'CKSODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  AKQODT,  'AKQODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  CKQODT,  'CKQODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  AKVODT,  'AKVODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  CKVODT,  'CKVODT',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,     ZLS,     'ZLS',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,    ZLES,    'ZLES',               RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, EDMF_PLUMES_W, 'EDMF_PLUMES_W',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, EDMF_PLUMES_QT, 'EDMF_PLUMES_QT', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, EDMF_PLUMES_THL, 'EDMF_PLUMES_THL', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dqrdt,  'EDMF_DQRDT', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dqsdt,  'EDMF_DQSDT', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_buoyf, 'EDMF_BUOYF',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_slqt, 'EDMF_SLQT', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_w2,   'EDMF_W2', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_w3,   'EDMF_W3', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_qt3,  'EDMF_QT3', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_sl3,  'EDMF_SL3', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  slqt,  'SLQT',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  w3,    'W3',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  w2,    'W2',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  sl3,   'SL3',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  sl2,   'SL2',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  lobukhov,  'LOBUKHOV', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_wqt,    'EDMF_WQT', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_wsl,    'EDMF_WSL', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_tke,    'EDMF_TKE', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_mfx,    'EDMF_MF', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dmfx,   'EDMF_DMF', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dry_a,  'EDMF_DRY_A',       RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_a, 'EDMF_MOIST_A',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  EDMF_FRC,    'EDMF_FRC', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dry_u,  'EDMF_DRY_U',      RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_u, 'EDMF_MOIST_U',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dry_v,   'EDMF_DRY_V',       RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_v, 'EDMF_MOIST_V',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dry_w,   'EDMF_DRY_W',      RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_w, 'EDMF_MOIST_W',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dry_qt,  'EDMF_DRY_QT',    RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_qt,  'EDMF_MOIST_QT',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_dry_thl,   'EDMF_DRY_THL',     RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_thl, 'EDMF_MOIST_THL',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_moist_qc,  'EDMF_MOIST_QC',    RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_entx,      'EDMF_ENTR', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT,  edmf_depth,     'EDMF_DEPTH', RC=STATUS)
     VERIFY_(STATUS)

!========== SHOC ===========
     call MAPL_StateGetPointer(EXPORT, TKEDISSx, 'TKEDISS',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, TKEBUOY, 'TKEBUOY',  RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, TKESHEAR,'TKESHEAR', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, ISOTROPY,'ISOTROPY', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, LSHOC,   'LSHOC',    RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, LSHOC1,  'LSHOC1',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, LSHOC2,  'LSHOC2',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, LSHOC3,  'LSHOC3',   RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, BRUNTSHOC, 'BRUNTSHOC', RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_StateGetPointer(EXPORT, SHOCPRNUM,'SHOCPRNUM', RC=STATUS)
     VERIFY_(STATUS)

     call MAPL_StateGetPointer(INTERNAL, TKEDISS, 'TKEDISS' , RC=STATUS); VERIFY_(STATUS)
! Initialize some arrays

      LWCRT = RADLW - RADLWC

      KH    = 0.0
      KM    = 0.0
      RI    = 0.0
      EKH   = 0.0
      EKM   = 0.0
      KHSFC = 0.0
      KHRAD = 0.0
      if(associated( ALH))  ALH = 0.0
      if(associated( ALM))  ALM = 0.0 
      if(associated(KHLS)) KHLS = 0.0
      if(associated(KMLS)) KMLS = 0.0

      ALLOC_ZPBL2 = .FALSE.
      CALC_ZPBL2 = .FALSE.
      if(associated(ZPBL2).OR.PBLHT_OPTION==1) CALC_ZPBL2 = .TRUE.
      if(.not.associated(ZPBL2 )) then
         allocate(ZPBL2(IM,JM))
         ALLOC_ZPBL2 = .TRUE.
      endif

      ALLOC_ZPBLRI = .FALSE.
      CALC_ZPBLRI = .FALSE.
      if(associated(ZPBLRI).OR.PBLHT_OPTION==5) CALC_ZPBLRI = .TRUE.
      if(.not.associated(ZPBLRI)) then
         allocate(ZPBLRI(IM,JM))
         ALLOC_ZPBLRI = .TRUE.
      endif

      ALLOC_ZPBL10p = .FALSE.
      CALC_ZPBL10p = .FALSE.
      if(associated(ZPBL10p).OR.PBLHT_OPTION==2.OR.PBLHT_OPTION==4) CALC_ZPBL10p = .TRUE.
      if(.not.associated(ZPBL10p )) then
         allocate(ZPBL10p(IM,JM))
         ALLOC_ZPBL10p = .TRUE.
      endif

      ALLOC_TCZPBL = .FALSE.
      CALC_TCZPBL = .FALSE.
      if(associated(TCZPBL).OR.PBLHT_OPTION==3.OR.PBLHT_OPTION==4) CALC_TCZPBL = .TRUE.
      if(.not.associated(TCZPBL)) then
                allocate(TCZPBL(IM,JM))
                   ALLOC_TCZPBL = .TRUE.
      endif

      do L=0,LM
         ZL0(:,:,L) = ZLE(:,:,L) - ZLE(:,:,LM) ! edge height above the surface 
      enddo

      if (SMTH_HGT > 0) then
         ! Use Thickness at the surface to determine index
         SMTH_LEV=LM
         do L=LM,1,-1
            where (ZL0(:,:,L) <= SMTH_HGT .and. ZL0(:,:,L-1) > SMTH_HGT)
               SMTH_LEV=L
            endwhere
         enddo
      else
         SMTH_LEV=LM-5
      end if

      call MAPL_GridCompTimerStart(GC, "---PRELIMS")

      ! Layer height, pressure, and virtual temperatures
      !-------------------------------------------------

      QL  = QLTOT+QRTOT
      QI  = QITOT+QSTOT+QGTOT
      QA  = FCLD
      Z   = 0.5*(ZL0(:,:,0:LM-1)+ZL0(:,:,1:LM)) ! layer height above surface
      PLO = 0.5*(PLE(:,:,0:LM-1)+PLE(:,:,1:LM))

      if (associated(ZLS))  ZLS = Z
      if (associated(ZLES)) ZLES = ZL0

      TV  = T *( 1.0 + MAPL_VIREPS * Q - QL - QI ) 
      THV = TV*(TH/T)

      TVE = (TV(:,:,1:LM-1) + TV(:,:,2:LM))*0.5

      ! Miscellaneous factors
      !----------------------

      RDZ = PLE(:,:,1:LM-1) / ( MAPL_RGAS * TVE )
      RDZ = RDZ(:,:,1:LM-1) / (Z(:,:,1:LM-1)-Z(:,:,2:LM))
      DMI = (MAPL_GRAV*DT)/(PLE(:,:,1:LM)-PLE(:,:,0:LM-1))

      TSM = THV
      USM = U
      VSM = V
      if (DO_SHOC == 0) then
      !===> Running 1-2-1 smooth of bottom levels of THV, U and V
      if (SMTH_HGT > 0) then
        TSM(:,:,LM) = THV(:,:,LM-1)*0.25 + THV(:,:,LM  )*0.75
        USM(:,:,LM) =   U(:,:,LM-1)*0.25 +   U(:,:,LM  )*0.75
        VSM(:,:,LM) =   V(:,:,LM-1)*0.25 +   V(:,:,LM  )*0.75
        do J=1,JM
         do I=1,IM
           do L=LM-1,SMTH_LEV(I,J),-1
              TSM(I,J,L) = THV(I,J,L-1)*0.25 + THV(I,J,L)*0.50 + THV(I,J,L+1)*0.25
              USM(I,J,L) =   U(I,J,L-1)*0.25 +   U(I,J,L)*0.50 +   U(I,J,L+1)*0.25
              VSM(I,J,L) =   V(I,J,L-1)*0.25 +   V(I,J,L)*0.50 +   V(I,J,L+1)*0.25
           end do
         end do
        end do
      else
        TSM(:,:,LM) = TSM(:,:,LM-1)*0.25 + TSM(:,:,LM  )*0.75
        USM(:,:,LM) = USM(:,:,LM-1)*0.25 + USM(:,:,LM  )*0.75
        VSM(:,:,LM) = VSM(:,:,LM-1)*0.25 + VSM(:,:,LM  )*0.75
        do J=1,JM
         do I=1,IM
           do L=LM-1,SMTH_LEV(I,J),-1
              TSM(I,J,L) = TSM(I,J,L-1)*0.25 + TSM(I,J,L)*0.50 + TSM(I,J,L+1)*0.25
              USM(I,J,L) = USM(I,J,L-1)*0.25 + USM(I,J,L)*0.50 + USM(I,J,L+1)*0.25
              VSM(I,J,L) = VSM(I,J,L-1)*0.25 + VSM(I,J,L)*0.50 + VSM(I,J,L+1)*0.25
           end do
         end do
        end do
      end if
      end if

      RHOE(:,:,1:LM-1)=PLE(:,:,1:LM-1)/(MAPL_RGAS*TVE) 
      RHOE(:,:,0)=PLE(:,:,0)/(MAPL_RGAS*TV(:,:,1))
      RHOE(:,:,LM)=PLE(:,:,LM)/(MAPL_RGAS*TV(:,:,LM))

      rho = plo/( MAPL_RGAS*tv )

      call MAPL_GridCompTimerStop(GC, "---PRELIMS")

   ! Calculate liquid water potential temperature (THL) and total water (QT)
    EXF=T/TH 
    THL=TH-(MAPL_ALHL*QLTOT+MAPL_ALHS*QITOT)/(MAPL_CP*EXF)
    QT=Q+QLTOT+QITOT

! get updraft constants
    call MAPL_GridCompGetResource(GC, "EDMF_DOMF", DOMF, default=0, RC=STATUS)
    
    if ( DOMF /= 0 ) then
      ! number of updrafts
      call MAPL_GridCompGetResource(GC, "EDMF_NUMUP", MFPARAMS%NUP, default=10, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_DOTRACERS", MFPARAMS%DOTRACERS, default=.true., RC=STATUS)
      ! boundaries for the updraft area (min/max sigma of w pdf)
      call MAPL_GridCompGetResource(GC, "EDMF_PWMIN", MFPARAMS%PWMIN, default=1.2, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_PWMAX", MFPARAMS%PWMAX, default=3., RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_ENTUFAC", MFPARAMS%ENTUFAC, default=2.0, RC=STATUS)  
      call MAPL_GridCompGetResource(GC, "EDMF_WA", MFPARAMS%WA, default=1.0, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_WB", MFPARAMS%WB, default=1.5, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_WC", MFPARAMS%WC, default=0.05, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_WCTHRESH", MFPARAMS%WCTHRESH, default=11., RC=STATUS)
      ! coefficients for surface forcing, appropriate for L137
      call MAPL_GridCompGetResource(GC, "EDMF_ALPHAW", MFPARAMS%AlphaW, default=0.05, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_ALPHAQT", MFPARAMS%AlphaQT, default=1.0, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_ALPHATH", MFPARAMS%AlphaTH, default=1.0, RC=STATUS) 
      ! Entrainment rate options
      call MAPL_GridCompGetResource(GC, "EDMF_ET", MFPARAMS%ET, default=2, RC=STATUS)
      ! constant entrainment rate   
      call MAPL_GridCompGetResource(GC, "EDMF_ENT0", MFPARAMS%ENT0, default=0.35, RC=STATUS)
      ! L0 if ET==1
      call MAPL_GridCompGetResource(GC, "EDMF_L0", MFPARAMS%L0, default=100., RC=STATUS)
      ! L0fac if ET==2
      call MAPL_GridCompGetResource(GC, "EDMF_L0FAC", MFPARAMS%L0fac, default=10., RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_MFLIMFAC", MFPARAMS%MFLIMFAC, default=2.0, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_ICE_RAMP", MFPARAMS%ICE_RAMP, default=-40.0, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_ENTRAIN", MFPARAMS%ENTRAIN, default=0, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_STOCHASTIC", MFPARAMS%STOCHFRAC, default=0.4, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_DISCRETE_TYPE", MFPARAMS%DISCRETE, default=1, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_IMPLICIT", MFPARAMS%IMPLICIT, default=1, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_PRCPCRIT", MFPARAMS%PRCPCRIT, default=-1., RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_UPABUOYDEP", MFPARAMS%UPABUOYDEP, default=1, RC=STATUS)
      call MAPL_GridCompGetResource(GC, "EDMF_TREFF", MFPARAMS%TREFF, default=1., RC=STATUS)
    else
       MFPARAMS%TREFF = 0.
       MFPARAMS%DOTRACERS = .false.
    end if

    call MAPL_GridCompGetResource(GC, 'SCM_SL', SCM_SL, DEFAULT=0)


if (SCM_SL /= 0) then
    call MAPL_GridCompGetResource(GC, 'SCM_SURF', SCM_SURF, DEFAULT=0)
    call MAPL_GridCompGetResource(GC, 'SCM_DTDT_SURF', SCM_DTDT_SURF, DEFAULT=0.)

    call MAPL_GridCompGetResource(GC, 'SCM_SL_FLUX', SCM_SL_FLUX, DEFAULT=0)
    call MAPL_GridCompGetResource(GC, 'SCM_SH', SCM_SH, DEFAULT=0.)
    call MAPL_GridCompGetResource(GC, 'SCM_LH', SCM_LH, DEFAULT=0.)
    call MAPL_GridCompGetResource(GC, 'SCM_Z0', SCM_Z0, DEFAULT=1.E-4)
    call MAPL_GridCompGetResource(GC, 'SCM_RH_SURF', SCM_RH_SURF, DEFAULT=0.98)
    call MAPL_GridCompGetResource(GC, 'SCM_TSURF', SCM_TSURF, DEFAULT=298.76) ! S6
!    call MAPL_GridCompGetResource(GC, 'SCM_TSURF', SCM_TSURF, DEFAULT=292.46) ! S11
!    call MAPL_GridCompGetResource(GC, 'SCM_TSURF', SCM_TSURF, DEFAULT=290.96) ! S12
    call MAPL_GridCompGetResource(GC, 'SCM_ZETA', SCM_ZETA, DEFAULT=-0.012957419628129) ! S6
!    call MAPL_GridCompGetResource(GC, 'SCM_ZETA', SCM_ZETA, DEFAULT=-0.013215659785478) ! S11
!    call MAPL_GridCompGetResource(GC, 'SCM_ZETA', SCM_ZETA, DEFAULT=-0.007700882024895) ! S12

       call MAPL_GridCompTimerStart(GC, "---SURFACE")

       call MAPL_StateGetPointer(IMPORT, SHOBS,'SHOBS', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_StateGetPointer(IMPORT, LHOBS,'LHOBS', RC=STATUS)
       VERIFY_(STATUS)


       if ( SCM_SL_FLUX == 1 ) then
          sh_scm(:,:)   = scm_sh
          evap_scm(:,:) = scm_lh/MAPL_ALHL
       elseif ( SCM_SL_FLUX == 2 ) then
          sh_scm(:,:)   = shobs
          evap_scm(:,:) = lhobs/MAPL_ALHL
       elseif ( SCM_SL_FLUX == 3 ) then
          zeta_scm(:,:) = scm_zeta
       end if

       call surface(IM, JM, LM, &                                         ! in
                    SCM_SURF, SCM_TSURF, SCM_RH_SURF, SCM_DTDT_SURF, & ! in
                    dt, ple, &                                            ! in
                    ssurf_scm, &                                          ! inout
                    qsurf_scm)                                            ! out

       call surface_layer(IM, JM, LM, &
                          SCM_SL_FLUX, SCM_Z0, &
                          zpbl, ssurf_scm, qsurf_scm, &
                          z, zl0, ple, rhoe, u, v, T, q, thv, &
                          sh_scm, evap_scm, zeta_scm, &
                          ustar_scm, cu_scm, ct_scm)

       cu_scm = 0.015
       cu => cu_scm
       ct => ct_scm
       cq => ct_scm

       sh    => sh_scm
       evap  => evap_scm
              
       ustar_scm = sqrt( CU*sqrt(U(:,:,LM)**2+V(:,:,LM)**2+0.01) / RHOE(:,:,LM) )

       bstar_scm = (MAPL_GRAV/(RHOE(:,:,LM)*ustar_scm)) *  &
                   (SH/MAPL_CP/THV(:,:,LM) + MAPL_VIREPS*EVAP)

       bstar => bstar_scm
       ustar => ustar_scm
       
!       print *,'bstar=',bstar_scm,'  ustar=',ustar_scm       

       call MAPL_GridCompTimerStop(GC, "---SURFACE")
end if




!===============================================================
!                      EDMF Mass Flux
!===============================================================
    call MAPL_GridCompTimerStart(GC, "---MASSFLUX")

! Initialize EDMF output variables
    mfsl2  = 0.0
    mfslqt = 0.0
    mfqt2  = 0.0
    mfw2   = 0.0
    mfw3   = 0.0
    mfqt3  = 0.0
    mfsl3  = 0.0
    mfwqt  = 0.0
    mfwsl  = 0.0
    mftke  = 0.0
    YS     = 0.0
    YQV    = 0.0
    YQL    = 0.0
    YQI    = 0.0
    YU     = 0.0
    YV     = 0.0
    awqv3  = 0.0
    awql3  = 0.0
    awqi3  = 0.0

    IF(DOMF /= 0) then

      call RUN_EDMF(1, IM, 1, JM, 1, LM, DT,  &
                    !== Inputs ==
                    PHIS,                     &
                    Z,                        &
                    ZL0,                      &
                    PLE,                      &
                    RHOE,                     &
                    TKESHOC,                  &
                    U,                        &
                    V,                        & 
                    T,                        & 
                    THL,                      & 
                    THV,                      & 
                    Q,                        & 
                    QLTOT,                    & 
                    QITOT,                    & 
                    SH,                       & 
                    EVAP,                     & 
                    FRLAND,                   &
                    ZPBL,                     &
                    !== Outputs for trisolver ==
                    ae3,                      &
                    aw3,                      &
                    aws3,                     &
                    awqv3,                    &
                    awql3,                    &
                    awqi3,                    &
                    awu3,                     &
                    awv3,                     &
                    YS,                       &
                    YQV,                      &
                    YQL,                      &
                    YQI,                      &
                    YU,                       &
                    YV,                       &
                    !== Outputs for ADG PDF ==
                    mfw2,                     &
                    mfw3,                     &
                    mfqt3,                    &
                    mfsl3,                    &
                    mfwqt,                    &
                    mfslqt,                   &
                    mfwsl,                    &
                    !== Outputs for SHOC ==
                    mftke,                    &
                    buoyf,                    &
                    edmf_mf,                  & ! needed for ADG PDF
                    edmfdrya, edmfmoista,     & ! outputs for ADG PDF
                    edmf_dqrdt, edmf_dqsdt,   & ! output for micro
                    edmf_dmfx,                &
                    !== Diagnostics, not used elsewhere ==
                    edmf_dry_w,               &
                    edmf_moist_w,             &
                    edmf_dry_qt,              &
                    edmf_moist_qt,            &
                    edmf_dry_thl,             &
                    edmf_moist_thl,           &
                    edmf_dry_u,               &
                    edmf_moist_u,             &
                    edmf_dry_v,               &
                    edmf_moist_v,             &
                    edmf_moist_qc,            &
                    edmf_entx,                &
                    edmf_depth,               &
                    EDMF_PLUMES_W,            &
                    EDMF_PLUMES_THL,          &
                    EDMF_PLUMES_QT )

      mfaw         = aw3
      slflxmf      = (aws3-awql3*mapl_alhl-awqi3*mapl_alhs)/mapl_cp
      qtflxmf      = awqv3+awql3+awqi3

      !=== Fill Exports ===
      if (associated(edmf_dry_a))     edmf_dry_a   = edmfdrya 
      if (associated(edmf_moist_a))   edmf_moist_a = edmfmoista 
      if (associated(edmf_buoyf))     edmf_buoyf   = buoyf 
      if (associated(edmf_mfx))       edmf_mfx     = edmf_mf
      if (associated(edmf_w2))        edmf_w2      = mfw2
      if (associated(edmf_w3))        edmf_w3      = mfw3
      if (associated(edmf_qt3))       edmf_qt3     = mfqt3
      if (associated(edmf_sl3))       edmf_sl3     = mfsl3
      if (associated(edmf_wqt))       edmf_wqt     = mfwqt
      if (associated(edmf_slqt))      edmf_slqt    = mfslqt
      if (associated(edmf_wsl))       edmf_wsl     = mfwsl
      if (associated(edmf_tke))       edmf_tke     = mftke
      if (associated(EDMF_FRC))       EDMF_FRC     = 0.5*(edmfdrya(:,:,0:LM-1)+edmfdrya(:,:,1:LM) &
                                                     + edmfmoista(:,:,0:LM-1)+edmfmoista(:,:,1:LM)) 
      do i = 1,IM
         do j = 1,JM
            k = LM
            do while (edmfdrya(i,j,k).gt.0.01 .and. edmfmoista(i,j,k).lt.1e-4 .and. k.gt.1)
               k = k-1
            end do
            drycblh(i,j) = ZL0(i,j,k)
         end do
      end do

    ELSE            ! if there is no mass-flux
      ae3   = 1.0
      aw3   = 0.0
      aws3  = 0.0
      awqv3 = 0.0
      awql3 = 0.0
      awqi3 = 0.0
      awu3  = 0.0
      awv3  = 0.0
      buoyf = 0.0  
      mfaw    = 0.0
      slflxmf = 0.0
      qtflxmf = 0.0

      if (associated(edmf_dry_a))     edmf_dry_a    = 0.0
      if (associated(edmf_moist_a))   edmf_moist_a  = 0.0
      if (associated(edmf_dry_w))     edmf_dry_w    = MAPL_UNDEF
      if (associated(edmf_moist_w))   edmf_moist_w  = MAPL_UNDEF 
      if (associated(edmf_dry_qt))    edmf_dry_qt   = MAPL_UNDEF
      if (associated(edmf_moist_qt))  edmf_moist_qt = MAPL_UNDEF 
      if (associated(edmf_dry_thl))   edmf_dry_thl  = MAPL_UNDEF 
      if (associated(edmf_moist_thl)) edmf_moist_thl= MAPL_UNDEF 
      if (associated(edmf_dry_u))     edmf_dry_u    = MAPL_UNDEF 
      if (associated(edmf_moist_u))   edmf_moist_u  = MAPL_UNDEF 
      if (associated(edmf_dry_v))     edmf_dry_v    = MAPL_UNDEF 
      if (associated(edmf_moist_v))   edmf_moist_v  = MAPL_UNDEF 
      if (associated(edmf_moist_qc))  edmf_moist_qc = MAPL_UNDEF 
      if (associated(edmf_buoyf))     edmf_buoyf    = 0.0
      if (associated(edmf_entx))      edmf_entx     = MAPL_UNDEF
      if (associated(edmf_mfx))       edmf_mfx      = 0.0 
      if (associated(edmf_w2))        edmf_w2       = mfw2
      if (associated(edmf_w3))        edmf_w3       = mfw3
      if (associated(edmf_qt3))       edmf_qt3      = mfqt3
      if (associated(edmf_sl3))       edmf_sl3      = mfsl3
      if (associated(edmf_wqt))       edmf_wqt      = mfwqt
      if (associated(edmf_slqt))      edmf_slqt     = mfslqt
      if (associated(edmf_wsl))       edmf_wsl      = mfwsl
      if (associated(edmf_tke))       edmf_tke      = mftke
      if (associated(EDMF_FRC))       EDMF_FRC = 0.

      drycblh = 0.     
   ENDIF



   
   call MAPL_GridCompTimerStop(GC, "---MASSFLUX")

   if (associated(lobukhov)) then
      do j = 1, JM
         do i = 1, IM
            obukhov_denom = 0.4*MAPL_GRAV*(sh(i,j)+mapl_epsilon*thv(i,j,LM)*evap(i,j)) / &
                            (MAPL_CP*rhoe(i,j,LM))
            if (obukhov_denom /= 0.0) then
               lobukhov(i,j) = -ustar(i,j)**3 * thv(i,j,LM) / obukhov_denom
            else
               ! SHOC uses zero Obukhov length as the neutral sentinel.
               lobukhov(i,j) = 0.0
            end if
         end do
      end do
   end if

!!!=================================================================
!!!===========================  SHOC  ==============================
!!!=================================================================

      if ( DO_SHOC /= 0 ) then

        LOCK_ON = 0
        ISOTROPY = 600.

        call MAPL_GridCompTimerStart(GC, "---SHOC", RC=STATUS)
        VERIFY_(STATUS)

        call RUN_SHOC( IM, JM, LM, LM+1, DT,  &
                       !== Inputs ==
                       SH(:,:),               &
                       LOBUKHOV(:,:),         &
                       PLO(:,:,1:LM),         &
                       ZL0(:,:,0:LM),         &
                       Z(:,:,1:LM),           &
                       U(:,:,1:LM),           &
                       V(:,:,1:LM),           &
                       OMEGA(:,:,1:LM),       &
                       T(:,:,1:LM),           &
                       Q(:,:,1:LM),           &
                       QITOT(:,:,1:LM),       &
                       QLTOT(:,:,1:LM),       &
                       QPI(:,:,1:LM),         &
                       QPL(:,:,1:LM),         &
                       QA(:,:,1:LM),          &
                       WTHV2(:,:,1:LM),       &
                       BUOYF(:,:,1:LM),       &
                       MFTKE(:,:,0:LM),       &
                       DRYCBLH(:,:),          &    
                       !== Input-Outputs ==
                       TKESHOC(:,:,1:LM),     &
                       TKH(:,:,1:LM),         &
                       !== Outputs ==
                       KM(:,:,1:LM),          &
                       ISOTROPY(:,:,1:LM),    &
                       TKEDISS,               &
                       !== Diagnostics ==  ! not used elsewhere
                       TKEBUOY,               &
                       TKESHEAR,              &
                       LSHOC,                 &
                       LSHOC1,                &
                       LSHOC2,                &
                       LSHOC3,                &
                       BRUNTSHOC,             &
                       RI,                    &
                       SHOCPRNUM,             &
                       !== Tuning params ==
                       SHOCPARAMS )

        if (associated(TKEDISSx)) TKEDISSx = TKEDISS
        KH(:,:,1:LM) = TKH(:,:,1:LM)

        call MAPL_GridCompTimerStop(GC, "---SHOC", RC=STATUS)
        VERIFY_(STATUS)

      end if  ! DOSHOC condition

!   Refresh diffusivities: First compute Louis...
!   ---------------------------------------------

      call MAPL_GridCompTimerStart(GC, "---LOUIS", RC=STATUS)
      VERIFY_(STATUS)

      if (DO_SHOC == 0) then
        if (JASON_LOUIS) then
        call LOUIS_KS( IM,JM,LM,            &
            Z,ZL0,TSM,USM,VSM,ZPBL,         &
            KH, KM, RI,                     &
            LOUIS_B_KH,LOUIS_B_KM,          &
            MINSHEAR, MINTHICK,             &
            LAMBDAM, LAMBDAM2,              & 
            LAMBDAH, LAMBDAH2,              & 
            ALHFAC, ALMFAC,                 &
            ZKHMENV, AKHMMAX,               &
            DU, ALH, KMLS, KHLS             )
        else
        call LOUIS_KS_OPTIMIZED( IM,JM,LM,MO_MAX_ITER,DT, &
            Z,ZL0,TSM,USM,VSM,ZPBL,Z0,Z0H,    &
            KH, KM, RI,                     &
            LOUIS_B_KH,LOUIS_B_KM,          & 
            LOUIS_C_KH,LOUIS_C_KM,          &
            LOUIS_D_KH,LOUIS_D_KM,          &
            MINSHEAR, MINTHICK,             &
            LAMBDAM, LAMBDAM2,              & 
            LAMBDAH, LAMBDAH2,              & 
            ALHFAC, ALMFAC,                 &
            ZKHMENV, AKHMMAX,               &
            DU, ALM, ALH, KMLS, KHLS             )
         endif
      end if


      call MAPL_GridCompTimerStop(GC, "---LOUIS", RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GridCompTimerStart(GC, "---LOCK", RC=STATUS)
      VERIFY_(STATUS)

      !   ...then add Lock.
      !--------------------

      DO_ENTRAIN: if (LOCK_ON==1) then

#ifdef _CUDA

         _ASSERT(LM <= GPU_MAXLEVS,'needs informative message') !If this is tripped, GNUmakefile
                                    !must be changed

         call MAPL_GridCompGetResource(GC, 'BLOCKSIZE_X', BLOCKSIZE_X, DEFAULT=16, __RC__)
         call MAPL_GridCompGetResource(GC, 'BLOCKSIZE_Y', BLOCKSIZE_Y, DEFAULT=8, __RC__)

         Block = dim3(blocksize_x,blocksize_y,1)
         Grid = dim3(ceiling(real(IM)/real(blocksize_x)),ceiling(real(JM)/real(blocksize_y)),1)

         call MAPL_GridCompTimerStart(GC, "----LOCK_ALLOC", __RC__)

         ! ----------------------
         ! Allocate device arrays
         ! ----------------------

         ! Inputs - Lock
         ! -------------
      
         ALLOCATE(TDTLW_IN_dev(IM,JM,LM), __STAT__)
         ALLOCATE(U_STAR_dev(IM,JM), __STAT__)
         ALLOCATE(B_STAR_dev(IM,JM), __STAT__)
         ALLOCATE(FRLAND_dev(IM,JM), __STAT__)
         ALLOCATE(T_dev(IM,JM,LM), __STAT__)
         ALLOCATE(QV_dev(IM,JM,LM), __STAT__)
         ALLOCATE(QL_dev(IM,JM,LM), __STAT__)
         ALLOCATE(QI_dev(IM,JM,LM), __STAT__)
         ALLOCATE(U_dev(IM,JM,LM), __STAT__)
         ALLOCATE(V_dev(IM,JM,LM), __STAT__)
         ALLOCATE(ZFULL_dev(IM,JM,LM), __STAT__)
         ALLOCATE(PFULL_dev(IM,JM,LM), __STAT__)
         ALLOCATE(ZHALF_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(PHALF_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(EIS_dev(IM,JM), __STAT__)                     

         ! Inoutputs - Lock
         ! ----------------
      
         ALLOCATE(DIFF_M_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(DIFF_T_dev(IM,JM,LM+1), __STAT__)

         ! Outputs - Lock
         ! --------------

         ALLOCATE(K_M_ENTR_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(K_T_ENTR_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(K_SFC_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(K_RAD_dev(IM,JM,LM+1), __STAT__)
         ALLOCATE(ZCLOUD_dev(IM,JM), __STAT__)
         ALLOCATE(ZRADML_dev(IM,JM), __STAT__)
         ALLOCATE(ZRADBASE_dev(IM,JM), __STAT__)
         ALLOCATE(ZSML_dev(IM,JM), __STAT__)

         ! Diagnostics - Lock
         ! ------------------

         ! MAT: Using device pointers on CUDA is a bit convoluted. First, we
         ! only allocate the actual working arrays on the device if the 
         ! EXPORT pointer is associated.

         IF (ASSOCIATED(ZCLDTOP))  ALLOCATE(ZCLDTOP_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(WESFC))    ALLOCATE(WENTR_SFC_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(WERAD))    ALLOCATE(WENTR_RAD_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(DBUOY))    ALLOCATE(DEL_BUOY_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(VSCSFC))   ALLOCATE(VSFC_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(VSCRAD))   ALLOCATE(VRAD_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(KERAD))    ALLOCATE(KENTRAD_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(VSCBRV))   ALLOCATE(VBRV_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(WEBRV))    ALLOCATE(WENTR_BRV_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(DSIEMS))   ALLOCATE(DSIEMS_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(CHIS))     ALLOCATE(CHIS_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(DELSINV))  ALLOCATE(DELSINV_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(SMIXT))    ALLOCATE(SLMIXTURE_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(CLDRF))    ALLOCATE(CLDRADF_DIAG_dev(IM,JM), __STAT__)
         IF (ASSOCIATED(RADRCODE)) ALLOCATE(RADRCODE_DIAG_dev(IM,JM), __STAT__)

         ! Then we associate the CUDA device pointer to the associated device
         ! array. That way CUDA knows what memory that pointer belongs to.
         ! We then pass in the pointer to the subroutine.

         IF (ASSOCIATED(ZCLDTOP))  ZCLDTOP_DIAG_dev_ptr   => ZCLDTOP_DIAG_dev
         IF (ASSOCIATED(WESFC))    WENTR_SFC_DIAG_dev_ptr => WENTR_SFC_DIAG_dev
         IF (ASSOCIATED(WERAD))    WENTR_RAD_DIAG_dev_ptr => WENTR_RAD_DIAG_dev
         IF (ASSOCIATED(DBUOY))    DEL_BUOY_DIAG_dev_ptr  => DEL_BUOY_DIAG_dev
         IF (ASSOCIATED(VSCSFC))   VSFC_DIAG_dev_ptr      => VSFC_DIAG_dev
         IF (ASSOCIATED(VSCRAD))   VRAD_DIAG_dev_ptr      => VRAD_DIAG_dev
         IF (ASSOCIATED(KERAD))    KENTRAD_DIAG_dev_ptr   => KENTRAD_DIAG_dev
         IF (ASSOCIATED(VSCBRV))   VBRV_DIAG_dev_ptr      => VBRV_DIAG_dev
         IF (ASSOCIATED(WEBRV))    WENTR_BRV_DIAG_dev_ptr => WENTR_BRV_DIAG_dev
         IF (ASSOCIATED(DSIEMS))   DSIEMS_DIAG_dev_ptr    => DSIEMS_DIAG_dev
         IF (ASSOCIATED(CHIS))     CHIS_DIAG_dev_ptr      => CHIS_DIAG_dev
         IF (ASSOCIATED(DELSINV))  DELSINV_DIAG_dev_ptr   => DELSINV_DIAG_dev
         IF (ASSOCIATED(SMIXT))    SLMIXTURE_DIAG_dev_ptr => SLMIXTURE_DIAG_dev
         IF (ASSOCIATED(CLDRF))    CLDRADF_DIAG_dev_ptr   => CLDRADF_DIAG_dev
         IF (ASSOCIATED(RADRCODE)) RADRCODE_DIAG_dev_ptr  => RADRCODE_DIAG_dev

         call MAPL_GridCompTimerStop(GC, "----LOCK_ALLOC", __RC__)

         call MAPL_GridCompTimerStart(GC, "----LOCK_DATA", __RC__)

         ! ---------------------
         ! Copy inputs to device
         ! ---------------------

         ! Inputs
         ! ------

                  TDTLW_IN_dev = RADLW
                  U_STAR_dev   = USTAR
                  B_STAR_dev   = BSTAR
                  FRLAND_dev   = FRLAND
                  EVAP_dev     = EVAP
                  SH_dev       = SH
                  T_dev        = T
                  QV_dev       = Q
                  QL_dev       = QLTOT
                  QI_dev       = QITOT
                  U_dev        = U
                  V_dev        = V
                  ZFULL_dev    = Z
                  PFULL_dev    = PLO
                  EIS_dev      = EIS
         ZHALF_dev(:,:,1:LM+1) = ZL0(:,:,0:LM)
         PHALF_dev(:,:,1:LM+1) = PLE(:,:,0:LM)

         ! Inoutputs - Lock
         ! ----------------
      
         DIFF_M_dev(:,:,1:LM+1) = KM(:,:,0:LM)
         DIFF_T_dev(:,:,1:LM+1) = KH(:,:,0:LM)

         call MAPL_GridCompTimerStop(GC, "----LOCK_DATA", __RC__)

         call MAPL_GridCompTimerStart(GC, "----LOCK_RUN", __RC__)

         call ENTRAIN<<<Grid,Block>>>(IM, JM, LM,   &
                                      ! Inputs
                                      TDTLW_IN_dev,   &
                                      U_STAR_dev,     &
                                      B_STAR_dev,     &
                                      FRLAND_dev,     &
                                      EVAP_dev,       &
                                      SH_dev,         &
                                      T_dev,          &
                                      QV_dev,         &
                                      QL_dev,         &
                                      QI_dev,         &
                                      U_dev,          &
                                      V_dev,          &
                                      ZFULL_dev,      &
                                      PFULL_dev,      &
                                      ZHALF_dev,      &
                                      PHALF_dev,      &
                                      EIS_dev,        &
                                      ! Inoutputs
                                      DIFF_M_dev,     &
                                      DIFF_T_dev,     &
                                      ! Outputs
                                      K_M_ENTR_dev,   &
                                      K_T_ENTR_dev,   &
                                      K_SFC_dev,      &
                                      K_RAD_dev,      &
                                      ZCLOUD_dev,     &
                                      ZRADML_dev,     &
                                      ZRADBASE_dev,   &
                                      ZSML_dev,       &
                                      ! Diagnostics
                                      ZCLDTOP_DIAG_dev_ptr, &
                                      WENTR_SFC_DIAG_dev_ptr, &
                                      WENTR_RAD_DIAG_dev_ptr, &
                                      DEL_BUOY_DIAG_dev_ptr, &
                                      VSFC_DIAG_dev_ptr, &
                                      VRAD_DIAG_dev_ptr, &
                                      KENTRAD_DIAG_dev_ptr, &
                                      VBRV_DIAG_dev_ptr, &
                                      WENTR_BRV_DIAG_dev_ptr, &
                                      DSIEMS_DIAG_dev_ptr, &
                                      CHIS_DIAG_dev_ptr, &
                                      DELSINV_DIAG_dev_ptr, &
                                      SLMIXTURE_DIAG_dev_ptr, &
                                      CLDRADF_DIAG_dev_ptr, &
                                      RADRCODE_DIAG_dev_ptr, &
                                      ! Input parameter constants
                                      USE_EIS, &
                                      PRANDTLSFC, PRANDTLRAD,   &
                                      BETA_SURF, BETA_RAD,      &
                                      TPFAC, ENTRATE_SURF, &
                                      PCEFF_SURF, VSCALE_SURF, KHRADFAC, KHSFCFAC_LND, KHSFCFAC_OCN )


         STATUS = cudaGetLastError()
         if (STATUS /= 0) then 
            write (*,*) "Error code from ENTRAIN kernel call: ", STATUS
            write (*,*) "Kernel call failed: ", cudaGetErrorString(STATUS)
            _ASSERT(.FALSE.,'needs informative message')
         end if

         ! --------------
         ! Kernel is done
         ! --------------

         call MAPL_GridCompTimerStop(GC, "----LOCK_RUN", __RC__)

         call MAPL_GridCompTimerStart(GC, "----LOCK_DATA", __RC__)

         ! ------------------------
         ! Copy outputs to the host
         ! ------------------------

         ! Inoutputs - Lock
         ! ----------------
      
            KM(:,:,0:LM) = DIFF_M_dev(:,:,1:LM+1)
            KH(:,:,0:LM) = DIFF_T_dev(:,:,1:LM+1)

         ! Outputs - Lock
         ! --------------
      
           EKM(:,:,0:LM) = K_M_ENTR_dev(:,:,1:LM+1)
           EKH(:,:,0:LM) = K_T_ENTR_dev(:,:,1:LM+1)
         KHSFC(:,:,0:LM) = K_SFC_dev(:,:,1:LM+1)
         KHRAD(:,:,0:LM) = K_RAD_dev(:,:,1:LM+1)
                  ZCLD   = ZCLOUD_dev
                  ZRADML = ZRADML_dev
                  ZRADBS = ZRADBASE_dev
                  ZSML   = ZSML_dev
      
         ! Diagnostics - Lock
         ! ------------------
      
         IF (ASSOCIATED(ZCLDTOP))  ZCLDTOP  = ZCLDTOP_DIAG_dev
         IF (ASSOCIATED(WESFC))    WESFC    = WENTR_SFC_DIAG_dev
         IF (ASSOCIATED(WERAD))    WERAD    = WENTR_RAD_DIAG_dev
         IF (ASSOCIATED(DBUOY))    DBUOY    = DEL_BUOY_DIAG_dev
         IF (ASSOCIATED(VSCSFC))   VSCSFC   = VSFC_DIAG_dev
         IF (ASSOCIATED(VSCRAD))   VSCRAD   = VRAD_DIAG_dev
         IF (ASSOCIATED(KERAD))    KERAD    = KENTRAD_DIAG_dev
         IF (ASSOCIATED(VSCBRV))   VSCBRV   = VBRV_DIAG_dev
         IF (ASSOCIATED(WEBRV))    WEBRV    = WENTR_BRV_DIAG_dev
         IF (ASSOCIATED(DSIEMS))   DSIEMS   = DSIEMS_DIAG_dev
         IF (ASSOCIATED(CHIS))     CHIS     = CHIS_DIAG_dev
         IF (ASSOCIATED(DELSINV))  DELSINV  = DELSINV_DIAG_dev
         IF (ASSOCIATED(SMIXT))    SMIXT    = SLMIXTURE_DIAG_dev
         IF (ASSOCIATED(CLDRF))    CLDRF    = CLDRADF_DIAG_dev
         IF (ASSOCIATED(RADRCODE)) RADRCODE = RADRCODE_DIAG_dev

         call MAPL_GridCompTimerStop(GC, "----LOCK_DATA", __RC__)

         call MAPL_GridCompTimerStart(GC, "----LOCK_DEALLOC", __RC__)

         ! ------------------------
         ! Deallocate device arrays
         ! ------------------------
   
         ! Inputs - Lock
         ! -------------
   
         DEALLOCATE(TDTLW_IN_dev)
         DEALLOCATE(U_STAR_dev)
         DEALLOCATE(B_STAR_dev)
         DEALLOCATE(FRLAND_dev)
         DEALLOCATE(EVAP_dev)
         DEALLOCATE(SH_dev)
         DEALLOCATE(T_dev)
         DEALLOCATE(QV_dev)
         DEALLOCATE(QL_dev)
         DEALLOCATE(QI_dev)
         DEALLOCATE(U_dev)
         DEALLOCATE(V_dev)
         DEALLOCATE(ZFULL_dev)
         DEALLOCATE(PFULL_dev)
         DEALLOCATE(ZHALF_dev)
         DEALLOCATE(PHALF_dev)
      
         ! Inoutputs - Lock
         ! ----------------
   
         DEALLOCATE(DIFF_M_dev)
         DEALLOCATE(DIFF_T_dev)
   
         ! Outputs - Lock
         ! --------------
      
         DEALLOCATE(K_M_ENTR_dev)
         DEALLOCATE(K_T_ENTR_dev)
         DEALLOCATE(K_SFC_dev)
         DEALLOCATE(K_RAD_dev)
         DEALLOCATE(ZCLOUD_dev)
         DEALLOCATE(ZRADML_dev)
         DEALLOCATE(ZRADBASE_dev)
         DEALLOCATE(ZSML_dev)
      
         ! Diagnostics - Lock
         ! ------------------

         ! MAT Again, we only deallocate a device array if the diagnostic
         ! was asked for.
      
         IF (ASSOCIATED(ZCLDTOP))  DEALLOCATE(ZCLDTOP_DIAG_dev)
         IF (ASSOCIATED(WESFC))    DEALLOCATE(WENTR_SFC_DIAG_dev)
         IF (ASSOCIATED(WERAD))    DEALLOCATE(WENTR_RAD_DIAG_dev)
         IF (ASSOCIATED(DBUOY))    DEALLOCATE(DEL_BUOY_DIAG_dev)
         IF (ASSOCIATED(VSCSFC))   DEALLOCATE(VSFC_DIAG_dev)
         IF (ASSOCIATED(VSCRAD))   DEALLOCATE(VRAD_DIAG_dev)
         IF (ASSOCIATED(KERAD))    DEALLOCATE(KENTRAD_DIAG_dev)
         IF (ASSOCIATED(VSCBRV))   DEALLOCATE(VBRV_DIAG_dev)
         IF (ASSOCIATED(WEBRV))    DEALLOCATE(WENTR_BRV_DIAG_dev)
         IF (ASSOCIATED(DSIEMS))   DEALLOCATE(DSIEMS_DIAG_dev)
         IF (ASSOCIATED(CHIS))     DEALLOCATE(CHIS_DIAG_dev)
         IF (ASSOCIATED(DELSINV))  DEALLOCATE(DELSINV_DIAG_dev)
         IF (ASSOCIATED(SMIXT))    DEALLOCATE(SLMIXTURE_DIAG_dev)
         IF (ASSOCIATED(CLDRF))    DEALLOCATE(CLDRADF_DIAG_dev)
         IF (ASSOCIATED(RADRCODE)) DEALLOCATE(RADRCODE_DIAG_dev)

         ! This step is probably unnecessary, but better safe than sorry
         ! as the lifetime of a device pointer is not really specified
         ! by NVIDIA
   
         IF (ASSOCIATED(ZCLDTOP))  NULLIFY(ZCLDTOP_DIAG_dev_ptr)
         IF (ASSOCIATED(WESFC))    NULLIFY(WENTR_SFC_DIAG_dev_ptr)
         IF (ASSOCIATED(WERAD))    NULLIFY(WENTR_RAD_DIAG_dev_ptr)
         IF (ASSOCIATED(DBUOY))    NULLIFY(DEL_BUOY_DIAG_dev_ptr)
         IF (ASSOCIATED(VSCSFC))   NULLIFY(VSFC_DIAG_dev_ptr)
         IF (ASSOCIATED(VSCRAD))   NULLIFY(VRAD_DIAG_dev_ptr)
         IF (ASSOCIATED(KERAD))    NULLIFY(KENTRAD_DIAG_dev_ptr)
         IF (ASSOCIATED(VSCBRV))   NULLIFY(VBRV_DIAG_dev_ptr)
         IF (ASSOCIATED(WEBRV))    NULLIFY(WENTR_BRV_DIAG_dev_ptr)
         IF (ASSOCIATED(DSIEMS))   NULLIFY(DSIEMS_DIAG_dev_ptr)
         IF (ASSOCIATED(CHIS))     NULLIFY(CHIS_DIAG_dev_ptr)
         IF (ASSOCIATED(DELSINV))  NULLIFY(DELSINV_DIAG_dev_ptr)
         IF (ASSOCIATED(SMIXT))    NULLIFY(SLMIXTURE_DIAG_dev_ptr)
         IF (ASSOCIATED(CLDRF))    NULLIFY(CLDRADF_DIAG_dev_ptr)
         IF (ASSOCIATED(RADRCODE)) NULLIFY(RADRCODE_DIAG_dev_ptr)

         call MAPL_GridCompTimerStop(GC, "----LOCK_DEALLOC", __RC__)

#else

!   ...then add Lock.
!--------------------

         CALL ENTRAIN(IM,JM,LM,                 &
                      ! Inputs
                      RADLW,                    &
                      USTAR,                    &
                      BSTAR,                    &
                      FRLAND,                   &
                      EVAP,                     &
                      SH,                       &
                      T,                        &
                      Q,                        &
                      QLTOT,                    &
                      QITOT,                    &
                      U,                        &
                      V,                        &
                      Z,                        &
                      PLO,                      &
                      ZL0,                      &
                      PLE,                      &
                      EIS,                      &
                      ! Inoutputs
                      KM,                       &
                      KH,                       &
                      ! Outputs
                      EKM,                      &
                      EKH,                      &
                      KHSFC,                    &
                      KHRAD,                    &
                      ZCLD,                     &
                      ZRADML,                   &
                      ZRADBS,                   &
                      ZSML,                     &
                      ! Diagnostics
                      ZCLDTOP,                  &
                      WESFC,                    &
                      WERAD,                    &
                      DBUOY,                    &
                      VSCSFC,                   &
                      VSCRAD,                   &
                      KERAD,                    &
                      VSCBRV,                   &
                      WEBRV,                    &
                      DSIEMS,                   &
                      CHIS,                     &
                      DELSINV,                  &
                      SMIXT,                    &
                      CLDRF,                    &
                      RADRCODE,                 &
                      ! Input parameter constants
                      USE_EIS,                  &
                      PRANDTLSFC, PRANDTLRAD,   &
                      BETA_SURF, BETA_RAD,      &
                      TPFAC, ENTRATE_SURF, &
                      PCEFF_SURF, VSCALE_SURF, KHRADFAC, KHSFCFAC_LND, KHSFCFAC_OCN )

#endif

      else ! Not running ENTRAIN...
         EKM   = 0.0
         EKH   = 0.0
         KHSFC = 0.0
         KHRAD = 0.0
      end if DO_ENTRAIN

      call MAPL_GridCompTimerStop(GC, "---LOCK", RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GridCompTimerStart(GC, "---POSTLOCK")



      ! TKE 
      if (associated(TKE)) then         ! Reminder: TKE is on model edges
        if (DO_SHOC /= 0) then          !           TKESHOC is not.
          TKE(:,:,1:LM-1) = 0.5*(TKESHOC(:,:,1:LM-1)+TKESHOC(:,:,2:LM))
          TKE(:,:,0) = 1e-6
          TKE(:,:,LM) = 1e-6
        else
          TKE = 1e-6 ! https://github.com/GEOS-ESM/GEOSgcm_GridComp/issues/594#issuecomment-1171360993
          do L = 1,LM-1
            TKE(:,:,L) = ( LAMBDADISS * &
            ( -1.*(KH(:,:,L)*MAPL_GRAV/((TSM(:,:,L) + TSM(:,:,L+1))*0.5)  *  ((TSM(:,:,L) - TSM(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))) +  &
            (KM(:,:,L)*((U(:,:,L) - U(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))*((U(:,:,L) - U(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1))))  +  &
            (KM(:,:,L)*((V(:,:,L) - V(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))*((V(:,:,L) - V(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))) )) ** 2
            TKE(:,:,L) = TKE(:,:,L) ** (1./3.)
          enddo
          TKE = max(1e-6, TKE) ! https://github.com/GEOS-ESM/GEOSgcm_GridComp/issues/594#issuecomment-1171360993

          ! If not running SHOC, estimate ISOTROPY from KH and TKE,
          ! based on Eq. 7 from Bogenschutz and Krueger (2013).
          ! This is a placeholder to allow use of the double-gaussian
          ! PDF without SHOC, but should be tested and revised!
          ISOTROPY(:,:,LM) = KH(:,:,LM-1) / max(0.01,0.1*TKE(:,:,LM-1))
          ISOTROPY(:,:,1) = KH(:,:,1) / max(0.01,0.1*TKE(:,:,1))
          do L = 2,LM-1
            ISOTROPY(:,:,L) = (KH(:,:,L)+KH(:,:,L-1)) / (0.1*(TKE(:,:,L)+TKE(:,:,L-1)))
          end do
          ISOTROPY = max(10.,min(2000.,ISOTROPY))

        end if
      end if ! TKE

      ! Update the higher order moments required for the ADG PDF
      if ( (PDFSHAPE.ge.5) ) then
      SL = T + (MAPL_GRAV*Z - MAPL_ALHL*QLTOT - MAPL_ALHS*QITOT)/MAPL_CP
      call update_moments(IM, JM, LM, DT, &
                          SH,             &  ! in
                          EVAP,           &
                          AREA,           &
                          ZPBL,           &
                          Z,              &
                          ZLE,            &
                          KH,             &
                          BRUNTSHOC,      &
                          TKESHOC,        &
                          ISOTROPY,       &
                          QT,             &
                          SL,             &
                          EDMF_FRC,       &
                          MFQT3,          &
                          MFSL3,          &
                          MFW2,           &
                          MFW3,           &
                          MFWQT,          &
                          MFWSL,          &
                          MFSLQT,         &
                          WQT_DC,         &
                          PDF_A,          &  ! inout
                          qt2,            &
                          qt3,            &
                          sl2,            &  ! out
                          sl3,            &
                          w2,             &
                          w3,             &
                          slqt,           &
                          doprogqt2,      &  ! tuning parameters
                          sl2tune,        &
                          qt2tune,        &
                          slqt2tune,      &
                          skew_tgen,      &
                          skew_tdis,      &
                          qt2_tdis,       &
                          free_atm_sigs )

       end if

      KPBLMIN  = count(PREF < 50000.)

                            ZPBL = MAPL_UNDEF
      if (associated(PPBL)) PPBL = MAPL_UNDEF

      if (CALC_TCZPBL) then
         TCZPBL = MAPL_UNDEF
         thetavs = T(:,:,LM)*(1.0+MAPL_VIREPS*Q(:,:,LM)/(1.0-Q(:,:,LM)))*(TH(:,:,LM)/T(:,:,LM))
         tcrib(:,:,LM) = 0.0
         do I = 1, IM
            do J = 1, JM
               do L=LM-1,1,-1
                  thetavh(I,J) = T(I,J,L)*(1.0+MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))*(TH(I,J,L)/T(I,J,L))
                  uv2h(I,J) = max(U(I,J,L)**2+V(I,J,L)**2,1.0E-8)
                  tcrib(I,J,L) = MAPL_GRAV*(thetavh(I,J)-thetavs(I,J))*Z(I,J,L)/(thetavs(I,J)*uv2h(I,J))
                  if (tcrib(I,J,L) >= tcri_crit) then
                     TCZPBL(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                     KPBLTC(I,J) = float(L)
                     exit
                  end if
               end do
            end do
         end do
         where (TCZPBL<0.)
            TCZPBL = Z(:,:,LM)
            KPBLTC = float(LM)
         end where
      end if ! CALC_TCZPBL

      if (CALC_ZPBL2) then
         ZPBL2 = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM
               do L=LM,2,-1
                  if ((KH(I,J,L-1) < 2.).and.(KH(I,J,L) >= 2.).and.(ZPBL2(I,J)==MAPL_UNDEF)) then
                     ZPBL2(I,J) = Z(I,J,L)
                     KPBL2(I,J) = float(L)
                  end if
               end do
            end do
         end do

         where ( ZPBL2 .eq. MAPL_UNDEF )
            ZPBL2 = Z(:,:,LM)
            KPBL2 = float(LM)
         end where
         ZPBL2 = MIN(ZPBL2,Z(:,:,KPBLMIN))
      end if ! CALC_ZPBL2

      if (CALC_ZPBL10p) then
         ZPBL10p = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM
               temparray(1:LM+1) = KH(I,J,0:LM)
               do L = LM,2,-1
                  locmax = maxloc(temparray,1)
                  minlval = max(0.001,0.0001*maxval(temparray))
                  if(temparray(locmax-1)<minlval.and.temparray(locmax+1)<minlval) temparray(locmax) = minlval
               enddo
               maxkh = temparray(LM)
               do L = LM-1,2,-1
                  if(temparray(L)>maxkh) maxkh = temparray(L)
                  if(temparray(L-1)<minlval) exit
               end do
               do L=LM-1,2,-1
                  if ( (temparray(L) < 0.1*maxkh) .and. (temparray(L+1) >= 0.1*maxkh)  &
                  .and. (ZPBL10p(I,J) == MAPL_UNDEF ) ) then
                     ZPBL10p(I,J) = ZL0(I,J,L)+ &
                  ((ZL0(I,J,L-1)-ZL0(I,J,L))/(temparray(L)-temparray(L+1))) * (0.1*maxkh-temparray(L+1))
                     KPBL10p(I,J) = float(L)
                  end if
               end do
               if (  ZPBL10p(I,J) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
                  ZPBL10p(I,J) = Z(I,J,LM)
                  KPBL10p(I,J) = float(LM)
               endif
            end do
         end do

         ZPBL10p = MIN(ZPBL10p,Z(:,:,KPBLMIN))
      end if ! CALC_ZPBL10p

      ! HTKE pbl height
      if (associated(ZPBLHTKE)) then
         ZPBLHTKE = MAPL_UNDEF
      end if ! ZPBLHTKE

      if (CALC_ZPBLRI) then
         ZPBLRI = MAPL_UNDEF
         
         thetavs = T(:,:,LM)*(1.0+MAPL_VIREPS*Q(:,:,LM)/(1.0-Q(:,:,LM)))*(TH(:,:,LM)/T(:,:,LM))
         tcrib(:,:,LM) = 0.0
         
         do I = 1, IM
            do J = 1, JM
               
               ! -----------------------------------------------------------------------
               ! Step 1: Identify boundary layer regime
               ! -----------------------------------------------------------------------
               is_marine_stratocumulus = .false.
               is_shallow_cumulus = .false.
               cloud_base_idx = -1
               cloud_top_idx = -1
               cloud_fraction_mean = 0.0
               cloud_thickness = 0.0

               ! Scan for low clouds (below 3 km only - ignore mid/high clouds)
               do L = LM-1, 1, -1
                  if (Z(I,J,L) > 3000.0) exit  ! Stop searching above 3 km
                  if (FCLD(I,J,L) > 0.05) then
                     if (cloud_base_idx < 0) cloud_base_idx = L  ! First cloud level from bottom
                     cloud_top_idx = L  ! Keep updating as we go up
                     cloud_fraction_mean = cloud_fraction_mean + FCLD(I,J,L)
                     cloud_thickness = cloud_thickness + 1.0
                  else
                     if (cloud_base_idx > 0) exit  ! Found cloud top, stop
                  end if
               end do

               ! Classify cloud regime (only if low clouds exist)
               if (cloud_base_idx > 0 .and. cloud_thickness > 0.0) then
                  cloud_fraction_mean = cloud_fraction_mean / cloud_thickness
                  cloud_thickness = Z(I,J,cloud_top_idx) - Z(I,J,cloud_base_idx)
                  ! Marine stratocumulus: solid deck, moderate thickness
                  if (cloud_fraction_mean > 0.6 .and. cloud_thickness > 200.0 .and. &
                      Z(I,J,cloud_base_idx) < 2000.0) then
                     is_marine_stratocumulus = .true.
                  ! Shallow cumulus: broken clouds, thin
                  else if (cloud_fraction_mean < 0.4 .and. cloud_thickness < 800.0 .and. &
                           Z(I,J,cloud_base_idx) < 1500.0) then
                     is_shallow_cumulus = .true.
                  end if
               end if

               ! -----------------------------------------------------------------------
               ! Step 2: Compute Richardson number profile
               ! -----------------------------------------------------------------------
               do L=LM-1,1,-1
                  thetavh(I,J) = T(I,J,L)*(1.0+MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))*(TH(I,J,L)/T(I,J,L))
                  uv2h(I,J) = max(U(I,J,L)**2+V(I,J,L)**2, 1.0E-8)
                  tcrib(I,J,L) = MAPL_GRAV*(thetavh(I,J)-thetavs(I,J))*Z(I,J,L)/(thetavs(I,J)*uv2h(I,J))
               end do

               ! -----------------------------------------------------------------------
               ! Step 3: Apply regime-specific PBL height diagnostic
               ! -----------------------------------------------------------------------
              
               ! REGIME 1: Marine Stratocumulus
               ! Place PBL top at cloud top (where strong inversion exists)
               if (is_marine_stratocumulus) then
                  ! Find the inversion at/above cloud top
                  do L = cloud_top_idx, 1, -1
                     if (tcrib(I,J,L) >= tcri_crit) then
                        ZPBLRI(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                        KPBLRI(I,J) = float(L)
                        exit
                     end if
                  end do
                  ! Fallback: if no inversion found above cloud, use cloud top
                  if (ZPBLRI(I,J) .eq. MAPL_UNDEF) then
                     ZPBLRI(I,J) = Z(I,J,cloud_top_idx)
                     KPBLRI(I,J) = float(cloud_top_idx)
                  end if
              
               ! REGIME 2: Shallow Cumulus
               ! Place PBL top at cloud base (sub-cloud mixed layer)
               else if (is_shallow_cumulus) then
                  ZPBLRI(I,J) = Z(I,J,cloud_base_idx)
                  KPBLRI(I,J) = float(cloud_base_idx)
              
               ! REGIME 3: Clear or Non-Boundary-Layer Clouds
               ! Use standard bulk Richardson criterion
               else
                  do L=LM-1,1,-1
                     if (tcrib(I,J,L) >= tcri_crit) then
                        ZPBLRI(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                        KPBLRI(I,J) = float(L)
                        exit
                     end if
                  end do
               end if

            end do
         end do 
         
         where (ZPBLRI<0.)
            ZPBLRI = Z(:,:,LM)
            KPBLRI = float(LM)
         end where
      end if

      ! thetav gradient based pbl height diagnostic
      if (associated(ZPBLTHV)) then
         ZPBLTHV = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM

               do L=LM,1,-1
                  thetav(L) = TH(I,J,L)*(1.0+MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))
               end do

               maxdthvdz = 0

               do L=LM-1,1,-1
                  if(Z(I,J,L)<=Z(I,J,KPBLMIN)) then
                     dthvdz = (thetav(L+1)-thetav(L))/(Z(I,J,L+1)-Z(I,J,L))
                     if(dthvdz>maxdthvdz) then
                        maxdthvdz = dthvdz
                        ZPBLTHV(I,J) = 0.5*(Z(I,J,L+1)+Z(I,J,L))
                     end if
                  end if
               end do

            end do 
         end do 
      end if ! ZPBLTHV

!=========================================================================                                      
!  ZPBL defined by minimum in vertical gradient of refractivity.                                                
!  As shown in Ao, et al, 2012: "Planetary boundary layer heights from                                          
!  GPS radio occultation refractivity and humidity profiles", Climate and                                       
!  Dynamics.  https://doi.org/10.1029/2012JD017598                                                              
!=========================================================================                                      
    if (associated(ZPBLRFRCT)) then

      a1 = 0.776    ! K/Pa                                                                                      
      a2 = 3.73e3   ! K2/Pa                                                                                     

      WVP = Q * PLO / (Q*(1.-0.622)+0.622)  ! water vapor partial pressure                                      

      ! Pressure gradient term                                                                                  
      dum3d(:,:,2:LM-1) = (PLO(:,:,1:LM-2)-PLO(:,:,3:LM)) / (Z(:,:,1:LM-2)-Z(:,:,3:LM))
      dum3d(:,:,1) = (PLO(:,:,1)-PLO(:,:,2)) / (Z(:,:,1)-Z(:,:,2))
      dum3d(:,:,LM) = (PLO(:,:,LM-1)-PLO(:,:,LM)) / (Z(:,:,LM-1)-Z(:,:,LM))
      tmp3d = a1 * dum3d / T

      ! Add Temperature gradient term                                                                           
      dum3d(:,:,2:LM-1) = (T(:,:,1:LM-2)-T(:,:,3:LM)) / (Z(:,:,1:LM-2)-Z(:,:,3:LM))
      dum3d(:,:,1) = (T(:,:,1)-T(:,:,2)) / (Z(:,:,1)-Z(:,:,2))
      dum3d(:,:,LM) = (T(:,:,LM-1)-T(:,:,LM)) / (Z(:,:,LM-1)-Z(:,:,LM))
      tmp3d = tmp3d - (a1*plo/T**2 + 2.*a2*WVP/T**3)*dum3d

      ! Add vapor pressure gradient term                                                                        
      dum3d(:,:,2:LM-1) = (WVP(:,:,1:LM-2)-WVP(:,:,3:LM)) / (Z(:,:,1:LM-2)-Z(:,:,3:LM))
      dum3d(:,:,1) = (WVP(:,:,1)-WVP(:,:,2)) / (Z(:,:,1)-Z(:,:,2))
      dum3d(:,:,LM) = (WVP(:,:,LM-1)-WVP(:,:,LM)) / (Z(:,:,LM-1)-Z(:,:,LM))
      tmp3d = tmp3d + (a2/T**2)*dum3d

      ! ZPBL is height of minimum in refractivity (tmp3d)                                                       
      do I = 1,IM
        do J = 1,JM
          K = MINLOC(tmp3d(I,J,:),DIM=1,BACK=.TRUE.)   ! return last index, if multiple                         
          ZPBLRFRCT(I,J) = Z(I,J,K)
        end do
      end do

    end if  ! ZPBLRFRCT 


      ! PBL height diagnostic based on specific humidity gradient
      ! PBLH defined as level with minimum QV gradient
      if (associated(ZPBLQV)) then
         ZPBLQV = MAPL_UNDEF

         do I = 1, IM
            do J = 1, JM

               maxdthvdz = 0  ! re-using variables from ZPBLTHV calc above

               do L=LM-1,1,-1
                  if(Z(I,J,L)<=Z(I,J,KPBLMIN)) then
                     dthvdz = -1.*(Q(I,J,L+1)-Q(I,J,L))/(Z(I,J,L+1)-Z(I,J,L))
                     if(dthvdz>maxdthvdz) then
                        maxdthvdz = dthvdz
                        ZPBLQV(I,J) = 0.5*(Z(I,J,L+1)+Z(I,J,L))
                     end if
                  end if
               end do

            end do 
         end do 
      end if ! ZPBLQV


     if (associated(SBITOP) .or. associated(SBIFRQ) ) then

        SBIFRQ = 0.
        SBITOP = MAPL_UNDEF

        do I = 1, IM
           do J = 1, JM
              if (T(I,J,LM-1).gt.T(I,J,LM)) then
                 SBIFRQ(I,J) = 1.
                 do L = LM-1,1,-1
                    if (T(I,J,L).gt.T(I,J,L+1)) then
                       SBITOP(I,J) = Z(I,J,L)
                    else
                       exit
                    end if
                 end do
              end if
           end do
        end do

     end if ! SBITOP, SBIFRQ

     ! Trade inversion base height
     if (associated(TRINVBS)) then
        TRINVBS = MAPL_UNDEF
        TRINVDELT = MAPL_UNDEF
        TRINVFRQ = 0.
        do I = 1,IM
           do J = 1,JM
              K = LM

              do while (PLO(I,J,K).gt.95000.)
                 K = K-1
              end do
              do L = K,1,-1    ! K is first level above 950mb
                 if (PLO(I,J,L).lt.60000.) exit
                 
                 if (T(I,J,L-1).ge.T(I,J,L)) then ! if next level is warmer...
                    LTOP = L                      ! L is index of minimum T so far
                    do while (T(I,J,LTOP).ge.T(I,J,L)) ! find depth of warm layer
                       LTOP = LTOP-1
                    end do
                    LTOP = LTOP+1   ! LTOP is index of highest level inside warm layer

                    if (  MAXVAL(T(I,J,LTOP:L))-T(I,J,L).ge.0.5 .or. &
                         (MAXVAL(T(I,J,LTOP:L))-T(I,J,L).gt.0.01 .and. PLO(I,J,L)-PLO(I,J,LTOP)>2500.) ) then

                       ! only save if DELTA-T exceeds any previous inversion
                       if ( TRINVFRQ(I,J).eq.0. .or. &
                            (TRINVFRQ(I,J).ne.0. .and. MAXVAL(T(I,J,LTOP:L))-T(I,J,L).gt.TRINVDELT(I,J)) ) then
                          TRINVBS(I,J)   = PLO(I,J,L)
                          TRINVDELT(I,J) = MAXVAL(T(I,J,LTOP:L))-T(I,J,L)
                          TRINVFRQ(I,J)  = 1.
                       end if

                    end if
                 end if ! next level warmer

              end do ! L vert loop

           end do
        end do
     end if

      SELECT CASE(PBLHT_OPTION)

      CASE( 1 )
         ZPBL = ZPBL2
         KPBL = KPBL2

      CASE( 2 )
         ZPBL = ZPBL10p
         KPBL = KPBL10P

      CASE( 3 )
         ZPBL = TCZPBL
         KPBL = KPBLTC

      CASE( 4 )
         WHERE (FRLAND(:,:)>0)
            ZPBL = TCZPBL
            KPBL = KPBLTC

         ELSEWHERE
            ZPBL = ZPBL10p
            KPBL = KPBL10P

         END WHERE

      CASE( 5 )
         ZPBL = ZPBLRI
         KPBL = KPBLRI

      END SELECT

      ZPBL = MIN(ZPBL,Z(:,:,KPBLMIN))
      KPBL = MAX(KPBL,float(KPBLMIN))
  
     ! Calc KPBL using surface turbulence, for use in UW shallow scheme
      if (associated(KPBL_SC)) then
        KPBL_SC = MAPL_UNDEF
        do I = 1, IM
          do J = 1, JM
            if ( JASON_PBL_SC ) then
                ! ----------------------------------------------------------------
                ! Use old 10% of HKHSFC
                ! ----------------------------------------------------------------
                temparray(1:LM+1) = KHSFC(I,J,0:LM)
                maxkh = maxval(temparray)
                kh_thresh = 0.1
                do L = LM-1, 2, -1
                  if ( (temparray(L) < kh_thresh*maxkh) .and. (temparray(L+1) >= kh_thresh*maxkh)  &
                  .and. (KPBL_SC(I,J) == MAPL_UNDEF ) ) then
                     KPBL_SC(I,J) = float(L)
                  end if
                end do
                ! =================================================================
                ! ROBUST FALLBACK FOR CALM / LAMINAR CONDITIONS
                ! =================================================================
                if ( KPBL_SC(I,J) == MAPL_UNDEF .or. maxkh < 1.0 ) then
                  KPBL_SC(I,J) = float(LM)
                endif
            else
                ! ----------------------------------------------------------------
                ! Use KPBL from PBLHT_OPTION
                ! ----------------------------------------------------------------
                KPBL_SC = KPBL
            endif
          end do
        end do
      endif
      ! =================================================================
      ! POST-PROCESSING LOOP (COHERENT LAYER SNAP)
      ! =================================================================
      if (associated(KPBL_SC) .and. associated(ZPBL_SC)) then
        do I = 1, IM
          do J = 1, JM
             ZPBL_SC(I,J) = Z(I,J, int(KPBL_SC(I,J)))
          end do
        end do
      endif

      if (associated(PPBL)) then
         do I = 1, IM
            do J = 1, JM
               PPBL(I,J) = PLO(I,J,nint(KPBL(I,J)))
            end do
         end do
         PPBL = MAX(PPBL,PLO(:,:,KPBLMIN))
      end if

      RHOAW3=RHOE*AW3  ! mass flux (positive)

      ! Second difference coefficients for scalars; RDZ is RHO/DZ, DMI is (G DT)/DP
      ! ---------------------------------------------------------------------------

      CKS(:,:,1:LM-1) = -(KH(:,:,1:LM-1)+MFPARAMS%TREFF*RHOAW3(:,:,1:LM-1)) * RDZ(:,:,1:LM-1)
      AKS(:,:,1     ) = 0.0
      AKS(:,:,2:LM  ) = CKS(:,:,1:LM-1) * DMI(:,:,2:LM  )
      CKS(:,:,1:LM-1) = CKS(:,:,1:LM-1) * DMI(:,:,1:LM-1)
      CKS(:,:,  LM  ) = -CT             * DMI(:,:,  LM  )

      ! Fill KH at level LM+1 with CT / RDZ for diagnostic output
      ! ---------------------------------------------------------

      KH(:,:,LM) = CT * Z(:,:,LM)*((MAPL_RGAS * TV(:,:,LM))/PLE(:,:,LM))
      TKH = KH

      ! Water vapor exchange coefficient can differ at the surface
      !-----------------------------------------------------------

      AKQ         = AKS
      CKQ         = CKS
      CKQ(:,:,LM) = -CQ * DMI(:,:,LM)

      ! Second difference coefficients for winds
      ! EKV is saved to use in the frictional heating calc.
      ! ---------------------------------------------------
      
      EKV(:,:,1:LM-1) = -KM(:,:,1:LM-1) * RDZ(:,:,1:LM-1)
      AKV(:,:,1     ) = 0.0
      AKV(:,:,2:LM  ) = EKV(:,:,1:LM-1) * DMI(:,:,2:LM  )
      CKV(:,:,1:LM-1) = EKV(:,:,1:LM-1) * DMI(:,:,1:LM-1)
      EKV(:,:,1:LM-1) = -MAPL_GRAV      * EKV(:,:,1:LM-1)

      CKV(:,:,  LM  ) = -  CU           * DMI(:,:,  LM  )
      EKV(:,:,  LM  ) =  MAPL_GRAV      * CU

      ! Fill KM at level LM with CU / RDZ for diagnostic output
      ! -------------------------------------------------------

      KM(:,:,LM) = CU * Z(:,:,LM)*(MAPL_RGAS * TV(:,:,LM)) / PLE(:,:,LM)   !Z/rho

      ! Setup the tridiagonal matrix
      ! ----------------------------

      BKS = 1.00 - (AKS+CKS)
      BKQ = 1.00 - (AKQ+CKQ)
      BKV = 1.00 - (AKV+CKV)

    !
    ! A,B,C,D-s for mass flux
    !
        
     AKSS(:,:,1)=0.0
     AKUU(:,:,1)=0.0

     RHOAW3=RHOE*AW3

     if (MFPARAMS%IMPLICIT == 1 .and. MFPARAMS%DISCRETE == 0) then
        AKSS(:,:,2:LM) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM) &
                         - 0.5*DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
        AKUU(:,:,2:LM) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM) &
                         - 0.5*DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
     else
        AKSS(:,:,2:LM) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM)
        AKUU(:,:,2:LM) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM)
     end if
     AKQQ = AKSS

     CKSS(:,:,LM)=-CT*DMI(:,:,LM)
     CKQQ(:,:,LM)=-CQ*DMI(:,:,LM)
     CKUU(:,:,LM)=-CU*DMI(:,:,LM)
  
     if (MFPARAMS%IMPLICIT == 1 .and. MFPARAMS%DISCRETE == 0) then
        CKSS(:,:,1:LM-1) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1) &
                           + 0.5*DMI(:,:,1:LM-1)*RHOAW3(:,:,1:LM-1)
        CKUU(:,:,1:LM-1) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1) &
                           + 0.5*DMI(:,:,1:LM-1)*RHOAW3(:,:,1:LM-1)
     else
        CKSS(:,:,1:LM-1) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1)
        CKUU(:,:,1:LM-1) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1)
     end if
     CKQQ(:,:,1:LM-1) = CKSS(:,:,1:LM-1)  
 
     BKSS = 1.0 - (CKSS+AKSS)
     BKQQ = 1.0 - (CKQQ+AKQQ)
     BKUU = 1.0 - (CKUU+AKUU)

! Add mass flux contribution
  
  if (MFPARAMS%IMPLICIT == 1) then
     if (MFPARAMS%DISCRETE == 0) then
        BKSS(:,:,LM) = BKSS(:,:,LM) - DMI(:,:,LM)*RHOAW3(:,:,LM-1)
        BKQQ(:,:,LM) = BKQQ(:,:,LM) - DMI(:,:,LM)*RHOAW3(:,:,LM-1)
        BKUU(:,:,LM) = BKUU(:,:,LM) - DMI(:,:,LM)*RHOAW3(:,:,LM-1)

        BKSS(:,:,1:LM-1) = BKSS(:,:,1:LM-1) + DMI(:,:,1:LM-1)*( RHOAW3(:,:,1:LM-1) - RHOAW3(:,:,0:LM-2) )
        BKQQ(:,:,1:LM-1) = BKQQ(:,:,1:LM-1) + DMI(:,:,1:LM-1)*( RHOAW3(:,:,1:LM-1) - RHOAW3(:,:,0:LM-2) )
        BKUU(:,:,1:LM-1) = BKUU(:,:,1:LM-1) + DMI(:,:,1:LM-1)*( RHOAW3(:,:,1:LM-1) - RHOAW3(:,:,0:LM-2) ) 
     else if (MFPARAMS%DISCRETE == 1) then
        AKSS(:,:,2:LM) = AKSS(:,:,2:LM) - DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
        AKQQ(:,:,2:LM) = AKQQ(:,:,2:LM) - DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
        AKUU(:,:,2:LM) = AKUU(:,:,2:LM) - DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)

        BKSS(:,:,2:LM-1) = BKSS(:,:,2:LM-1) + DMI(:,:,2:LM-1)*RHOAW3(:,:,2:LM-1)
        BKQQ(:,:,2:LM-1) = BKQQ(:,:,2:LM-1) + DMI(:,:,2:LM-1)*RHOAW3(:,:,2:LM-1)
        BKUU(:,:,2:LM-1) = BKUU(:,:,2:LM-1) + DMI(:,:,2:LM-1)*RHOAW3(:,:,2:LM-1)
     end if
  end if

! Y-s ... these are rhs - mean value - surface flux 
! (these are added in the diffuse and vrtisolve)
   
   ! Add prescribed surface fluxes
   if ( SCM_SL /= 0 .and. (SCM_SL_FLUX == 1 .or. SCM_SL_FLUX == 2) ) then
      YS(:,:,LM)  = YS(:,:,LM)  + DMI(:,:,LM)*SH(:,:)!/RHOE(:,:,LM)
      YQV(:,:,LM) = YQV(:,:,LM) + DMI(:,:,LM)*EVAP(:,:)!/RHOE(:,:,LM)
   end if

      ! Add the topographic roughness term
      ! ----------------------------------

      if (associated(AKSODT)) then
         AKSODT = -AKS/DT
         AKSODT(:,:,1) = 0.0
      end if

      if (associated(CKSODT)) then
         CKSODT = -CKS/DT
         CKSODT(:,:,LM) = 0.0
      end if

      if (associated(AKQODT)) then
         AKQODT = -AKQ/DT
         AKQODT(:,:,1) = 0.0
      end if

      if (associated(CKQODT)) then
         CKQODT = -CKQ/DT
         CKQODT(:,:,LM) = 0.0
      end if

      if (associated(AKVODT)) AKVODT = -AKV/DT
      if (associated(CKVODT)) CKVODT = -CKV/DT

      call MAPL_GridCompTimerStop(GC, "---POSTLOCK")

!BOP
!
!   Orograpghic drag follows  Beljaars (2003):
!   $$
!   \frac{\partial}{\partial z}\frac{\tau}{\rho} = \frac{C_B}{\lambda_B} |U(z)| U(z) 
!          e^{-\tilde{z}^\frac{3}{2}}\tilde{z}^{-1.2},
!   $$
!   where $z$ is the height above the surface in meters, 
!   $\tilde{z}=\frac{z}{\lambda_B}$, $\tau$ is the orographic stress at $z$,
!   $\rho$ is the air density, $U(z)$ is the wind velocity, and $\lambda_B$ is a vertical length scale.
!   Beljaars uses $\lambda_B = 1500$m, for which the non-dimensional parameter $C_B = 2.5101471 \times 10^{-8}$.
!   These are the default values, but both can be modified from the configuration. To avoid underflow.
!   the tendency is set to zero once $\tilde{z}$ exceeds 4 (i.e., 6 km from the surface for default values). 
!
!EOP

      call MAPL_GridCompTimerStart(GC, "---BELJAARS")
      if (C_B /= 0.0) then
      call BELJAARS(IM, JM, LM, DT, &
                    LAMBDA_B, C_B,  &
                    KPBL,           &
                    U, V, Z, AREA,  &
                    VARFLT, PLE,    &
                    BKV, BKUU, FKV  )
      endif
      call MAPL_GridCompTimerStop(GC, "---BELJAARS")

      call MAPL_GridCompTimerStart(GC, "---DECOMP")

! Do LU decomposition; C is not modified.
! On exit, B is the main diagonals of the LU
! decomposition, and A is the r.h.s multiplier.
!----------------------------------------------

      AKX = AKS
      BKX = BKS
      call VTRILU(AKX,BKX,CKS)
      AKS = AKX
      BKS = BKX

      AKX = AKQ
      BKX = BKQ
      call VTRILU(AKX,BKX,CKQ)
      AKQ = AKX
      BKQ = BKX

      AKX = AKV
      BKX = BKV
      call VTRILU(AKX,BKX,CKV)
      AKV = AKX
      BKV = BKX

! Get the sensitivity of solution to a unit
! change in the surface value. B and C are
! not modified.
!------------------------------------------

      call VTRISOLVESURF(BKS,CKS,DKS)
      call VTRISOLVESURF(BKQ,CKQ,DKQ)
      call VTRISOLVESURF(BKV,CKV,DKV)

 !
 ! LU decomposition for the mass-flux variables
 !     
      AKX=AKSS
      BKX=BKSS
      call VTRILU(AKX,BKX,CKSS)
      BKSS=BKX
      AKSS=AKX
     
      AKX=AKQQ
      BKX=BKQQ
      call VTRILU(AKX,BKX,CKQQ)
      BKQQ=BKX
      AKQQ=AKX  

      AKX=AKUU
      BKX=BKUU
      call VTRILU(AKX,BKX,CKUU)
      BKUU=BKX
      AKUU=AKX  

      call MAPL_GridCompTimerStop(GC, "---DECOMP")

      if(ALLOC_TCZPBL) deallocate(TCZPBL)
      if(ALLOC_ZPBL2) deallocate(ZPBL2)
      if(ALLOC_ZPBL10p) deallocate(ZPBL10p)
      if(ALLOC_ZPBLRI) deallocate(ZPBLRI)

      RETURN_(ESMF_SUCCESS)
     end subroutine REFRESH

!=============================================================================
!=============================================================================

!BOP

! !CROUTINE: DIFFUSE -- Solves for semi-implicit diffusive tendencies assuming fixed surface conditions.  

! !INTERFACE:

  subroutine DIFFUSE(IM,JM,LM,RC)

! !ARGUMENTS:

    integer,           intent(IN)       :: IM,JM,LM
    integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: {\tt DIFFUSE} computes semi-implicit tendencies of all fields in
!  the TR bundle. Each field is examined for three attributes: {\tt DiffuseLike},
!  {\tt FriendlyToTURBULENCE}, and {\tt WeightedTendency}. These determine the behavior of
!  {\tt DIFFUSE} for that field. {\tt DiffuseLike} can be either 'U', 'Q', or 'S'; the default is 'Q'.
!  {\tt FriendlyToTURBULENCE}, and {\tt WeightedTendency} are ESMF logicals.
!  If {\tt FriendlyToTURBULENCE} is true, the field in TR is updated directly; otherwise
!  it is left untouched. In either case, If the corresponding pointer TRI bundle is associated, the
!  tendencies are returned there. If {\tt WeightedTendency} is true, the tendency in TRI, if any,
!  is pressure weighted.

!EOP

    character(len=ESMF_MAXSTR)          :: IAm='Diffuse'
    integer                             :: STATUS

    character(len=ESMF_MAXSTR)          :: TYPE
    character(len=ESMF_MAXSTR)          :: NAME
    type (ESMF_Field)                   :: FIELD
    type (ESMF_Array)                   :: ARRAY
    type (ESMF_FieldBundle)             :: TR
    type (ESMF_FieldBundle)             :: TRI
    real, dimension(:,:,:), pointer     :: S, SOI, SOD
    real, dimension(:,:),   pointer     :: SG, SF, SDF, CX, SRG
    real, dimension(:,:,:), pointer     :: DX
    real, dimension(:,:,:), pointer     :: AK, BK, CK
    type(ESMF_Info) :: infoh

    integer                             :: KM, K,L
    logical                             :: FRIENDLY
    logical                             :: WEIGHTED
    logical                             :: OPT      ! selects algorithm in VTRISOLVE

    real,               dimension(IM,JM,LM) :: DP,DZ
    real(kind=MAPL_R8), dimension(IM,JM,LM) :: SX

    real :: DOMF

    integer :: i, j, ll

    ! Parameters for idealized SCM surface layer
    integer :: SCM_SL, SCM_SL_FLUX
    real    :: SCM_SH, SCM_LH

    ! EDMF transport
    real    :: EntExp, EntDyn
    real, dimension(LM) :: UPSX
    real, dimension(:,:,:), pointer :: edmf_mf, edmf_entx
    
    ! pointers to exports after diffuse
    real, dimension(:,:,:), pointer     :: UAFDIFFUSE, VAFDIFFUSE, SAFDIFFUSE, QAFDIFFUSE

    real, dimension(:,:),   pointer     :: SHOBS, LHOBS

! Sea Spray
    real, dimension(:,:), pointer       :: SH_SPRAY_ => NULL()
    real, dimension(:,:), pointer       :: LH_SPRAY_ => NULL()
    real, dimension(IM,JM)              :: SH_SPRAY
    real, dimension(IM,JM)              :: LH_SPRAY

    real, parameter :: SH_SPRAY_MIN = -500.0
    real, parameter :: SH_SPRAY_MAX =  500.0
    real, parameter :: LH_SPRAY_MIN = -500.0
    real, parameter :: LH_SPRAY_MAX =  500.0


    ! Get info for idealized SCM surface layer
    call MAPL_GridCompGetResource(GC, 'SCM_SL', SCM_SL, default=0, RC=STATUS)
    VERIFY_(STATUS)

    ! Prescribed surface exchange coefficients
    if ( SCM_SL /= 0 ) then
       call MAPL_GridCompGetResource(GC, 'SCM_SL_FLUX', SCM_SL_FLUX, default=0, RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, 'SCM_SH', SCM_SH, default=0., RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GridCompGetResource(GC, 'SCM_LH', SCM_LH, default=0., RC=STATUS)
       VERIFY_(STATUS)

       CU => cu_scm
       CT => ct_scm
       CQ => ct_scm

       call MAPL_StateGetPointer(IMPORT, SHOBS,'SHOBS', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_StateGetPointer(IMPORT, LHOBS,'LHOBS', RC=STATUS)
       VERIFY_(STATUS)
    end if



! Get the bundles containing the quantities to be diffused, 
!     their tendencies, their surface values, their surface
!     fluxes, and the derivatives of their surface fluxes
!     wrt the surface values. 
!----------------------------------------------------------

    call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)

    if (DO_WAVES/=0 .and. DO_SEA_SPRAY/=0) then
       call MAPL_StateGetPointer(IMPORT, SH_SPRAY_, 'SHFX_SPRAY',   RC=STATUS)
       VERIFY_(STATUS)

       call MAPL_StateGetPointer(IMPORT, LH_SPRAY_, 'LHFX_SPRAY',   RC=STATUS)
       VERIFY_(STATUS)

       SH_SPRAY = SH_SPRAY_
       LH_SPRAY = LH_SPRAY_

       where (SH_SPRAY < SH_SPRAY_MIN)  SH_SPRAY = SH_SPRAY_MIN
       where (SH_SPRAY > SH_SPRAY_MAX)  SH_SPRAY = SH_SPRAY_MAX

       where (LH_SPRAY < LH_SPRAY_MIN)  LH_SPRAY = LH_SPRAY_MIN
       where (LH_SPRAY > LH_SPRAY_MAX)  LH_SPRAY = LH_SPRAY_MAX
    end if

    call ESMF_StateGet(EXPORT, 'TRI',    TRI,    RC=STATUS); VERIFY_(STATUS)

! Get pointers to exports of U,V and S that diffuse sees
!  Required for SYNCTQ (ALLOC=.TRUE.)
    call MAPL_StateGetPointer(EXPORT, UAFDIFFUSE ,  'UAFDIFFUSE' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_StateGetPointer(EXPORT, VAFDIFFUSE ,  'VAFDIFFUSE' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_StateGetPointer(EXPORT, SAFDIFFUSE ,  'SAFDIFFUSE' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_StateGetPointer(EXPORT, QAFDIFFUSE ,  'QAFDIFFUSE' , RC=STATUS); VERIFY_(STATUS)


    call MAPL_StateGetPointer(EXPORT,  edmf_mf,    'EDMF_MF', RC=STATUS); VERIFY_(STATUS)
    call MAPL_StateGetPointer(EXPORT,  edmf_entx,  'EDMF_ENTR', RC=STATUS); VERIFY_(STATUS)

    
! Count the firlds in TR...
!--------------------------

    call ESMF_FieldBundleGet(TR, fieldCOUNT=KM, RC=STATUS)
    VERIFY_(STATUS)

! ...and make sure the other bundles are the same.
!-------------------------------------------------

    call ESMF_FieldBundleGet(TRI,    FieldCount=K , RC=STATUS)
    VERIFY_(STATUS)
    _ASSERT(KM==K,'needs informative message')

! Pressure thickness of layers
!-----------------------------

    DP = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)
    DZ = ZLE(:,:,0:LM-1)-ZLE(:,:,1:LM)

! Loop over all quantities to be diffused.
!----------------------------------------

    do K=1,KM

! Get the Kth Field and its name from tracer bundle
!--------------------------------------------------

       call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
       VERIFY_(STATUS)

       call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
       VERIFY_(STATUS)

! Get item's diffusion type (U, S or Q; default is Q)
!----------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="DiffuseLike",         &
            VALUE=TYPE, DEFAULTVALUE=dflt_q,    RC=STATUS)
       VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,key='FriendlyToTURBULENCE',value=FRIENDLY,default=dflt_false,RC=STATUS)
       VERIFY_(STATUS)

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

       call ESMF_AttributeGet(FIELD, NAME="WeightedTendency",   &
            VALUE=WEIGHTED, DEFAULTVALUE=dflt_false,    RC=STATUS)
       VERIFY_(STATUS)

! Get pointer to the quantity, its tendency, its surface value,
!   the surface flux, and the sensitivity of the surface flux.
! -------------------------------------------------------------

       call TURB_BundleGetPointerToData(TR    ,      NAME,         S  , RC=STATUS)
       VERIFY_(STATUS)
       call TURB_BundleGetPointerToData(TRI   , trim(NAME)//'IT' , SOI, RC=STATUS)
       VERIFY_(STATUS)

! Surface values, fluxes and sensitivities exist only for S, U, V and Q;
!  all other quantities are left unassociated (zero flux).
!------------------------------------------------------------------
       nullify(SRG, SF, SDF)
       select case (trim(NAME))
       case ('S','U','V','Q')
          call MAPL_StateGetPointer(IMPORT, SRG, trim(NAME)//'HAT', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_StateGetPointer(EXPORT, SF , trim(NAME)//'FLX', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_StateGetPointer(EXPORT, SDF, trim(NAME)//'DFL', RC=STATUS)
          VERIFY_(STATUS)
       end select

! The quantity must exist; others are optional.
!----------------------------------------------

       _ASSERT(associated(S ),'needs informative message')

! If the surface values does not exists, we assume zero flux.
!------------------------------------------------------------
       
       if(associated(SRG)) then
          SG => SRG
       else
          allocate (SG(0,0), stat=STATUS)
          VERIFY_(STATUS)
       end if

       ! Add presribed fluxes
       if ( SCM_SL /= 0 .and. (SCM_SL_FLUX /= 1 .and. SCM_SL_FLUX /= 2) ) then
          if ( trim(name) == 'S' ) then
             SG => ssurf_scm
          end if
          if ( trim(name) == 'Q' ) then
             SG => qsurf_scm
          end if
       end if


       
! Pick the right exchange coefficients
!-------------------------------------

OPT = .TRUE.

if ( (trim(name) /= 'S'   ) .and. (trim(name) /= 'Q'   ) .and. &
     (trim(name) /= 'QLLS') .and. (trim(name) /= 'QILS') .and. &
     (trim(name) /= 'U'   ) .and. (trim(name) /= 'V'   )) then
    

       if     ( TYPE=='U' ) then ! Momentum
          CX => CU
          DX => DKV
          AK => AKV; BK => BKV; CK => CKV
       else if( TYPE=='Q' ) then ! Water Vapor or other tracers
          CX => CQ
          DX => DKQ
          AK => AKQ; BK => BKQ; CK => CKQ
       else if( TYPE=='S' ) then ! Heat
          CX => CT
          DX => DKS
          AK => AKS; BK => BKS; CK => CKS
       else
          RETURN_(ESMF_FAILURE)
       endif

! Copy diffused quantity to temp buffer
! ------------------------------------------
       
       SX = S

       ! Calculate EDMF tracer transport
       if (MFPARAMS%DOTRACERS) then
         do I=1,IM
           do J=1,JM
             if (edmf_mf(I,J,LM-1).gt.1e-8) then              
               UPSX(:)   = 0.
               UPSX(LM-1) = SX(I,J,LM)
               L = LM-2
               do while (edmf_mf(I,J,L).gt.1e-8 .and. L.gt.1)
                 entdyn  = max(0.,edmf_mf(I,J,L)-edmf_mf(I,J,L+1))/(edmf_mf(I,J,L+1)*DZ(I,J,L+1)) ! dynamical entrainment
                 entexp  = exp(-(entdyn+EDMF_ENTX(I,J,L+1))*DZ(I,J,L+1)) 

                 ! Effect of mixing on tracers in updraft
                 UPSX(L)  = SX(I,J,L+1)*(1.-entexp)+UPSX(L+1)*entexp
                 L = L-1
               end do
               do ll = 1,2    ! substep
                 SX(I,J,2:LM) = SX(I,J,2:LM) + 0.5*(DT*MAPL_GRAV/DP(I,J,2:LM))* &
                                      ( edmf_mf(I,J,2:LM)   * (UPSX(2:LM)   - SX(I,J,2:LM) )  &
                                       -edmf_mf(I,J,1:LM-1) * (UPSX(1:LM-1) - SX(I,J,1:LM-1) ) )
               end do
             end if
           end do ! JM
        end do ! IM
        SX = max( 0., SX )  ! prevent negative values from roundoff
       end if
       
 elseif (trim(name) =='S') then
          CX => CT
          DX => DKS
          AK => AKSS; BK => BKSS; CK => CKSS
          SX=S+YS      
 elseif (trim(name)=='Q') then
          CX => CQ
          DX => DKQ
          AK => AKQQ; BK => BKQQ; CK => CKQQ
          SX=S+YQV
 elseif (trim(name)=='QLLS') then
          CX => CQ
          DX => DKQ
          AK => AKQQ; BK => BKQQ; CK => CKQQ
          SX=S+YQL
          if (DO_SHOC/=0.) OPT = .FALSE.
 elseif (trim(name)=='QILS') then
          CX => CQ
          DX => DKQ
          AK => AKQQ; BK => BKQQ; CK => CKQQ
          SX=S+YQI
          if (DO_SHOC/=0.) OPT = .FALSE.
 elseif (trim(name)=='U') then       
         CX => CU
         DX => DKV
         AK => AKUU; BK => BKUU; CK => CKUU
         SX=S+YU
 elseif (trim(name)=='V') then       
         CX => CU
         DX => DKV
         AK => AKUU; BK => BKUU; CK => CKUU
         SX=S+YV        
 end if


! Solve for semi-implicit changes. This modifies SX
! -------------------------------------------------

       call VTRISOLVE(AK,BK,CK,SX,SG,OPT)

! Compute the surface fluxes
!---------------------------

       if(associated(SF)) then
          if ( SCM_SL /= 0 .and. SCM_SL_FLUX == 1 ) then
             if ( trim(name) == 'S' ) then
                SF(:,:) = scm_sh
             elseif ( trim(name) == 'Q' ) then
                SF(:,:) = scm_lh/mapl_alhl
             end if
          else if ( SCM_SL /= 0 .and. SCM_SL_FLUX ==2 ) then
             if ( trim(name) == 'S' ) then
                SF(:,:) = SHOBS 
             elseif ( trim(name) == 'Q' ) then
                SF(:,:) = LHOBS/MAPL_ALHL 
             end if
          else
             if(size(SG)>0) then
                SF = CX*(SG - SX(:,:,LM))
             else
                SF = 0.0
             end if
          end if
       end if

       if (DO_WAVES /= 0 .and. DO_SEA_SPRAY /= 0) then
          if (trim(name) == 'S') then
             SF = SF + SH_SPRAY
          end if

          if (trim(name) == 'Q') then 
             SF = SF + LH_SPRAY/MAPL_ALHL
          end if
       end if

! Create tendencies
!------------------

       if(associated(SOI)) then
          if( WEIGHTED ) then
             SOI = ( (SX - S)/DT )*DP
          else
             SOI = ( (SX - S)/DT )
          endif
       end if

       if (DO_WAVES /= 0 .and. DO_SEA_SPRAY /= 0) then
          if (trim(name) == 'S') then
             SX(:,:,LM) = SX(:,:,LM) + (SH_SPRAY/(DP(:,:,LM)/MAPL_GRAV))*DT
          end if

          if (trim(name) == 'Q') then
             SX(:,:,LM) = SX(:,:,LM) + (LH_SPRAY/(MAPL_ALHL*DP(:,:,LM)/MAPL_GRAV))*DT
          end if
       end if

       if( trim(name)=='S' ) then
          SINC = ( (SX - S)/DT )
       end if

! Update friendlies
!------------------

       if(FRIENDLY) then
          S = SX
       end if

! Fill exports of U,V and S after diffusion
      if( trim(name) == 'U' ) then
          if(associated(UAFDIFFUSE)) UAFDIFFUSE = SX
       endif
      if( trim(name) == 'V' ) then
          if(associated(VAFDIFFUSE)) VAFDIFFUSE = SX
       endif
       if( trim(name) == 'S' ) then 
          if(associated(SAFDIFFUSE)) SAFDIFFUSE = SX
       endif
       if( trim(name) == 'Q' ) then
          if(associated(QAFDIFFUSE)) QAFDIFFUSE = SX
       endif

! Compute the derivative of the surface flux wrt the surface value
!-----------------------------------------------------------------

       if(associated(SDF)) then
          SDF = CX * (1.0-DX(:,:,LM))
       endif

       if(.not.associated(SRG)) then
          deallocate (SG)
       end if

    enddo ! End loop over all quantities to be diffused
! -----------------------------------------------------

    RETURN_(ESMF_SUCCESS)
  end subroutine DIFFUSE

end subroutine RUN1


!*********************************************************************
!*********************************************************************
!*********************************************************************


!BOP

! !IROUTINE: RUN2 -- The second run stage for the TURBULENCE component

! !INTERFACE:

   subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp)                :: GC     ! Gridded component 
    type(ESMF_State)                   :: IMPORT ! Import state
    type(ESMF_State)                   :: EXPORT ! Export state
    type(ESMF_Clock)                   :: CLOCK  ! The clock
    integer,              intent(out)  :: RC     ! Error code:

! !DESCRIPTION: Second run stage of {\tt GEOS\_TurbulenceGridComp} performs
!    the updates due to changes in surface quantities. Its input are the changes in 
!    surface quantities during the time step. It can also compute the frictional
!    dissipation terms as exports, but these are not added to the temperatures.


!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (ESMF_State       )            :: INTERNAL 
    type (ESMF_Grid        )            :: esmfgrid

! Local variables

    integer                             :: IM, JM, LM
    real                                :: DT

    real, pointer, dimension(:,:)       :: VARFLT
    real, allocatable, dimension(:,:)   :: LATS, LONS

! Begin... 
!---------

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'Run2'

    call MAPL_GridCompTimerStart(GC, "TOTAL")
    call MAPL_GridCompTimerStart(GC, "-RUN2")

! Get parameters from the grid/internal state.
! ---------------------------------------------

    call MAPL_GridCompGet(GC, grid=esmfgrid, num_levels=LM, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridGet(esmfgrid, IM=IM, JM=JM, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridGetCoordinates(esmfgrid, longitudes=LONS, latitudes=LATS, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GridCompGetInternalState(GC, INTERNAL, RC=STATUS)
    VERIFY_(STATUS)

! Get application's timestep from configuration
!----------------------------------------------

    call MAPL_GridCompGetResource(GC, "RUN_DT", DT, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_StateGetPointer(IMPORT,VARFLT,  'VARFLT', RC=STATUS)
    VERIFY_(STATUS)

! Solve the free atmosphere problem
! ---------------------------------

    call MAPL_GridCompTimerStart(GC, "--UPDATE")
      call UPDATE(IM,JM,LM,LATS,RC=STATUS)
      VERIFY_(STATUS)
    call MAPL_GridCompTimerStop(GC, "--UPDATE")

!  All done with RUN
!-------------------

    call MAPL_GridCompTimerStop(GC, "-RUN2")
    call MAPL_GridCompTimerStop(GC, "TOTAL")
    RETURN_(ESMF_SUCCESS)

  contains

!BOP

! !CROUTINE: UPDATE -- Updates diffusive effects for changes at surface.

! !INTERFACE:

    subroutine UPDATE(IM,JM,LM,LATS,RC)

! !ARGUMENTS:

      integer,           intent(IN)       :: IM,JM,LM
      integer, optional, intent(OUT)      :: RC

! !DESCRIPTION: 
!    Some description

!EOP
  
 
      character(len=ESMF_MAXSTR)          :: IAm='Update'
      integer                             :: STATUS

      character(len=ESMF_MAXSTR)          :: TYPE
      character(len=ESMF_MAXSTR)          :: NAME
      type (ESMF_Field)                   :: FIELD
      type (ESMF_FieldBundle)             :: TR
      type (ESMF_FieldBundle)             :: TRI
      real, dimension(:,:,:), pointer     :: PLE
      real, dimension(:,:,:), pointer     :: ZLE
      real, dimension(:,:,:), pointer     :: S, SOI, SINC, INTDIS, TOPDIS
      real, dimension(:,:  ), pointer     :: DSG, SF, SDF, SRFDIS
      real, dimension(:,:  ), pointer     :: HGTLM5, LM50M
      real, dimension(:,:  ), pointer     :: KETRB, KESRF, KETOP, KEINT
      real, dimension(:,:,:), pointer     :: DKS, DKV, DKQ, DKX, EKV, FKV
      real, dimension(:,:,:), pointer     :: DPDTTRB
      real, dimension(:,:,:), pointer     :: QTFLXTRB, SLFLXTRB, WSL, WQT, &
                                             TKH, UFLXTRB, VFLXTRB, QTX, SLX, &
                                             SLFLXMF, QTFLXMF, MFAW, TKEDISS

      integer                             :: KM, K, L, I, J
      logical                             :: FRIENDLY
      logical                             :: WEIGHTED
      real, dimension(IM,JM,LM)           :: DZ, DP, SX
      real, dimension(IM,JM,LM)           :: DF
      real, dimension(IM,JM,LM)           :: QT,SL,U,V,ZLO
      real, dimension(IM,JM,0:LM)         :: ZL0
      real, allocatable                   :: tmp3d(:,:,:)
      integer, allocatable                :: KK(:)
      !  pointers to export of S after update
      real, dimension(:,:,:), pointer     :: SAFUPDATE

! The following variables are for SHVC parameterization

      real, dimension(IM,JM,LM)           :: SOIOFS, XINC
      real,    dimension(IM,JM)           :: z500, z1500, z7000, STDV
      integer, dimension(IM,JM)           :: L300, L500, L1500, L7000, LSURF
      integer, dimension(IM,JM)           :: LTOPS,LBOT,LTOPQ
      logical, dimension(IM,JM)           :: DidSHVC
      real                                :: REDUFAC, SUMSOI
      real                                :: SHVC_CRIT
      real                                :: SHVC_1500, SHVC_ZDEPTH
      real                                :: lat_in_degrees, lat_effect
      real,  dimension(IM,JM)             :: LATS
      real                                :: SHVC_ALPHA, SHVC_EFFECT, SHVC_SCALING 
      logical                             :: DO_SHVC
      logical                             :: ALLOC_TMP
      integer                             :: KS, DO_SHOC
      real                                :: HGT_SURFACE, WGTSUM, CAP_INTDIS, CAP_TOPDIS
      type(ESMF_Info) :: infoh

      ! For idealized SCM surface layer
      integer :: SCM_SL

      character(len=ESMF_MAXSTR) :: GRIDNAME
      character(len=4)           :: imchar
      character(len=2)           :: dateline
      integer                    :: imsize,nn

! Pressure-weighted dissipation heating rates
!--------------------------------------------

      ALLOC_TMP = .FALSE.

      call MAPL_StateGetPointer(INTERNAL, TKH , 'TKH' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, TKEDISS, 'TKEDISS' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_StateGetPointer(EXPORT, QTX      , 'QT'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, SLX      , 'SL'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, QTFLXTRB , 'QTFLXTRB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, SLFLXTRB , 'SLFLXTRB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, UFLXTRB  , 'UFLXTRB'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, VFLXTRB  , 'VFLXTRB'  , RC=STATUS); VERIFY_(STATUS)

      ! MF contribution, used to calculate TRB fluxes above
      call MAPL_StateGetPointer(INTERNAL, SLFLXMF  , 'SLFLXMF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, QTFLXMF  , 'QTFLXMF'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, MFAW     , 'MFAW'     , RC=STATUS); VERIFY_(STATUS)
      
      ! Used in update_moments for ADG PDF (requires all of above)
      call MAPL_StateGetPointer(EXPORT, WSL,     'WSL'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, WQT,     'WQT'   , RC=STATUS); VERIFY_(STATUS)

      call MAPL_StateGetPointer(EXPORT, KETRB ,  'KETRB' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, KESRF ,  'KESRF' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, KETOP ,  'KETOP' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, KEINT ,  'KEINT' , RC=STATUS); VERIFY_(STATUS)

      call MAPL_StateGetPointer(EXPORT, DPDTTRB,  'DPDTTRB', RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_StateGetPointer(EXPORT, SRFDIS,  'SRFDIS', RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, INTDIS,  'INTDIS', RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(EXPORT, TOPDIS,  'TOPDIS', RC=STATUS)
      VERIFY_(STATUS)

                      HGT_SURFACE = 50.0
      if (LM .eq. 72) HGT_SURFACE =  0.0
      call MAPL_GridCompGetResource(GC, 'HGT_SURFACE', HGT_SURFACE, default=HGT_SURFACE, RC=STATUS)
      VERIFY_(STATUS)

                      CAP_INTDIS = 1.0 ! Kelvin [per time step done when applied]                 
      if (LM .eq. 72) CAP_INTDIS = 0.0
      call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_CAP_INTDIS", CAP_INTDIS, default=CAP_INTDIS, RC=STATUS); VERIFY_(STATUS)

                      CAP_TOPDIS = 1.0 ! Kelvin [per time step done when applied]                 
      if (LM .eq. 72) CAP_TOPDIS = 0.0
      call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_CAP_TOPDIS", CAP_TOPDIS, default=CAP_TOPDIS, RC=STATUS); VERIFY_(STATUS)

      call MAPL_GridCompGetResource(GC, trim(COMP_NAME)//"_DO_SHOC", DO_SHOC, default=0, RC=STATUS); VERIFY_(STATUS)


! SHVC Resource parameters. SHVC_EFFECT can be set to zero to turn-off SHVC.
! SHVC_EFFECT = 1. is the tuned value for  2 degree horizontal resolution.
! It should be set to a lower number at higher resolution.

      call MAPL_GridCompGetResource(GC, 'SHVC_EFFECT', SHVC_EFFECT, default=0., RC=STATUS)
      VERIFY_(STATUS)

      DO_SHVC = SHVC_EFFECT > 0.0

      if(DO_SHVC) then
         call MAPL_GridCompGetResource(GC, 'SHVC_CRIT', SHVC_CRIT, default=300., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GridCompGetResource(GC, 'SHVC_ALPHA', SHVC_ALPHA, default=1., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GridCompGetResource(GC, 'SHVC_1500', SHVC_1500, default=2100., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GridCompGetResource(GC, 'SHVC_ZDEPTH', SHVC_ZDEPTH, default=3500., RC=STATUS)
         VERIFY_(STATUS)
         call MAPL_GridCompGetResource(GC, 'SHVC_SCALING', SHVC_SCALING, default=1.0, RC=STATUS)
      end if

! Determine whether running idealized SCM surface layer
!------------------------------------------------------

      call MAPL_GridCompGetResource(GC, 'SCM_SL', SCM_SL, DEFAULT=0)

! Get imports
!------------

      call MAPL_StateGetPointer(IMPORT,    PLE,     'PLE', RC=STATUS); VERIFY_(STATUS)
      call MAPL_StateGetPointer(IMPORT,    ZLE,     'ZLE', RC=STATUS); VERIFY_(STATUS)

! Get the tendecy sensitivities computed in RUN1
!-----------------------------------------------

      call MAPL_StateGetPointer(INTERNAL, DKS,   'DKS',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, DKV,   'DKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, DKQ,   'DKQ',   RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_StateGetPointer(INTERNAL, EKV,   'EKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, FKV,   'FKV',   RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_StateGetPointer(INTERNAL, SINC,  'SINC',    RC=STATUS)
      VERIFY_(STATUS)

! Get the bundles containing the quantities to be diffused, 
!     their tendencies, their surface values, their surface
!     fluxes, and the derivatives of their surface fluxes
!     wrt the surface values. 
!----------------------------------------------------------

      call ESMF_StateGet(IMPORT, 'TR' ,    TR,     RC=STATUS); VERIFY_(STATUS)

      call ESMF_StateGet(EXPORT, 'TRI',    TRI,    RC=STATUS); VERIFY_(STATUS)

! Count them...
!--------------

      call ESMF_FieldBundleGet(TR , FieldCount=KM, RC=STATUS)
      VERIFY_(STATUS)

! and make sure the other bundles are the same.
!----------------------------------------------


! KK gives the order in which quantities will be process.
!--------------------------------------------------------

      allocate(KK(KM), stat=STATUS)
      VERIFY_(STATUS)

      do K = 1,KM
         KK(K) = K
      end do

! Clear the accumulators for the dissipation.
!--------------------------------------------

      if(associated(SRFDIS)) SRFDIS = 0.0
      if(associated(INTDIS)) INTDIS = 0.0
      if(associated(TOPDIS)) TOPDIS = 0.0
      if(associated(KETRB )) KETRB  = 0.0
      if(associated(KESRF )) KESRF  = 0.0
      if(associated(KETOP )) KETOP  = 0.0
      if(associated(KEINT )) KEINT  = 0.0

! Pressure thickness of layers
!-----------------------------

      DP = PLE(:,:,1:LM)-PLE(:,:,0:LM-1)

      do L=0,LM
         ZL0(:,:,L) = ZLE(:,:,L) - ZLE(:,:,LM) ! Edge heights above the surface 
      enddo
      ZLO = 0.5*(ZL0(:,:,1:LM)+ZL0(:,:,0:LM-1)) ! Layer heights above the surface

      DZ = ZLE(:,:,0:LM-1) - ZLE(:,:,1:LM) ! Layer thickness (positive m)

! Diagnostics
      call MAPL_StateGetPointer(EXPORT, HGTLM5 ,  'HGTLM5' , RC=STATUS); VERIFY_(STATUS)
      if(associated(HGTLM5)) then
         HGTLM5 = ZL0(:,:,LM-5)
      end if
      call MAPL_StateGetPointer(EXPORT, LM50M ,  'LM50M' , RC=STATUS); VERIFY_(STATUS)
      if(associated(LM50M)) then
         LM50M = LM
         do L=LM,2,-1
            where (ZL0(:,:,L) <= 50. .and. ZL0(:,:,L-1) > 50.)
               LM50M=L-1
            endwhere
         enddo
      end if

      LSURF=LM
      do L=LM,1,-1
         where (ZL0(:,:,L) <= HGT_SURFACE .and. ZL0(:,:,L-1) > HGT_SURFACE)
            LSURF=L
         endwhere
      enddo

      L300=1.
      do L=LM,1,-1
         where (ZL0(:,:,L) <= 300. .and. ZL0(:,:,L-1) > 300.)
            L300=L
         endwhere
      enddo

      L500=1.
      do L=LM,1,-1
         where (ZL0(:,:,L) <= 500. .and. ZL0(:,:,L-1) > 500.)
            L500=L
         endwhere
      enddo

      if (associated(QTFLXTRB).or.associated(QTX).or.associated(WQT)) then
        QT = 0.0
        ALLOC_TMP = .TRUE.
      end if
      if (associated(SLFLXTRB).or.associated(SLX).or.associated(WSL)) then
        SL = 0.
        ALLOC_TMP = .TRUE.
      end if

      if (associated(UFLXTRB))  U = 0.0
      if (associated(VFLXTRB))  V = 0.0

! Section 1 of 2. SHVC parameterization (W. Chao, J. Atmos. Sci., May 2012, P.1547) 
!  Defining the top and bottom levels of the heat and moisture redistribution layer
!----------------------------------------------------------------------------------

      SHVC_INIT: if(DO_SHVC) then

! Ensure that S is processed first. This only matters for SHVC
!-------------------------------------------------------------

         KS = 0

         do K = 1,KM
            call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
            VERIFY_(STATUS)

            call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
            VERIFY_(STATUS)

            if    (NAME == 'S') then
               KS=KK(1); KK(1)=K; KK(K)=KS
            end if
         end do

         _ASSERT(KS /= 0 ,'needs informative message')

! SHVC super-layers
!------------------

         z500  =  500.
         z1500 = 1500.
         z7000 = 7000.

         STDV = sqrt(varflt*SHVC_SCALING)   ! Scaling VARFLT based on resolution

         where (STDV >=700.)
            z1500 = SHVC_1500                   
         endwhere

         where ( (STDV >300.) .and. (STDV <700.) )
            z1500 = 1500.+ (SHVC_1500-1500.)* (STDV - 300.)/400.
         endwhere  

         z7000 = z1500 + SHVC_ZDEPTH

         L500=1.
         do L=LM,2,-1
            where (ZL0(:,:,L) <= z500 .and. ZL0(:,:,L-1) > z500)     
               L500=L-1    
            endwhere
         enddo

         L1500=1.
         do L=LM,2,-1
            where (ZL0(:,:,L) <= z1500 .and. ZL0(:,:,L-1) > z1500)    
               L1500=L-1
            endwhere
         enddo

         L7000=1.
         do L=LM,2,-1
            where (ZL0(:,:,L) <= z7000 .and. ZL0(:,:,L-1) > z7000)    
               L7000=L-1
            endwhere
         enddo

         LBOT  = L1500-1         
         LTOPS = L7000
         LTOPQ = L1500-(LM-L500)*2

         SOIOFS = 0.0

      end if SHVC_INIT

! Get pointer to export S after update required for SYNCTQ (ALLOC=.TRUE.)
!----------------------------------------------------
    call MAPL_StateGetPointer(EXPORT, SAFUPDATE ,  'SAFUPDATE' , RC=STATUS); VERIFY_(STATUS)

! Loop over all quantities to be diffused.
!-----------------------------------------

       if (associated(INTDIS) .and. DO_SHOC /= 0) then
          INTDIS(:,:,1:LM) = -1.*TKEDISS(:,:,1:LM)*DP(:,:,1:LM)/MAPL_CP
       endif

       TRACERS: do KS=1,KM

         K = KK(KS)

! Get Kth field from bundle
!--------------------------
          
         call ESMF_FieldBundleGet(TR, K, FIELD, RC=STATUS)
         VERIFY_(STATUS)

         call ESMF_FieldGet(FIELD, name=NAME, RC=STATUS)
         VERIFY_(STATUS)

! Get item's diffusion type (U, S or Q; default is Q)
!----------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="DiffuseLike",         &
              VALUE=TYPE, DEFAULTVALUE=dflt_Q,    RC=STATUS)
         VERIFY_(STATUS)

! Get item's friendly status (default is not friendly)
!-----------------------------------------------------

         call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
         call ESMF_InfoGet(infoh,key='FriendlyToTURBULENCE',value=Friendly,default=dflt_false,RC=STATUS)
         VERIFY_(STATUS)

! Get item's weighting (default is unweighted tendencies)
!--------------------------------------------------------

         call ESMF_AttributeGet(FIELD, NAME="WeightedTendency",    &
              VALUE=WEIGHTED, DEFAULTVALUE=dflt_false,    RC=STATUS)
         VERIFY_(STATUS)

! Get pointers to the quantity, its tendency, its surface increment,
!   the preliminary surface flux, and the sensitivity of the surface
!   flux to the surface value.
! ------------------------------------------------------------------

         call TURB_BundleGetPointerToData(TR    ,      NAME,         S  , RC=STATUS)
         VERIFY_(STATUS)
         call TURB_BundleGetPointerToData(TRI   , trim(NAME)//'IT' , SOI, RC=STATUS)
         VERIFY_(STATUS)
         nullify(DSG, SF, SDF)
         select case (trim(NAME))
         case ('S','U','V','Q')
            call MAPL_StateGetPointer(IMPORT, DSG, trim(NAME)//'DEL', RC=STATUS)
            VERIFY_(STATUS)
            call MAPL_StateGetPointer(EXPORT, SF , trim(NAME)//'FLX', RC=STATUS)
            VERIFY_(STATUS)
            call MAPL_StateGetPointer(EXPORT, SDF, trim(NAME)//'DFL', RC=STATUS)
            VERIFY_(STATUS)
         end select

!  Point to the appropriate sensitivity
!--------------------------------------

         if      ( TYPE=='U' ) then
            DKX => DKV
         else if ( TYPE=='Q' ) then
            DKX => DKQ
         else if ( TYPE=='S' ) then
            DKX => DKS
         else
            RETURN_(ESMF_FAILURE)
         end if
         if( trim(NAME)=='QV' ) then
            DKX => DKQ
         end if
         if( trim(NAME)=='S') then
            DKX => DKS
         end if
         if( trim(NAME)=='U' .or. trim(NAME)=='V' ) then
            DKX => DKV
         end if

! Update diffused quantity
!-------------------------

         SX = S

         if( associated(DSG) .and. SCM_SL == 0 ) then
            do L=1,LM 
               SX(:,:,L) = SX(:,:,L) + DKX(:,:,L)*DSG 
            end do
         end if

! Increment the dissipation
!-------------------------- 

         if( TYPE=='U' ) then
            if(associated(INTDIS)) then
               if (DO_SHOC==0) then
                 DF(:,:,1:LM-1) = (0.5/MAPL_CP)*EKV(:,:,1:LM-1)*(SX(:,:,1:LM-1)-SX(:,:,2:LM))**2 ! Shear
                 DF(:,:,  LM  ) =  0.0 ! no shear at the surface, surface friction added later
                 INTDIS(:,:,1:LM-1) = INTDIS(:,:,1:LM-1) + DF(:,:,1:LM-1)
                 INTDIS(:,:,2:LM  ) = INTDIS(:,:,2:LM  ) + DF(:,:,1:LM-1)
               endif
               ! Add surface dissipation to lower levels
               do J=1,JM
                  do I=1,IM
                     WGTSUM = 0.0
                     do L=L300(I,J),LM
                        WGTSUM = WGTSUM + DZ(I,J,L)*(1.0-ZL0(I,J,L)/ZL0(I,J,L300(I,J)))**2
                     end do
                    !  weighted by the layer thickness
                     DF(I,J,LM) = (1.0/MAPL_CP)*EKV(I,J,LM)*SX(I,J,LSURF(I,J))**2 ! Use Surface Winds
                     DF(I,J,LM) = DF(I,J,LM)/WGTSUM
                     do L=L300(I,J),LM
                        INTDIS(I,J,L) = INTDIS(I,J,L) + DF(I,J,LM)*DZ(I,J,L)*(1.0-ZL0(I,J,L)/ZL0(I,J,L300(I,J)))**2
                     end do 
                  end do
               end do
               if (CAP_INTDIS > 0.0) then 
                  ! limit frictional heating from INTDIS by CAP_INTDIS/DT [K/s]
                  do L=1,LM
                     do J=1,JM
                        do I=1,IM
                           INTDIS(I,J,L) = SIGN(MIN(ABS(INTDIS(I,J,L)/DP(I,J,L)),CAP_INTDIS/DT)*DP(I,J,L),INTDIS(I,J,L))
                        end do 
                     end do 
                  end do 
               endif
            endif
            if(associated(TOPDIS)) then
               TOPDIS = TOPDIS + (1.0/MAPL_CP)*FKV*SX**2
               if (CAP_TOPDIS > 0.0) then
                  ! limit frictional heating from TOPDIS by CAP_TOPDIS/DT [K/s]
                  do L=1,LM                                   
                     do J=1,JM
                        do I=1,IM
                           TOPDIS(I,J,L) = SIGN(MIN(ABS(TOPDIS(I,J,L)/DP(I,J,L)),CAP_TOPDIS/DT)*DP(I,J,L),TOPDIS(I,J,L))
                        end do 
                     end do 
                  end do
               endif
             endif
            if(associated(SRFDIS)) then
               SRFDIS = SRFDIS + (1.0/MAPL_CP)*EKV(:,:,LM)*SX(:,:,LM)**2
            endif
         end if

! Update tendencies
! -----------------

         if( associated(SOI) .and. associated(DSG) .and. SCM_SL == 0 ) then
            if( WEIGHTED ) then
               do L=1,LM
                  SOI(:,:,L) = SOI(:,:,L) +  (DKX(:,:,L)*DSG/DT)*DP(:,:,L)
               end do
            else
               do L=1,LM
                  SOI(:,:,L) = SOI(:,:,L) +  (DKX(:,:,L)*DSG/DT)
               end do
            endif
         end if

! Section 2 of 2. SHVC parameterization   (W. Chao,  J. Atmos. Sci., 2012, p1547)   
!  To use SHVC set SHVC_EFFECT in AGCM.rc to > 0.0.
!--------------------------------------------------------------------------------

         RUN_SHVC: if (DO_SHVC) then

            XINC = 0.0

            S_or_Q: if (NAME=='S') then

               if( associated(DSG) .and. SCM_SL == 0 ) then
                  do L=1,LM
                     SINC(:,:,L) = SINC(:,:,L) +  (DKX(:,:,L)*DSG/DT)
                  end do
               end if

               do I=1,IM
                  do J=1,JM
                lat_effect = 1.
                lat_in_degrees= ABS(LATS(I,J)/(3.14159/2.)*90.)
                if (lat_in_degrees >=42.) lat_effect=0.
                if (lat_in_degrees >37. .and. lat_in_degrees < 42.)   &
                lat_effect = 1.0 - (lat_in_degrees-37.)/(42.-37.)
                     if (STDV(I,J) > SHVC_CRIT) then

                        SUMSOI = sum(SINC(I,J,L500(I,J):LM)*DP(I,J,L500(I,J):LM))
                        DidSHVC(I,J) = SUMSOI >= 0.0

                        if (DidSHVC(I,J)) then
                           if (STDV(I,J) >= 800.) then
                              REDUFAC = 1.0
                           elseif (STDV(i,j) >700. .and. STDV(I,J) <800.) then
                              REDUFAC = 0.95 + 0.05*(STDV(I,J)-700.)/100.
                           else
                              REDUFAC = max(min((STDV(I,J)-SHVC_CRIT)/100.,0.95),0.0)
                           end if

                           REDUFAC = REDUFAC * SHVC_EFFECT  *lat_effect       

                           SUMSOI = 0.
                           do L=L500(i,j),LM
                              SUMSOI        = SUMSOI + SINC(I,J,L)*REDUFAC*DP(I,J,L)
                              XINC  (I,J,L) = -SINC(I,J,L) * REDUFAC
                              SOIOFS(I,J,L) =  XINC(I,J,L) / SX(I,J,L)
                           enddo   !do L

                           XINC(I,J,LTOPS(I,J):LBOT(I,J)) = SUMSOI/SUM(DP(I,J,LTOPS(I,J):LBOT(I,J)))
                        endif
                     else
                        DidSHVC(I,J) = .false.
                     endif   ! end of if (STDV>SHVC_CRIT)
                  enddo   !do J
               enddo   !do I

            elseif (NAME == 'Q')  then

! SHVC_ALPHA below is the alpha factor mentioned on page 1552 of Chao (2012, cited above)
!----------------------------------------------------------------------------------------

               do J=1,JM
                  do I=1,IM
                     if (DidSHVC(I,J)) then
                        SUMSOI = 0.
                        do L=L500(I,J),LM
                           XINC(I,J,L) = SHVC_ALPHA*(SOIOFS(I,J,L)*SX(I,J,L))
                           SUMSOI      = SUMSOI +  XINC(I,J,L)*DP(I,J,L)
                        enddo

                        XINC(I,J,LTOPQ(I,J):LBOT(I,J)) = - SUMSOI/SUM(DP(I,J,LTOPQ(I,J):LBOT(I,J)))
                     endif
                  enddo
               enddo

            endif S_or_Q

            if (name == 'S' .or. name == 'Q') then
               SX  = SX + XINC * DT

               if(associated(SOI)) then
                  if(WEIGHTED) then
                     SOI = SOI + XINC*DP
                  else
                     SOI = SOI + XINC
                  end if
               end if
            end if


         end if RUN_SHVC

! Replace friendly
!-----------------

         if(FRIENDLY) then
            S = SX
         end if

! Fill export uf S after update
       if( name=='S' ) then 
          if(associated(SAFUPDATE)) SAFUPDATE = SX
       endif

! Update surface fluxes
! ---------------------

         if( associated(SF) .and. associated(DSG) .and. SCM_SL == 0 ) then
            SF = SF + DSG*SDF
         end if

         if(associated(DPDTTRB)) then
            if( name=='Q' ) then
               DPDTTRB(:,:,1:LM-1) = 0.0
               DPDTTRB(:,:,LM)     = MAPL_GRAV*SF
            end if
         end if

       if( name=='Q' .or. name=='QLLS'  .or. name=='QLCN'  .or. &
                          name=='QILS'  .or. name=='QICN' ) then
          if(associated(QTFLXTRB).or.associated(QTX)) QT = QT + SX
       endif

       if( name=='S' ) then
           if(associated(SLFLXTRB).or.associated(SLX).or.associated(WSL)) SL = SL + SX
       end if

       if( name=='QLLS' .or. name=='QLCN' ) then
          if(associated(SLFLXTRB).or.associated(SLX).or.associated(WSL)) SL = SL - MAPL_ALHL*SX
       endif

       if( name=='QILS' .or. name=='QICN' ) then
          if(associated(SLFLXTRB).or.associated(SLX).or.associated(WSL)) SL = SL - MAPL_ALHS*SX
       endif

       if( name=='U' ) then
           if(associated(UFLXTRB)) U = U + SX
       end if

       if( name=='V' ) then
           if(associated(VFLXTRB)) V = V + SX
       end if

      enddo TRACERS

! End loop over all quantities to be diffused
!--------------------------------------------

      deallocate(KK)

      if (ALLOC_TMP) allocate(tmp3d(IM,JM,0:LM))

      if (associated(QTX)) QTX = QT
      if (associated(SLX)) SLX = SL

! Calculate diagnostic fluxes due to ED and MF (edges)
! and total flux for ADG PDF (centers)
!--------------------------------------------
      if (associated(QTFLXTRB).or.associated(WQT)) then
         tmp3d(:,:,1:LM-1) = (QT(:,:,1:LM-1)-QT(:,:,2:LM))/(ZLO(:,:,1:LM-1)-ZLO(:,:,2:LM))
         tmp3d(:,:,1:LM-1) = -1.*TKH(:,:,1:LM-1)*tmp3d(:,:,1:LM-1)
         tmp3d(:,:,LM) = tmp3d(:,:,LM-1)
         tmp3d(:,:,0) = 0.0
         if (associated(QTFLXMF).and.MFPARAMS%IMPLICIT.eq.1) then
            QTFLXMF(:,:,1:LM-1) = QTFLXMF(:,:,1:LM-1)-MFAW(:,:,1:LM-1)*QT(:,:,1:LM-1)
            QTFLXMF(:,:,LM) = QTFLXMF(:,:,LM-1)
            QTFLXMF(:,:,0) = 0.
         end if
         if (associated(QTFLXTRB)) QTFLXTRB = tmp3d + QTFLXMF
         if (associated(WQT)) WQT = 0.5*( tmp3d(:,:,1:LM)+tmp3d(:,:,0:LM-1) + QTFLXMF(:,:,1:LM)+QTFLXMF(:,:,0:LM-1) ) 
      end if

      if (associated(SLFLXTRB).or.associated(WSL)) then
         tmp3d(:,:,1:LM-1) = (SL(:,:,1:LM-1)-SL(:,:,2:LM))/(ZLO(:,:,1:LM-1)-ZLO(:,:,2:LM))
         tmp3d(:,:,1:LM-1) = -1.*TKH(:,:,1:LM-1)*tmp3d(:,:,1:LM-1)
         tmp3d(:,:,LM) = tmp3d(:,:,LM-1)
         tmp3d(:,:,0) = 0.0
         if (associated(SLFLXMF).and.MFPARAMS%IMPLICIT.eq.1) then
            SLFLXMF(:,:,1:LM-1) = SLFLXMF(:,:,1:LM-1)-MFAW(:,:,1:LM-1)*SL(:,:,1:LM-1)/MAPL_CP
            SLFLXMF(:,:,LM) = 0.
            SLFLXMF(:,:,0) = 0.
         end if
         if (associated(SLFLXTRB)) SLFLXTRB = tmp3d/MAPL_CP + SLFLXMF
         if (associated(WSL)) WSL = 0.5*( (tmp3d(:,:,1:LM)+tmp3d(:,:,0:LM-1))/MAPL_CP + SLFLXMF(:,:,1:LM)+SLFLXMF(:,:,0:LM-1) )         
      end if
      if (ALLOC_TMP) deallocate(tmp3d)
      if (associated(UFLXTRB)) then
         UFLXTRB(:,:,1:LM-1) = (U(:,:,1:LM-1)-U(:,:,2:LM))/(ZLO(:,:,1:LM-1)-ZLO(:,:,2:LM))
         UFLXTRB(:,:,1:LM-1) = -1.*TKH(:,:,1:LM-1)*UFLXTRB(:,:,1:LM-1)
         UFLXTRB(:,:,LM) = UFLXTRB(:,:,LM-1)
         UFLXTRB(:,:,0) = 0.0
      end if
      if (associated(VFLXTRB)) then
         VFLXTRB(:,:,1:LM-1) = (V(:,:,1:LM-1)-V(:,:,2:LM))/(ZLO(:,:,1:LM-1)-ZLO(:,:,2:LM))
         VFLXTRB(:,:,1:LM-1) = -1.*TKH(:,:,1:LM-1)*VFLXTRB(:,:,1:LM-1)
         VFLXTRB(:,:,LM) = VFLXTRB(:,:,LM-1)
         VFLXTRB(:,:,0) = 0.0
      end if

      if(associated(INTDIS)) then
         if(associated(KETRB)) then
            do L=1,LM
               KETRB = KETRB - INTDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
            end do
         end if
         if(associated(KEINT)) then
            do L=1,LM
               KEINT = KEINT - INTDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
            end do
            if(associated(SRFDIS)) KEINT = KEINT + SRFDIS* (MAPL_CP/MAPL_GRAV) ! avoid double-counting SRF in INT
         end if
      end if

      if(associated(TOPDIS)) then
         if(associated(KETRB)) then
            do L=1,LM
               KETRB = KETRB - TOPDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
            end do
         end if
         if(associated(KETOP)) then
            do L=1,LM
               KETOP = KETOP - TOPDIS(:,:,L)* (MAPL_CP/MAPL_GRAV)
            end do
         end if
      endif
      if(associated(SRFDIS)) then
         if(associated(KESRF)) KESRF = KESRF - SRFDIS* (MAPL_CP/MAPL_GRAV)
      endif

      RETURN_(ESMF_SUCCESS)
    end subroutine UPDATE

  end subroutine RUN2


!*********************************************************************
!*********************************************************************
!*********************************************************************

! TURB_MaxMin: DEBUG_TRB-only diagnostic print, replaces the old
! MAPL2-era `MAPL_MaxMin(label, array)` (print a label plus the array's
! max/min to the log). MAPL3's real `MAPL_MaxMin` (mp_utils/
! ArrayReductions.F90) is a completely different API now - a pure
! parallel-reduction FUNCTION returning a [max,min] 2-element array and
! requiring an explicit MPI communicator, with no label/print behavior
! at all. Since every call site here is purely for interactive
! DEBUG_TRB=.true. debugging (not correctness-relevant), a plain local
! (per-PET, non-globally-reduced) min/max print is a reasonable
! like-for-like stand-in rather than wiring an MPI communicator through
! for a debug aid.

   subroutine TURB_MaxMin2(label, array)
     character(len=*),      intent(in) :: label
     real, dimension(:,:),  intent(in) :: array

     print *, trim(label), ' min=', minval(array), ' max=', maxval(array)
   end subroutine TURB_MaxMin2

   subroutine TURB_MaxMin3(label, array)
     character(len=*),        intent(in) :: label
     real, dimension(:,:,:),  intent(in) :: array

     print *, trim(label), ' min=', minval(array), ' max=', maxval(array)
   end subroutine TURB_MaxMin3

!*********************************************************************
!*********************************************************************
!*********************************************************************

! TURB_BundleGetPointerToData: replaces the deleted MAPL2 ESMFL_Mod
! routine ESMFL_BundleGetPointerToData (base/ESMFL_Mod.F90, removed
! upstream in MAPL commit 20019e0e "Delete ESMFL_Mod.F90 and migrate all
! internal callers" - no MAPL3 replacement was ever added since nothing
! inside MAPL itself still called it). Its old body was plain ESMF (no
! MAPL-specific logic): fetch the named field out of the bundle, then
! its data pointer if-and-only-if the field is actually allocated
! (ESMF_FIELDSTATUS_COMPLETE), else leave the pointer null. Reproduced
! here verbatim (see mapl3-porting-notes.md), for the rank-3 TR/TRI
! bundle members used by DIFFUSE/UPDATE. (The surface-only quantities that
! used to live in the TRG/DTG/FSTAR/DFSTAR bundles are now plain
! imports/exports.)

   subroutine TURB_BundleGetPointerToData3(BUNDLE, NAME, PTR, RC)
     type(ESMF_FieldBundle), intent(inout) :: BUNDLE
     character(len=*),       intent(in)    :: NAME
     real, dimension(:,:,:), pointer       :: PTR
     integer, optional,      intent(out)   :: RC

     type(ESMF_Field)            :: FIELD
     type(ESMF_FieldStatus_Flag) :: fieldStatus
     integer                     :: STATUS

     call ESMF_FieldBundleGet(BUNDLE, NAME, FIELD=FIELD, RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_FieldGet(FIELD, status=fieldStatus, RC=STATUS)
     VERIFY_(STATUS)
     if (fieldStatus == ESMF_FIELDSTATUS_COMPLETE) then
        call ESMF_FieldGet(FIELD, 0, PTR, RC=STATUS)
        VERIFY_(STATUS)
     else
        nullify(PTR)
     end if

     RETURN_(ESMF_SUCCESS)
   end subroutine TURB_BundleGetPointerToData3

!*********************************************************************
!*********************************************************************
!*********************************************************************

!*********************************************************************

!*********************************************************************

!BOP

! !IROUTINE:  LOUIS_KS -- Computes atmospheric diffusivities at interior levels

! !INTERFACE:

   subroutine LOUIS_KS( IM,JM,LM,    &
         ZZ,ZE,PV,UU,VV,ZPBL,        &
         KH,KM,RI,LOUISKH,LOUISKM,   & 
         MINSHEAR, MINTHICK,         &
         LAMBDAM, LAMBDAM2,          &
         LAMBDAH, LAMBDAH2,          &
         ALHFAC, ALMFAC,             &
         ZKHMENV, AKHMMAX,           &
         DU_DIAG, ALH_DIAG,KMLS_DIAG,KHLS_DIAG)

! !ARGUMENTS:

      ! Inputs
      integer, intent(IN   ) :: IM,JM,LM
      real,    intent(IN   ) ::   ZZ(IM,JM,  LM) ! Height of layer center above the surface (m).
      real,    intent(IN   ) ::   PV(IM,JM,  LM) ! Virtual potential temperature at layer center (K).
      real,    intent(IN   ) ::   UU(IM,JM,  LM) ! Eastward velocity at layer center (m s-1).
      real,    intent(IN   ) ::   VV(IM,JM,  LM) ! Northward velocity at layer center (m s-1).
      real,    intent(IN   ) ::   ZE(IM,JM,0:LM) ! Height of layer base above the surface (m).
      real,    intent(IN   ) :: ZPBL(IM,JM     ) ! PBL Depth (m)

      ! Outputs
      real,    intent(  OUT) ::   KM(IM,JM,0:LM) ! Momentum diffusivity at base of each layer (m+2 s-1).
      real,    intent(  OUT) ::   KH(IM,JM,0:LM) ! Heat diffusivity at base of each layer  (m+2 s-1).
      real,    intent(  OUT) ::   RI(IM,JM,0:LM) ! Richardson number
   
      ! Diagnostic outputs
      real,    pointer       ::   DU_DIAG(:,:,:) ! Magnitude of wind shear (s-1).
      real,    pointer       ::  ALH_DIAG(:,:,:) ! Blackadar Length Scale diagnostic (m) [Optional] 
      real,    pointer       :: KMLS_DIAG(:,:,:) ! Momentum diffusivity at base of each layer (m+2 s-1).
      real,    pointer       :: KHLS_DIAG(:,:,:) ! Heat diffusivity at base of each layer  (m+2 s-1).

      ! These are constants
      real,    intent(IN   ) :: LOUISKH, LOUISKM ! Louis scheme parameters (usually 5).
      real,    intent(IN   ) :: MINSHEAR     ! Min shear allowed in Ri calculation (s-1).
      real,    intent(IN   ) :: MINTHICK     ! Min layer thickness (m).
      real,    intent(IN   ) :: LAMBDAM      ! Blackadar(1962) length scale parameter for momentum (m).
      real,    intent(IN   ) :: LAMBDAM2     ! Second Blackadar parameter for momentum (m).
      real,    intent(IN   ) :: LAMBDAH      ! Blackadar(1962) length scale parameter for heat (m).
      real,    intent(IN   ) :: LAMBDAH2     ! Second Blackadar parameter for heat (m).
      real,    intent(IN   ) :: ALHFAC
      real,    intent(IN   ) :: ALMFAC
      real,    intent(IN   ) :: ZKHMENV      ! Transition height for reduction of diffusivity in the free atm (m)
      real,    intent(IN   ) :: AKHMMAX      ! Maximum allowe diffusivity (m+2 s-1).

! !DESCRIPTION: Computes Louis et al.(1979) Richardson-number-based diffusivites,
!                as well as an additional ``entrainment'' diffusivity.
!                The Louis diffusivities for momentum, $K_m$, and for heat
!   and moisture, $K_h$, are defined at the interior layer edges. For LM layers,
!   we define diffusivities at the base of the top LM-1 layers. All indexing
!   is from top to bottom of the atmosphere. 
!
!
!  The Richardson number, Ri, is defined at the same edges as the diffusivities. 
!  $$
!  {\rm Ri}_l = \frac{ \frac{g}{\left(\overline{\theta_v}\right)_l}\left(\frac{\delta \theta_v}{\delta z}\right)_l }
!                    { \left(\frac{\delta {\bf |V|}}{\delta z}\right)^2_l             }, \, \,  l=1,LM-1
!  $$
!  where $\theta_v=\theta(1+\epsilon q)$ is the virtual potential temperature,
!  $\epsilon=\frac{M_a}{M_w}-1$, $M_a$ and $M_w$ are the molecular weights of
!  dry air and water, and $q$ is the specific humidity.
!  $\delta \theta_v$ is the difference of $\theta_v$ in the layers above and below the edge 
!  at which Ri$_l$ is defined; $\overline{\theta_v}$ is their average.
!
!  The diffusivities at the layer edges have the form:
!  $$
!  K^m_l = (\ell^2_m)_l \left(\frac{\delta {\bf |V|}}{\delta z}\right)_l f_m({\rm Ri}_l)
!  $$
!  and
!  $$
!  K^h_l = (\ell^2_h)_l \left(\frac{\delta {\bf |V|}}{\delta z}\right)_l f_h({\rm Ri}_l),
!  $$
!  where $k$ is the Von Karman constant, and $\ell$ is the 
!  Blackdar(1962) length scale, also defined at the layer edges.
!
!  Different turbulent length scales can be used for heat and momentum. 
!  in both cases, we use the  traditional formulation:
!  $$
!  (\ell_{(m,h)})_l  = \frac{kz_l}{1 + \frac{kz_l}{\lambda_{(m,h)}}},
!  $$
!  where, near the surface, the scale is proportional to $z_l$, the height above 
!  the surface of edge level $l$, and far from the surface it approaches $\lambda$.
!  The length scale $\lambda$ is usually taken to be a constant (order 150 m), assuming
!  the same scale for the outre boundary layer and the free atmosphere. We make it
!  a function of height, reducing its value in the free atmosphere. The momentum
!  length scale written as:
!  $$
!  \lambda_m = \max(\lambda_1 e^{\left(\frac{z_l}{z_T}\right)^2}, \lambda_2)
!  $$
!  where $\lambda_2 \le \lambda_1$ and $z_T$ is the top of the boundary layer.
!  The length scale for heat and other scalers is taken as: $\lambda_h =  \sqrt\frac{3d}{2} \lambda_m$,
!  following the scheme used at ECMWF.
!
!  The two universal functions of the Richardson number,  $f_m$ and $f_h$,
!  are taken from Louis et al (1982). For unstable conditions (Ri$\le 0$),
!  they are:
!  $$
!  f_m = (1 - 2b \psi)
!  $$
!  and
!  $$
!  f_h = (1 - 3b \psi),
!  $$
!  where
!  $$
!  \psi = \frac{ {\rm Ri} }{ 1 + 3bC(z)\sqrt{-{\rm Ri}} },
!  $$
!  and
!  $$
!  C(z)=
!  $$

!  For stable condition (Ri$\ge 0$), they are
!  $$
!  f_m = \frac{1}{1.0 + \frac{2b{\rm Ri}}{\psi}}
!  $$
!  and
!  $$
!  f_h = \frac{1}{1.0 + 3b{\rm Ri}\psi},
!  $$
!  where
!  $$
!  \psi =  \sqrt{1+d{\rm Ri}}.
!  $$
!  As in Louis et al (1982), the parameters appearing in these are taken  
!  as $b = c = d = 5$. 


!EOP

! Locals

      real, dimension(IM,JM,LM-1) :: ALH, ALM, DV, DZ, DT, TM, LAMBDAM_X, LAMBDAH_X, RLS
      real, dimension(IM,JM     ) :: pbllocal

      integer :: I,J,L
      real :: PS
      real, parameter :: r13 = 1.0/3.0
      real, parameter :: r32 = 3.0/2.0

! Begin...

!===>   Initialize output arrays

      KH = 0.0
      KM = 0.0
      RI = 0.0

!===>   Initialize pbllocal

      pbllocal = ZPBL
      where ( pbllocal .LE. ZZ(:,:,LM) ) pbllocal = ZZ(:,:,LM)

!===>   Quantities needed for Richardson number (all layer edges above the surface layer)

      DZ(:,:,:) = (ZZ(:,:,1:LM-1) - ZZ(:,:,2:LM))
      TM(:,:,:) = (PV(:,:,1:LM-1) + PV(:,:,2:LM))*0.5
      DT(:,:,:) = (PV(:,:,1:LM-1) - PV(:,:,2:LM))
      DV(:,:,:) = (UU(:,:,1:LM-1) - UU(:,:,2:LM))**2 + &
                  (VV(:,:,1:LM-1) - VV(:,:,2:LM))**2

!===>   Limits on distance between layer centers and vertical shear at edges.

      DZ = max(DZ, MINTHICK)
      DV = SQRT(DV)/DZ
      DT = DT/DZ

!===>   Richardson number  ( RI = G*(DTheta_v/DZ) / (Theta_v*|DV/DZ|^2) )

      RI(:,:,1:LM-1) = MAPL_GRAV*DT/(TM*(MAX(DV, MINSHEAR)**2))

!===>   Blackadar(1962) length scale: $1/l = 1/(kz) + 1/\lambda$

     ! 0.1 * local PBL includes diurnal variability on the Blackadar length scale
      pbllocal = 0.1*pbllocal
      do L = 1, LM-1
         LAMBDAM_X(:,:,L) = MAX( pbllocal , LAMBDAM2)
         LAMBDAH_X(:,:,L) = MAX( pbllocal , LAMBDAH2)
      end do

     ! cap the Blackadar length scales
      LAMBDAM_X = MIN(LAMBDAM_X,LAMBDAM)
      LAMBDAH_X = MIN(LAMBDAH_X,LAMBDAH)

      ALM = ( MAPL_KARMAN*ZE(:,:,1:LM-1)/( 1.0 + MAPL_KARMAN*(ZE(:,:,1:LM-1)/LAMBDAM_X) ) )**2
      ALH = ( MAPL_KARMAN*ZE(:,:,1:LM-1)/( 1.0 + MAPL_KARMAN*(ZE(:,:,1:LM-1)/LAMBDAH_X) ) )**2

      if (associated(ALH_DIAG)) then
         ALH_DIAG(:,:,0) = 0.0
         ALH_DIAG(:,:,1:LM-1) = SQRT( ALH )
         ALH_DIAG(:,:,LM) = 0.0
      endif

      do L=1,LM-1
        do J=1,JM
          do I=1,IM
            if ( RI(I,J,L) < 0.0 ) then
             !===> UnStable case PS from eqs 14 and 24 of Louis, 1979
              PS = ( (ZZ(I,J,L)/ZZ(I,J,L+1))**r13 - 1.0 )**3
              PS = SQRT( (PS/(ZE(I,J,L)*(DZ(I,J,L)**3))) * ABS(RI(I,J,L)) )

              KM(I,J,L) = 1.0 - 10.0*RI(I,J,L) / &
                               (1.0 + 15.0*LOUISKM*ALM(I,J,L)*PS)
              KH(I,J,L) = 1.0 - 15.0*RI(I,J,L) / &
                               (1.0 + 15.0*LOUISKH*ALH(I,J,L)*PS)
            else
             !===>   Stable case
              PS = sqrt(1.0 + 5.0*RI(I,J,L))

              KM(I,J,L) = 1.0 / (1.0 + 10.0*RI(I,J,L)/PS)
              KH(I,J,L) = 1.0 / (1.0 + 15.0*RI(I,J,L)*PS)
            end if
          end do
        end do
      end do

!===>   Reduction length in the free atmosphere eq 3.12 (IFS Documentation Cycle CY25r1)

      RLS = (0.2 + (0.8)/(1.0 + (ZE(:,:,1:LM-1)/ZKHMENV)**2))**2

!===>   DIMENSIONALIZE Kz and LIMIT DIFFUSIVITY

      KM(:,:,1:LM-1) = ALMFAC*ALM*KM(:,:,1:LM-1)*DV*RLS
      KH(:,:,1:LM-1) = ALHFAC*ALH*KH(:,:,1:LM-1)*DV*RLS

   if (DEBUG_TRB) then
      call TURB_MaxMin('LOUIS: RI', RI)
      call TURB_MaxMin('LOUIS: KM', KM)
      call TURB_MaxMin('LOUIS: KH', KH)
   endif 

      KM  = min(KM, AKHMMAX)
      KH  = min(KH, AKHMMAX)

      if (associated(  DU_DIAG)) then
         DU_DIAG(:,:,0) = 0.0
         DU_DIAG(:,:,1:LM-1) = DV
         DU_DIAG(:,:,LM) = 0.0
      endif
      if (associated(KMLS_DIAG)) KMLS_DIAG = KM
      if (associated(KHLS_DIAG)) KHLS_DIAG = KH

  end subroutine LOUIS_KS


! !IROUTINE:  LOUIS_KS_OPTIMIZED -- Optimized Computes atmospheric diffusivities at interior levels

subroutine LOUIS_KS_OPTIMIZED( IM,JM,LM,MO_MAX_ITER,DTIME, &
      ZZ,ZE,PV,UU,VV,ZPBL,Z0,Z0H,    &
      KH,KM,RI,                   &
      LOUIS_B_KH,LOUIS_B_KM,      & 
      LOUIS_C_KH,LOUIS_C_KM,      &
      LOUIS_D_KH,LOUIS_D_KM,      &
      MINSHEAR, MINTHICK,         &
      LAMBDAM, LAMBDAM2,          &
      LAMBDAH, LAMBDAH2,          &
      ALHFAC, ALMFAC,             &
      ZKHMENV, AKHMMAX,           &
      DU_DIAG,ALM_DIAG,ALH_DIAG,KMLS_DIAG,KHLS_DIAG)

   implicit none

   ! Physical constants
   real, parameter :: BLACKADAR_SCALE = 0.1     ! PBL local scaling factor
   real, parameter :: R13 = 1.0/3.0             ! 1/3 power
   real, parameter :: MIN_DIFFUSIVITY = 0.01    ! Minimum diffusivity
   real, parameter :: STABILITY_EPS = 1.e-10    ! Small number for stability
   real, parameter :: MAX_PS_DIVISOR = 0.1      ! Maximum PS divisor
   
   ! Arguments (same as original)
   integer, intent(IN) :: IM,JM,LM,MO_MAX_ITER
   real, intent(IN) :: DTIME
   real, intent(IN) :: ZZ(IM,JM,LM), PV(IM,JM,LM), UU(IM,JM,LM), VV(IM,JM,LM)
   real, intent(IN) :: ZE(IM,JM,0:LM), ZPBL(IM,JM), Z0(IM,JM), Z0H(IM,JM)
   real, intent(OUT) :: KM(IM,JM,0:LM), KH(IM,JM,0:LM), RI(IM,JM,0:LM)
   real, intent(IN) :: LOUIS_B_KH,LOUIS_B_KM
   real, intent(IN) :: LOUIS_C_KH,LOUIS_C_KM
   real, intent(IN) :: LOUIS_D_KH,LOUIS_D_KM
   real, intent(IN) :: MINSHEAR, MINTHICK
   real, intent(IN) :: LAMBDAM, LAMBDAM2, LAMBDAH, LAMBDAH2
   real, intent(IN) :: ALHFAC, ALMFAC, ZKHMENV, AKHMMAX
   real, pointer :: DU_DIAG(:,:,:), ALM_DIAG(:,:,:), ALH_DIAG(:,:,:), KMLS_DIAG(:,:,:), KHLS_DIAG(:,:,:)

   ! Local variables - optimized memory layout
   real :: z_local, dz, dz_inv, th_avg, dth_local, shear, shear_sq
   real :: ri_local, ps_local, km_local, kh_local
   real :: lm_local, lh_local, am_local, ah_local, z_exp, zm_factor, zh_factor, cfl_limit
   real :: pbl_local, z_ratio, karman_z
   real, allocatable :: S_M(:,:,:), S_H(:,:,:)

   integer :: i, j, l

   if (associated(ALM_DIAG)) ALM_DIAG = 0.0
   if (associated(ALH_DIAG)) ALH_DIAG = 0.0

   allocate(S_M(IM,JM,LM))
   allocate(S_H(IM,JM,LM))

   if (MO_MAX_ITER > 0) then
       call ComputeMOScalingInterlevel(MO_MAX_ITER, IM, JM, LM, ZE, Z0, Z0H, UU, VV, PV, &
                                       S_M, S_H)
   else
       S_M = 1.0
       S_H = 1.0
   endif

   ! Initialize arrays
   KM = MIN_DIFFUSIVITY
   KH = MIN_DIFFUSIVITY
   RI = 0.0

   ! Optimized main computation loop
   !$OMP PARALLEL DO PRIVATE(i,j,l,dz_inv,dth_local,shear,shear_sq, &
   !$OMP                     z_local,dz,th_avg,ri_local,ps_local, &
   !$OMP                     km_local,kh_local,lm_local,lh_local,am_local,ah_local,z_exp,zm_factor,zh_factor, &
   !$OMP                     cfl_limit,pbl_local,z_ratio, &
   !$OMP                     karman_z) &
   !$OMP             COLLAPSE(2) SCHEDULE(STATIC)
   do l = 1, LM-1
      do j = 1, JM
         ! Vectorizable inner loop with minimal branching
         do i = 1, IM
            
            ! Pre-compute frequently used values

            ! layer thickness
            z_local = ZZ(i,j,l)
            dz = max(ZE(i,j,l) - ZE(i,j,l+1), MINTHICK)
            dz_inv = 1.0 / dz

            ! delta-theta and theta-layer
            dth_local =  PV(i,j,l) - PV(i,j,l+1)
            th_avg    = (PV(i,j,l) + PV(i,j,l+1))*0.5

            ! Wind shear squared
            shear_sq = ((UU(i,j,l) - UU(i,j,l+1)) * dz_inv)**2 + &
                       ((VV(i,j,l) - VV(i,j,l+1)) * dz_inv)**2
            shear = sqrt(shear_sq)
            if (associated(DU_DIAG)) DU_DIAG(i,j,l) = shear
            
            ! Richardson number
            shear_sq = max(shear_sq, MINSHEAR**2) ! Limit SHEAR^2 in RI calculation
            ri_local = (MAPL_GRAV/th_avg) * dth_local * dz_inv / shear_sq
            RI(i,j,l) = ri_local

            ! Cap asymptotic mixing lengths to boundary layer height if provided
            pbl_local = BLACKADAR_SCALE * max(ZPBL(i,j),ZZ(i,j,LM))
            lm_local = min(max(pbl_local,LAMBDAM2),LAMBDAM)
            lh_local = min(max(pbl_local,LAMBDAH2),LAMBDAH)
            ! Blackadar mixing Lengths
            karman_z = MAPL_KARMAN * ZE(i,j,l)
            am_local = ( karman_z / ( 1.0 + karman_z/lm_local ) )**2
            ah_local = ( karman_z / ( 1.0 + karman_z/lh_local ) )**2

            if (associated(ALM_DIAG)) ALM_DIAG(i,j,l) = SQRT(am_local)
            if (associated(ALH_DIAG)) ALH_DIAG(i,j,l) = SQRT(ah_local)

            ! Stability functions - optimized with minimal branching
            if (ri_local < 0.0) then
               ps_local = ( (ZZ(i,j,l)/ZZ(i,j,l+1))**r13 - 1.0 )**3
               ps_local = SQRT( (ps_local/(ZE(i,j,l)*(dz**3))) * abs(ri_local) )
               ! Unstable case
                ! Momentum
               km_local = am_local * ( 1.0 - 2.0*LOUIS_B_KM*ri_local / &
                                      (1.0 + 2.0*LOUIS_B_KM*am_local*LOUIS_C_KM*ps_local) )
                ! Heat
               kh_local = ah_local * ( 1.0 - 2.0*LOUIS_B_KH*ri_local / &
                                      (1.0 + 3.0*LOUIS_B_KH*ah_local*LOUIS_C_KH*ps_local) )
            else
               ! Stable case
                ! Momentum 
               ps_local = max(sqrt(1.0 + LOUIS_D_KM*ri_local), STABILITY_EPS)
               km_local = am_local * ( 1.0 / (1.0 + 2.0*LOUIS_B_KM*ri_local/ps_local) )
                ! Heat
               ps_local = max(sqrt(1.0 + LOUIS_D_KH*ri_local), STABILITY_EPS)
               kh_local = ah_local * ( 1.0 / (1.0 + 3.0*LOUIS_B_KH*ri_local*ps_local) )
            endif

            ! Free atmosphere reduction
            z_exp = exp(-1.0 * ZE(i,j,l)/ZKHMENV)
            zm_factor = S_M(i,j,l) * z_exp
            zh_factor = S_H(i,j,l) * z_exp

            ! Eddy diffusivities
            KM(i,j,l) = ALMFAC * km_local * shear * zm_factor
            KH(i,j,l) = ALHFAC * kh_local * shear * zh_factor

            ! CFL limiting at 1.9xCFL
            if (AKHMMAX < 0.0) then
               cfl_limit = 1.9 * dz*dz / max(DTIME, 1.0e-12)
               KM(i,j,l) = max(MIN_DIFFUSIVITY, min(KM(i,j,l), cfl_limit, ABS(AKHMMAX)))
               KH(i,j,l) = max(MIN_DIFFUSIVITY, min(KH(i,j,l), cfl_limit, ABS(AKHMMAX)))
            else
               KM(i,j,l) = max(MIN_DIFFUSIVITY, min(KM(i,j,l), AKHMMAX))
               KH(i,j,l) = max(MIN_DIFFUSIVITY, min(KH(i,j,l), AKHMMAX))
            endif
         end do
      end do
   end do
   !$OMP END PARALLEL DO

   if (DEBUG_TRB) then
      call TURB_MaxMin('LOUIS: RI', RI)
      call TURB_MaxMin('LOUIS: KM', KM)
      call TURB_MaxMin('LOUIS: KH', KH)
   endif

   if (associated(KMLS_DIAG)) KMLS_DIAG = KM
   if (associated(KHLS_DIAG)) KHLS_DIAG = KH

end subroutine LOUIS_KS_OPTIMIZED

subroutine ComputeMOScalingInterlevel(MO_MAX_ITER, IM, JM, LM, ZE, Z0, Z0H, UU, VV, PV, &
                                      S_M, S_H)
   use MAPL_Constants, only: MAPL_GRAV, MAPL_KARMAN
   implicit none

   integer, intent(IN) :: MO_MAX_ITER, IM, JM, LM
   real, intent(IN) :: ZE(IM,JM,0:LM)          ! interfaces: 0=top ... LM=surface
   real, intent(IN) :: Z0(IM,JM), Z0H(IM,JM)
   real, intent(IN) :: UU(IM,JM,LM), VV(IM,JM,LM), PV(IM,JM,LM)  ! layer centers 1..LM
   real, intent(OUT):: S_M(IM,JM,LM), S_H(IM,JM,LM)              ! scalings per layer center

   ! Parameters / safeties
   real, parameter :: MO_TOL = 1.0e-4
   real, parameter :: ASSUMED_DTH = 0.1
   real, parameter :: EPS = 1.0e-12
   real, parameter :: UMIN = 1.0e-4
   real, parameter :: USTAR_MIN = 1.0e-6
   real, parameter :: ZMIN = 1.0e-6
   real, parameter :: LOG_MIN = 1.0e-6
   real, parameter :: LOG_ARG_FLOOR = 1.01

   integer :: i, j, k, iter
   real :: zc_k, zc_km1, dz_cent         ! midpoint heights and center spacing
   real :: Uref, tref, ustar, Cd, Ch, Lmo
   real :: log_m, log_h, psi_m, psi_h, denom
   real :: dtheta_ref, tmpval

   !$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
   !$OMP   PRIVATE(i,j,k,iter,zc_k,zc_km1,dz_cent,Uref,tref,ustar,Cd,Ch,Lmo,log_m,log_h,psi_m,psi_h,denom,dtheta_ref,tmpval) &
   !$OMP   SHARED(MO_MAX_ITER, IM, JM, LM, ZE, Z0, Z0H, UU, VV, PV, S_M, S_H)
   do k = 1, LM                ! layer centers 1..LM (top -> bottom)
      do j = 1, JM
         do i = 1, IM

            ! default outputs
            S_M(i,j,k) = 1.0
            S_H(i,j,k) = 1.0

            ! compute layer center height zc(k) = 0.5*(ZE(k-1)+ZE(k))
            zc_k = 0.5*( ZE(i,j,k-1) + ZE(i,j,k) )
            if (zc_k < ZMIN) zc_k = ZMIN

            if (k == LM) then
               ! ------------------------------------------------------------
               ! SURFACE-LAYER CASE: apply MO relative to surface roughness
               !  - use center height zc_k as the reference height above ground
               !  - use absolute wind at level LM as reference velocity
               ! ------------------------------------------------------------
               Uref = sqrt( UU(i,j,k)**2 + VV(i,j,k)**2 )
               if (Uref < UMIN) Uref = UMIN
               tref = PV(i,j,k)

               log_m = log( max( zc_k / max(Z0(i,j), EPS), LOG_ARG_FLOOR ) )
               log_h = log( max( zc_k / max(Z0H(i,j), EPS), LOG_ARG_FLOOR ) )

            else
               ! ------------------------------------------------------------
               ! INTER-LEVEL CASE: friction/stability BETWEEN layer centers k-1 and k
               !  - compute midpoint heights zc(k) and zc(k-1)
               !  - use spacing between centers dz_cent = zc(k) - zc(k-1)
               !  - use shear between levels as velocity scale: Uref = |U(k)-U(k-1)|
               !  - use theta difference tref = PV(k) - PV(k-1)
               ! ------------------------------------------------------------
               if (k > 1) then  ! Add this check
                  zc_km1 = 0.5*( ZE(i,j,k-2) + ZE(i,j,k-1) )
               else
                  zc_km1 = ZE(i,j,0)  ! or appropriate boundary value
               endif
               ! defensive floor (zc_km1 is well-defined for k>=2 because ZE index k-2 exists)
               if (zc_km1 < ZMIN) zc_km1 = ZMIN

               dz_cent = zc_k - zc_km1
               if (dz_cent < ZMIN) dz_cent = ZMIN

               ! velocity scale = shear across the interface between centers
               Uref = sqrt( (UU(i,j,k) - UU(i,j,k-1))**2 + (VV(i,j,k) - VV(i,j,k-1))**2 )
               if (Uref < UMIN) Uref = UMIN

               tref = PV(i,j,k) - PV(i,j,k-1)

               ! use dz_cent for log argument (analogue of z/z0 for inter-level)
               log_m = log( max( dz_cent / max(EPS, ZMIN), LOG_ARG_FLOOR ) )
               log_h = log_m
            end if

            ! initial neutral estimates
            Cd = (MAPL_KARMAN / max(log_m, LOG_MIN))**2
            Ch = (MAPL_KARMAN / max(log_h, LOG_MIN))**2
            ustar = sqrt( max(0.0, Cd) ) * Uref
            if (ustar < USTAR_MIN) ustar = USTAR_MIN

            dtheta_ref = ASSUMED_DTH
            Lmo = 1.0e6

            ! iterative MO correction (applies for both surface and inter-level)
            do iter = 1, MO_MAX_ITER

               if (abs(Lmo) < EPS) then
                  psi_m = 0.0
                  psi_h = 0.0
               else if (Lmo > 0.0) then
                  psi_m = -5.0 * zc_k / Lmo
                  psi_h = psi_m
               else
                  tmpval = 1.0 - 16.0 * zc_k / Lmo
                  if (tmpval <= EPS) then
                     psi_m = 0.0
                     psi_h = 0.0
                  else
                     psi_m = 2.0 * log( (1.0 + tmpval**0.25) / 2.0 )
                     psi_h = 2.0 * log( (1.0 + tmpval**0.5) / 2.0 )
                  end if
               end if

               denom = log_m - psi_m
               if (denom <= EPS) denom = EPS
               Cd = (MAPL_KARMAN / denom)**2

               denom = (log_m - psi_m) * (log_h - psi_h)
               if (denom <= EPS) denom = EPS
               Ch = MAPL_KARMAN**2 / denom

               ustar = sqrt( max(0.0, Cd) ) * Uref
               if (ustar < USTAR_MIN) ustar = USTAR_MIN

               Lmo = - (ustar**3 * tref) / ( MAPL_KARMAN * MAPL_GRAV * Ch * dtheta_ref + EPS )

               if (abs(zc_k / Lmo) < MO_TOL) exit
            end do  ! iter

            ! final scalings (safe divide)
            S_M(i,j,k) = Cd / max(Cd, EPS)
            S_H(i,j,k) = Ch / max(Ch, EPS)

         end do  ! i
      end do     ! j
   end do        ! k
   !$OMP END PARALLEL DO

end subroutine ComputeMOScalingInterlevel

!BOP
! !IROUTINE: BELJAARS -- orographic drag following Beljaars (2003)
!
! !INTERFACE:
!
   subroutine BELJAARS(IM, JM, LM, DT, &
                       LAMBDA_B, C_B,  &
                       KPBL,           &
                       U, V, Z, AREA,  &
                       VARFLT, PLE,    &
                       BKV, BKVV, FKV  )

!BOP
!
!   Orographic drag follows  Beljaars (2003):
!   $$
!   \frac{\partial}{\partial z}\frac{\tau}{\rho} = \frac{C_B}{\lambda_B} |U(z)| U(z) 
!          e^{-\tilde{z}^\frac{3}{2}}\tilde{z}^{-1.2},
!   $$
!   where $z$ is the height above the surface in meters, 
!   $\tilde{z}=\frac{z}{\lambda_B}$, $\tau$ is the orographic stress at $z$,
!   $\rho$ is the air density, $U(z)$ is the wind velocity, and $\lambda_B$ is a vertical length scale.
!   Beljaars uses $\lambda_B = 1500$m, for which the non-dimensional parameter $C_B = 2.5101471 \times 10^{-8}$.
!   These are the default values, but both can be modified from the configuration. To avoid underflow.
!   the tendency is set to zero once $\tilde{z}$ exceeds 4 (i.e., 6 km from the surface for default values). 
!
!EOP

      integer, intent(IN   )                    :: IM,JM,LM
      real,    intent(IN   )                    :: DT
      real,    intent(IN   )                    :: LAMBDA_B
      real,    intent(IN   )                    :: C_B

      real,    intent(IN   ), dimension(:,:,: ) :: U
      real,    intent(IN   ), dimension(:,:,: ) :: V
      real,    intent(IN   ), dimension(:,:,: ) :: Z
      real,    intent(IN   ), dimension(:,:   ) :: KPBL, AREA, VARFLT
      real,    intent(IN   ), dimension(:,:,0:) :: PLE

      real,    intent(INOUT), dimension(:,:,: ) :: BKV,BKVV

      real,    intent(  OUT), dimension(:,:,: ) :: FKV

      integer :: I,J,L
      real    :: CBl, wsp, FKV_temp
      real, parameter :: C_TOFD = 9.031E-09 * 12.0

      if (C_B > 0.0) then
      do I = 1, IM
         do J = 1, JM
            CBl = C_B*1.e-7*VARFLT(I,J)
            do L = LM, 1, -1
               FKV(I,J,L) = 0.0
               if (CBl > 0.0 .AND. Z(I,J,L) < 4.0*LAMBDA_B ) then
                  FKV_temp = Z(I,J,L)/LAMBDA_B
                  FKV_temp = exp(-FKV_temp*sqrt(FKV_temp))*(FKV_temp**(-1.2))
                  FKV_temp = CBl*(FKV_temp/LAMBDA_B)*min(5.0,sqrt(U(I,J,L)**2+V(I,J,L)**2))

                  BKV(I,J,L)  = BKV(I,J,L)  + DT*FKV_temp
                  BKVV(I,J,L) = BKVV(I,J,L) + DT*FKV_temp
                  FKV(I,J,L)  = FKV_temp * (PLE(I,J,L)-PLE(I,J,L-1))
               end if
            end do
         end do 
      end do 
      else
    ! C_TOFD is the end product of all coeficients in eq 16 of Beljaars, 2003 (doi: 10.1256/qj.03.73)
    ! C_B is a factor used to amplify the variance of the filtered topography
      CBl = C_TOFD * C_B**2
      do L = LM, 1, -1
        do J = 1, JM
          do I = 1, IM
            FKV(I,J,L) = 0.0
            if (VARFLT(i,j) > 0.0 .AND. Z(I,J,L) < 4.0*LAMBDA_B) then
                wsp = SQRT(U(I,J,L)**2+V(I,J,L)**2)
                FKV_temp = exp(-1*(Z(I,J,L)/LAMBDA_B)**1.5) * Z(I,J,L)**(-1.2)
                FKV_temp = CBl * VARFLT(i,j) * FKV_temp * wsp
                FKV(I,J,L)  = MIN(20.0,FKV_temp * (PLE(I,J,L)-PLE(I,J,L-1))) ! include limit on this forcing for stability
                FKV_temp = FKV(I,J,L)/(PLE(I,J,L)-PLE(I,J,L-1))
                BKV(I,J,L)  = BKV(I,J,L)  + DT*FKV_temp
                BKVV(I,J,L) = BKVV(I,J,L) + DT*FKV_temp
            end if
          end do
        end do
      end do

      if (DEBUG_TRB) call TURB_MaxMin('TOFD: BKV ', BKV)
      if (DEBUG_TRB) call TURB_MaxMin('TOFD: FKV ', FKV*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1)))
      if (DEBUG_TRB) call TURB_MaxMin('TOFD: FKVP', FKV)

      endif

   end subroutine BELJAARS

!*********************************************************************

!BOP

! !IROUTINE:  VTRILU --  Does LU decomposition of tridiagonal matrix.

! !INTERFACE:

   subroutine VTRILU(A,B,C)

! !ARGUMENTS:

      real,               dimension(:,:,:), intent(IN   ) :: C
      real(kind=MAPL_R8), dimension(:,:,:), intent(INOUT) :: A, B

! !DESCRIPTION: {\tt VTRILU} performs an $LU$ decomposition on
! a tridiagonal matrix $M=LU$.
!
! $$
! M = \left( \begin{array}{ccccccc}
!      b_1 & c_1 & & & & & \\
!      a_2 & b_2 & c_2 & & & &  \\
!      &  \cdot& \cdot & \cdot & & &  \\
!      & & \cdot& \cdot & \cdot & &  \\
!      &&  & \cdot& \cdot & \cdot &  \\
!      &&&& a_{K-1} & b_{K-1} & c_{K-1}   \\
!      &&&&& a_{K} & b_{K}
!    \end{array} \right)
! $$
!
!
! $$
! \begin{array}{lr}
! L = \left( \begin{array}{ccccccc}
!      1 &&&&&& \\
!      \hat{a}_2 & 1 & &&&&  \\
!      &  \cdot& \cdot &  & & &  \\
!      & & \cdot& \cdot &  &&  \\
!      &&  & \cdot& \cdot &  &  \\
!      &&&& \hat{a}_{K-1} & 1 &   \\
!      &&&&& \hat{a}_{K} & 1
!    \end{array} \right)
! &
! U = \left( \begin{array}{ccccccc}
!      \hat{b}_1 & c_1 &&&&& \\
!       & \hat{b}_2 & c_2 &&&&  \\
!      &  & \cdot & \cdot & & &  \\
!      & & & \cdot & \cdot &&  \\
!      &&  & & \cdot & \cdot &  \\
!      &&&&  & \hat{b}_{K-1} & c_{K-1}   \\
!      &&&&&  & \hat{b}_{K}
!    \end{array} \right)
! \end{array}
! $$
!
!
! On input, A, B, and C contain, $a_k$, $b_k$, and $c_k$
! the lower, main, and upper diagonals of the matrix, respectively.
! On output, B contains $1/\hat{b}_k$, the inverse of the main diagonal of $U$,
! and A contains $\hat{a}_k$,
! the lower diagonal of $L$. C contains the upper diagonal of the original matrix and of $U$.
!
! The new diagonals $\hat{a}_k$ and $\hat{b}_k$ are:
! $$
! \begin{array}{rcl}
! \hat{b}_1 & = & b_1, \\
! \hat{a}_k & = & \makebox[2 in][l]{$a_k / \hat{b}_{k-1}$,}  k=2, K, \\
! \hat{b}_k & = & \makebox[2 in][l]{$b_k - c_{k-1} \hat{a}_k$,} k=2, K. 
! \end{array}
! $$
!EOP

      integer :: LM, L

      LM = size(C,3)

      B(:,:,1) = 1. / B(:,:,1)

      do L = 2,LM
         A(:,:,L) = A(:,:,L) * B(:,:,L-1)
         B(:,:,L) = 1. / ( B(:,:,L) - C(:,:,L-1) * A(:,:,L) )
      end do

   end subroutine VTRILU

!*********************************************************************

!BOP

! !IROUTINE:  VTRISOLVESURF -- Solves for sensitivity to surface value


! !INTERFACE:

   subroutine VTRISOLVESURF(B,C,Y)

! !ARGUMENTS:

      real,    dimension(:,:,:), intent(IN   ) :: B, C
      real,    dimension(:,:,:), intent(  OUT) :: Y

! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
!   for the special case
!   where the surface Y (YG) is 1 and the rest of the input Ys are 0.
!   Everything else is as in {\tt VTRISOLVE}. This gives the sensitivity of the
!   solution to a unit change in the surface values.

!EOP

      integer :: LM, L

      LM = size(B,3)

      Y(:,:,LM) = -C(:,:,LM) * B(:,:,LM)

      do L = LM-1,1,-1
         Y(:,:,L) = -C(:,:,L) * Y(:,:,L+1) * B(:,:,L)
      end do

   end subroutine VTRISOLVESURF

!BOP

! !IROUTINE:  VTRISOLVE -- Solves for tridiagonal system that has been decomposed by VTRILU


! !INTERFACE:

  subroutine VTRISOLVE ( A,B,C,Y,YG,OPT )

! !ARGUMENTS:

    real,               dimension(:,:,:),  intent(IN   ) ::  A, B, C
    real(kind=MAPL_R8), dimension(:,:,:),  intent(INOUT) ::  Y
    real,               dimension(:,:),    intent(IN)    ::  YG
    logical,                               intent(IN)    :: OPT

! !DESCRIPTION: Solves tridiagonal system that has been LU decomposed
!   $LU x = f$. This is done by first solving $L g = f$ for $g$, and 
!   then solving $U x = g$ for $x$. The solutions are:
! $$
! \begin{array}{rcl}
! g_1 & = & f_1, \\
! g_k & = & \makebox[2 in][l]{$f_k - g_{k-1} \hat{a}_{k}$,}  k=2, K, \\
! \end{array}
! $$
! and  
! $$
! \begin{array}{rcl}
! x_K & = & g_K /\hat{b}_K, \\
! x_k & = & \makebox[2 in][l]{($g_k - c_k g_{k+1}) / \hat{b}_{k}$,}  k=K-1, 1 \\
! \end{array}
! $$
!  
!  On input A contains the $\hat{a}_k$, the lower diagonal of $L$,
!   B contains the $1/\hat{b}_k$, inverse of the  main diagonal of $U$,
!   C contains the $c_k$, the upper diagonal of $U$. The forcing, $f_k$ is
!   
!   It returns the
!   solution in the r.h.s input vector, Y. A has the multiplier from the
!   decomposition, B the 
!   matrix (U), and C the upper diagonal of the original matrix and of U.
!   YG is the LM+1 (Ground) value of Y.

!EOP

    integer :: LM, L

    LM = size(Y,3)

! Sweep down, modifying rhs with multiplier A

    do L = 2,LM
       Y(:,:,L) = Y(:,:,L) - Y(:,:,L-1) * A(:,:,L)
    enddo

! Sweep up, solving for updated value. Note B has the inverse of the main diagonal

    if(size(YG)>0) then
       Y(:,:,LM)   = (Y(:,:,LM) - C(:,:,LM) * YG        )*B(:,:,LM)
    else if (OPT) then
       Y(:,:,LM)   =  Y(:,:,LM)*B(:,:,LM-1)/(B(:,:,LM-1) - A(:,:,LM)*(1.0+C(:,:,LM-1)*B(:,:,LM-1) ))
    else
       Y(:,:,LM)   =  Y(:,:,LM)*B(:,:,LM)/( 1.0+C(:,:,LM)*B(:,:,LM) )     ! Alternate formulation
    endif

    do L = LM-1,1,-1
       Y(:,:,L) = (Y(:,:,L ) - C(:,:,L ) * Y(:,:,L+1))*B(:,:,L )
    enddo

    return
  end subroutine VTRISOLVE


end module GEOS_TurbulenceGridCompMod


