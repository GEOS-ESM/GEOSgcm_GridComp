!  $Id$

#include "MAPL_Generic.h"

module GEOS_GwdGridCompMod

   !BOP
   !MODULE: GEOS_Gwd -- A Module to compute the forcing due to parameterized gravity wave drag
   !DESCRIPTION:
   !
   ! {\tt GWD} is a light-weight gridded component to compute the forcing
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
   ! (Sassi et al. 2003; http://acd.ucar.edu/models/WACCM). Both parameteriz-
   ! ations are based on Lindzen's [1981] formulation. The interested reader
   ! is referred to those publications for details of the mathematical
   ! derivations.

   !USES:

   use esmf
   use MAPL, only: MAPL_Verify, MAPL_Assert, MAPL_Return
   use MAPL, only: MAPL_get_current_thread, MAPL_get_num_threads
   use MAPL, only: MAPL_find_bounds, MAPL_Interval
   use MAPL, only: MAPL_AM_I_ROOT, MAPL_ArrayGather
   use MAPL_Constants, only: MAPL_RADIUS, MAPL_RGAS, MAPL_GRAV, MAPL_VIREPS, MAPL_PI, MAPL_P00, MAPL_CP
   use MAPL, only: MAPL_GridGet, MAPL_GridGetCoordinates, mapl_GridGetGlobalCellCountPerDim
   use MAPL, only: MAPL_GridCompSetEntryPoint
   use MAPL, only: MAPL_GridCompGet, MAPL_GridCompGetResource
   use MAPL, only: MAPL_GridCompGetInternalState
   use MAPL, only: MAPL_GridCompAddSpec, MAPL_GridCompTimerStart, MAPL_GridCompTimerStop
   use MAPL, only: MAPL_StateGetPointer, MAPL_ClockGet
   use MAPL, only: MAPL_RESTART_SKIP
   use MAPL, only: MAPL_VERTICAL_STAGGER_NONE, MAPL_VERTICAL_STAGGER_CENTER, MAPL_VERTICAL_STAGGER_EDGE
   use MAPL, only: MAPL_UngriddedDims, MAPL_UngriddedDim

   use gw_rdg, only : gw_rdg_init
   use gw_oro, only : gw_oro_init
   use gw_convect, only : gw_beres_init, BeresSourceDesc
   use gw_common, only: GWBand, gw_common_init, gw_newtonian_set
   use gw_drag_ncar, only: gw_intr_ncar
   use gw_drag, only: gw_intr

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !EOP
   ! config params
   type :: ThreadWorkspace
      type(GWBand) :: beres_band
      type(BeresSourceDesc) :: beres_dc_desc
      type(GWBand) :: oro_band
      type(GWBand) :: rdg_band
   end type ThreadWorkspace

   logical :: use_threads
   integer :: num_threads

   type :: GEOS_GwdGridComp
      real :: GEOS_BGSTRESS
      real :: GEOS_EFFGWBKG
      real :: GEOS_EFFGWORO
      integer :: GEOS_PGWV
      real :: NCAR_EFFGWBKG
      real :: NCAR_EFFGWORO
      integer :: NCAR_NRDG
      real :: Z1
      real :: TAU1
      real :: H0
      real :: HH
      real, allocatable :: alpha(:)
      type(ThreadWorkspace), allocatable :: workspaces(:)
   end type GEOS_GwdGridComp

   character(*), parameter :: PRIVATE_STATE = "GWD_PRIVATE_STATE"

   logical :: DEBUG_TQ_ERRORS
   logical :: DEBUG_GWD

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !INTERFACE:
   subroutine SetServices(gc, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, optional :: rc ! return code

      !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
      !                the Initialize and Finalize services, as well as allocating
      !   our instance of a generic state and putting it in the
      !   gridded component (GC). Here we only need to set the run method and
      !   add the state variable specifications (also generic) to our instance
      !   of the generic state. This is the way our true state variables get into
      !   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

      type(GEOS_GwdGridComp), pointer :: self
      type(MAPL_UngriddedDim) :: ungrd_16
      integer :: status

      _SET_NAMED_PRIVATE_STATE(gc, GEOS_GwdGridComp, PRIVATE_STATE)
      _GET_NAMED_PRIVATE_STATE(gc, GEOS_GwdGridComp, PRIVATE_STATE, self)

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="run", _RC)

      num_threads = 1
      if (use_threads) then
         num_threads = MAPL_get_num_threads()
      end if
      allocate(self%workspaces(0:num_threads - 1), _STAT)

      ! We need to get NCAR_NRDG because this is used in the auto-generated
      ! code GWD_Internal___.h via ACG
      call MAPL_GridCompGetResource(gc, "NCAR_NRDG", self%NCAR_NRDG, default=0, _RC)

      ungrd_16 = MAPL_UngriddedDim(16, NAME="sixteen", units="1")
#include "GWD_Import___.h"
#include "GWD_Export___.h"
#include "GWD_Internal___.h"

      _RETURN(_SUCCESS)
   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code

      !DESCRIPTION: The Initialize method of the GWD Physics Gridded Component first
      !   calls the Initialize method of the children.  Then, if using the NCAR GWD
      !   scheme, calls the initialization routines.
      !EOP

      type(ESMF_Grid) :: grid
      integer :: im, jm, jm_thread, LM
      real, allocatable :: lats(:, :), lons(:, :)

      character(len=:), allocatable :: gridname
      character(len=4) :: imchar
      character(len=2) :: dateline
      integer :: imsize, nn
      real :: sigma, stretch_factor

      real, pointer, dimension(:) :: PREF

      ! NCAR GWD variables
      character(len=:), allocatable :: BERES_FILE_NAME
      character(len=ESMF_MAXSTR) :: ERRstring
      real :: NCAR_QBO_HDEPTH_SCALING
      real :: NCAR_TAU_TOP_ZERO
      real :: NCAR_PRNDL
      real :: NCAR_ORO_GW_DC, NCAR_BKG_GW_DC
      real :: NCAR_ORO_FCRIT2, NCAR_BKG_FCRIT2
      real :: NCAR_ORO_WAVELENGTH, NCAR_BKG_WAVELENGTH
      real :: NCAR_ORO_SOUTH_FAC
      real :: NCAR_ORO_TNDMAX, NCAR_BKG_TNDMAX
      real :: NCAR_HR_CF ! Grid cell convective conversion factor
      real :: NCAR_TR_EFF ! Convective region efficiency factor
      real :: NCAR_ET_EFF ! Frontal region efficiency factor
      real :: NCAR_ET_TAUBGND ! Extratropical background frontal forcing
      real :: NCAR_EFFGWBKG
      real :: NCAR_DC_BERES_SRC_LEVEL
      integer :: NCAR_ORO_PGWV, NCAR_BKG_PGWV
      integer :: GEOS_PGWV
      logical :: JASON_BKG, JASON_ORO
      logical :: NCAR_ET_USE_DQCDT
      logical :: NCAR_DC_BERES

      type(GEOS_GwdGridComp), pointer :: self
      integer :: thread
      type(MAPL_Interval), allocatable :: bounds(:)
      integer :: status

      _GET_NAMED_PRIVATE_STATE(gc, GEOS_GwdGridComp, PRIVATE_STATE, self)

      ! Grid info
      call MAPL_GridCompGet(gc, grid=grid, num_levels=LM, _RC)
      call MAPL_GridGet(grid, im=im, jm=jm, _RC)
      call MAPL_GridGetCoordinates(grid, longitudes=lons, latitudes=lats, _RC)

      ! Get grid name to determine IMSIZE
      call MAPL_GridCompGetResource(gc, 'AGCM.GRIDNAME', gridname, _RC)
      gridname = AdjustL(gridname)
      nn = len_trim(gridname)
      dateline = gridname(nn - 1:nn)
      imchar = gridname(3:index(gridname, 'x') - 1)
      read(imchar, *) imsize
      if (dateline == 'CF') imsize = imsize * 4
      call MAPL_GridCompGetResource(gc, 'AGCM.STRETCH_FACTOR', stretch_factor, default=1.0, _RC)
      imsize = imsize * CEILING(stretch_factor)
      sigma = 1.0 - 0.9839 * exp(-0.09835 * 4.e7 * 0.9 / imsize / 1000.) ! Based on Arakawa 2011 sigma used in GF2020

      ! Background Gravity wave drag
      call MAPL_GridCompGetResource(gc, 'JASON_BKG', JASON_BKG, default=(LM == 72), _RC)
      if (JASON_BKG) then
         GEOS_PGWV = 4
         call MAPL_GridCompGetResource(gc, "GEOS_PGWV", self%GEOS_PGWV, default=GEOS_PGWV, _RC)
         call MAPL_GridCompGetResource(gc, "GEOS_BGSTRESS", self%GEOS_BGSTRESS, default=0.900, _RC)
         call MAPL_GridCompGetResource(gc, "GEOS_EFFGWBKG", self%GEOS_EFFGWBKG, default=0.125, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_EFFGWBKG", self%NCAR_EFFGWBKG, default=0.000, _RC)
         !! call MAPL_GridCompGetResource(gc, "RAYLEIGH_TAU1", self%TAU1,          default=172800., _RC)
         call MAPL_GridCompGetResource(gc, "RAYLEIGH_TAU1", self%TAU1, default=0.000, _RC)
      else
         GEOS_PGWV = NINT(32 * LM / 181.0)
         call MAPL_GridCompGetResource(gc, "GEOS_PGWV", self%GEOS_PGWV, default=GEOS_PGWV, _RC)
         call MAPL_GridCompGetResource(gc, "GEOS_BGSTRESS", self%GEOS_BGSTRESS, default=0.000, _RC)
         call MAPL_GridCompGetResource(gc, "GEOS_EFFGWBKG", self%GEOS_EFFGWBKG, default=0.000, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_EFFGWBKG", self%NCAR_EFFGWBKG, default=0.375, _RC)
         call MAPL_GridCompGetResource(gc, "RAYLEIGH_TAU1", self%TAU1, default=0.00, _RC)
      end if

      ! Orographic Gravity wave drag
      call MAPL_GridCompGetResource(gc, 'JASON_ORO', JASON_ORO, default=(LM == 72), _RC)
      if (JASON_ORO) then
         call MAPL_GridCompGetResource(gc, "GEOS_EFFGWORO", self%GEOS_EFFGWORO, default=0.250, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_EFFGWORO", self%NCAR_EFFGWORO, default=0.000, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_NRDG", self%NCAR_NRDG, default=0, _RC)
      else
         call MAPL_GridCompGetResource(gc, "GEOS_EFFGWORO", self%GEOS_EFFGWORO, default=0.000, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_NRDG", self%NCAR_NRDG, default=0, _RC)
         if (self%NCAR_NRDG == 16) then
            call MAPL_GridCompGetResource(gc, "NCAR_EFFGWORO", self%NCAR_EFFGWORO, default=1.000, _RC)
         else
            call MAPL_GridCompGetResource(gc, "NCAR_EFFGWORO", self%NCAR_EFFGWORO, default=0.250, _RC)
         end if
      end if

      ! Rayleigh friction
      if (self%TAU1 > 0.0) then
         call MAPL_GridCompGetResource(gc, "RAYLEIGH_Z1", self%Z1, default=75000., _RC)
         call MAPL_GridCompGetResource(gc, "RAYLEIGH_H0", self%H0, default=7000., _RC)
         call MAPL_GridCompGetResource(gc, "RAYLEIGH_HH", self%HH, default=7500., _RC)
      end if

      ! NCAR GWD settings
      call MAPL_GridCompGetResource(gc, "NCAR_TAU_TOP_ZERO", NCAR_TAU_TOP_ZERO, default=50.0, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_PRNDL", NCAR_PRNDL, default=0.50, _RC)
      NCAR_QBO_HDEPTH_SCALING = 1.0 - 0.75 * sigma
      call MAPL_GridCompGetResource( &
           gc, &
           "NCAR_QBO_HDEPTH_SCALING", &
           NCAR_QBO_HDEPTH_SCALING, &
           default=NCAR_QBO_HDEPTH_SCALING, &
           _RC)
      NCAR_HR_CF = CEILING(20.0 * sigma)
      call MAPL_GridCompGetResource(gc, "NCAR_HR_CF", NCAR_HR_CF, default=NCAR_HR_CF, _RC)

      call gw_common_init(NCAR_TAU_TOP_ZERO, 1, &
           MAPL_GRAV, &
           MAPL_RGAS, &
           MAPL_CP, &
           NCAR_PRNDL, NCAR_QBO_HDEPTH_SCALING, NCAR_HR_CF, ERRstring)

      ! Beres Scheme File
      call MAPL_GridCompGetResource(gc, &
           "BERES_FILE_NAME", BERES_FILE_NAME, &
           default='ExtData/g5gcm/gwd/newmfspectra40_dc25.nc', _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_BKG_PGWV", NCAR_BKG_PGWV, default=32, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_BKG_GW_DC", NCAR_BKG_GW_DC, default=2.5, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_BKG_FCRIT2", NCAR_BKG_FCRIT2, default=1.0, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_BKG_WAVELENGTH", NCAR_BKG_WAVELENGTH, default=1.e5, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_TR_EFF", NCAR_TR_EFF, default=1.0, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_ET_EFF", NCAR_ET_EFF, default=1.0, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_ET_TAUBGND", NCAR_ET_TAUBGND, default=6.4, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_ET_USE_DQCDT", NCAR_ET_USE_DQCDT, default=.false., _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_BKG_TNDMAX", NCAR_BKG_TNDMAX, default=250.0, _RC)
      NCAR_BKG_TNDMAX = NCAR_BKG_TNDMAX / 86400.0
      ! Beres DeepCu
      call MAPL_GridCompGetResource(gc, "NCAR_DC_BERES_SRC_LEVEL", NCAR_DC_BERES_SRC_LEVEL, default=70000.0, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_DC_BERES", NCAR_DC_BERES, default=.true., _RC)
      if (use_threads) then
         bounds = MAPL_find_bounds(jm, num_threads)
         do thread = 0, num_threads - 1
            jm_thread = bounds(thread + 1)%max - bounds(thread + 1)%min + 1
            call gw_beres_init(BERES_FILE_NAME, &
                 self%workspaces(thread)%beres_band, &
                 self%workspaces(thread)%beres_dc_desc, &
                 NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2, &
                 NCAR_BKG_WAVELENGTH, NCAR_DC_BERES_SRC_LEVEL, &
                 1000.0, .true., NCAR_TR_EFF, NCAR_ET_EFF, NCAR_ET_TAUBGND, NCAR_ET_USE_DQCDT, &
                 NCAR_BKG_TNDMAX, NCAR_DC_BERES, &
                 im * jm_thread, lats(:, bounds(thread + 1)%min:bounds(thread + 1)%max))
         end do
      else
         call gw_beres_init(BERES_FILE_NAME, &
              self%workspaces(0)%beres_band, &
              self%workspaces(0)%beres_dc_desc, &
              NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2, &
              NCAR_BKG_WAVELENGTH, NCAR_DC_BERES_SRC_LEVEL, &
              1000.0, .true., NCAR_TR_EFF, NCAR_ET_EFF, NCAR_ET_TAUBGND, NCAR_ET_USE_DQCDT, &
              NCAR_BKG_TNDMAX, NCAR_DC_BERES, &
              im * jm, lats)
      end if

      ! Orographic Scheme
      call MAPL_GridCompGetResource(gc, "NCAR_ORO_PGWV", NCAR_ORO_PGWV, default=0, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_ORO_GW_DC", NCAR_ORO_GW_DC, default=2.5, _RC)
      call MAPL_GridCompGetResource(gc, "NCAR_ORO_WAVELENGTH", NCAR_ORO_WAVELENGTH, default=1.e5, _RC)

      if (self%NCAR_NRDG > 0) then
         call MAPL_GridCompGetResource(gc, "NCAR_ORO_FCRIT2", NCAR_ORO_FCRIT2, default=1.0, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_ORO_TNDMAX", NCAR_ORO_TNDMAX, default=400.0, _RC)
         NCAR_ORO_TNDMAX = NCAR_ORO_TNDMAX / 86400.0
         ! Ridge Scheme
         do thread = 0, num_threads - 1
            call gw_rdg_init(self%workspaces(thread)%rdg_band, NCAR_ORO_GW_DC, NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, &
                 NCAR_ORO_TNDMAX, NCAR_ORO_PGWV)
         end do
      else
         ! Old Scheme
         call MAPL_GridCompGetResource(gc, "NCAR_ORO_FCRIT2", NCAR_ORO_FCRIT2, default=0.5, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_ORO_SOUTH_FAC", NCAR_ORO_SOUTH_FAC, default=1.0, _RC)
         call MAPL_GridCompGetResource(gc, "NCAR_ORO_TNDMAX", NCAR_ORO_TNDMAX, default=400.0, _RC)
         NCAR_ORO_TNDMAX = NCAR_ORO_TNDMAX / 86400.0
         do thread = 0, num_threads - 1
            call gw_oro_init(self%workspaces(thread)%oro_band, NCAR_ORO_GW_DC, &
                 NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_PGWV, &
                 NCAR_ORO_SOUTH_FAC, NCAR_ORO_TNDMAX)
         end do
      end if

      call MAPL_GridCompGetResource(gc, "DEBUG_GWD", DEBUG_GWD, default=.false., _RC)
      call MAPL_GridCompGetResource(gc, "DEBUG_TQ_ERRORS", DEBUG_TQ_ERRORS, default=.false., _RC)

      allocate(self%alpha(LM + 1), _STAT)
      call MAPL_StateGetPointer(import, PREF, 'PREF', _RC)
      call gw_newtonian_set(LM, PREF, self%alpha)

      _RETURN(_SUCCESS)
   end subroutine Initialize

   !BOP
   !IROUTINE: RUN -- Run method for the GWD component
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
      !                the Initialize and Finalize services, as well as allocating
      !EOP

      type(ESMF_Alarm) :: alarm
      type(ESMF_Grid) :: grid
      integer :: i, j, L, im, jm, LM
      real :: tcrib
      real, allocatable, dimension(:, :) :: lats, lons
      ! Rayleigh friction parameters
      real :: H0, HH, Z1, TAU1
      type(GEOS_GwdGridComp), pointer :: self
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread, status

      _GET_NAMED_PRIVATE_STATE(gc, GEOS_GwdGridComp, PRIVATE_STATE, self)

      H0 = self%H0
      HH = self%HH
      Z1 = self%Z1
      TAU1 = self%TAU1

      ! Local aliases to the state, grid, and configuration
      ! Grid info
      call MAPL_GridCompGet(gc, grid=grid, num_levels=LM, _RC)
      call MAPL_GridGetCoordinates(grid, longitudes=lons, latitudes=lats, _RC)
      call MAPL_GridGet(grid, im=im, jm=jm, _RC)

      ! If its time, recalculate the GWD tendency
      ! if ( ESMF_AlarmIsRinging( ALARM ) ) then
      ! call ESMF_AlarmRingerOff(ALARM, _RC)
      call MAPL_GridCompTimerStart(gc, "gwd_driver", _RC)
      call Gwd_Driver(_RC)
      call MAPL_GridCompTimerStop(gc, "gwd_driver", _RC)
      ! endif

      _RETURN(_SUCCESS)

   contains

      subroutine Gwd_Driver(rc)
         integer, optional, intent(out) :: rc

         type(ESMF_TimeInterval) :: TINT
         integer :: status

#include "GWD_DeclarePointer___.h"
         real, pointer, dimension(:, :, :) :: PTR3D
         real, pointer, dimension(:, :) :: PTR2D
         real, dimension(im, jm, LM) :: TMP3D
         real, dimension(im, jm, LM) :: ZM, PMID, PDEL, RPDEL, PMLN
         real, dimension(im, jm) :: a2, Hefold
         real, dimension(im, jm, LM) :: DUDT_ORG, DVDT_ORG, DTDT_ORG
         real, dimension(im, jm, LM) :: DUDT_GWD, DVDT_GWD, DTDT_GWD
         real, dimension(im, jm, LM) :: DUDT_RAH, DVDT_RAH, DTDT_RAH
         real, dimension(im, jm, LM) :: DUDT_TOT, DVDT_TOT, DTDT_TOT
         real, dimension(im, jm, LM + 1) :: PILN, ZI
         real, dimension(LM) :: ZREF, KRAY
         real, dimension(im, jm) :: GBXAR_TMP
         real, dimension(im, jm) :: TAUXO_TMP, TAUYO_TMP
         real, dimension(im, jm) :: TAUXB_TMP, TAUYB_TMP
         real, dimension(im, jm, LM + 1) :: TAUXO_3D, TAUYO_3D, FEO_3D, FEPO_3D
         real, dimension(im, jm, LM + 1) :: TAUXB_3D, TAUYB_3D, FEB_3D, FEPB_3D
         real, dimension(im, jm, LM) :: DUBKGSRC, DVBKGSRC, DTBKGSRC
         real, dimension(im, jm) :: KEGWD_X, KEORO_X, KERAY_X, KEBKG_X, KERES_X
         real, dimension(im, jm) :: PEGWD_X, PEORO_X, PERAY_X, PEBKG_X, BKGERR_X
         real, dimension(im, jm, LM) :: DUDT_GWD_GEOS, DVDT_GWD_GEOS, DTDT_GWD_GEOS
         real, dimension(im, jm, LM) :: DUDT_ORG_GEOS, DVDT_ORG_GEOS, DTDT_ORG_GEOS
         real, dimension(im, jm) :: TAUXB_TMP_GEOS, TAUYB_TMP_GEOS
         real, dimension(im, jm) :: TAUXO_TMP_GEOS, TAUYO_TMP_GEOS
         real, dimension(im, jm, LM) :: DUDT_GWD_NCAR, DVDT_GWD_NCAR, DTDT_GWD_NCAR
         real, dimension(im, jm, LM) :: DUDT_ORG_NCAR, DVDT_ORG_NCAR, DTDT_ORG_NCAR
         real, dimension(im, jm) :: TAUXB_TMP_NCAR, TAUYB_TMP_NCAR
         real, dimension(im, jm) :: TAUXO_TMP_NCAR, TAUYO_TMP_NCAR
         real, allocatable, target, dimension(:, :, :) :: scratch_ridge

         integer :: j, K, L, nrdg, ikpbl
         real :: DT ! time interval in sec
         real :: a1, wsp, var_temp
         integer :: i, irun
         type(ESMF_State) :: internal

         call MAPL_ClockGet(clock, dt=DT, _RC)

         ! Pointers to import, export and internal variables
         call MAPL_GridCompGetInternalState(gc, internal)
#include "GWD_GetPointer___.h"

         call PREGEO(im * jm, LM, PLE, lats, PMID, PDEL, RPDEL, PILN, PMLN)

         ! Compute ZM
         call GEOPOTENTIAL(im * jm, LM, PILN, PMLN, PLE, PMID, PDEL, RPDEL, T, Q, ZI, ZM)

         ! Do gravity wave drag calculations on a list of soundings
         if (self%NCAR_NRDG /= 0.0) then
            GBXAR_TMP = GBXAR * (MAPL_RADIUS / 1000.)**2 ! transform to km^2
            where (ANGLL < -180)
               ANGLL = 0.0
            end where
            do nrdg = 1, self%NCAR_NRDG
               KWVRDG(:, :, nrdg) = 0.001 / (HWDTH(:, :, nrdg) + 0.001)
               EFFRDG(:, :, nrdg) = self%NCAR_EFFGWORO * (HWDTH(:, :, nrdg) * CLNGT(:, :, nrdg)) / GBXAR_TMP
            end do
         else
            allocate(scratch_ridge(im, jm, 16))
            scratch_ridge = 0.0
            MXDIS => scratch_ridge
            HWDTH => scratch_ridge
            CLNGT => scratch_ridge
            ANGLL => scratch_ridge
            ANIXY => scratch_ridge
            KWVRDG => scratch_ridge
            EFFRDG => scratch_ridge
            GBXAR_TMP = 0.0
         end if

         ! Use new NCAR code convective+oro (excludes extratropical bkg sources)
         DUDT_GWD_NCAR = 0.0
         DVDT_GWD_NCAR = 0.0
         DTDT_GWD_NCAR = 0.0
         TAUXB_TMP_NCAR = 0.0
         TAUYB_TMP_NCAR = 0.0
         DUDT_ORG_NCAR = 0.0
         DVDT_ORG_NCAR = 0.0
         DTDT_ORG_NCAR = 0.0
         TAUXO_TMP_NCAR = 0.0
         TAUYO_TMP_NCAR = 0.0
         call MAPL_GridCompTimerStart(gc, "gw_intr_ncar", _RC)
         if ((self%NCAR_EFFGWORO /= 0.0) .or. (self%NCAR_EFFGWBKG /= 0.0)) then
            do L = 1, LM
               TMP3D(:, :, L) = (1.0 - CNV_FRC) * (DQLDT(:, :, L) + DQIDT(:, :, L))
            end do
            if (associated(DQCDT_LS)) DQCDT_LS = TMP3D
            thread = MAPL_get_current_thread()
            workspace => self%workspaces(thread)
            call gw_intr_ncar(im * jm, LM, DT, self%NCAR_NRDG, &
                 workspace%beres_dc_desc, &
                 workspace%beres_band, workspace%oro_band, workspace%rdg_band, &
                 PLE, T, U, V, &
                 HT_dc, TMP3D, &
                 SGH, MXDIS, HWDTH, CLNGT, ANGLL, &
                 ANIXY, GBXAR_TMP, KWVRDG, EFFRDG, PREF, &
                 PMID, PDEL, RPDEL, PILN, ZM, lats, &
                 PHIS, &
                 DUDT_GWD_NCAR, DVDT_GWD_NCAR, DTDT_GWD_NCAR, &
                 DUDT_ORG_NCAR, DVDT_ORG_NCAR, DTDT_ORG_NCAR, &
                 TAUXO_TMP_NCAR, TAUYO_TMP_NCAR, &
                 TAUXB_TMP_NCAR, TAUYB_TMP_NCAR, &
                 self%NCAR_EFFGWORO, &
                 self%NCAR_EFFGWBKG, self%alpha, &
                 _RC)
         end if
         call MAPL_GridCompTimerStop(gc, "gw_intr_ncar", _RC)

         ! Use GEOS GWD only for Extratropical background sources...
         DUDT_GWD_GEOS = 0.0
         DVDT_GWD_GEOS = 0.0
         DTDT_GWD_GEOS = 0.0
         TAUXB_TMP_GEOS = 0.0
         TAUYB_TMP_GEOS = 0.0
         DUDT_ORG_GEOS = 0.0
         DVDT_ORG_GEOS = 0.0
         DTDT_ORG_GEOS = 0.0
         TAUXO_TMP_GEOS = 0.0
         TAUYO_TMP_GEOS = 0.0
         ! call MAPL_GridCompTimerStart(gc, "gw_intr_geos", _RC)
         if ((self%GEOS_EFFGWORO /= 0.0) .or. (self%GEOS_EFFGWBKG /= 0.0)) then
            call gw_intr(im * jm, LM, DT, &
                 self%GEOS_PGWV, &
                 PLE, T, U, V, SGH, PREF, &
                 PMID, PDEL, RPDEL, PILN, ZM, lats, &
                 DUDT_GWD_GEOS, DVDT_GWD_GEOS, DTDT_GWD_GEOS, &
                 DUDT_ORG_GEOS, DVDT_ORG_GEOS, DTDT_ORG_GEOS, &
                 TAUXO_TMP_GEOS, TAUYO_TMP_GEOS, TAUXO_3D, TAUYO_3D, FEO_3D, &
                 TAUXB_TMP_GEOS, TAUYB_TMP_GEOS, TAUXB_3D, TAUYB_3D, FEB_3D, &
                 FEPO_3D, FEPB_3D, DUBKGSRC, DVBKGSRC, DTBKGSRC, &
                 self%GEOS_BGSTRESS, &
                 self%GEOS_EFFGWORO, &
                 self%GEOS_EFFGWBKG, &
                 _RC)
         end if
         ! call MAPL_GridCompTimerStop(gc, "gw_intr_geos", _RC)

         ! Total
         DUDT_GWD = DUDT_GWD_GEOS + DUDT_GWD_NCAR
         DVDT_GWD = DVDT_GWD_GEOS + DVDT_GWD_NCAR
         DTDT_GWD = DTDT_GWD_GEOS + DTDT_GWD_NCAR
         ! Background
         TAUXB_TMP = TAUXB_TMP_GEOS + TAUXB_TMP_NCAR
         TAUYB_TMP = TAUYB_TMP_GEOS + TAUYB_TMP_NCAR
         ! Orographic
         DUDT_ORG = DUDT_ORG_GEOS + DUDT_ORG_NCAR
         DVDT_ORG = DVDT_ORG_GEOS + DVDT_ORG_NCAR
         DTDT_ORG = DTDT_ORG_GEOS + DTDT_ORG_NCAR
         TAUXO_TMP = TAUXO_TMP_GEOS + TAUXO_TMP_NCAR
         TAUYO_TMP = TAUYO_TMP_GEOS + TAUYO_TMP_NCAR
         !call MAPL_TimerOff(MAPL,"-INTR")

         call POSTINTR(im * jm, LM, DT, H0, HH, Z1, TAU1, &
              PREF, &
              PDEL, &
              U, &
              V, &
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
              PEGWD_X, &
              PEORO_X, &
              PERAY_X, &
              PEBKG_X, &
              KEGWD_X, &
              KEORO_X, &
              KERAY_X, &
              KEBKG_X, &
              KERES_X, &
              BKGERR_X)

         !! Tendency diagnostics
         if (associated(DUDT)) DUDT = DUDT_TOT
         if (associated(DVDT)) DVDT = DVDT_TOT
         if (associated(DTDT)) DTDT = DTDT_TOT * PDEL ! DTDT has to be pressure weighted for dynamics

         if (associated(DUDT_RAY)) DUDT_RAY = DUDT_RAH
         if (associated(DVDT_RAY)) DVDT_RAY = DVDT_RAH
         if (associated(DTDT_RAY)) DTDT_RAY = DTDT_RAH

         !! KE diagnostics
         if (associated(PEGWD)) PEGWD = PEGWD_X
         if (associated(PEORO)) PEORO = PEORO_X
         if (associated(PERAY)) PERAY = PERAY_X
         if (associated(PEBKG)) PEBKG = PEBKG_X
         if (associated(KEGWD)) KEGWD = KEGWD_X
         if (associated(KEORO)) KEORO = KEORO_X
         if (associated(KERAY)) KERAY = KERAY_X
         if (associated(KEBKG)) KEBKG = KEBKG_X
         if (associated(KERES)) KERES = KERES_X
         if (associated(BKGERR)) BKGERR = BKGERR_X

         !! Tendency diagnostics
         if (associated(DUDT_ORO)) DUDT_ORO = DUDT_ORG
         if (associated(DVDT_ORO)) DVDT_ORO = DVDT_ORG
         if (associated(DTDT_ORO)) DTDT_ORO = DTDT_ORG

         if (associated(DUDT_BKG)) DUDT_BKG = DUDT_GWD - DUDT_ORG
         if (associated(DVDT_BKG)) DVDT_BKG = DVDT_GWD - DVDT_ORG
         if (associated(DTDT_BKG)) DTDT_BKG = DTDT_GWD - DTDT_ORG

         ! Orographic stress
         if (associated(TAUGWX)) TAUGWX = TAUXO_TMP + TAUXB_TMP
         if (associated(TAUGWY)) TAUGWY = TAUYO_TMP + TAUYB_TMP
         if (associated(TAUOROX)) TAUOROX = TAUXO_TMP
         if (associated(TAUOROY)) TAUOROY = TAUYO_TMP
         if (associated(TAUBKGX)) TAUBKGX = TAUXB_TMP
         if (associated(TAUBKGY)) TAUBKGY = TAUYB_TMP

         ! Export unweighted T Tendency
         if (associated(TTMGW)) TTMGW = DTDT_TOT

         ! Fill additional exports
         if (associated(Q_EXP)) Q_EXP = Q
         if (associated(U_EXP)) U_EXP = U + DUDT_TOT * DT
         if (associated(V_EXP)) V_EXP = V + DVDT_TOT * DT
         if (associated(T_EXP)) T_EXP = T + DTDT_TOT * DT
         if (associated(PREF_EXP)) PREF_EXP = PREF
         if (associated(SGH_EXP)) SGH_EXP = SGH
         if (associated(PLE_EXP)) PLE_EXP = PLE

         if (DEBUG_GWD) then
            ! [TODO] Purnendu: The MAPL_MaxMin calls are now different in MAPL3
            !if(associated( T_EXP )) call MAPL_MaxMin('GWD: T_AF_GWD ', T_EXP)
            !if(associated( U_EXP )) call MAPL_MaxMin('GWD: U_AF_GWD ', U_EXP)
            !if(associated( V_EXP )) call MAPL_MaxMin('GWD: V_AF_GWD ', V_EXP)
         end if

         if (allocated(scratch_ridge)) deallocate(scratch_ridge)

         if (associated(T_EXP) .and. DEBUG_TQ_ERRORS) then
            do L = 1, LM
               do j = 1, jm
                  do i = 1, im
                     if (T_EXP(i, j, L) > 333.0) then
                        print *, "Temperature spike detected : ", T_EXP(i, j, L)
                        print *, "    GWD TOT Temp Increment : ", DTDT_GWD(i, j, L) * DT
                        print *, "    GWD ORO Temp Increment : ", DTDT_ORG(i, j, L) * DT
                        print *, "    GWD BKG Temp Increment : ", (DTDT_GWD(i, j, L) - DTDT_ORG(i, j, L)) * DT
                        print *, "    AFTER GWD Parameterization"
                        print *, "  Latitude       =", lats(i, j) * 180.0 / MAPL_PI
                        print *, "  Longitude      =", lons(i, j) * 180.0 / MAPL_PI
                        print *, "  Pressure (mb)  =", PMID(i, j, L) / 100.0
                        if (associated(U_EXP) .and. associated(V_EXP)) then
                           print *, "            UWND =", U_EXP(i, j, L)
                           print *, "            VWND =", V_EXP(i, j, L)
                        end if
                     end if
                  end do ! im loop
               end do ! jm loop
            end do ! LM loop
         end if

         ! All done
         _RETURN(_SUCCESS)
      end subroutine Gwd_Driver

   end subroutine Run

   subroutine GEOPOTENTIAL(pcols, pver, PILN, PMLN, pint, PMID, PDEL, RPDEL, T, Q, ZI, ZM)

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
      integer, intent(in) :: pcols ! Number of longitudes
      integer, intent(in) :: pver ! Number of vertical layers

      real, intent(in) :: PILN(pcols, pver + 1) ! Log interface pressures
      real, intent(in) :: PMLN(pcols, pver) ! Log midpoint pressures
      real, intent(in) :: pint(pcols, pver + 1) ! Interface pressures
      real, intent(in) :: PMID(pcols, pver) ! Midpoint pressures
      real, intent(in) :: PDEL(pcols, pver) ! layer thickness
      real, intent(in) :: RPDEL(pcols, pver) ! inverse of layer thickness
      real, intent(in) :: T(pcols, pver) ! temperature
      real, intent(in) :: Q(pcols, pver) ! specific humidity

      ! Output arguments

      real, intent(out) :: ZI(pcols, pver + 1) ! Height above surface at interfaces
      real, intent(out) :: ZM(pcols, pver) ! Geopotential height at mid level
      !
      !---------------------------Local variables-----------------------------
      !
      logical :: fvdyn ! finite volume dynamics
      integer :: i, K ! Lon, level indices
      real :: hkk ! diagonal element of hydrostatic matrix
      real :: hkl ! off-diagonal element
      real :: tv ! virtual temperature
      real :: tvfac ! Tv/T

      real, parameter :: ROG = MAPL_RGAS / MAPL_GRAV
      !
      !-----------------------------------------------------------------------
      !

      ! Set dynamics flag

      fvdyn = .true.

      ! The surface height is zero by definition.

      I_LOOP: do i = 1, pcols

         ZI(i, pver + 1) = 0.0

         ! Compute zi, zm from bottom up.
         ! Note, zi(i,k) is the interface above zm(i,k)

         do K = pver, 1, -1

            ! First set hydrostatic elements consistent with dynamics

            if (fvdyn) then
               hkl = PILN(i, K + 1) - PILN(i, K)
               hkk = PILN(i, K + 1) - PMLN(i, K)
            else
               hkl = PDEL(i, K) / PMID(i, K)
               hkk = 0.5 * hkl
            end if

            ! Now compute tv, zm, zi

            tvfac = 1. + MAPL_VIREPS * Q(i, K)
            tv = T(i, K) * tvfac

            ZM(i, K) = ZI(i, K + 1) + ROG * tv * hkk
            ZI(i, K) = ZI(i, K + 1) + ROG * tv * hkl
         end do
      end do I_LOOP

      return
   end subroutine GEOPOTENTIAL

   !-----------------------------------------------------------------------

   subroutine PREGEO(pcols, pver, PLE, lats, PMID, PDEL, RPDEL, PILN, PMLN)

      implicit none

      !------------------------------Arguments--------------------------------
      !
      ! Input arguments
      !

      integer, intent(in) :: pcols ! Number of longitudes
      integer, intent(in) :: pver ! Number of vertical layers

      real, intent(in) :: PLE(pcols, pver + 1) ! Interface pressures
      real, intent(in) :: lats(pcols) ! latitude in radian

      ! Output arguments

      real, intent(out) :: PMID(pcols, pver) ! Midpoint pressures
      real, intent(out) :: PDEL(pcols, pver) ! layer thickness
      real, intent(out) :: RPDEL(pcols, pver) ! inverse of layer thickness
      real, intent(out) :: PILN(pcols, pver + 1) ! Log interface pressures
      real, intent(out) :: PMLN(pcols, pver) ! Log midpoint pressures

      !
      !---------------------------Local variables-----------------------------
      !
      integer :: i, K

      real :: hvsd ! Efficiency factor

      real, parameter :: PI_GWD = 4.0 * atan(1.0) ! This is *not* MAPL_PI

      !
      !-----------------------------------------------------------------------
      !

      ! Form pressure factors
      !----------------------

      I_LOOP: do i = 1, pcols

         do K = 1, pver
            PMID(i, K) = 0.5 * (PLE(i, K) + PLE(i, K + 1))
            PDEL(i, K) = PLE(i, K + 1) - PLE(i, K)
            RPDEL(i, K) = 1.0 / PDEL(i, K)
            PILN(i, K) = log(PLE(i, K))
            PMLN(i, K) = log(PMID(i, K)) !
         end do
         PILN(i, pver + 1) = log(PLE(i, pver + 1))
      end do I_LOOP

   end subroutine PREGEO

   subroutine POSTINTR(pcols, pver, DT, H0, HH, Z1, TAU1, &
        PREF, &
        PDEL, &
        U, &
        V, &
        DUDT_GWD, &
        DVDT_GWD, &
        DTDT_GWD, &
        DUDT_ORG, &
        DVDT_ORG, &
        DTDT_ORG, &

   ! Outputs
        DUDT_TOT, &
        DVDT_TOT, &
        DTDT_TOT, &
        DUDT_RAH, &
        DVDT_RAH, &
        DTDT_RAH, &
        PEGWD, &
        PEORO, &
        PERAY, &
        PEBKG, &
        KEGWD, &
        KEORO, &
        KERAY, &
        KEBKG, &
        KERES, &
        BKGERR)

      implicit none

      !------------------------------Arguments--------------------------------
      !
      ! Input arguments
      !

      integer, intent(in) :: pcols ! Number of longitudes
      integer, intent(in) :: pver ! Number of vertical layers
      real, intent(in) :: DT ! Time step
      real, intent(in) :: H0, HH, Z1, TAU1 ! Rayleigh friction parameters

      real, intent(in) :: PREF(pver + 1)
      real, intent(in) :: PDEL(pcols, pver)
      real, intent(in) :: U(pcols, pver)
      real, intent(in) :: V(pcols, pver)

      real, intent(in) :: DUDT_GWD(pcols, pver)
      real, intent(in) :: DVDT_GWD(pcols, pver)
      real, intent(in) :: DTDT_GWD(pcols, pver)
      real, intent(in) :: DUDT_ORG(pcols, pver)
      real, intent(in) :: DVDT_ORG(pcols, pver)
      real, intent(in) :: DTDT_ORG(pcols, pver)

      real, intent(out) :: DUDT_TOT(pcols, pver)
      real, intent(out) :: DVDT_TOT(pcols, pver)
      real, intent(out) :: DTDT_TOT(pcols, pver)
      real, intent(out) :: DUDT_RAH(pcols, pver)
      real, intent(out) :: DVDT_RAH(pcols, pver)
      real, intent(out) :: DTDT_RAH(pcols, pver)
      real, intent(out) :: PEGWD(pcols)
      real, intent(out) :: PEORO(pcols)
      real, intent(out) :: PERAY(pcols)
      real, intent(out) :: PEBKG(pcols)
      real, intent(out) :: KEGWD(pcols)
      real, intent(out) :: KEORO(pcols)
      real, intent(out) :: KERAY(pcols)
      real, intent(out) :: KEBKG(pcols)
      real, intent(out) :: KERES(pcols)
      real, intent(out) :: BKGERR(pcols)

      !
      !---------------------------Local variables-----------------------------
      !
      integer :: i, K
      real :: ZREF, KRAY
      !
      !-----------------------------------------------------------------------
      !

      I_LOOP: do i = 1, pcols

         PEGWD(i) = 0.0
         PEORO(i) = 0.0
         PERAY(i) = 0.0
         PEBKG(i) = 0.0
         KEGWD(i) = 0.0
         KEORO(i) = 0.0
         KERAY(i) = 0.0
         KEBKG(i) = 0.0
         KERES(i) = 0.0
         BKGERR(i) = 0.0

         do K = 1, pver

            ! Rayleigh friction
            !------------------
            if (TAU1 > 0.0) then
               ZREF = H0 * log(MAPL_P00 / (0.5 * (PREF(K) + PREF(K + 1))))
               KRAY = (1.0 / TAU1) * (1.0 - TANH((Z1 - ZREF) / HH))
               KRAY = KRAY / (1 + DT * KRAY)
               DUDT_RAH(i, K) = -U(i, K) * KRAY
               DVDT_RAH(i, K) = -V(i, K) * KRAY
               DTDT_RAH(i, K) = -((U(i, K) + (0.5 * DT) * DUDT_RAH(i, K)) * DUDT_RAH(i, K) + &
                    (V(i, K) + (0.5 * DT) * DVDT_RAH(i, K)) * DVDT_RAH(i, K)) * (1.0 / MAPL_CP)
            else
               DUDT_RAH(i, K) = 0.0
               DVDT_RAH(i, K) = 0.0
               DTDT_RAH(i, K) = 0.0
            end if

            DUDT_TOT(i, K) = DUDT_RAH(i, K) + DUDT_GWD(i, K)
            DVDT_TOT(i, K) = DVDT_RAH(i, K) + DVDT_GWD(i, K)
            DTDT_TOT(i, K) = DTDT_RAH(i, K) + DTDT_GWD(i, K)

            ! KE dIagnostics
            !----------------

            PEGWD(i) = PEGWD(i) + DTDT_TOT(i, K) * PDEL(i, K) * (MAPL_CP / MAPL_GRAV)
            PEORO(i) = PEORO(i) + DTDT_ORG(i, K) * PDEL(i, K) * (MAPL_CP / MAPL_GRAV)
            PERAY(i) = PERAY(i) + DTDT_RAH(i, K) * PDEL(i, K) * (MAPL_CP / MAPL_GRAV)
            PEBKG(i) = PEBKG(i) + (DTDT_GWD(i, K) - DTDT_ORG(i, K)) * PDEL(i, K) * (MAPL_CP / MAPL_GRAV)

            KEGWD(i) = KEGWD(i) + ((U(i, K) + (0.5 * DT) * DUDT_TOT(i, K)) * DUDT_TOT(i, K) + &
                 (V(i, K) + (0.5 * DT) * DVDT_TOT(i, K)) * DVDT_TOT(i, K)) * PDEL(i, K) * (1.0 / MAPL_GRAV)

            KEORO(i) = KEORO(i) + ((U(i, K) + (0.5 * DT) * DUDT_ORG(i, K)) * DUDT_ORG(i, K) + &
                 (V(i, K) + (0.5 * DT) * DVDT_ORG(i, K)) * DVDT_ORG(i, K)) * PDEL(i, K) * (1.0 / MAPL_GRAV)

            KERAY(i) = KERAY(i) + ((U(i, K) + (0.5 * DT) * DUDT_RAH(i, K)) * DUDT_RAH(i, K) + &
                 (V(i, K) + (0.5 * DT) * DVDT_RAH(i, K)) * DVDT_RAH(i, K)) * PDEL(i, K) * (1.0 / MAPL_GRAV)

            KEBKG(i) = KEBKG(i) &
                 + ((U(i, K) + (0.5 * DT) * (DUDT_GWD(i, K) - DUDT_ORG(i, K))) * (DUDT_GWD(i, K) - DUDT_ORG(i, K)) + &
                 (V(i, K) + (0.5 * DT) * (DVDT_GWD(i, K) - DVDT_ORG(i, K))) * (DVDT_GWD(i, K) - DVDT_ORG(i, K))) * &
                 PDEL(i, K) * (1.0 / MAPL_GRAV)
         end do

         BKGERR(i) = -(PEBKG(i) + KEBKG(i))
         KERES(i) = PEGWD(i) + KEGWD(i) + BKGERR(i)

      end do I_LOOP

   end subroutine POSTINTR

   subroutine Write_Profile(avar, area, grid, NAME)
      type(ESMF_Grid), intent(in) :: grid
      real, intent(in) :: avar(:, :)
      real, intent(in) :: area(:, :)
      character(len=*), intent(in) :: NAME

      real(kind=ESMF_KIND_R8), allocatable :: locArr(:, :)
      real(kind=ESMF_KIND_R8), allocatable :: glbArr(:, :)
      real, allocatable :: area_global(:, :)
      real, allocatable :: avar_global(:, :)
      real :: rng(3)
      integer :: DIMS(3), im, jm, status, rc
      integer, allocatable :: global_dims(:)

      call MAPL_GridGet(grid, im=im, jm=jm, _RC)
      DIMS(1:2) = [im, jm]
      allocate(locArr(DIMS(1), DIMS(2)))

      call mapl_GridGetGlobalCellCountPerDim(grid, globalCellCountPerDim=global_dims, _RC)
      DIMS(1:size(global_dims)) = global_dims
      allocate(glbArr(DIMS(1), DIMS(2)))
      allocate(area_global(DIMS(1), DIMS(2)))
      allocate(avar_global(DIMS(1), DIMS(2)))

#if 1
      locArr = avar
      call MAPL_ArrayGather(locArr, glbArr, grid)
      avar_global = glbArr

      locArr = area
      call MAPL_ArrayGather(locArr, glbArr, grid)
      area_global = glbArr

      if (MAPL_AM_I_ROOT()) then
         rng(1) = MINVAL(MINVAL(avar_global, DIM=1), DIM=1)
         rng(2) = MAXVAL(MAXVAL(avar_global, DIM=1), DIM=1)
         rng(3) = SUM(SUM(avar_global * area_global, DIM=1), DIM=1) / &
              SUM(SUM(area_global, DIM=1), DIM=1)
         write(*, '(A," ",3(f21.9,1x))') trim(NAME), rng(:)
      end if
#else
      rng(1) = MINVAL(MINVAL(avar, DIM=1), DIM=1)
      rng(2) = MAXVAL(MAXVAL(avar, DIM=1), DIM=1)
      rng(3) = SUM(SUM(avar * area, DIM=1), DIM=1) / &
           SUM(SUM(area, DIM=1), DIM=1)
      write(*, '(A," ",3(f21.9,1x))'), trim(NAME), rng(:)
#endif

      deallocate(locArr)
      deallocate(glbArr)
      deallocate(area_global)
      deallocate(avar_global)

   end subroutine Write_Profile

end module GEOS_GwdGridCompMod

subroutine GWD_SetServices(gc, rc)
   use esmf
   use GEOS_GwdGridCompMod, only : mySetservices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetservices(gc, rc=rc)
end subroutine GWD_SetServices
