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

   !USES:

   use ESMF
   use mapl3g_Generic, only: MAPL_GridCompSetEntryPoint
   use mapl3g_Generic, only: MAPL_GridCompGet, MAPL_GridCompGetInternalState, MAPL_GridCompGetResource
   use mapl3g_Geom_API, only: MAPL_GridGet, MAPL_GridGetCoordinates
   use mapl3g_State_API, only: MAPL_StateGetPointer

   ! use MAPL, only: MAPL_AddImportSpec, MAPL_AddExportSpec, MAPL_AddInternalSpec
   ! use MAPL, only: MAPL_DimsHorzVert, MAPL_DimsHorzOnly
   ! use MAPL, only: MAPL_VLocationCenter, MAPL_VLocationEdge, MAPL_VLocationNone
   ! use MAPL, only: MAPL_RestartSkip
   ! use MAPL, only: MAPL_GridCompSetEntryPoint, MAPL_GenericSetServices, MAPL_GenericInitialize
   ! use MAPL, only: MAPL_MetaComp
   ! use MAPL, only: MAPL_Get, MAPL_GetObjectFromGC, MAPL_GetPointer, MAPL_GetResource
   ! use MAPL, only: MAPL_VerifyFriendly
   ! use MAPL, only: MAPL_TimerAdd, MAPL_TimerOn, MAPL_TimerOff
   use MAPL, only: MAPL_Verify, MAPL_Return, MAPL_ASRT
   use MAPL_Constants, only: MAPL_PI, MAPL_GRAV, MAPL_P00, MAPL_KAPPA, MAPL_RGAS, MAPL_CP

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:

   public SetServices

   !EOP

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !DESCRIPTION:  This version uses the MAPL\_GenericSetServices, which sets
   !              the Initialize and Finalize services, as well as allocating
   !   our instance of the MAPL\_MetaComp and putting it in the
   !   gridded component (GC). The MAPL\_MetaComp contains an ESMF state to use as
   !   an internal state. This routine describes the contents of this Internal state,
   !   as well as of the conventional Import and Export states by making call
   !   to MAPL.
   !INTERFACE:
   subroutine SetServices(gc, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      integer, optional, intent(out) :: rc
      !EOP

      integer :: status
      character(len=ESMF_MAXSTR) :: comp_name

      !IROUTINE: State Descriptions
      !DESCRIPTION: The component uses all three states (Import, Export
      !  and Internal). There is no Private (non-ESMF) Internal state. All
      !  three are managed by MAPL. The Internal state contains only invariant,
      !  horizontally dependent quantities that are set in Initialize, therefore
      !  it never needs checkpointing.

      ! Set state variable specs
#include "HS_Import___.h"
#include "HS_Export___.h"
#include "HS_Internal___.h"

      ! Set the Initialize and Run entry points
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize
   !DESCRIPTION: Here we initialize the internal state, which in HS contains
   !   two-dimensional invariant arrays used in the forcing calculations. They
   !   are done here simply for economy. Initialize calls MAPL\_GenericInitialize.
   !   If the import state needs to be restarted, MAPL will do it
   !   if a restart file is provided.

   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, optional, intent(out) :: rc
      !EOP

      type(ESMF_Grid) :: grid
      type(ESMF_State) :: internal
      real, pointer :: lats(:, :), lons(:, :)
      real :: dx, dy, x0, y0, afac, phi0, qmax
      ! Pointers to internals
#include "HS_DeclarePointer___.h"
      integer :: status

      ! Get coordinate information
      call MAPL_GridCompGet(gc, grid=grid, _RC)
      call MAPL_GridGetCoordinates(grid, longitudes=lons, latitudes=lats, _RC)

      ! Get internal state
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get pointers to internal variables
#include "HS_GetPointer___.h"

      ! Initialize geometric factors
      SPHI2 = sin(lats)**2
      CPHI2 = cos(lats)**2

      ! Precompute Local heating distribution

      ! K day-1 :: Local heating
      call MAPL_GridCompGetResource(gc, 'QMAX', qmax, default=.000, _RC)

      if (qmax /= 0.0) then

         ! Get local heating parameters from the configuration
         call MAPL_GridCompGetResource(gc, 'X0', x0, default=0.0, _RC)
         call MAPL_GridCompGetResource(gc, 'Y0', y0, default=0.0, _RC)
         call MAPL_GridCompGetResource(gc, 'DX', dx, default=30.0, _RC)
         call MAPL_GridCompGetResource(gc, 'DY', dy, default=5.0, _RC)

         ! Horizontal structure of optional stationary heat source
         x0 = x0 * (MAPL_PI / 180.0)
         y0 = y0 * (MAPL_PI / 180.0)
         dx = dx * (MAPL_PI / 180.0)
         dy = dy * (MAPL_PI / 180.0)

         HFCN = ((lons - x0) / dx)**2
         HFCN = min(HFCN, ((lons + 360. - x0) / dx)**2)
         HFCN = min(HFCN, ((lons - 360. - x0) / dx)**2)

         HFCN = HFCN + ((lats - y0) / dy)**2

         where (HFCN < 10.)
            HFCN = exp(-0.5 * HFCN)
            elsewhere
            HFCN = 0.0
         end where
      else
         HFCN = 0.0
      end if

      ! Get Williamson parameters (alpha, phi0) from the configuration
      call MAPL_GridCompGetResource(gc, 'AFAC', afac, default=2.65 / 15., _RC)
      call MAPL_GridCompGetResource(gc, 'PHI0', phi0, default=60., _RC)

      ! Horizontal structure of Willianson parameterization
      P_I = (1.0 + tanh(afac * (abs(lats) - phi0)))

      _RETURN(_SUCCESS)
   end subroutine Initialize

   !BOP
   !IROUTINE: Run
   !DESCRIPTION: The Run method of the HS Gridded Component. It computes
   !  tendencies for temperature and winds as described above,
   !  as well as a number of
   !  diagnostics. These are all optional Exports. The input winds
   !  and/or temperatures may be Friendly. All Friendlies are updated
   !  using a time step
   !  obtained from 'RUN\_DT:' in the configuration. Both U and V must be
   !  Friendly for the wind to be updated. It is an error for one to be Friendly
   !  and not the other.

   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)
      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, optional, intent(out) :: rc
      !EOP

      type(ESMF_State) :: internal
      type(ESMF_Grid) :: grid
      type(ESMF_Logical) :: friendly

      ! Pointers to imports/internals/exports
#include "HS_DeclarePointer___.h"

      ! Scratch arrays and working pointers
      real, allocatable, dimension(:, :) :: pii, dp, pl, uu, vv, vr, te, f1, rr, ds, dm, pk
      real, pointer, dimension(:, :) :: ps, pt
      integer :: level, im, jm, lm, fricq, status
      logical :: friendly_temp, friendly_wind
      real :: dt, ka, ks, kf

      ! 8 Held-Suarez parameters
      real :: sig1, taua, taus, tauf, delv1, delh, t0, tstrt

      ! Localized heating parameters
      real :: p_1, qmax

      ! 3 Williamson parameters
      real :: p_d, gam_i, gam_d

      real, parameter :: DAYLEN = 86400.

      ! Get esmf internal state from generic state.
      call MAPL_GridCompGet(gc, grid=grid, num_levels=lm, _RC)
      call MAPL_GridGet(grid, im=im, jm=jm, _RC)

      ! Pointers to internals
      call MAPL_GridCompGetInternalState(gc, internal, _RC)
#include "HS_GetPointer___.h"

      ! Get parameters from the configuration

      ! 8 Held-Suarez parameters
      call MAPL_GridCompGetResource(gc, 'TAUA', taua, default=40., _RC)
      call MAPL_GridCompGetResource(gc, 'TAUS', taus, default=4.0, _RC)
      call MAPL_GridCompGetResource(gc, 'TAUF', tauf, default=1.0, _RC)
      call MAPL_GridCompGetResource(gc, 'SIG1', sig1, default=0.7, _RC)
      call MAPL_GridCompGetResource(gc, 'T0', t0, default=315., _RC)
      call MAPL_GridCompGetResource(gc, 'DELV1', delv1, default=10., _RC)
      call MAPL_GridCompGetResource(gc, 'DELH', delh, default=60., _RC)
      call MAPL_GridCompGetResource(gc, 'TSTRT', tstrt, default=200., _RC)
      call MAPL_GridCompGetResource(gc, 'FRICQ', fricq, default=1, _RC)

      ! Localized heating parameters
      call MAPL_GridCompGetResource(gc, 'P_1', p_1, default=0.0, _RC)
      call MAPL_GridCompGetResource(gc, 'QMAX', qmax, default=0.0, _RC)

      ! 3 Williamson parameters
      call MAPL_GridCompGetResource(gc, 'GAM_I', gam_i, default=-3.345e-3, _RC)
      call MAPL_GridCompGetResource(gc, 'GAM_D', gam_d, default=.002, _RC)
      call MAPL_GridCompGetResource(gc, 'P_D', p_d, default=1.E4, _RC)
      !EOR

      ! Check for Friendliness
      friendly = ESMF_TRUE ! MAPL_VerifyFriendly(import, 'V', trim(comp_name), _RC)
      friendly_wind = friendly == ESMF_TRUE

      friendly = ESMF_TRUE ! MAPL_VerifyFriendly(import, 'U', trim(comp_name), _RC)
      ASSERT_((friendly == ESMF_TRUE) .eqv. friendly_wind)

      friendly = ESMF_TRUE ! MAPL_VerifyFriendly(import, 'TEMP', trim(comp_name), _RC)
      friendly_temp = friendly == ESMF_TRUE

      ! If we have a friendly, we need a time step
      if (friendly_temp .or. friendly_wind) then
         call MAPL_GridCompGetResource(gc, 'RUN_DT', dt, default=0.0, _RC)
      end if

      ! Allocate 10 2D scratch arrays
      allocate(dp(im, jm), _STAT)
      allocate(pl(im, jm), _STAT)
      allocate(uu(im, jm), _STAT)
      allocate(vv(im, jm), _STAT)
      allocate(vr(im, jm), _STAT)
      allocate(te(im, jm), _STAT)
      allocate(ds(im, jm), _STAT)
      allocate(pii(im, jm), _STAT)
      allocate(f1(im, jm), _STAT)
      allocate(rr(im, jm), _STAT)
      allocate(dm(im, jm), _STAT)
      allocate(pk(im, jm), _STAT)

      ! Begin calculations.
      ps => PLE(:, :, lm)
      pt => PLE(:, :, 0)

      pii = p_d - (p_d - pt) * 0.5 * P_I

      ! Initialize vertically integrated diagnostics
      if (associated(DISS)) DISS = 0.0
      if (associated(TAUX)) TAUX = 0.0
      if (associated(TAUY)) TAUY = 0.0

      ! Loop invariants
      ka = 1.0 / (DAYLEN * taua)
      ks = 1.0 / (DAYLEN * taus)
      kf = 1.0 / (DAYLEN * tauf)

      LEVELS: do level = 1, lm

         dp = (PLE(:, :, level) - PLE(:, :, level - 1))
         pl = (PLE(:, :, level) + PLE(:, :, level - 1)) * 0.5
         dm = dp / MAPL_GRAV
         pk = (pl / MAPL_P00)**MAPL_KAPPA

         ! H&S equilibrium temperature
         te = pk * (t0 - delh * SPHI2 - delv1 * CPHI2 * log(pl / MAPL_P00))
         te = max(te, tstrt)

         ! Williamson Stratospheric modifications to equilibrium temperature
         where (pl < p_d)
            te = tstrt * (min(1.0, pl / p_d)**(MAPL_RGAS * gam_d / MAPL_GRAV) &
                 + min(1.0, pl / pii)**(MAPL_RGAS * gam_i / MAPL_GRAV) - 1)
         end where

         !  Exports of equilibrium T and Theta
         if (associated(T_EQ)) T_EQ(:, :, level) = te
         if (associated(THEQ)) THEQ(:, :, level) = te / pk

         ! Vertical structure of timescales in H&S.
         f1 = max(0.0, ((pl / ps) - sig1) / (1.0 - sig1))

         ! Atmospheric heating from H&S
         rr = (ka + (ks - ka) * f1 * CPHI2**2) * (te - TEMP(:, :, level))

         if (associated(DTDT)) DTDT(:, :, level) = dp * rr
         if (friendly_temp) TEMP(:, :, level) = TEMP(:, :, level) + dt * rr

         ! Wind tendencies
         uu = -U(:, :, level) * (f1 * kf)
         vv = -V(:, :, level) * (f1 * kf)

         if (associated(DUDT)) DUDT(:, :, level) = uu
         if (associated(DVDT)) DVDT(:, :, level) = vv

         if (friendly_wind) then
            U(:, :, level) = U(:, :, level) + dt * uu
            V(:, :, level) = V(:, :, level) + dt * vv
         end if

         !  Frictional heating from H&S drag
         ds = U(:, :, level) * uu + V(:, :, level) * vv

         if (associated(DISS)) DISS = DISS - ds * dm
         if (fricq /= 0) then
            if (associated(DTDT)) DTDT(:, :, level) = DTDT(:, :, level) - ds * (dp / MAPL_CP)
            if (friendly_temp) TEMP(:, :, level) = TEMP(:, :, level) - ds * (dt / MAPL_CP)
         end if

         !  Surface stresses from vertically integrated H&S surface drag
         if (associated(TAUX)) TAUX = TAUX - uu * dm
         if (associated(TAUY)) TAUY = TAUY - vv * dm

         ! Localized heat source, if any
         if ((associated(DTDT) .or. FriendlyTemp) .and. qmax/=0.0) then
            where (pl > p_1)
               vr = HFCN * (qmax / DAYLEN) * sin(MAPL_PI * (ps - pl) / (ps - p_1))
               elsewhere
               vr = 0.
            end where

            if (associated(DTDT)) DTDT(:, :, level) = DTDT(:, :, level) + dp * vr
            if (FriendlyTemp) TEMP(:, :, level) = TEMP(:, :, level) + dt * vr
         end if

      end do LEVELS

      _RETURN(_SUCCESS)
   end subroutine Run

end module GEOS_HSGridCompMod
