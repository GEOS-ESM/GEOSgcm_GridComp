!  $Id$

!------------------------------------------------------------------------------
!BOP
! !ROUTINE:  FVdycore_wrapper --- Wrapper for NASA finite-volume dynamical core
!
! !INTERFACE:
   subroutine FVdycore_wrapper( phisxy, txy,  qqq, STATE,               &
                                pkxy, pelnxz, oma_xy,                   &
                                convcpt, convthv, epvxyz, cxxyz, cyxyz, &
                                ptfxxyz, ptfyxyz,                       &
                                mfxxyz_ur, mfyxyz_ur,                   &
                                mfxxyz, mfyxyz, mfzxyz, convt,          &
                                kenrg1, penrg1, tenrg1,                 &
                                kenrg2, penrg2, tenrg2,                 &
                                dkedtad, dkedtpg, dkedtdp, dkedtho,     &
                                dthdtremap, dthdtconsv, dtmp,           &
                                P00, PI, CP, KAPPA, OMEGA,              &
                                RADIUS, GRAV, RGAS, EPS,                &
                                NKE, NPHI, CONSV, FILL                  )


! !USES:
   use shr_kind_mod,  only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use mod_comm, only: commglobal, mp_swapirr, mp_swapirr_r4,           &
                       mp_barrier, mp_send4d_ns, mp_recv4d_ns
   use dynamics_vars, only : T_TRACERS, T_FVDYCORE_VARS,        &
                             T_FVDYCORE_GRID, T_FVDYCORE_STATE, &
                             d2a3d
   use FVperf_module, only : FVstartclock, FVstopclock
   use diag_module,   only : compute_vdot_gradp
   use sw_core,       only : d2a2c_winds

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   type (T_FVDYCORE_STATE), intent(INOUT), target :: STATE

! !INPUT PARAMETERS:
! Surface geopotential
   real(r8), target     :: phisxy(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                  STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
! Dry air temperature
   real(r8), intent(in) :: txy(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                               STATE%GRID%jfirstxy:STATE%GRID%jlastxy, &
                               STATE%GRID%km )

   type(T_TRACERS)             :: qqq2    ! Specific Humidity (for benergy)
   type(T_TRACERS), intent(in) :: qqq     ! Specific Humidity (for benergy)
   logical, intent(in)         :: convt   ! true: output pt, the virtual temperature
                                          ! false: pt is updated

   real(r8), intent(in)        :: P00     ! Reference Surface Pressure
   real(r8), intent(in)        :: PI      ! Constants, passed as arguments for portability:
   real(r8), intent(in)        :: CP      ! Specific heat
   real(r8), intent(in)        :: KAPPA   ! Kappa (=2/7)
   real(r8), intent(in)        :: OMEGA   ! Angular velocity of earth
   real(r8), intent(in)        :: RADIUS  ! Radius of earth
   real(r8), intent(in)        :: GRAV    ! Gravitational acceleration 
   real(r8), intent(in)        :: RGAS    ! Dry air gas constant
   real(r8), intent(in)        :: EPS     ! Constant for the virtual effect

   integer,  intent(in)        :: NKE, NPHI ! OMEGA-ALPHA Tracer Indices (default: -999)


!
!  Other interesting stuff which the dynamical core can provide
!

! pe**kappa
   real(r8), intent(INOUT) :: pkxy(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                   STATE%GRID%km+1)

! log pressure (pe) at layer edges
   real(r8), intent(INOUT) :: pelnxz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                     STATE%GRID%km+1,&
                                     STATE%GRID%jfirstxy:STATE%GRID%jlastxy)
! ertel's potential vorticity ( K*m^2 / (kg*sec) )
   real(r8), intent(INOUT) :: epvxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                     STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                     STATE%GRID%km)
! zonal accumulated C-grid winds
   real(r8), intent(INOUT) :: cxxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                    STATE%GRID%km)
! meridonal accumulated C-grid winds
   real(r8), intent(INOUT) :: cyxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                    STATE%GRID%km)
! zonal total mass flux ( K * pa * m^2/s )
   real(r8), intent(INOUT) :: ptfxxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                      STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                      STATE%GRID%km)
! meridonal total mass flux ( K * pa * m^2/s )
   real(r8), intent(INOUT) :: ptfyxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                      STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                      STATE%GRID%km)

! zonal total mass flux ( pa * m^2/s )
   real(r8), intent(INOUT) :: mfxxyz_ur(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                        STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                        STATE%GRID%km)
! meridonal total mass flux ( pa * m^2/s )
   real(r8), intent(INOUT) :: mfyxyz_ur(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                        STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                        STATE%GRID%km)
! Remapped zonal total mass flux ( pa * m^2/s )
   real(r8), intent(INOUT) :: mfxxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                     STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                     STATE%GRID%km)
! Remapped meridonal total mass flux ( pa * m^2/s )
   real(r8), intent(INOUT) :: mfyxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                     STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                     STATE%GRID%km)
! Remapped vertical total mass flux ( pa * m^2/s )
   real(r8), intent(INOUT) :: mfzxyz(STATE%GRID%ifirstxy:STATE%GRID%ilastxy,&
                                     STATE%GRID%jfirstxy:STATE%GRID%jlastxy,&
                                     STATE%GRID%km+1)

! !OUTPUT PARAMETERS:
   real(r8), intent(OUT) ::  oma_xy(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy, &
                                    STATE%GRID%km                           )
   real(r8), intent(OUT) :: convcpt(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy, &
                                    STATE%GRID%km                           )
   real(r8), intent(OUT) :: convthv(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy, &
                                    STATE%GRID%km                           )

   real(r8), intent(OUT) :: kenrg1(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: penrg1(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: tenrg1(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: kenrg2(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: penrg2(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: tenrg2(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                   STATE%GRID%jfirstxy:STATE%GRID%jlastxy )

   real(r8), intent(OUT) :: dkedtad(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: dkedtpg(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: dkedtdp(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: dkedtho(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                    STATE%GRID%jfirstxy:STATE%GRID%jlastxy )

   real(r8), intent(OUT) :: dthdtremap(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                       STATE%GRID%jfirstxy:STATE%GRID%jlastxy )
   real(r8), intent(OUT) :: dthdtconsv(STATE%GRID%ifirstxy:STATE%GRID%ilastxy, &
                                       STATE%GRID%jfirstxy:STATE%GRID%jlastxy )

   real(r8), intent(OUT) :: dtmp ! Temperature Change due to CONSV=TRUE

! Developer: Shian-Jiann Lin, NASA/GSFC; email: lin@dao.gsfc.nasa.gov
!
! Top view of D-grid prognostatic variables: u, v, and delp (and other scalars)
!
!               u(i,j+1)
!                 |
!      v(i,j)---delp(i,j)---v(i+1,j)
!                 |
!               u(i,j)
!
! External routine required: the user needs to supply a subroutine to set up
!                            "Eulerian vertical coordinate" for remapping purpose.
!                             Currently this routine is named as set_eta()
!                             In principle any terrian following vertical
!                             coordinate can be used. The input to fvcore
!                             need not be on the same vertical coordinate
!                             as the output.
!                             If SPMD is defined the Pilgrim communication
!                             library developed by Will Sawyer will be needed.
!
! Remarks: values at poles for both u and v need not be defined; but values for
!          all other scalars needed to be defined at both poles (as polar cap mean
!          quantities). Tracer advection is done "off-line" using the
!          large time step. Consistency is maintained by using the time accumulated
!          Courant numbers and horizontal mass fluxes for the FFSL algorithm.
!          The input "pt" can be either dry potential temperature
!          defined as T/pkz (adiabatic case) or virtual potential temperature
!          defined as T*/pkz (full phys case). IF convt is true, pt is not updated.
!          Instead, virtual temperature is ouput.
!          ipt is updated if convt is false.
!          The user may set the value of nx to optimize the SMP performance
!          The optimal valuse of nx depends on the total number of available
!          shared memory CPUs per node (NS). Assuming the maximm MPI 
!          decomposition is used in the y-direction, set nx=1 if the
!          NS <=4; nx=4 if NS=16.
!
! !REVISION HISTORY:
!   WS  03.07.16:  From dynpkg
!   WS  03.08.06   Improvements w.r.t. newest FVdycore from CAM
!   WS  03.08.13   Added logic for 2d domain decomposition
!   WS  03.11.19   Integrated 1D decomposition case (data copies only)
!   WS  04.09.20   Transition away from dynamics_vars (code now reentrant)
!   WS  05.05.17   Split off into a separate file for CAM unification
!   WP  06.01.18   Added horizontal/vertical mass fluxes and mfz_comp
!   WS  06.05.17   Added compute_vdot_gradp for proper OMGA calculation
!   WS  06.06.26   Revised for newest benergy routine (conservation mode)
!   LT  07.06.04   Modified OMGA calculation
!   WS  09.04.01   Upgraded to PILGRIM from cam3_6_33
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Local variables

   type (T_FVDYCORE_GRID), pointer :: grid
   type (T_FVDYCORE_VARS), pointer :: vars

   integer i, j, k, iq          ! Loop indices
   integer            :: rc     ! return code
   integer            :: klastp
   integer            :: kord   ! parameter controlling monotonicity in mapping
                                ! recommendation: kord=4
   integer            :: te_method   ! parameter controlling total energy remapping
                                     ! recommendation: kord=1 (cubic interpolation)
   integer            :: ntotq  ! declared dimension of q3
   integer            :: ndt    ! the large time step in seconds
                                ! Also the mapping time step in this setup
!
! Variables which used to come from dynamics_vars
!
   integer            :: icd, jcd            ! C-grid algorithm order in X and Y
   integer            :: iord, jord          ! D-grid algorithm order in X and Y
   integer            :: im, jm, km, nq      ! global domain dimensions
   integer            :: ng_s, ng_c, ng_d    ! halo region sizes

   integer            :: jfirst, jlast, kfirst, klast  ! YZ decomp. ranges
   integer            :: ifirstxy, ilastxy, jfirstxy, jlastxy ! XY decomp. ranges
   integer            :: twod_decomp, myid_z, npr_z, npr_y

   real(r8) rcap      ! Radius of polar cap
   real(r8) ptop      ! Pressure at the top of atmosphere
   real(r8) sump      ! Summation for Pole Values
   real(r8) sumpg     ! Summation for Pole Values
   real(r8) sumdp     ! Summation for Pole Values
   real(r8) sumho     ! Summation for Pole Values
   real(r8) sumad     ! Summation for Pole Values

   integer            :: nsplit

   real(r8)   umax              ! Maximum winds, m/s
   parameter (umax = 300.0)
   real(r8), parameter        :: HALF                 = 0.5_r8
   real(r8), parameter        :: ONE                  = 1.0_r8

   logical    consv             ! Flag to force conservation of total energy
   logical    fill              ! Flag to use fill algorithm
   logical    BUDGETS           ! Logical to Calculate Exact Budet Terms

   integer    nx          ! # of split pieces in x-direction; for performance, the
   parameter (nx = 1)     ! user may set nx=1 if there is NO shared memory multitasking
   integer ipe, it
   integer n
   integer incount, outcount
   integer iqa, iqb, iqc, iqd, mq  ! used for tracer transpose grouping
   integer im1, ip1, js2gd, jn2gd, jpole, jstar, msgn

! Geometric arrays
! column integrated Total Energy
   real(r8) :: tte(state%grid%jm)


! The control variables (with YZ decomposition)

   real(r8), allocatable  :: u(:,:,:)    ! u wind velocities, staggered grid
   real(r8), allocatable  :: v(:,:,:)    ! v wind velocities, staggered grid
   real(r8), allocatable  :: pt(:,:,:)   ! scaled (virtual) potential temperature
   real(r8), pointer      :: pe(:,:,:)   ! Pressure at layer edges 
   type(T_TRACERS), allocatable        :: q_internal(:)   ! WAS: q3

! Other arrays
   real(r8), pointer     :: phis (:,:)   ! Surface geopotential
   real(r8), pointer     :: delp (:,:,:) ! Pressure thickness  (must be calculated)
   real(r8), pointer     :: delp0(:,:,:) ! Pressure thickness  (must be calculated)
   real(r8), pointer     :: delpd(:,:,:) ! Pressure thickness  (must be calculated)

!----------------------------------------------------------------------------
! The three arrays PE, PK, PKZ must be pre-computed as input to benergy(). 
! They are NOT needed if consv=.F.; updated on output (to be used by physdrv)
! Please refer to routine pkez on the algorithm for computing pkz
! from pe and pk
!----------------------------------------------------------------------------

   real(r8), allocatable :: cosp (:,:)  ! Cosine of Area Grid-Box
   real(r8), allocatable ::  ua(:,:,:)  ! u-wind on the A-Grid
   real(r8), allocatable ::  va(:,:,:)  ! v-wind on the A-Grid
   real(r8), allocatable ::  pk(:,:,:)  ! pe to the kappa
   real(r8), allocatable :: pke(:,:,:)  ! pe to the kappa
   real(r8), allocatable :: dum(:,:,:)  ! layer-mean pk for converting t to pt

   real(r8), allocatable ::      ptfx(:,:,:),      ptfy(:,:,:)
   real(r8), allocatable ::ptfx_accum(:,:,:),ptfy_accum(:,:,:)
   real(r8), allocatable ::       mfx(:,:,:),       mfy(:,:,:)
   real(r8), allocatable :: mfx_accum(:,:,:), mfy_accum(:,:,:)
   real(r8), allocatable ::        cx(:,:,:),        cy(:,:,:)
   real(r8), allocatable ::  cx_accum(:,:,:),  cy_accum(:,:,:)
   real(r8), allocatable ::       tv (:,:,:)

   real(r8), allocatable ::     dthconsv   (:,:,:)
   real(r8), allocatable ::     dthdtremap1(:,:)
   real(r8), allocatable ::     dthdtremap2(:,:)
   real(r8), allocatable ::     dkedt      (:,:)
   real(r8), allocatable ::     dkedt_xy   (:,:,:)
   real(r8), allocatable ::     dkedt_pg_yz(:,:,:)
   real(r8), allocatable ::     dkedt_dp_yz(:,:,:)
   real(r8), allocatable ::     dkedt_ad_yz(:,:,:)
   real(r8), allocatable ::     dkedt_ho_yz(:,:,:)

   real(r8), allocatable ::   thvx(:,:,:)
   real(r8), allocatable ::   thvy(:,:,:)
   real(r8), allocatable ::   difx(:,:,:)
   real(r8), allocatable ::   dify(:,:,:)
   real(r8), allocatable ::   pbrx(:,:,:)
   real(r8), allocatable ::   pbry(:,:,:)
   real(r8), allocatable ::   facx(:,:,:)
   real(r8), allocatable ::   facy(:,:,:)
   real(r8), allocatable ::   dpdx(:,:,:)
   real(r8), allocatable ::   dpdy(:,:,:)
   real(r8), allocatable :: oma_yz(:,:,:)
   real(r8), allocatable :: cpt_yz(:,:,:)
   real(r8), allocatable :: thv_yz(:,:,:)
   real(r8), allocatable ::   pkz (:,:,:)
   real(r8), allocatable ::   pkz0(:,:,:)
   real(r8), allocatable ::   pkzb(:,:,:)
   real(r8), allocatable ::  thdp (:,:,:)
   real(r8), allocatable ::  thdp0(:,:,:)
   real(r8), allocatable ::  thdpb(:,:,:)

   real(r8), allocatable ::       gze(:,:,:)
   real(r8), allocatable ::       phi(:,:,:)
   real(r8), allocatable ::   del_KE (:,:,:), tot_KE (:,:,:)
   real(r8), allocatable ::   del_PHI(:,:,:), tot_PHI(:,:,:)

   real(r8), allocatable :: worka(:,:,:),dp0(:,:,:)
   real(r8), allocatable :: delpf(:,:,:)

   real(r8), allocatable :: u0   (:,:,:), v0   (:,:,:)
   real(r8), allocatable :: ubar (:,:,:), vbar (:,:,:)
   real(r8), allocatable :: dupg (:,:,:), dvpg (:,:,:)
   real(r8), allocatable :: duad (:,:,:), dvad (:,:,:)

   real(r8), allocatable :: ud_yz(:,:,:), vd_yz(:,:,:)
   real(r8), allocatable :: ua_yz(:,:,:), va_yz(:,:,:)
   real(r8), allocatable :: uc_yz(:,:,:), vc_yz(:,:,:)

   real(r8), allocatable :: dwz(:,:,:), pkc(:,:,:), wz(:,:,:)
   real(r8), allocatable :: dpt(:,:,:)
   real(r8), allocatable :: pkcc(:,:,:), wzc(:,:,:)
   real(r8), allocatable :: tmp2(:,:,:)

! The following variables are work arrays for xy=>yz transpose
   real(r8), allocatable :: pkkp(:,:,:), wzkp(:,:,:)

! The following variables are xy instanciations
   real(r8), allocatable :: tmpxy(:,:,:), dp0xy(:,:,:), wzxy(:,:,:)
   real(r8), allocatable :: tmp3dfor2d(:,:,:)
   real(r8), pointer     :: pexy(:,:,:), psxy(:,:)
   real(r8), pointer     :: delpxy (:,:,:)
   real(r8), pointer     :: delpxy0(:,:,:)

   real(r8) dt
   real(r8) bdt
   real(r8) te0,te1

   GRID => STATE%GRID     ! For convenience
   VARS => STATE%VARS     ! For convenience

! Initialize Budget Diagnostics
! -----------------------------
   convcpt    = 0.0
   convthv    = 0.0
   kenrg1     = 0.0
   penrg1     = 0.0
   tenrg1     = 0.0
   kenrg2     = 0.0
   penrg2     = 0.0
   tenrg2     = 0.0
   dkedtad    = 0.0
   dkedtpg    = 0.0
   dkedtdp    = 0.0
   dkedtho    = 0.0
   dthdtremap = 0.0
   dthdtconsv = 0.0

!
! Set the local variables which used to be imported from dynamics_vars
!
   nsplit   = STATE%NSPLIT
   icd      = STATE%icd
   jcd      = STATE%jcd
   iord     = STATE%iord
   jord     = STATE%jord
   kord     = STATE%kord
   te_method= STATE%te_method
   im       = GRID%im
   jm       = GRID%jm
   km       = GRID%km
   nq       = GRID%nq
   ntotq    = GRID%nq
   GRID%ntotq = ntotq
   ng_s     = GRID%ng_s
   ng_c     = GRID%ng_c
   ng_d     = GRID%ng_d
   jfirst   = GRID%jfirst
   kfirst   = GRID%kfirst
   jlast    = GRID%jlast
   klast    = GRID%klast
   ifirstxy = GRID%ifirstxy
   jfirstxy = GRID%jfirstxy
   ilastxy  = GRID%ilastxy
   jlastxy  = GRID%jlastxy

   ptop     = GRID%ptop
   rcap     = GRID%rcap

   twod_decomp = GRID%twod_decomp
   myid_z      = GRID%myid_z
   npr_y       = GRID%npr_y
   npr_z       = GRID%npr_z

   js2gd = max(   2,jfirst-ng_d)   ! NG latitudes on S (starting at 2   )
   jn2gd = min(jm-1,jlast +ng_d)   ! NG latitudes on N (ending   at jm-1)

!
! Allocate the YZ decomposed control variables 
! (adapted from prognostics.F90)
!
   klastp = klast+1
   ndt    = nint( state%dt )  !  IS THIS CORRECT???

   allocate( delpxy0(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
   allocate( delpxy (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
   allocate( pexy(ifirstxy:ilastxy,km+1,jfirstxy:jlastxy) )
   allocate( psxy(ifirstxy:ilastxy,jfirstxy:jlastxy) )

   allocate(   ua(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
   allocate(   va(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

   allocate(   u(im,jfirst-ng_d:jlast+ng_s,kfirst:klast  ) )
   allocate(   v(im,jfirst-ng_s:jlast+ng_d,kfirst:klast  ) )
   allocate(  pt(im,jfirst-ng_d:jlast+ng_d,kfirst:klast  ) )
   allocate( dum(im,jfirst     :jlast     ,kfirst:klast  ) )
   allocate(  pk(im,jfirst     :jlast     ,kfirst:klast+1) )
   allocate( pke(im,kfirst:klast+1,jfirst:jlast) )

   allocate(  cosp(im,jfirst:jlast ) )
   allocate(oma_yz(im,jfirst:jlast,kfirst:klast ) )
   allocate(cpt_yz(im,jfirst:jlast,kfirst:klast ) )
   allocate(thv_yz(im,jfirst:jlast,kfirst:klast ) )

   allocate(  thvx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  thvy(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  difx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  dify(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  pbrx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  pbry(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  facx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  facy(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  dpdx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  dpdy(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  pkz (im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  pkz0(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate(  pkzb(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate( thdp (im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate( thdp0(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
   allocate( thdpb(im,jfirst-ng_d:jlast+ng_d,kfirst:klast ) )
!
! Allocation of tracers
!
   allocate( q_internal(nq) )    ! WAS: q3

!
! Other arrays
!
   allocate( delp (im, jfirst:jlast, kfirst:klast) )
   allocate( delp0(im, jfirst:jlast, kfirst:klast) )
   allocate( delpd(im, jfirst:jlast, kfirst:klast) )

!
!  Determine surface pressure (PS) and pressure difference (DELP)
!
!$omp parallel do private(i,j)
   do j=jfirstxy,jlastxy
     do i=ifirstxy,ilastxy
       psxy(i,j) = vars%pe(i,j,km+1)
     enddo
   enddo

!
! Define local arrays derived from the edge pressure
!    

!$omp parallel do private(i,j,k)
   do k=1,km
     do j=jfirstxy,jlastxy
       do i=ifirstxy,ilastxy
         delpxy (i,j,k) = vars%pe(i,j,k+1) - vars%pe(i,j,k)
         delpxy0(i,j,k) = delpxy (i,j,k)
         pexy   (i,k,j) = vars%pe(i,j,k) 
         pelnxz (i,k,j) = log(vars%pe(i,j,k))
       enddo
     enddo
   enddo

! Clean up loop for level km+1 

!$omp parallel do private(i,j)
   do j=jfirstxy,jlastxy
     do i=ifirstxy,ilastxy
       pexy  (i,km+1,j) = vars%pe(i,j,km+1) 
       pelnxz(i,km+1,j) = log(vars%pe(i,j,km+1))
     enddo
   enddo

   if (twod_decomp == 1) then   ! true 2D decomposition --> transpose

     allocate( phis(im, jfirst:jlast ) )
     allocate( pe(im,kfirst:klastp,jfirst:jlast) )

!
! Transpose phisxy to phis (special 2D case)
!
     allocate( tmp3dfor2d(ifirstxy:ilastxy,jfirstxy:jlastxy,npr_z) )

     do k=1,npr_z
       do j=jfirstxy,jlastxy
         do i=ifirstxy,ilastxy
           tmp3dfor2d(i,j,k) = phisxy(i,j)
         enddo
       enddo
     enddo
     call mp_swapirr( commglobal, grid%xy2d_to_yz2d%SendDesc, &
                      grid%xy2d_to_yz2d%RecvDesc, tmp3dfor2d, phis )

     do k=1,npr_z
         do j=MAX(2,jfirstxy),MIN(jlastxy,jm-1)
            tmp3dfor2d(:,j,k) = grid%cosp(j)
         enddo
         if ( jfirstxy == 1  ) tmp3dfor2d(:,jfirstxy,k) = 0.0
         if ( jlastxy  == jm ) tmp3dfor2d(:,jlastxy, k) = 0.0
     enddo
     call mp_swapirr( commglobal, grid%xy2d_to_yz2d%SendDesc, &
                      grid%xy2d_to_yz2d%RecvDesc, tmp3dfor2d, cosp )

     deallocate(tmp3dfor2d)

     call mp_swapirr( commglobal, grid%ijk_xy_to_yz%SendDesc, &
                      grid%ijk_xy_to_yz%RecvDesc, delpxy, delp )

!
! State control variables first have to be transposed from XY to YZ
!
     call mp_swapirr( commglobal, grid%uxy_to_u%SendDesc,     &
                      grid%uxy_to_u%RecvDesc, vars%u, u )

     call mp_swapirr( commglobal, grid%vxy_to_v%SendDesc,     &
                      grid%vxy_to_v%RecvDesc, vars%v, v )

     call mp_swapirr( commglobal, grid%ptxy_to_pt%SendDesc,   &
                      grid%ptxy_to_pt%RecvDesc, vars%pt, pt,  &
                      a2in = vars%pkz, a2out = pkz )

     call mp_swapirr( commglobal, grid%pexy_to_pe%SendDesc,   &
                      grid%pexy_to_pe%RecvDesc, pexy, pe )

!
! Allocate internal tracers
!
     do mq = 1, nq
       q_internal(mq)%is_r4 = vars%tracer(mq)%is_r4
       if (  vars%tracer(mq)%is_r4 ) then
         allocate( q_internal(mq)%content_r4(im,jfirst:jlast,kfirst:klast) )
! mp_swapirr_r4 not yet available
         call mp_swapirr_r4( commglobal, grid%r4_xy_to_yz%SendDesc, &
                             grid%r4_xy_to_yz%RecvDesc,             &
                             vars%tracer(mq)%content_r4,            &
                             q_internal(mq)%content_r4 )
       else
         allocate( q_internal(mq)%content(im,jfirst:jlast,kfirst:klast) )
         call mp_swapirr( commglobal, grid%ijk_xy_to_yz%SendDesc,   &
                          grid%ijk_xy_to_yz%RecvDesc,               &
                          vars%tracer(mq)%content,                  &
                          q_internal(mq)%content )
       endif
     enddo

   else    ! 1D decomposition --> arrays can be copied

!
! In this case:   kfirst=1, klast=km,  ifirstxy=1 ilastxy=im
!                 jfirst=jfirstxy, jlast=jlastxy

!$omp parallel do private(i,j,k)
     do k=1,km
       do j=jfirst,jlast
         do i=1,im
              u(i,j,k) = vars%u  (i,j,k)
              v(i,j,k) = vars%v  (i,j,k)
             pt(i,j,k) = vars%pt (i,j,k)
            pkz(i,j,k) = vars%pkz(i,j,k)
           delp(i,j,k) =   delpxy(i,j,k)
         enddo
       enddo
     enddo

!$omp parallel do private(i,j,k,mq)
     do mq = 1, nq
       q_internal(mq)%is_r4 = vars%tracer(mq)%is_r4
       if ( vars%tracer(mq)%is_r4 ) then
         q_internal(mq)%content_r4 => vars%tracer(mq)%content_r4
       else
         q_internal(mq)%content => vars%tracer(mq)%content
       endif
     enddo

! PHIS, PE can point directly to their XY equivalents in this case
     phis => phisxy
     pe   => pexy

     do j=MAX(2,jfirstxy),MIN(jlastxy,jm-1)
        cosp(:,j) = grid%cosp(j)
     enddo
     if ( jfirstxy == 1  ) cosp(:,jfirstxy) = 0.0
     if ( jlastxy  == jm ) cosp(:,jlastxy ) = 0.0

   endif


! Allocate temporary work arrays
! Change later to use pointers for SMP performance???
! (prime candidates: uc, vc, delpf)

      allocate( worka(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   dp0(im,jfirst:     jlast,     kfirst:klast) )
      allocate(  ptfx(im,jfirst:     jlast,     kfirst:klast) )
      allocate(  ptfy(im,jfirst:     jlast+1,   kfirst:klast) )
      allocate(ptfx_accum(im,jfirst: jlast,     kfirst:klast) )
      allocate(ptfy_accum(im,jfirst: jlast+1,   kfirst:klast) )
      allocate(   mfx(im,jfirst:     jlast,     kfirst:klast) )
      allocate(   mfy(im,jfirst:     jlast+1,   kfirst:klast) )
      allocate( mfx_accum(im,jfirst: jlast,     kfirst:klast) )
      allocate( mfy_accum(im,jfirst: jlast+1,   kfirst:klast) )
      allocate(    cx(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(    cy(im,jfirst:     jlast+1,   kfirst:klast) )
      allocate( cx_accum(im,jfirst:  jlast     ,kfirst:klast) )
      allocate( cy_accum(im,jfirst:  jlast     ,kfirst:klast) )
      allocate( delpf(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(   dpt(im,jfirst-1:   jlast+1,   kfirst:klast) )
      allocate(   dwz(im,jfirst-1:    jlast,    kfirst:klast+1) )
      allocate(   pkc(im,jfirst-1:   jlast+1,   kfirst:klast+1) ) 
      allocate(    wz(im,jfirst-1:   jlast+1,   kfirst:klast+1) )
      allocate(  pkcc(im,jfirst  :   jlast  ,   kfirst:klast+1) ) 
      allocate(   wzc(im,jfirst  :   jlast  ,   kfirst:klast+1) ) 
      allocate(  pkkp(im,jfirst:jlast,kfirst:klast+1))
      allocate(  wzkp(im,jfirst:jlast,kfirst:klast+1))

      allocate(  wzxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      allocate( tmpxy(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )
      allocate( dp0xy(ifirstxy:ilastxy,jfirstxy:jlastxy,km  ) )

      allocate( gze      (ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      allocate( phi      (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( del_PHI  (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( del_KE   (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( tot_PHI  (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( tot_KE   (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( tv       (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

      allocate(  u0   (im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate(  v0   (im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
      allocate(  ubar (im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate(  vbar (im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
      allocate(  dupg (im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate(  dvpg (im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
      allocate(  duad (im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate(  dvad (im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )

      allocate(  ud_yz(im,jfirst-ng_d:jlast+ng_s,kfirst:klast) )
      allocate(  vd_yz(im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
      allocate(  ua_yz(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(  va_yz(im,jfirst-ng_s:jlast+ng_d,kfirst:klast) )
      allocate(  uc_yz(im,jfirst-ng_d:jlast+ng_d,kfirst:klast) )
      allocate(  vc_yz(im,jfirst-2:   jlast+2,   kfirst:klast) )

      allocate( dkedt_ho_yz(im,jfirst:jlast,kfirst:klast)         )
      allocate( dkedt_dp_yz(im,jfirst:jlast,kfirst:klast)         )
      allocate( dkedt_pg_yz(im,jfirst:jlast,kfirst:klast)         )
      allocate( dkedt_ad_yz(im,jfirst:jlast,kfirst:klast)         )
      allocate( dkedt_xy   (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( dkedt      (ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      allocate( dthdtremap1(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      allocate( dthdtremap2(ifirstxy:ilastxy,jfirstxy:jlastxy)    )
      allocate( dthconsv   (ifirstxy:ilastxy,jfirstxy:jlastxy,km) )

! Initialize Budget Tracer Variables
! ----------------------------------
            BUDGETS = ( NPHI.ne.-999 ) .and. ( NKE.ne.-999 )
       if ( BUDGETS ) then
            if ( vars%tracer(NPHI)%is_r4 ) then
                 tot_PHI  = vars%tracer(NPHI)%content_r4 * delpxy
            else
                 tot_PHI  = vars%tracer(NPHI)%content    * delpxy
            endif
                 del_PHI = tot_PHI

            if ( vars%tracer(NKE )%is_r4 ) then
                 tot_KE   = vars%tracer(NKE )%content_r4 * delpxy
            else
                 tot_KE   = vars%tracer(NKE )%content    * delpxy
            endif
                 del_KE  = tot_KE 
       endif

       do k=kfirst,klast
       do j=jfirst,jlast
          pkz0(:,j,k) =  pkz(:,j,k)
         thdp0(:,j,k) =   pt(:,j,k)*delp(:,j,k)
         delp0(:,j,k) = delp(:,j,k)
            u0(:,j,k) =    u(:,j,k)
            v0(:,j,k) =    v(:,j,k)
       enddo
       enddo

! ---------------------------------------------------

! First touch pkc and wz??? (bufferpack is multitask in vertical but geop
! computations are parallel in j-loop)

   if ( km > 1 ) then         ! not shallow water equations

      if( consv ) then
! Compute globally integrated Total Energy (te0)

! WS: 2006.05.05   This section now in synch with last GFDL code
!                  NOTE: new benergy takes temperature as input
         call FVstartclock(grid,"--BENERGY")
         qqq2 = qqq
         call benergy(grid, vars%u,  vars%v, txy,     delpxy,           &
                      qqq,  vars%pe, pelnxz, phisxy,                    &
                      eps,  cp,      rgas,   tte,     te0 )
         call FVstopclock(grid,"--BENERGY")
      endif

   endif

! Construct 1 Single Outer Loop for Dynamics and Advection
! --------------------------------------------------------
      bdt = ndt
       dt = bdt / float(nsplit)

   if( nq > 0 ) then

!$omp parallel do private(i, j, k)
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=1,im
              ptfx_accum(i,j,k) = 0.
               mfx_accum(i,j,k) = 0.
                cx_accum(i,j,k) = 0.
                cy_accum(i,j,k) = 0.
            enddo
         enddo
         do j=jfirst,jlast+1
            do i=1,im
              ptfy_accum(i,j,k) = 0.
               mfy_accum(i,j,k) = 0.
            enddo
         enddo
      enddo

   endif

  oma_yz  = 0.0
  cpt_yz  = 0.0
  thv_yz  = 0.0

  do 2000 n=1, nsplit

   if( nq > 0 ) then

!$omp parallel do private(i, j, k)
      do k=kfirst,klast
         do j=jfirst,jlast
            do i=1,im
! Save initial delp field before the small-time-step
! Initialize the CFL number accumulators: (cx, cy)
! Initialize total mass fluxes: (mfx, mfy)
               dp0(i,j,k) = delp(i,j,k)
                cx(i,j,k) = 0.
                cy(i,j,k) = 0.
               mfx(i,j,k) = 0.
            enddo
         enddo
         do j=jfirst,jlast+1
            do i=1,im
               mfy(i,j,k) = 0.
            enddo
         enddo
      enddo

   endif

      ptfx = 0.
      ptfy = 0.

! Force Updates of pexy every timestep
! ------------------------------------
      ipe = 1

! Call the Lagrangian dynamical core using small time step
! --------------------------------------------------------
      call FVstartclock(grid,"--CDCORE")

      call cd_core(grid,    nx,      u,       v,      pt,            &
                   delp,    pe,      pk,      ipe,    dt,            &
                   ptop,    umax,    pi,      radius, cp,            &
                   kappa,   icd,     jcd,     iord,   jord,          &
                   ipe,     omega,   phis,    cx,     cy,            &
                   mfx,     mfy,     ptfx,    ptfy,                  &
                   delpf,   uc_yz,   vc_yz,                          &
                   ubar,    vbar,    duad,    dvad,   dupg,  dvpg,   &
                   dum,     dpt,     worka,   dwz,    pkc,           &
                   wz,      phisxy,  vars%pt, pkxy,   pexy,          &
                   pkcc,    wzc,     wzxy,    delpxy, pkkp, wzkp     )

     call FVstopclock(grid,"--CDCORE")

! ---------------------------------------------------------
!                    Compute Omega*Alpha
! ---------------------------------------------------------

! Remove TimeStep from TH*DelP fluxes
! -----------------------------------
     call FVstartclock(grid,"--OMEGA")
     ptfx = ptfx/dt
     ptfy = ptfy/dt

! Load TH and PE**Kappa (A-Grid)
! ------------------------------
     if (twod_decomp == 1) then
         call mp_swapirr( commglobal, grid%ptxy_to_pt%SendDesc, &
                          grid%ptxy_to_pt%RecvDesc, vars%pt, pt  )
         call mp_swapirr( commglobal, grid%pexy_to_pe%SendDesc, &
                          grid%pexy_to_pe%RecvDesc, pexy, pe  )
     else
         pe => pexy
         do k=kfirst,klast
         do j=jfirst,jlast
         do i=1,im
            pt(i,j,k) = vars%pt(i,j,k)
         enddo
         enddo
         enddo
     endif
         pke = pe**kappa

! Compute TH*DelP and P**Kappa (A-Grid)
! -------------------------------------
     do k=kfirst,klast
     do j=jfirst,jlast
     do i=1,im
        pkz(i,j,k) = ( pke(i,k+1,j)-pke(i,k,j) ) / ( kappa * log(pe(i,k+1,j)/pe(i,k,j)) )
       thdp(i,j,k) = pt(i,j,k)*delp(i,j,k)
     enddo
     enddo
     enddo

! Compute Time-Averaged TH*DelP and P**Kappa (A-Grid)
! ---------------------------------------------------
     do k=kfirst,klast
     do j=jfirst,jlast
     do i=1,im
      pkzb(i,j,k) = 0.5_r8*(  pkz0(i,j,k) + pkz (i,j,k) )
     thdpb(i,j,k) = 0.5_r8*( thdp0(i,j,k) + thdp(i,j,k) )
     enddo
     enddo
     enddo

! Compute PTFX*D(PK)/Dx and PTFY*D(PK)/Dy at C-Grid Locations
! -----------------------------------------------------------
     call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                        kfirst, klast, ng_d, ng_d, pkzb )
     call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                        kfirst, klast, ng_d, ng_d, pkzb )
     do k=kfirst,klast
     do j=jfirst,jlast
        im1=im
     do i=1,im
        dpdx(i,j,k) = ptfx(i,j,k)*( pkzb(i,j,k)-pkzb(im1,j,k) )
        pbrx(i,j,k) = ptfx(i,j,k)*( pkzb(i,j,k)+pkzb(im1,j,k) )*0.5_r8
        im1=i
     enddo
     enddo
     do j=max(2,jfirst),jlast
     do i=1,im
        dpdy(i,j,k) = ptfy(i,j,k)*( pkzb(i,j,k)-pkzb(i,j-1,k) )
        pbry(i,j,k) = ptfy(i,j,k)*( pkzb(i,j,k)+pkzb(i,j-1,k) )*0.5_r8
     enddo
     enddo
     enddo

! Average Above C-Grid Quantities back to A-Grid Locations (Away from Poles)
! --------------------------------------------------------------------------
     call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                        kfirst, klast, ng_d, ng_d, dpdy )
     call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                        kfirst, klast, ng_d, ng_d, dpdy )
     call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                        kfirst, klast, ng_d, ng_d, pbry )
     call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                        kfirst, klast, ng_d, ng_d, pbry )

     do k=kfirst,klast
     do j=jfirst,jlast
        i=im
     do ip1=1,im
        facx(i,j,k) = 0.5_r8*( dpdx(ip1,j,k)+dpdx(i,j,k) )
        i=ip1
     enddo
     enddo
     do j=max(2,jfirst),min(jm-1,jlast)
     do i=1,im
        facy(i,j,k) = 0.5_r8*( dpdy(i,j+1,k)+dpdy(i,j,k) )
     enddo
     enddo
     enddo

! Compute Omega*Alpha*DelP at A-Grid Locations (Away from Poles)
! --------------------------------------------------------------
     if ( jfirst ==  1 ) then
          jpole  =   1
          jstar  =   2
          msgn   =   1
     endif
     if ( jlast  == jm ) then
          jpole  =  jm
          jstar  =  jm
          msgn   =  -1
     endif

     do k=kfirst,klast
        do j=max(2,jfirst),min(jm-1,jlast)
        do i=1,im
        oma_yz(i,j,k) = oma_yz(i,j,k) + cp*( thdpb(i,j,k)*( pkz(i,j,k)-pkz0(i,j,k) )/dt  +  facx(i,j,k) + facy(i,j,k)/cosp(i,j) )
        enddo
        enddo

! Poles
! -----
        if ( jfirst == 1 .or. jlast == jm ) then
          sump = 0.0_r8
          do i=1,im
          sump = sump + thdpb(i,jpole,k)*( pkz(i,jpole,k)-pkz0(i,jpole,k) )
          enddo
          sump = sump / (dt*im)
          do i=1,im
          oma_yz(i,jpole,k) = oma_yz(i,jpole,k) + cp*sump
          enddo
          sump = 0.0_r8
          do i=1,im
          sump = sump + pbry(i,jstar,k) - pkzb(i,jpole,k)*ptfy(i,jstar,k)
          enddo
          sump = msgn*sump*rcap
          do i=1,im
          oma_yz(i,jpole,k) = oma_yz(i,jpole,k) + cp*sump
          enddo
        endif

     enddo  ! End K-Loop
     call FVstopclock(grid,"--OMEGA")

! ---------------------------------------------------------
!               Compute Energy Budgets
! ---------------------------------------------------------

     call FVstartclock(grid,"--BUDGETS")
     if( BUDGETS ) then

     do k=kfirst,klast
     do j=jfirst,jlast
        i=im
     do ip1=1,im
        difx(i,j,k) = ( pbrx(ip1,j,k)-pbrx(i,j,k) )
        thvx(i,j,k) = ( ptfx(ip1,j,k)-ptfx(i,j,k) )
        i=ip1
     enddo
     enddo
     do j=max(2,jfirst),min(jm-1,jlast)
     do i=1,im
        dify(i,j,k) = ( pbry(i,j+1,k)-pbry(i,j,k) )
        thvy(i,j,k) = ( ptfy(i,j+1,k)-ptfy(i,j,k) )
     enddo
     enddo
     enddo

! Compute Budget Diagnostics at A-Grid Locations (Away from Poles)
! ----------------------------------------------------------------
     if ( jfirst ==  1 ) then
          jpole  =   1
          jstar  =   2
          msgn   =   1
     endif
     if ( jlast  == jm ) then
          jpole  =  jm
          jstar  =  jm
          msgn   =  -1
     endif

     do k=kfirst,klast
        do j=max(2,jfirst),min(jm-1,jlast)
        do i=1,im
        cpt_yz(i,j,k) = cpt_yz(i,j,k) - cp*( difx(i,j,k) + dify(i,j,k)/cosp(i,j) )
        thv_yz(i,j,k) = thv_yz(i,j,k) -    ( thvx(i,j,k) + thvy(i,j,k)/cosp(i,j) )
        enddo
        enddo

! Poles
! -----
        if ( jfirst == 1 .or. jlast == jm ) then
          sump = 0.0_r8
          do i=1,im
          sump = sump + pbry(i,jstar,k)
          enddo
          sump = msgn*sump*rcap
          do i=1,im
          cpt_yz(i,jpole,k) = cpt_yz(i,jpole,k) - cp*sump
          enddo

          sump = 0.0_r8
          do i=1,im
          sump = sump + ptfy(i,jstar,k)
          enddo
          sump = msgn*sump*rcap
          do i=1,im
          thv_yz(i,jpole,k) = thv_yz(i,jpole,k) - sump
          enddo
        endif

     enddo  ! End K-Loop

! ---------------------------------------------------------
!  Compute Component Energetics Tendencies Across CD_CORE
! ---------------------------------------------------------

     delpd = delp - delp0

! Pressure-Gradient Component
! ---------------------------
      do j=max(2,jfirst),jlast           ; ud_yz(:,j,:) = u0(:,j,:)*dupg(:,j,:) ; enddo
      do j=max(2,jfirst),min(jm-1,jlast) ; vd_yz(:,j,:) = v0(:,j,:)*dvpg(:,j,:) ; enddo

      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      do k=kfirst,klast
         call d2a2c_winds(grid, ud_yz(1,jfirst-ng_d,k), vd_yz(1,jfirst-ng_s,k),        &
                                ua_yz(1,jfirst-ng_d,k), va_yz(1,jfirst-ng_s,k),        &
                                uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k),        &
                       .false., uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k))
      enddo
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jstar,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jpole,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      dkedt_pg_yz(:,jfirst:jlast,:) = 0.5_r8 *( ua_yz(:,jfirst:jlast,:) + va_yz(:,jfirst:jlast,:) )*delp(:,jfirst:jlast,:)

! Advective (Inertial) Component
! ------------------------------
      do j=max(2,jfirst),jlast           ; ud_yz(:,j,:) = u0(:,j,:)*duad(:,j,:) ; enddo
      do j=max(2,jfirst),min(jm-1,jlast) ; vd_yz(:,j,:) = v0(:,j,:)*dvad(:,j,:) ; enddo

      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      do k=kfirst,klast
         call d2a2c_winds(grid, ud_yz(1,jfirst-ng_d,k), vd_yz(1,jfirst-ng_s,k),        &
                                ua_yz(1,jfirst-ng_d,k), va_yz(1,jfirst-ng_s,k),        &
                                uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k),        &
                       .false., uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k))
      enddo
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jstar,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jpole,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      dkedt_ad_yz(:,jfirst:jlast,:) = 0.5_r8 *( ua_yz(:,jfirst:jlast,:) + va_yz(:,jfirst:jlast,:) )*delp(:,jfirst:jlast,:)

! Mean Wind Component
! -------------------
      do j=max(2,jfirst),jlast           ; ud_yz(:,j,:) = 0.5_r8*u0(:,j,:)*u0(:,j,:) ; enddo
      do j=max(2,jfirst),min(jm-1,jlast) ; vd_yz(:,j,:) = 0.5_r8*v0(:,j,:)*v0(:,j,:) ; enddo
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      do k=kfirst,klast
         call d2a2c_winds(grid, ud_yz(1,jfirst-ng_d,k), vd_yz(1,jfirst-ng_s,k),        &
                                ua_yz(1,jfirst-ng_d,k), va_yz(1,jfirst-ng_s,k),        &
                                uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k),        &
                       .false., uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k))
      enddo
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jstar,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jpole,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      dkedt_dp_yz(:,jfirst:jlast,:) = 0.5_r8 *( ua_yz(:,jfirst:jlast,:) + va_yz(:,jfirst:jlast,:) )*delpd(:,jfirst:jlast,:)

! Higher-Order Terms
! ------------------
      do j=max(2,jfirst),jlast           ; ud_yz(:,j,:) = 0.5_r8*( dupg(:,j,:)+duad(:,j,:) )**2 ; enddo
      do j=max(2,jfirst),min(jm-1,jlast) ; vd_yz(:,j,:) = 0.5_r8*( dvpg(:,j,:)+dvad(:,j,:) )**2 ; enddo

      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      do k=kfirst,klast
         call d2a2c_winds(grid, ud_yz(1,jfirst-ng_d,k), vd_yz(1,jfirst-ng_s,k),        &
                                ua_yz(1,jfirst-ng_d,k), va_yz(1,jfirst-ng_s,k),        &
                                uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k),        &
                       .false., uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k))
      enddo
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jstar,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jpole,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      dkedt_ho_yz(:,jfirst:jlast,:) = 0.5_r8 *( ua_yz(:,jfirst:jlast,:) + va_yz(:,jfirst:jlast,:) )*delp(:,jfirst:jlast,:)

! Pole Values
! -----------
      if ( jfirst ==  1 .or. jlast == jm ) then
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
      endif
!              j  = jstar
               j  = jpole
      do k=kfirst,klast
           sumho  = 0.0_r8
           sumdp  = 0.0_r8
           sumpg  = 0.0_r8
           sumad  = 0.0_r8
           do i=1,im
           sumho  = sumho  + dkedt_ho_yz(i,j,k)
           sumdp  = sumdp  + dkedt_dp_yz(i,j,k)
           sumpg  = sumpg  + dkedt_pg_yz(i,j,k)
           sumad  = sumad  + dkedt_ad_yz(i,j,k)
           enddo
           sumho  = sumho /im
           sumdp  = sumdp /im
           sumpg  = sumpg /im
           sumad  = sumad /im
           do i=1,im
           dkedt_ho_yz(i,j,k) = sumho
           dkedt_dp_yz(i,j,k) = sumdp
           dkedt_pg_yz(i,j,k) = sumpg
           dkedt_ad_yz(i,j,k) = sumad
           enddo
      enddo
      endif

      endif  ! End BUDGETS Test
      call FVstopclock(grid,"--BUDGETS")


! Reset Initial Values
! --------------------
     do k=kfirst,klast
     do j=jfirst,jlast
        u0(:,j,k) =    u(:,j,k)
        v0(:,j,k) =    v(:,j,k)
      pkz0(:,j,k) =  pkz(:,j,k)
     thdp0(:,j,k) = thdp(:,j,k)
     delp0(:,j,k) = delp(:,j,k)
     enddo
     enddo

! Accumulate Fluxes
! -----------------
  ptfx_accum(:,:,:) =ptfx_accum(:,:,:) +ptfx(:,:,:)
  ptfy_accum(:,:,:) =ptfy_accum(:,:,:) +ptfy(:,:,:)
   mfx_accum(:,:,:) = mfx_accum(:,:,:) + mfx(:,:,:)
   mfy_accum(:,:,:) = mfy_accum(:,:,:) + mfy(:,:,:)
    cx_accum(:,:,:) =  cx_accum(:,:,:) +  cx(:,jfirst:jlast,:)
    cy_accum(:,:,:) =  cy_accum(:,:,:) +  cy(:,jfirst:jlast,:)


! Perform small-time-step scalar transport using instantaneous CFL and mass fluxes
! --------------------------------------------------------------------------------
   if( nq .ne. 0 ) then
      call FVstartclock(grid,"--TRAC2D")
      call trac2d( grid,   dp0, q_internal,         cx,    cy,       &
                   mfx,          mfy, STATE%iord, STATE%jord,  fill, &
                   dum,        worka  )
      call FVstopclock(grid,"--TRAC2D")
   endif


! Recover BUDGET Tracer Variables into xy-Decomposition
! -----------------------------------------------------
   call FVstartclock(grid,"--BUDGETS")
   if ( BUDGETS ) then
   if ( twod_decomp == 1 ) then
        do mq = 1, nq
           if( mq.eq.NPHI .or. mq.eq.NKE ) then 
               if ( vars%tracer(mq)%is_r4 ) then
!                        mp_swapirr_r4 not yet available
                    call mp_swapirr_r4(commglobal, grid%r4_yz_to_xy%SendDesc, &
                                       grid%r4_yz_to_xy%RecvDesc,             &
				       q_internal(mq)%content_r4,             &
				       vars%tracer(mq)%content_r4  )
               else
                    call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc,  &
                                     grid%ijk_yz_to_xy%RecvDesc,              &
                                     q_internal(mq)%content,                  &
                                     vars%tracer(mq)%content )
               endif
           endif
        enddo
   else
        do mq = 1, nq
           if( mq.eq.NPHI .or. mq.eq.NKE ) then 
               if ( vars%tracer(mq)%is_r4 ) then
                    vars%tracer(mq)%content_r4 => q_internal(mq)%content_r4
               else
                    vars%tracer(mq)%content    => q_internal(mq)%content
               endif
           endif
        enddo
   endif

! Compute and Accumulate Budget Tracer Variable Tendencies
! --------------------------------------------------------
        if (               vars%tracer(NPHI)%is_r4 ) then
             del_PHI  =  ( vars%tracer(NPHI)%content_r4 * delpxy - del_PHI )
        else
             del_PHI  =  ( vars%tracer(NPHI)%content    * delpxy - del_PHI )
        endif
             tot_PHI = tot_PHI + del_PHI

        if (               vars%tracer(NKE )%is_r4 ) then
             del_KE   =  ( vars%tracer(NKE )%content_r4 * delpxy - del_KE  )
        else
             del_KE   =  ( vars%tracer(NKE )%content    * delpxy - del_KE  )
        endif
             tot_KE  = tot_KE  + del_KE 

! Accumulate Kinetic Energy Tendencies across CD_CORE
! ---------------------------------------------------
   if ( twod_decomp == 1 ) then
        call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc,             &
                         grid%ijk_yz_to_xy%RecvDesc, dkedt_ad_yz, dkedt_xy )
   else
        dkedt_xy = dkedt_ad_yz
   endif
   dkedt = 0.0
   do k=1,km
   dkedt = dkedt + dkedt_xy(:,:,k)
   enddo
   dkedt   = dkedt / (grav*dt)
   dkedtad = dkedtad + dkedt


   if ( twod_decomp == 1 ) then
        call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc, &
                         grid%ijk_yz_to_xy%RecvDesc, dkedt_pg_yz, dkedt_xy )
   else
        dkedt_xy = dkedt_pg_yz
   endif
   dkedt = 0.0
   do k=1,km
   dkedt = dkedt + dkedt_xy(:,:,k)
   enddo
   dkedt   = dkedt / (grav*dt)
   dkedtpg = dkedtpg  + dkedt
 

   if ( twod_decomp == 1 ) then
        call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc,            &
                         grid%ijk_yz_to_xy%RecvDesc, dkedt_dp_yz, dkedt_xy )
   else
        dkedt_xy = dkedt_dp_yz
   endif
   dkedt = 0.0
   do k=1,km
   dkedt = dkedt + dkedt_xy(:,:,k)
   enddo
   dkedt   = dkedt / (grav*dt)
   dkedtdp = dkedtdp + dkedt

   if ( twod_decomp == 1 ) then
        call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc,            &
                         grid%ijk_yz_to_xy%RecvDesc, dkedt_ho_yz, dkedt_xy )
   else
        dkedt_xy = dkedt_ho_yz
   endif
   dkedt = 0.0
   do k=1,km
   dkedt = dkedt + dkedt_xy(:,:,k)
   enddo
   dkedt   = dkedt / (grav*dt)
   dkedtho = dkedtho + dkedt

! Update D-Grid Kinetic Energy for HOT
! ------------------------------------

      do j=max(2,jfirst),jlast           ; ud_yz(:,j,:) = u(:,j,:)*u(:,j,:) ; enddo
      do j=max(2,jfirst),min(jm-1,jlast) ; vd_yz(:,j,:) = v(:,j,:)*v(:,j,:) ; enddo
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_send4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_s, ng_d, ud_yz )
      call mp_recv4d_ns( commglobal, im, jm, km, 1, jfirst, jlast, &
                         kfirst, klast, ng_d, ng_s, vd_yz )
      do k=kfirst,klast
         call d2a2c_winds(grid, ud_yz(1,jfirst-ng_d,k), vd_yz(1,jfirst-ng_s,k),        &
                                ua_yz(1,jfirst-ng_d,k), va_yz(1,jfirst-ng_s,k),        &
                                uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k),        &
                       .false., uc_yz(1,jfirst-ng_d,k), vc_yz(1,jfirst-2   ,k))
      enddo
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jstar,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
           ua_yz(:,jpole,:) = 2.0_r8 * ud_yz(:,jpole,:)
           va_yz(:,jpole,:) = 2.0_r8 * vd_yz(:,jstar,:)
      endif
      dkedt_pg_yz(:,jfirst:jlast,:) = 0.25_r8 *( ua_yz(:,jfirst:jlast,:) + va_yz(:,jfirst:jlast,:) )

      if ( jfirst ==  1 .or. jlast == jm ) then
      if ( jfirst ==  1 ) then
           jpole  =   1
           jstar  =   2
      endif
      if ( jlast  == jm ) then
           jpole  =  jm
           jstar  =  jm-1
      endif
      do k=kfirst,klast
           sump = 0.0_r8
           do i=1,im
!          sump = sump + dkedt_pg_yz(i,jstar,k)
           sump = sump + dkedt_pg_yz(i,jpole,k)
           enddo
           sump = sump/im
           do i=1,im
           dkedt_pg_yz(i,jpole,k) = sump
           enddo
      enddo
      endif

      if ( twod_decomp == 1 ) then
           call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc, &
                            grid%ijk_yz_to_xy%RecvDesc, dkedt_pg_yz, dkedt_xy )
      else
           dkedt_xy = dkedt_pg_yz
      endif
 
! Update Latest Values of BUDGET Tracer Variables
! -----------------------------------------------
      gze(:,:,km+1) = phisxy
      do k=km,1,-1
      do j=jfirstxy,jlastxy
      do i=ifirstxy,ilastxy
           gze(i,j,k) = gze(i,j,k+1) + cp*vars%pt(i,j,k)*( pkxy(i,j,k+1)-pkxy(i,j,k) )
            tv(i,j,k) = vars%pt(i,j,k) * ( ( pkxy(i,j,k+1)-pkxy(i,j,k) )/( kappa * log(pexy(i,k+1,j)/pexy(i,k,j)) ) )
           phi(i,j,k) = ( gze(i,j,k+1)*pexy(i,k+1,j)-gze(i,j,k)*pexy(i,k,j) )/delpxy(i,j,k) + (1+kappa)*cp*tv(i,j,k)
      enddo
      enddo
      enddo

! Update Budget Tracer Variables with Latest Values
! -------------------------------------------------
      if ( vars%tracer(NPHI)%is_r4 ) then
           vars%tracer(NPHI)%content_r4 =   phi
      else
           vars%tracer(NPHI)%content    =   phi
      endif
           del_PHI =      phi * delpxy

      if ( vars%tracer(NKE)%is_r4 ) then
           vars%tracer(NKE)%content_r4  = dkedt_xy
      else
           vars%tracer(NKE)%content     = dkedt_xy
      endif
           del_KE  = dkedt_xy * delpxy
 
! Transpose Updated BUDGET Tracer Variables into yz-Decomposition
! ---------------------------------------------------------------
   if ( twod_decomp == 1 ) then
        do mq = 1, nq
           if( mq.eq.NPHI .or. mq.eq.NKE ) then 
               if (  vars%tracer(mq)%is_r4 ) then
!                         mp_swapirr_r4 not yet available
                     call mp_swapirr_r4( commglobal, grid%r4_xy_to_yz%SendDesc, &
                                         grid%r4_xy_to_yz%RecvDesc,  &
                                         vars%tracer(mq)%content_r4, &
                                         q_internal(mq)%content_r4   )
               else
                     call mp_swapirr( commglobal, grid%ijk_xy_to_yz%SendDesc, &
                                      grid%ijk_xy_to_yz%RecvDesc, &
                                      vars%tracer(mq)%content, &
                                      q_internal(mq)%content   )
               endif
           endif
        enddo
   else
        do mq = 1, nq
           if( mq.eq.NPHI .or. mq.eq.NKE ) then 
               if ( vars%tracer(mq)%is_r4 ) then
                    q_internal (mq)%content_r4 => vars%tracer(mq)%content_r4
               else
                    q_internal (mq)%content    => vars%tracer(mq)%content
               endif
           endif
        enddo
   endif
   endif ! End BUDGET Test
   call FVstopclock(grid,"--BUDGETS")

2000  continue

! Compute Average of Accumulated Variables
! ----------------------------------------
    cx_accum = cx_accum / nsplit
    cy_accum = cy_accum / nsplit
    oma_yz   = oma_yz   / nsplit

   call FVstartclock(grid,"--BUDGETS")
   if( BUDGETS ) then
    dkedtad  = dkedtad  / nsplit
    dkedtpg  = dkedtpg  / nsplit
    dkedtdp  = dkedtdp  / nsplit
    dkedtho  = dkedtho  / nsplit
    cpt_yz   = cpt_yz   / nsplit
    thv_yz   = thv_yz   / nsplit
   endif
   call FVstopclock(grid,"--BUDGETS")

! Transpose ps, u, v, and q from yz to xy decomposition
! -----------------------------------------------------
!
! Note: pt, pe and pk will have already been transposed through
! call to geopk in cd_core. geopk does not actually require
! secondary xy decomposition; direct 16-byte technique works just
! as well, perhaps better. However, transpose method is used on last
! call to avoid having to compute these three transposes now.
!

    call FVstartclock(grid,"--TRANSPOSE_FWD")
    if ( twod_decomp == 1 ) then

! Transpose Omega*Alpha
          call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc, &
                           grid%ijk_yz_to_xy%RecvDesc, oma_yz, oma_xy )

! Transpose Convergence of CpT, THv
       if( BUDGETS ) then
          call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc,      &
                           grid%ijk_yz_to_xy%RecvDesc, cpt_yz, convcpt, &
                           a2in=thv_yz, a2out=convthv )
       endif

! Transpose u, v
      call mp_swapirr( commglobal, grid%u_to_uxy%SendDesc,             &
                       grid%u_to_uxy%RecvDesc, u, vars%u )
      call mp_swapirr( commglobal, grid%v_to_vxy%SendDesc,             &
                       grid%v_to_vxy%RecvDesc, v, vars%v )

! Transpose accumulated courant numbers cx, cy
!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               cx_accum(i,j,k) = cx_accum(i,j,k) / grid%dtdx(j)
            enddo
         enddo
      enddo
!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               cy_accum(i,j,k) = cy_accum(i,j,k) / grid%dtdy
            enddo
         enddo
      enddo
      call mp_swapirr(commglobal, grid%ijk_yz_to_xy%SendDesc,       &
                      grid%ijk_yz_to_xy%RecvDesc, cx_accum, cxxyz,  &
                      a2in=cy_accum, a2out=cyxyz )

! Transpose Tracers 
      do mq = 1, nq
        if ( q_internal(mq)%is_r4 ) THEN
!              mp_swapirr_r4 not yet available
          call mp_swapirr_r4(commglobal, grid%r4_yz_to_xy%SendDesc, &
                             grid%r4_yz_to_xy%RecvDesc,             &
                             q_internal(mq)%content_r4,             &
                             vars%tracer(mq)%content_r4 )
          deallocate( q_internal(mq)%content_r4 )
        else
          call mp_swapirr( commglobal, grid%ijk_yz_to_xy%SendDesc,  &
                           grid%ijk_yz_to_xy%RecvDesc,              &
                           q_internal(mq)%content,                  &
                           vars%tracer(mq)%content )
          deallocate( q_internal(mq)%content )
        endif
      enddo

! Horizontal mass fluxes
!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               worka(i,j,k) = mfx_accum(i,j,k) * (grid%dl*radius*grid%cosp(j)) * (radius*grid%dp) / ndt ! Pa m^2 / s
            enddo
         enddo
      enddo
      call mp_swapirr(commglobal, grid%ijk_yz_to_xy%SendDesc,           &
                      grid%ijk_yz_to_xy%RecvDesc, worka, mfxxyz)

!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               worka(i,j,k) = mfy_accum(i,j,k) * (radius*grid%dp) * (grid%dl*radius) / ndt ! Pa m^2 / s  cos(lat) factor is removed
            enddo
         enddo
      enddo
      call mp_swapirr(commglobal, grid%ijk_yz_to_xy%SendDesc,           &
                      grid%ijk_yz_to_xy%RecvDesc, worka, mfyxyz)

! Horizontal PT fluxes
!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               worka(i,j,k) = ptfx_accum(i,j,k) * (grid%dl*radius*grid%cosp(j)) * (radius*grid%dp) / ndt ! K Pa m^2 / s
            enddo
         enddo
      enddo
      call mp_swapirr(commglobal, grid%ijk_yz_to_xy%SendDesc,           &
                      grid%ijk_yz_to_xy%RecvDesc, worka, ptfxxyz)

!$omp parallel do private(i,j,k)
      do k = kfirst,klast
         do j = jfirst,jlast
            do i = 1,im
               worka(i,j,k) = ptfy_accum(i,j,k) * (radius*grid%dp) * (grid%dl*radius) / ndt ! K Pa m^2 / s  cos(lat) factor is removed
            enddo
         enddo
      enddo
      call mp_swapirr(commglobal, grid%ijk_yz_to_xy%SendDesc,           &
                      grid%ijk_yz_to_xy%RecvDesc, worka, ptfyxyz)


    else    ! 1D decomposition --> arrays can be copied

!
! In this case:   kfirst=1, klast=km,  ifirstxy=1 ilastxy=im
!                 jfirst=jfirstxy, jlast=jlastxy

      do k=1,km
        do j=jfirst,jlast
          do i=1,im
            oma_xy(i,j,k)  = oma_yz(i,j,k)
          enddo
        enddo
      enddo

      if( BUDGETS ) then
      do k=1,km
        do j=jfirst,jlast
          do i=1,im
           convcpt(i,j,k)  = cpt_yz(i,j,k)
           convthv(i,j,k)  = thv_yz(i,j,k)
          enddo
        enddo
      enddo
      endif

!$omp parallel do private(i,j,k)
      do k=1,km
        do j=jfirst,jlast
          do i=1,im
            vars%u(i,j,k)  = u(i,j,k)
            vars%v(i,j,k)  = v(i,j,k)
          enddo
        enddo
      enddo

! Horizontal accumulate C-grid winds
!$omp parallel do private(i,j,k)
      do k=1,km
        do j=jfirst,jlast
          do i=1,im
            cxxyz(i,j,k) = cx_accum(i,j,k) / grid%dtdx(j)
          enddo
        enddo
        do j=jfirst,jlast
          do i=1,im
            cyxyz(i,j,k) = cy_accum(i,j,k) / grid%dtdy
          enddo
        enddo
      enddo

!$omp parallel do private(i,j,k,mq)
      do mq = 1, nq
! Copying pointers back is not strictly necessary
        if ( vars%tracer(mq)%is_r4 ) then
          vars%tracer(mq)%content_r4 => q_internal(mq)%content_r4
        else
          vars%tracer(mq)%content    => q_internal(mq)%content
        endif
      enddo

! Horizontal mass and PT fluxes
!$omp parallel do private(i,j,k)
      do k=1,km
        do j=jfirst,jlast
          do i=1,im
            mfxxyz(i,j,k) = mfx_accum(i,j,k) * (grid%dl*radius*grid%cosp(j)) * (radius*grid%dp) / ndt !   Pa m^2 / s
           ptfxxyz(i,j,k) =ptfx_accum(i,j,k) * (grid%dl*radius*grid%cosp(j)) * (radius*grid%dp) / ndt ! K Pa m^2 / s
          enddo
        enddo
        do j=jfirst,jlast
          do i=1,im
            mfyxyz(i,j,k) = mfy_accum(i,j,k) * (radius*grid%dp) * (grid%dl*radius) / ndt !   Pa m^2 / s  cos(lat) factor is removed
           ptfyxyz(i,j,k) =ptfy_accum(i,j,k) * (radius*grid%dp) * (grid%dl*radius) / ndt ! K Pa m^2 / s  cos(lat) factor is removed
          enddo
        enddo
      enddo

    endif   ! twod_decomp == 0/1
    call FVstopclock(grid,"--TRANSPOSE_FWD")

! Recover Total Update for Omega*Alpha and BUDGET Tracers
! -------------------------------------------------------

    oma_xy = oma_xy / delpxy

    call FVstartclock(grid,"--BUDGETS")
    if ( BUDGETS ) then
         if ( vars%tracer(NPHI)%is_r4 ) then
              vars%tracer(NPHI)%content_r4 = tot_PHI / delpxy
         else
              vars%tracer(NPHI)%content    = tot_PHI / delpxy
         endif

         if ( vars%tracer(NKE )%is_r4 ) then
              vars%tracer(NKE )%content_r4 = tot_KE  / delpxy
         else
              vars%tracer(NKE )%content    = tot_KE  / delpxy
         endif
    endif
    call FVstopclock(grid,"--BUDGETS")


! Perform vertical remapping for non-shallow water model
! ------------------------------------------------------
    if ( km > 1 ) then

! Perform vertical remapping from Lagrangian control-volume to
! the Eulerian coordinate as specified by the routine set_eta.
! Note that this finite-volume dycore is otherwise independent of the vertical
! Eulerian coordinate.

!
! te_map requires uxy, vxy, psxy, pexy, pkxy, phisxy, q3xy, and ptxy
!
! Compute Energetics After CD-CORE and Before TE-MAP
! --------------------------------------------------
      call FVstartclock(grid,"--BUDGETS")
      if( BUDGETS ) then
          call Energetics (state,vars%u,vars%v,vars%pt,pexy,phisxy,kenrg1,penrg1,tenrg1)
          dthdtremap1 = 0.0
          do k=1,km
          dthdtremap1 = dthdtremap1 + vars%pt(:,:,k)*(pexy(:,k+1,:)-pexy(:,k,:))
          enddo
      endif
      call FVstopclock(grid,"--BUDGETS")

      mfxxyz_ur = mfxxyz
      mfyxyz_ur = mfyxyz
      call FVstartclock(grid,"--REMAP")
      call te_map(grid, consv,  convt,  psxy,  oma_xy,                             &
                  pexy,   delpxy, vars%pkz,  pkxy,   ndt,                          &
                  nx,     vars%u, vars%v,   vars%pt, vars%tracer,                  &
                  phisxy, cp,     kappa,     kord,   pelnxz,                       &
                  te0,    tmpxy,  dp0xy, mfxxyz, mfyxyz, te_method, dthconsv, dtmp )
      call FVstopclock(grid,"--REMAP")

! Compute Energetics After TE-MAP
! -------------------------------
      call FVstartclock(grid,"--BUDGETS")
      if( BUDGETS ) then
          dthdtremap2 = 0.0
          do k=1,km
          dthdtremap2 = dthdtremap2 + ( vars%pt(:,:,k)-dthconsv(:,:,k) )*(pexy(:,k+1,:)-pexy(:,k,:))
          enddo
          dthdtremap  = (dthdtremap2-dthdtremap1) * P00**KAPPA / (GRAV*NDT)

          dthdtconsv  = 0.0
          do k=1,km
          dthdtconsv  = dthdtconsv + dthconsv(:,:,k)*(pexy(:,k+1,:)-pexy(:,k,:))
          enddo
          dthdtconsv  = dthdtconsv * P00**KAPPA / (GRAV*NDT)
          call Energetics (state,vars%u,vars%v,vars%pt,pexy,phisxy,kenrg2,penrg2,tenrg2)
      endif
      call FVstopclock(grid,"--BUDGETS")

    endif


! Compute Ertel's potential vorticity
! -----------------------------------
   call FVstartclock(grid,"--EPVD")
   call epvd( grid, vars%u, vars%v, vars%pt, delpxy, grav, radius, omega, epvxyz )

! Compute vertical mass flux
! -----------------------------------
   call mfz_comp( grid, radius, grav, ndt, mfxxyz, mfyxyz, mfzxyz, delpxy0, delpxy )
   call FVstopclock(grid,"--EPVD")

!
! Finally, store pexy (IKJ indexing) back into vars%pe (IJK)
!
    do k=1,km+1
      do j=jfirstxy,jlastxy
        do i=ifirstxy,ilastxy
          vars%pe(i,j,k) = pexy(i,k,j)
        enddo
      enddo
    enddo

#if 0
   allocate( tmp2(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
   if( vars%tracer(1)%is_r4 ) then
         qqq2%content_r4 = vars%tracer(1)%content_r4
         tmp2 = vars%pt*vars%pkz/(1.0+eps*qqq2%content_r4)
       else
         qqq2%content    = vars%tracer(1)%content
         tmp2 = vars%pt*vars%pkz/(1.0+eps*qqq2%content)
   endif
   do k=1,km+1
   do j=jfirstxy,jlastxy
   do i=ifirstxy,ilastxy
          vars%pe(i,j,k) = pexy(i,k,j)
          pelnxz (i,k,j) = log(vars%pe(i,j,k))
   enddo
   enddo
   enddo
   call benergy(grid, vars%u,  vars%v, tmp2,    delpxy,           &
                qqq2, vars%pe, pelnxz, phisxy,                    &
                eps,  cp,      rgas,   tte,     te1 )
   if(grid%iam==0) write(6,*) 'Inside FVwrap, TE_beg: ',te0,' TE_end: ',te1,' Diff: ',te0-te1  
   deallocate( tmp2 )
#endif

    if ( twod_decomp == 1 ) then
      deallocate( phis )
      deallocate( pe )
    endif

    deallocate( cosp )
    deallocate( pke  )
    deallocate( pkz  )
    deallocate( pkz0 )
    deallocate( pkzb )
    deallocate( thvx )
    deallocate( thvy )
    deallocate( difx )
    deallocate( dify )
    deallocate( pbrx )
    deallocate( pbry )
    deallocate( facx )
    deallocate( facy )
    deallocate( dpdx )
    deallocate( dpdy )
    deallocate( oma_yz )
    deallocate( cpt_yz )
    deallocate( thv_yz )

    deallocate( thdp  )
    deallocate( thdp0 )
    deallocate( thdpb )

    deallocate( delp0       )
    deallocate( delpd       )
    deallocate( dkedt       )
    deallocate( dkedt_xy    )
    deallocate( dkedt_ho_yz )
    deallocate( dkedt_dp_yz )
    deallocate( dkedt_pg_yz )
    deallocate( dkedt_ad_yz )
    deallocate( dthdtremap1 )
    deallocate( dthdtremap2 )
    deallocate( dthconsv    )

    deallocate(  u0    )
    deallocate(  v0    )
    deallocate(  ubar  )
    deallocate(  vbar  )
    deallocate(  dupg  )
    deallocate(  dvpg  )
    deallocate(  duad  )
    deallocate(  dvad  )

    deallocate( ptfy )
    deallocate( ptfx )
    deallocate( ptfy_accum )
    deallocate( ptfx_accum )
    deallocate( mfy )
    deallocate( mfx )
    deallocate( mfy_accum )
    deallocate( mfx_accum )
    deallocate(  cy )
    deallocate(  cx )
    deallocate(  cy_accum )
    deallocate(  cx_accum )
    deallocate( dp0 )
    deallocate( delpf )
    deallocate( ua_yz )
    deallocate( va_yz )
    deallocate( uc_yz )
    deallocate( vc_yz )
    deallocate( ud_yz )
    deallocate( vd_yz )
    deallocate( dpt   )
    deallocate( pkc   )
    deallocate( dwz   )
    deallocate(  wz   )
    deallocate( worka )
    deallocate( pkcc )
    deallocate( wzc )
    deallocate( pkkp )
    deallocate( wzkp )

! Other arrays
    deallocate( delp )
    deallocate( dum )
    deallocate( pk )

!
! Deallocate the control variables

    deallocate( q_internal )
    deallocate( pt )
    deallocate( v  )
    deallocate( u  )

!
! Deallocate XY decomposition stuff
!
    deallocate( delpxy0 )
    deallocate( delpxy  )
    deallocate( pexy )
    deallocate( psxy )
    deallocate( wzxy )
    deallocate( tmpxy  )
    deallocate( dp0xy  )

    deallocate( gze    )
    deallocate( phi    )
    deallocate( del_KE , tot_KE  )
    deallocate( del_PHI, tot_PHI )
    deallocate( ua )
    deallocate( va )
    deallocate( tv )

    return

contains

  subroutine Energetics (state,ud,vd,thv,plexzy,phis,ke,pe,te)
  use dynamics_vars, only : d2a3d

  type (T_FVDYCORE_STATE) :: STATE
  real(8)   ud(:,:,:)
  real(8)   vd(:,:,:)
  real(8)  thv(:,:,:)
  real(8)   ke(:,:)
  real(8)   pe(:,:)
  real(8)   te(:,:)
  real(8) phis(:,:)

  real(8) plexzy(:,:,:)
  real(8) kinetic, potential, sump
  integer i,ifirst,ilast
  integer j,jfirst,jlast
  integer k,im,jm,km

  real(8), allocatable ::   ud2(:,:,:)
  real(8), allocatable ::   vd2(:,:,:)
  real(8), allocatable ::   ua2(:,:,:)
  real(8), allocatable ::   va2(:,:,:)

  real(8), allocatable ::  delp(:,:,:)
  real(8), allocatable ::    pk(:,:,:)
  real(8), allocatable ::   ple(:,:,:)
  real(8), allocatable ::   pke(:,:,:)
  real(8), allocatable :: gztop(:,:)

  ifirst = lbound( ud,1 )
  ilast  = ubound( ud,1 )
  jfirst = lbound( ud,2 )
  jlast  = ubound( ud,2 )
  km     = ubound( ud,3 )
  im     = state%grid%im
  jm     = state%grid%jm

  allocate( gztop( ifirst:ilast, jfirst:jlast          ) )
  allocate( delp ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( pk   ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( ple  ( ifirst:ilast, jfirst:jlast , 1:km+1 ) )
  allocate( pke  ( ifirst:ilast, jfirst:jlast , 1:km+1 ) )
  allocate( ua2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( va2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( ud2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )
  allocate( vd2  ( ifirst:ilast, jfirst:jlast , 1:km   ) )

  do k=1,km+1
  do j=jfirst,jlast
  do i=ifirst,ilast
     ple(i,j,k) = plexzy(i,k,j)
  enddo
  enddo
  enddo

! Compute Model Top Height
! ------------------------
      pke   = ple**kappa
      gztop = phis
      do k=km,1,-1
        gztop(:,:)   = gztop(:,:) + cp*thv(:,:,k)*( pke(:,:,k+1)-pke(:,:,k) )
         delp(:,:,k) = ple(:,:,k+1)-ple(:,:,k)
           pk(:,:,k) = ( pke(:,:,k+1)-pke(:,:,k) )/( kappa*log(ple(:,:,k+1)/ple(:,:,k)) )
      enddo

! Compute D-Grid Kinetic Energy
! -----------------------------
    ud2 = ud*ud
    vd2 = vd*vd
    call d2a3d( state%grid, ud2, vd2, ua2, va2 )

    if( state%grid%jfirstxy.eq.1 ) then
        ua2(:,1,:) = ud2(:,2,:)
        va2(:,1,:) = vd2(:,2,:)
    endif
    if( state%grid%jlastxy.eq.jm ) then
        ua2(:,jlast,:) = ud2(:,jlast  ,:)
        va2(:,jlast,:) = vd2(:,jlast-1,:)
    endif

! Compute Energetics:  Cp*Tv + K + PHI
! ------------------------------------
       ke = 0.0
       pe = 0.0
  do k=1,km
  do j=jfirst,jlast
  do i=ifirst,ilast
       kinetic   = 0.5_r8*( ua2(i,j,k) + va2(i,j,k) )
       potential = cp*thv(i,j,k)*pk(i,j,k)
       ke(i,j)   = ke(i,j) +   kinetic*delp(i,j,k)
       pe(i,j)   = pe(i,j) + potential*delp(i,j,k)
  enddo
  enddo
  enddo
       ke(:,:) =  ke(:,:)/grav
       pe(:,:) =  pe(:,:)/grav
       te(:,:) = (phis(:,:)*ple(:,:,km+1)-gztop(:,:)*ple(:,:,1))/grav

  if( state%grid%jfirstxy.eq.1 ) then
      call par_xsum ( state%grid, ke(ifirst:ilast,1), 1, sump )
!     call par_xsum ( state%grid, ke(ifirst:ilast,2), 1, sump )
      sump = sump/im
      do i=ifirst,ilast
         ke(i,1) = sump
      enddo
  endif
  if( state%grid%jlastxy.eq.jm ) then
      call par_xsum ( state%grid, ke(ifirst:ilast,jlast  ), 1, sump )
!     call par_xsum ( state%grid, ke(ifirst:ilast,jlast-1), 1, sump )
      sump = sump/im
      do i=ifirst,ilast
         ke(i,jlast) = sump
      enddo
  endif

  deallocate ( gztop )
  deallocate ( delp  )
  deallocate ( pk    )
  deallocate ( ple   )
  deallocate ( pke   )
  deallocate ( ua2   )
  deallocate ( va2   )
  deallocate ( ud2   )
  deallocate ( vd2   )

  return
  end subroutine Energetics

!EOC
end subroutine FVdycore_wrapper
!------------------------------------------------------------------------------
