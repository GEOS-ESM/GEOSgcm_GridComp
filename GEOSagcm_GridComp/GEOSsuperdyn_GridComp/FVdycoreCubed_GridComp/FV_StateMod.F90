#include "MAPL_Generic.h"

module FV_StateMod
!BOP
!
! !MODULE: FV_StateMod --- GEOS5/CAM Cubed-Sphere fvcore state variables/init/run/finalize
!
! !USES:
#if defined( MAPL_MODE )
   use ESMF                ! ESMF base class
   use MAPL_Mod            ! MAPL base class
#endif

   use MAPL_ConstantsMod, only: MAPL_CP, MAPL_RGAS, MAPL_RVAP, MAPL_GRAV, MAPL_RADIUS, &
                                MAPL_KAPPA, MAPL_PI_R8, MAPL_ALHL, MAPL_PSDRY

   use fms_mod, only: fms_init, set_domain, nullify_domain
   use mpp_domains_mod, only: mpp_update_domains, CGRID_NE, DGRID_NE, mpp_get_boundary
   use mpp_parameter_mod, only: AGRID_PARAM=>AGRID, CORNER
   use fv_timing_mod,    only: timing_on, timing_off, timing_init, timing_prt
   use mpp_mod, only: mpp_pe, mpp_root_pe
   use mpp_domains_mod,  only: domain2d

   use fv_grid_utils_mod,  only: inner_prod, mid_pt_sphere, cubed_to_latlon, &
                           ptop_min, g_sum
   use fv_grid_tools_mod,  only: get_unit_vector
   use fv_control_mod, only: fv_init1, fv_init2
   use fv_arrays_mod , only: fv_atmos_type, FVPRC, REAL4, REAL8
   use init_hydro_mod, only: p_var
   use fv_dynamics_mod, only: fv_dynamics
   use fv_update_phys_mod, only: fv_update_phys
   use sw_core_mod, only: d2a2c_vect
   use fv_sg_mod, only: fv_subgrid_z
   use gfdl_cloud_microphys_mod, only: gfdl_cloud_microphys_init

   use fv_diagnostics_mod, only: prt_maxmin, prt_minmax

implicit none
private

#include "mpif.h"

  real(REAL8), save :: elapsed_time = 0
  real(REAL8), save :: massD0
  real(REAL8), save :: massADDED = 0.0

  real(FVPRC), parameter :: tiny_number=1.e-15
  real(FVPRC), parameter :: huge_number=1.e+15

! !PUBLIC DATA MEMBERS:

  logical :: FV_OFF = .false.
  integer :: INT_FV_OFF = 0
  logical :: DEBUG = .false.
  logical :: COLDSTART = .false.
  logical :: SW_DYNAMICS = .false.
  integer :: INT_ADIABATIC = 0
  logical :: ADIABATIC = .false.
  logical :: FV_HYDROSTATIC = .true.
  integer :: INT_check_mass = 0
  logical :: check_mass = .false.
  integer :: INT_fix_mass = 1
  logical :: fix_mass = .true.
  integer :: CASE_ID = 11
  integer :: AdvCore_Advection = 0

  public FV_Atm
  public FV_Setup, FV_InitState, FV_Run, FV_Finalize, FV_DA_Incs
  public FV_HYDROSTATIC, ADIABATIC, DEBUG, COLDSTART, CASE_ID, SW_DYNAMICS, AdvCore_Advection
  public FV_RESET_CONSTANTS
  public FV_To_State, State_To_FV
  public T_TRACERS, T_FVDYCORE_VARS, T_FVDYCORE_GRID, T_FVDYCORE_STATE
  public fv_fillMassFluxes
  public fv_computeMassFluxes
  public fv_getVerticalMassFlux
  public fv_getPK
  public fv_getOmega
  public fv_getVorticity
  public fv_getDivergence
  public fv_getUpdraftHelicity
  public fv_getEPV
  public fv_getDELZ
  public fv_getPKZ
  public fv_getQ

  public debug_fv_state

  public INTERP_DGRID_TO_AGRID
  public INTERP_AGRID_TO_DGRID

  interface fv_computeMassFluxes
     module procedure fv_computeMassFluxes_r4
     module procedure fv_computeMassFluxes_r8
  end interface  

  INTERFACE INTERP_DGRID_TO_AGRID

   MODULE PROCEDURE fv_getAgridWinds_3D
   MODULE PROCEDURE fv_getAgridWinds_2D

  END INTERFACE

  INTERFACE INTERP_AGRID_TO_DGRID

   MODULE PROCEDURE a2d3d    ! 2d to 3d
   MODULE PROCEDURE a2d2d    ! 2d to 2d

  END INTERFACE

  logical, save :: Init_FV_Domain = .true.

  type(fv_atmos_type), allocatable, save :: FV_Atm(:)
  logical, allocatable, save             :: grids_on_this_pe(:)

  type T_TRACERS
       logical                                   :: is_r4
       real(REAL8), dimension(:,:,:  ), pointer     :: content
       real(REAL4), dimension(:,:,:  ), pointer     :: content_r4
       character(LEN=ESMF_MAXSTR)                          :: tname
  end type T_TRACERS

! T_FVDYCORE_VARS contains the prognostic variables for FVdycore
  type T_FVDYCORE_VARS
       real(REAL8), dimension(:,:,:  ), pointer     :: U      => NULL() ! U winds (D-grid)
       real(REAL8), dimension(:,:,:  ), pointer     :: V      => NULL() ! V winds (D-grid)
       real(REAL8), dimension(:,:,:  ), pointer     :: PT     => NULL() ! scaled virtual pot. temp.
       real(REAL8), dimension(:,:,:  ), pointer     :: PE     => NULL() ! Pressure at layer edges
       real(REAL8), dimension(:,:,:  ), pointer     :: PKZ    => NULL() ! P^kappa mean
       real(REAL8), dimension(:,:,:  ), pointer     :: DZ     => NULL() ! Height Thickness
       real(REAL8), dimension(:,:,:  ), pointer     :: W      => NULL() ! Vertical Velocity
       type(T_TRACERS), dimension(:), pointer    :: tracer => NULL() ! Tracers
       integer :: nwat ! Number of water species
  end type T_FVDYCORE_VARS

! T_FVDYCORE_GRID contains information about the horizontal and vertical
! discretization, unlike in ARIES where these data are split into HORZ_GRID
! and VERT_GRID.  The reason for this: currently all of this information is
! initialized in one call to FVCAM dynamics_init.

  type T_FVDYCORE_GRID
!
#if defined( MAPL_MODE )
    type (MAPL_MetaComp),   pointer :: FVgenstate
    type (ESMF_Grid)                :: GRID           ! The 'horizontal' grid (2D decomp only)
#endif

    integer                         :: NG               ! Ghosting
!
    integer                         :: IS               ! Start X-index (exclusive, unghosted)
    integer                         :: IE               ! End X-index (exclusive, unghosted)
    integer                         :: JS               ! Start Y-index (exclusive, unghosted)
    integer                         :: JE               ! End Y-index (exclusive, unghosted)
!
    integer                         :: ISD              ! Start X-index (exclusive, ghosted)
    integer                         :: IED              ! End X-index (exclusive, ghosted)
    integer                         :: JSD              ! Start Y-index (exclusive, ghosted)
    integer                         :: JED              ! End Y-index (exclusive, ghosted)
!
    integer                         :: NPX             ! Full X- dim
    integer                         :: NPY             ! Full Y- dim
    integer                         :: NPZ             ! Numer of levels
    integer                         :: NPZ_P1          ! NPZ+1 (?)

    integer                         :: NTILES          ! How many log-rectangular tiles does my grid Have
                                                       ! lat-lon      = 1
                                                       ! cubed-sphere = 6

    real(REAL8), allocatable           :: DXC(:,:)     ! local C-Gird DeltaX
    real(REAL8), allocatable           :: DYC(:,:)     ! local C-Gird DeltaY 

    real(REAL8), allocatable           :: AREA(:,:)    ! local cell area
    real(REAL8)                        :: GLOBALAREA   ! global area

    integer                         :: KS              ! Number of true pressure levels (out of NPZ+1)
    real(REAL8)                        :: PTOP            ! pressure at top (ak(1))
    real(REAL8)                        :: PINT            ! initial pressure (ak(npz+1))
    real(REAL8), dimension(:), pointer :: AK => NULL()    ! Sigma mapping
    real(REAL8), dimension(:), pointer :: BK => NULL()    ! Sigma mapping
    real(REAL8)                        :: f_coriolis_angle = 0
!
! Tracers
!
    integer                         :: NQ              ! Number of advected tracers
  end type T_FVDYCORE_GRID

  integer, parameter :: NUM_FVDYCORE_ALARMS        = 3
  integer, parameter :: NUM_TIMES      = 8
  integer, parameter :: TIME_TO_RUN  = 1

  type T_FVDYCORE_STATE
!!!    private
    type (T_FVDYCORE_VARS)               :: VARS
    type (T_FVDYCORE_GRID )              :: GRID
#if defined( MAPL_MODE )
    type (ESMF_Clock), pointer           :: CLOCK
    type (ESMF_Alarm)                    :: ALARMS(NUM_FVDYCORE_ALARMS)

#endif
    integer(kind=8)                      :: RUN_TIMES(4,NUM_TIMES)
    logical                              :: DOTIME, DODYN
    real(REAL8)                          :: DT          ! Large time step
    integer                              :: NSPLIT
    integer                              :: NUM_CALLS
  end type T_FVDYCORE_STATE

! Constants used by fvcore
    real(REAL8)                             :: pi
    real(REAL8)                             :: omega    ! angular velocity of earth's rotation
    real(FVPRC)                             :: cp       ! heat capacity of air at constant pressure
    real(REAL8)                             :: radius   ! radius of the earth (m)
    real(REAL8)                             :: rgas     ! Gas constant of the air
    real(REAL8)                             :: rvap     ! Gas constant of vapor
    real(FVPRC)                             :: kappa    ! kappa
    real(REAL8)                             :: grav     ! Gravity
    real(REAL8)                             :: hlv      ! latent heat of evaporation
    real(FVPRC)                             :: zvir     ! RWV/RAIR-1

  real(kind=4), pointer             :: phis(:,:)

  logical :: fv_first_run = .true.

  integer :: ntracers=11

  type (ESMF_Alarm)                    :: MASSALARM
  type (ESMF_TimeInterval)             :: MassAlarmInt
  logical :: first_mass_fix = .true.

!
! !DESCRIPTION:
!
!      This module provides variables which are specific to the Lin-Rood
!      dynamical core.  Most of them were previously SAVE variables in
!      different routines and were set with an "if (first)" statement.
!
!      \begin{tabular}{|l|l|} \hline \hline
!        lr\_init    &  Initialize the Lin-Rood variables  \\ \hline
!        lr\_clean   &  Deallocate all internal data structures \\ \hline
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   2007.07.17   Putman     Created from lat-lon core
!
!EOP
!-----------------------------------------------------------------------
   real(REAL8), parameter ::  D0_0                    =   0.0
   real(REAL8), parameter ::  D0_5                    =   0.5
   real(REAL8), parameter ::  D1_0                    =   1.0
   real(REAL8), parameter ::  D2_0                    =   2.0
   real(REAL8), parameter ::  D4_0                    =   4.0
   real(REAL8), parameter ::  D180_0                  = 180.0
   real(REAL8), parameter ::  ratmax                  =  0.81

contains

!-----------------------------------------------------------------------
!
! !INTERFACE:
 subroutine FV_RESET_CONSTANTS(FV_PI, FV_OMEGA, FV_CP, FV_RADIUS, FV_RGAS, &
                               FV_RVAP, FV_KAPPA, FV_GRAV, FV_HLV, FV_ZVIR)
   real (REAL8), optional, intent(IN) :: FV_PI
   real (REAL8), optional, intent(IN) :: FV_OMEGA
   real (REAL8), optional, intent(IN) :: FV_CP
   real (REAL8), optional, intent(IN) :: FV_RADIUS
   real (REAL8), optional, intent(IN) :: FV_RGAS
   real (REAL8), optional, intent(IN) :: FV_RVAP
   real (REAL8), optional, intent(IN) :: FV_KAPPA
   real (REAL8), optional, intent(IN) :: FV_GRAV
   real (REAL8), optional, intent(IN) :: FV_HLV
   real (REAL8), optional, intent(IN) :: FV_ZVIR

   if (present(FV_PI)) then
      pi = FV_PI
   endif
   if (present(FV_OMEGA)) then
      omega = FV_OMEGA
   endif
   if (present(FV_CP)) then
      cp = FV_CP
   endif
   if (present(FV_RADIUS)) then
      radius = FV_RADIUS
   endif
   if (present(FV_RGAS)) then
      rgas = FV_RGAS
   endif
   if (present(FV_RVAP)) then
      rvap = FV_RVAP
   endif
   if (present(FV_KAPPA)) then
      kappa = FV_KAPPA
   endif
   if (present(FV_GRAV)) then
      grav   = FV_GRAV
   endif
   if (present(FV_HLV)) then
      hlv = FV_HLV
   endif
   if (present(FV_ZVIR)) then
      zvir = FV_ZVIR
   endif

 end subroutine FV_RESET_CONSTANTS
!
!-----------------------------------------------------------------------
 subroutine FV_Setup(GC,LAYOUT_FILE, RC)

  use test_cases_mod, only : test_case

  type (ESMF_GridComp)         , intent(INOUT) :: GC
  character(LEN=*)             , intent(IN   ) :: LAYOUT_FILE
  integer, optional            , intent(OUT  ) :: RC
! Local
   character(len=ESMF_MAXSTR)       :: IAm='FV_StateMod:FV_Setup'
! Local variables

  type (ESMF_Config)           :: cf
  type (ESMF_VM)               :: VM
  integer              :: status
  real(FVPRC) :: DT

  integer   :: ks                 !  True # press. levs
  integer   :: ndt,nx,ny

  type (MAPL_MetaComp),          pointer :: MAPL  => NULL()

  integer :: comm
  integer :: p_split=1

! BEGIN

  call ESMF_VMGetCurrent(VM, rc=STATUS)
  VERIFY_(STATUS)

    call MAPL_MemUtilsWrite(VM, trim(IAm), RC=STATUS )
    VERIFY_(STATUS)

! Retrieve the pointer to the state
! ---------------------------------

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"--FMS_INIT")
    call ESMF_VMGet(VM,mpiCommunicator=comm,rc=status)
    VERIFY_(STATUS)
    call fms_init(comm)
    call MAPL_TimerOff(MAPL,"--FMS_INIT")
    call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FMS_INIT', RC=STATUS )
    VERIFY_(STATUS)
! Start up FV                   
    call MAPL_TimerOn(MAPL,"--FV_INIT")
    call fv_init1(FV_Atm, DT, grids_on_this_pe, p_split)
    call MAPL_TimerOff(MAPL,"--FV_INIT")
    call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FV_INIT', RC=STATUS )
    VERIFY_(STATUS)

  if (FV_Atm(1)%flagstruct%npz == 1) SW_DYNAMICS = .true.


! FV grid dimensions setup from MAPL
      call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%npx, 'AGCM_IM:', default= 32, RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%npy, 'AGCM_JM:', default=192, RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%npz, 'AGCM_LM:', default= 72, RC=STATUS )
      VERIFY_(STATUS)
! FV likes npx;npy in terms of cell vertices
      if (FV_Atm(1)%flagstruct%npy == 6*FV_Atm(1)%flagstruct%npx) then
         FV_Atm(1)%flagstruct%ntiles = 6
         FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npx+1
         FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
      else
         FV_Atm(1)%flagstruct%ntiles = 1
         FV_Atm(1)%flagstruct%npy    = FV_Atm(1)%flagstruct%npy+1
         FV_Atm(1)%flagstruct%npx    = FV_Atm(1)%flagstruct%npx+1
      endif
! Check for Doubly Periodic Domain Info
      call MAPL_GetResource( MAPL, FV_Atm(1)%flagstruct%deglat, label='FIXED_LATS:', default=FV_Atm(1)%flagstruct%deglat, rc=status )
      VERIFY_(STATUS)
! MPI decomp setup
      call MAPL_GetResource( MAPL, nx, 'NX:', default=0, RC=STATUS )
      VERIFY_(STATUS)
      FV_Atm(1)%layout(1) = nx
      call MAPL_GetResource( MAPL, ny, 'NY:', default=0, RC=STATUS )
      VERIFY_(STATUS)
      if (FV_Atm(1)%flagstruct%grid_type == 4) then
         FV_Atm(1)%layout(2) = ny 
      else
         FV_Atm(1)%layout(2) = ny / 6
      end if

! Get other scalars
! -----------------

  call MAPL_GetResource( MAPL, ndt, 'RUN_DT:', default=0, RC=STATUS )
  VERIFY_(STATUS)
  DT = ndt

! Advect tracers within DynCore(AdvCore_Advection=.false.)
!             or within AdvCore(AdvCore_Advection=.true.)
  call MAPL_GetResource( MAPL, AdvCore_Advection, label='AdvCore_Advection:', default=AdvCore_Advection, rc=status )
  VERIFY_(STATUS)

  call MAPL_GetResource( MAPL, INT_fix_mass,    label='fix_mass:'    , default=INT_fix_mass, rc=status )
  call MAPL_GetResource( MAPL, INT_check_mass,  label='check_mass:'  , default=INT_check_mass, rc=status )
  call MAPL_GetResource( MAPL, INT_ADIABATIC,   label='ADIABATIC:'   , default=INT_adiabatic, rc=status )
  call MAPL_GetResource( MAPL, INT_FV_OFF,      label='FV_OFF:'      , default=INT_FV_OFF, rc=status )

  ! MAT The Fortran Standard, and thus gfortran, *does not allow* the use
  !     of if (integer). So, we must convert integer resources to logicals

  if (INT_fix_mass == 0) then
     fix_mass = .FALSE.
  else
     fix_mass = .TRUE.
  end if

  if (INT_check_mass == 0) then
     check_mass = .FALSE.
  else
     check_mass = .TRUE.
  end if

  if (INT_ADIABATIC == 0) then
     ADIABATIC = .FALSE.
  else
     ADIABATIC = .TRUE.
  end if

  if (INT_FV_OFF == 0) then
     FV_OFF = .FALSE.
  else
     FV_OFF = .TRUE.
  end if

! Constants
!
    pi     = MAPL_PI_R8
    omega  = MAPL_OMEGA    ! angular velocity of earth's rotation
    cp     = MAPL_CP       ! heat capacity of air at constant pressure
    radius = MAPL_RADIUS   ! radius of the earth (m)
    rgas   = MAPL_RGAS     ! Gas constant of the air
    rvap   = MAPL_RVAP     ! Gas constant of vapor
    kappa  = MAPL_KAPPA    ! kappa
    grav   = MAPL_GRAV     ! Gravity
    hlv    = MAPL_ALHL     ! latent heat of evaporation
    zvir   = MAPL_RVAP/MAPL_RGAS - 1.   ! RWV/RAIR-1

!
! Create some resolution dependent defaults for FV3 in GEOS...
!    These can be overrided in fv_core_nml in fvcore_layout.rc linked to input.nml
!-------------------------------------------------------------
  ! Number of water species for FV3 determined later
  ! when reading the tracer bundle in fv_first_run 
   FV_Atm(1)%flagstruct%nwat = 0
  ! Veritical resolution dependencies
   if (FV_Atm(1)%flagstruct%npz == 72) then
     FV_Atm(1)%flagstruct%n_sponge = 9 ! ~0.2mb
     FV_Atm(1)%flagstruct%n_zfilter = 25 ! ~10mb
   endif
   if (FV_Atm(1)%flagstruct%npz == 132) then
     FV_Atm(1)%flagstruct%n_sponge = 9 ! ~0.2mb
     FV_Atm(1)%flagstruct%n_zfilter = 30 ! ~10mb
   endif
   FV_Atm(1)%flagstruct%tau = 0.
   FV_Atm(1)%flagstruct%rf_cutoff = 7.5e2
   FV_Atm(1)%flagstruct%d2_bg_k1 = 0.20
   FV_Atm(1)%flagstruct%d2_bg_k2 = 0.06
   FV_Atm(1)%flagstruct%remap_option = 0
   FV_Atm(1)%flagstruct%fv_sg_adj = DT
   FV_Atm(1)%flagstruct%kord_tm =  9
   FV_Atm(1)%flagstruct%kord_mt =  9
   FV_Atm(1)%flagstruct%kord_wz =  9
   FV_Atm(1)%flagstruct%kord_tr =  9
   FV_Atm(1)%flagstruct%z_tracer = .true.
  ! Some default horizontal flags
   FV_Atm(1)%flagstruct%adjust_dry_mass = fix_mass
   FV_Atm(1)%flagstruct%consv_te = 1.
   FV_Atm(1)%flagstruct%consv_am = .false.
   FV_Atm(1)%flagstruct%fill = .true.
   FV_Atm(1)%flagstruct%dwind_2d = .false.
   FV_Atm(1)%flagstruct%delt_max = 0.002
   FV_Atm(1)%flagstruct%ke_bg = 0.0
  ! Some default damping options
   FV_Atm(1)%flagstruct%nord = 2
   FV_Atm(1)%flagstruct%dddmp = 0.2
   FV_Atm(1)%flagstruct%d4_bg = 0.12
   FV_Atm(1)%flagstruct%d2_bg = 0.0
   FV_Atm(1)%flagstruct%d_ext = 0.0
  ! Some default time-splitting options
   FV_Atm(1)%flagstruct%n_split = 0
   FV_Atm(1)%flagstruct%k_split = 1
   if (FV_Atm(1)%flagstruct%ntiles == 6) then
     ! Cubed-sphere grid resolution and DT dependence 
     !              based on ideal remapping DT
      if (FV_Atm(1)%flagstruct%npx >= 48) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/1800.0  )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 90) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 900.0   )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 180) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 450.0   )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 360) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 225.0   )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 720) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/ 112.5   )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 1440) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/  56.25  )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 2880) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/  28.125 )
      endif
      if (FV_Atm(1)%flagstruct%npx >= 5760) then
         FV_Atm(1)%flagstruct%k_split = CEILING(DT/  14.0625)
      endif
   endif
  ! default NonHydrostatic settings (irrelavent to Hydrostatic)
   FV_Atm(1)%flagstruct%beta = 0.0
   FV_Atm(1)%flagstruct%a_imp = 1.0
   FV_Atm(1)%flagstruct%p_fac = 0.1
  ! Resolution dependence (only on full cube, not doubly periodic)
   if (FV_Atm(1)%flagstruct%ntiles == 6) then
     ! Monotonic Hydrostatic defaults
       FV_Atm(1)%flagstruct%hydrostatic = .true.
       FV_Atm(1)%flagstruct%make_nh = .false.
       FV_Atm(1)%flagstruct%vtdm4 = 0.0
       FV_Atm(1)%flagstruct%do_vort_damp = .false.
       FV_Atm(1)%flagstruct%d_con = 0.
       FV_Atm(1)%flagstruct%hord_mt =  10
       FV_Atm(1)%flagstruct%hord_vt =  10
       FV_Atm(1)%flagstruct%hord_tm =  10
       FV_Atm(1)%flagstruct%hord_dp =  10
       FV_Atm(1)%flagstruct%hord_tr =  8
     ! NonMonotonic defaults for c360 (~25km) and finer
       if (FV_Atm(1)%flagstruct%npx >= 360) then
       ! This combination of horizontal advection schemes is critical 
       ! for anomaly correlation NWP skill. 
       ! Using all = 5 (like GFS) produces a substantial degredation in skill
         FV_Atm(1)%flagstruct%hord_mt =  5
         FV_Atm(1)%flagstruct%hord_vt =  6
         FV_Atm(1)%flagstruct%hord_tm =  6
         FV_Atm(1)%flagstruct%hord_dp = -6
       ! This is the best/fastest option for tracers
         FV_Atm(1)%flagstruct%hord_tr =  8
       ! Must now include explicit vorticity damping
         FV_Atm(1)%flagstruct%d_con = 1.
         FV_Atm(1)%flagstruct%do_vort_damp = .true.
         FV_Atm(1)%flagstruct%vtdm4 = 0.02
       endif
     ! NonHydrostatuic defaults for c1440 (~6km) and finer
     !    and continue to adjust vorticity damping with
     !    increasing resolution
       if (FV_Atm(1)%flagstruct%npx >= 1440) then
         FV_Atm(1)%flagstruct%hydrostatic = .false.
         FV_Atm(1)%flagstruct%make_nh = .false.
         FV_Atm(1)%flagstruct%vtdm4 = 0.04
       endif
       if (FV_Atm(1)%flagstruct%npx >= 2880) then
         FV_Atm(1)%flagstruct%vtdm4 = 0.06
       endif
       if (FV_Atm(1)%flagstruct%npx >= 5760) then
         FV_Atm(1)%flagstruct%vtdm4 = 0.08
       endif

   endif

!! Start up FV                   
    call MAPL_TimerOn(MAPL,"--FV_INIT")
    call fv_init2(FV_Atm, DT, grids_on_this_pe, p_split)
    call MAPL_TimerOff(MAPL,"--FV_INIT")
    call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FV_INIT', RC=STATUS )
    VERIFY_(STATUS)

!! Setup GFDL microphysics module
    call gfdl_cloud_microphys_init()

 ASSERT_(DT > 0.0)

  call WRITE_PARALLEL("Dynamics PE Layout ")
  call WRITE_PARALLEL(FV_Atm(1)%layout(1)    ,format='("NPES_X  : ",(   I3))')
  call WRITE_PARALLEL(FV_Atm(1)%layout(2)    ,format='("NPES_Y  : ",(   I3))')

  call WRITE_PARALLEL((/FV_Atm(1)%flagstruct%npx,FV_Atm(1)%flagstruct%npy,FV_Atm(1)%flagstruct%npz/)       , &
    format='("Resolution of dynamics restart     =",3I5)'  )

  ks = FV_Atm(1)%ks ! ALT: this was the value when we read "old" style FV_internal restart
                    !      if needed, we could compute, ks by count(BK==0.0)
                    !      then FV will try to run slightly more efficient code
                    !      So far, GEOS-5 has used ks = 0
  ASSERT_(ks <= FV_Atm(1)%flagstruct%NPZ+1)
  call WRITE_PARALLEL(ks                          , &
     format='("Number of true pressure levels =", I5)'   )

  FV_HYDROSTATIC = FV_Atm(1)%flagstruct%hydrostatic
  DEBUG          = FV_Atm(1)%flagstruct%fv_debug
  prt_minmax     = FV_Atm(1)%flagstruct%fv_debug

  call MAPL_MemUtilsWrite(VM, trim(Iam), RC=STATUS )
  VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)

contains

!-----------------------------------------------------------------------
! BOP
! !IROUTINE:  init_nsplit --- find proper value for nsplit if not specified
!
! !INTERFACE:
  integer function INIT_NSPLIT(dtime,npx,npy)
!
    implicit none

! !INPUT PARAMETERS:
    real (REAL8), intent(in) :: dtime      !  time step
    integer, intent(in)   :: npx,npy    !  Global horizontal resolution

! !DESCRIPTION:
! 
!    If nsplit=0 (module variable) then determine a good value
!    for ns (used in fvdycore) based on resolution and the large-time-step
!    (dtime). The user may have to set this manually if instability occurs.
! 
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!   03.12.04   Sawyer  Moved here from dynamics_vars.  Now a function
!   07.16.07   Putman  Modified for cubed-sphere
!
! EOP
!-----------------------------------------------------------------------
! !LOCAL VARIABLES:
    real (REAL8)   umax
    real (REAL8)   dimx
    real (REAL8)   dim0                      ! base dimension
    real (REAL8)   dt0                       ! base time step              
    real (REAL8)   ns0                       ! base nsplit for base dimension
    integer     ns                        ! final value to be returned
                     
    parameter ( dim0 = 180.  )
    parameter ( dt0  = 1800. )
    parameter ( umax = 350.  )
 
    ns0  = 7.
    dimx = 4.0*npx
    if (FV_Atm(1)%flagstruct%grid_type < 4) then
       ns = nint ( ns0*abs(dtime)*dimx/(dt0*dim0) + 0.49 )
       if (.not. FV_Atm(1)%flagstruct%hydrostatic) ns = ns*1.5 ! time-step needs to be shortened for NH stabilitiy
       ns = max ( 1, ns )
    else
      !ns = nint ( 2.*umax*dtime/sqrt(dx_const**2 + dy_const**2) + 0.49 )
       ns = nint ( ns0*dtime/sqrt(FV_Atm(1)%flagstruct%dx_const**2 + FV_Atm(1)%flagstruct%dy_const**2) + 0.49 )
    endif

    init_nsplit = ns/FV_Atm(1)%flagstruct%k_split

    return
  end function INIT_NSPLIT
!---------------------------------------------------------------------

 end subroutine FV_Setup

 subroutine FV_InitState (STATE, CLOCK, INTERNAL, IMPORT, GC, RC)

  use test_cases_mod, only : test_case, init_double_periodic 

  type (T_FVDYCORE_STATE),pointer              :: STATE

  type (ESMF_Clock), target,     intent(INOUT) :: CLOCK
  type (ESMF_GridComp)         , intent(INOUT) :: GC
  type (ESMF_State)            , intent(INOUT) :: INTERNAL
  type (ESMF_State)            , intent(INOUT) :: IMPORT
  integer, optional            , intent(OUT  ) :: RC

! Local variables

! Pointers to geography info in the MAPL MetaComp

  real,                 pointer :: LATS (:,:)
  real,                 pointer :: LONS (:,:)

  type (ESMF_TimeInterval)     :: Time2Run
  type (ESMF_VM)               :: VM
  type (T_FVDYCORE_GRID) , pointer :: GRID
  integer              :: status
  real(REAL8) :: DT

  integer   :: is ,ie , js ,je    !  Local dims
  integer   :: isc,iec, jsc,jec   !  Local dims
  integer   :: isd,ied, jsd,jed   !  Local dims
  integer   :: k                  !  Vertical loop index
  integer   :: ng
  integer   :: ndt

  integer   :: i,j

  type (ESMF_Time) :: fv_time
  integer :: days, seconds

  character(len=ESMF_MAXSTR)       :: IAm='FV:FV_InitState'

  real(REAL8), pointer                   :: AK(:) => NULL()
  real(REAL8), pointer                   :: BK(:) => NULL()
  real(REAL8), dimension(:,:,:), pointer :: U     => NULL()
  real(REAL8), dimension(:,:,:), pointer :: V     => NULL()
  real(REAL8), dimension(:,:,:), pointer :: PT    => NULL()
  real(REAL8), dimension(:,:,:), pointer :: PE    => NULL()
  real(REAL8), dimension(:,:,:), pointer :: PKZ   => NULL()
  real(REAL8), dimension(:,:,:), pointer :: DZ    => NULL()
  real(REAL8), dimension(:,:,:), pointer :: W     => NULL()
  type (MAPL_MetaComp),          pointer :: mapl  => NULL()

  real(REAL8), ALLOCATABLE :: UA(:,:,:)
  real(REAL8), ALLOCATABLE :: VA(:,:,:)
  real(REAL8), ALLOCATABLE :: UD(:,:,:)
  real(REAL8), ALLOCATABLE :: VD(:,:,:)

  logical    :: hybrid
  integer    :: tile_in
  integer    :: gid, masterproc

! BEGIN

! Retrieve the pointer to the state
! ---------------------------------

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

  call MAPL_GetResource( MAPL, ndt, 'RUN_DT:', default=0, RC=STATUS )
  VERIFY_(STATUS)
  DT = ndt

  STATE%GRID%FVgenstate => MAPL
  GRID => STATE%GRID     ! For convenience
  STATE%DOTIME= .TRUE.
  STATE%DT        = DT
  STATE%NSPLIT    = FV_Atm(1)%flagstruct%N_SPLIT
  GRID%NG     = 3 ; ng = 3
  GRID%NPX    = FV_Atm(1)%flagstruct%NPX-1
  GRID%NPY    = FV_Atm(1)%flagstruct%NPY-1
  GRID%NPZ    = FV_Atm(1)%flagstruct%NPZ
  GRID%NPZ_P1 = FV_Atm(1)%flagstruct%NPZ+1
  GRID%NTILES = 6
  GRID%NQ     = MAX(1,FV_Atm(1)%flagstruct%ncnst)

  masterproc = mpp_root_pe()
  gid = mpp_pe()

  call ESMF_GridCompGet(gc, grid=GRID%GRID, VM=VM, rc=STATUS)
    VERIFY_(STATUS)

  call WRITE_PARALLEL(' ')
  call WRITE_PARALLEL(STATE%DT,format='("Dynamics time step : ",(F10.4))')
  call WRITE_PARALLEL(' ')

! Get pointers to internal state vars
  call MAPL_GetPointer(internal, ak, "AK",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, bk, "BK",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, u, "U",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, v, "V",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pt, "PT",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pe, "PE",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, pkz, "PKZ",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, dz, "DZ",rc=status)
  VERIFY_(STATUS)
  call MAPL_GetPointer(internal, w, "W",rc=status)
  VERIFY_(STATUS)

  call CREATE_VARS ( FV_Atm(1)%bd%isc, FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc, FV_Atm(1)%bd%jec,     &
                     1, FV_Atm(1)%flagstruct%npz, FV_Atm(1)%flagstruct%npz+1,            &
                     U, V, PT, PE, PKZ, DZ, W, &
                     STATE%VARS )
  call MAPL_MemUtilsWrite(VM, 'FV_StateMod: CREATE_VARS', RC=STATUS )
  VERIFY_(STATUS)

  GRID%IS     = FV_Atm(1)%bd%isc
  GRID%IE     = FV_Atm(1)%bd%iec
  GRID%JS     = FV_Atm(1)%bd%jsc
  GRID%JE     = FV_Atm(1)%bd%jec
  GRID%isd    = FV_Atm(1)%bd%isd
  GRID%ied    = FV_Atm(1)%bd%ied
  GRID%jsd    = FV_Atm(1)%bd%jsd
  GRID%jed    = FV_Atm(1)%bd%jed
  if(.not.associated(GRID%AK)) allocate(GRID%AK(size(ak)))
  if(.not.associated(GRID%BK)) allocate(GRID%BK(size(bk)))
  GRID%AK     = ak
  GRID%BK     = bk

! Local Copy of dimensions

  IS     = FV_Atm(1)%bd%isc
  IE     = FV_Atm(1)%bd%iec
  JS     = FV_Atm(1)%bd%jsc
  JE     = FV_Atm(1)%bd%jec
  ISC    = FV_Atm(1)%bd%isc         
  IEC    = FV_Atm(1)%bd%iec 
  JSC    = FV_Atm(1)%bd%jsc 
  JEC    = FV_Atm(1)%bd%jec
  ISD    = FV_Atm(1)%bd%isd
  IED    = FV_Atm(1)%bd%ied
  JSD    = FV_Atm(1)%bd%jsd
  JED    = FV_Atm(1)%bd%jed

  allocate( GRID%DXC(IS:IE,JS:JE) )
  GRID%DXC = fv_atm(1)%gridstruct%dxc(IS:IE,JS:JE)

  allocate( GRID%DYC(IS:IE,JS:JE) )
  GRID%DYC = fv_atm(1)%gridstruct%dyc(IS:IE,JS:JE)

  allocate( GRID%AREA(IS:IE,JS:JE) )
  GRID%AREA = fv_atm(1)%gridstruct%area(IS:IE,JS:JE)
  GRID%GLOBALAREA = fv_atm(1)%gridstruct%globalarea

  if (FV_Atm(1)%flagstruct%grid_type == 4) then
     fv_atm(1)%gridstruct%fC(:,:) = 2.*MAPL_OMEGA*sin(FV_Atm(1)%flagstruct%deglat/180.*MAPL_PI_R8)
     fv_atm(1)%gridstruct%f0(:,:) = 2.*MAPL_OMEGA*sin(FV_Atm(1)%flagstruct%deglat/180.*MAPL_PI_R8)
  else
   if (GRID%f_coriolis_angle == -999) then
     fv_atm(1)%gridstruct%fC(:,:) = 0.0
     fv_atm(1)%gridstruct%f0(:,:) = 0.0
   else
     do j=jsd,jed+1
        do i=isd,ied+1
           fv_atm(1)%gridstruct%fC(i,j) = 2.*MAPL_OMEGA*( -COS(FV_Atm(1)%gridstruct%grid(i,j,1))*COS(FV_Atm(1)%gridstruct%grid(i,j,2))*SIN(GRID%f_coriolis_angle) + &
                                      SIN(FV_Atm(1)%gridstruct%grid(i,j,2))*COS(GRID%f_coriolis_angle) )
        enddo
     enddo
     do j=jsd,jed
        do i=isd,ied
           fv_atm(1)%gridstruct%f0(i,j) = 2.*MAPL_OMEGA*( -COS(FV_Atm(1)%gridstruct%agrid(i,j,1))*COS(FV_Atm(1)%gridstruct%agrid(i,j,2))*SIN(GRID%f_coriolis_angle) + &
                                      SIN(FV_Atm(1)%gridstruct%agrid(i,j,2))*COS(GRID%f_coriolis_angle) )
        enddo
     enddo
   endif
  endif

! Check coordinate information from MAPL_MetaComp
!--------------------------------------------
    call MAPL_Get(MAPL,                &
       LATS          = LATS,           & ! These are in radians
       LONS          = LONS,           & ! These are in radians
       INTERNAL_ESMF_STATE=INTERNAL,   &
                             RC=STATUS )
    VERIFY_(STATUS)

  STATE%CLOCK => CLOCK
  call ESMF_TimeIntervalSet(Time2Run, &
                            S=nint(STATE%DT), rc=status)
  VERIFY_(status)

  STATE%ALARMS(TIME_TO_RUN) = ESMF_AlarmCreate(name="Time2Run", clock=clock, &
                              ringInterval=Time2Run, &
                              Enabled=.TRUE., rc=status) ; VERIFY_(status)
  call ESMF_AlarmEnable(STATE%ALARMS(TIME_TO_RUN), rc=status); VERIFY_(status)
  call ESMF_AlarmRingerOn(STATE%ALARMS(TIME_TO_RUN), rc=status); VERIFY_(status)

  call WRITE_PARALLEL(' ')
  call WRITE_PARALLEL(STATE%DT, &
    format='("INITIALIZED ALARM: DYN_TIME_TO_RUN EVERY ",F9.1," secs.")')

!  Clear wall clock time clocks and global budgets

  STATE%RUN_TIMES = 0
  STATE%NUM_CALLS = 0

  call ESMF_ClockGet( CLOCK, currTime=fv_time, rc=STATUS )
  VERIFY_(STATUS)
  call ESMF_TimeGet( fv_time, dayOfYear=days, s=seconds, rc=STATUS )
  VERIFY_(STATUS)

  ! ---------------------------------------
  ! Create alarm for dry mass fix reporting
  ! ---------------------------------------

  ! Set an interval for printing. Currently hard-coded to
  ! six hours, but could be a MAPL_GetResource value 
  call ESMF_TimeIntervalSet(MassAlarmInt, H=6, rc=STATUS)
  VERIFY_(STATUS)

  ! Create the alarm with the above interval
  MASSALARM = ESMF_AlarmCreate(CLOCK = CLOCK, &
     name         = trim(Iam)//"_MassAlarm",  &
     RingInterval = MassAlarmInt,             &
     RingTime     = fv_time,                  &
     RefTime      = fv_time,                  &
     Enabled      = .true.,                   &
     sticky       = .false.,                  &
     RC           = STATUS                    )
  VERIFY_(STATUS) 

  call MAPL_GetPointer ( import, phis, 'PHIS', RC=STATUS )
  VERIFY_(STATUS)

 ! Set FV3 surface geopotential
  FV_Atm(1)%phis(isc:iec,jsc:jec) = real(phis,kind=REAL8)
  call mpp_update_domains(FV_Atm(1)%phis, FV_Atm(1)%domain, complete=.true.)

  FV_Atm(1)%ak = ak
  FV_Atm(1)%bk = bk
  FV_Atm(1)%ptop = FV_Atm(1)%ak(1)
  FV_Atm(1)%q(:,:,:,:) = 0.0 ! We Don't Have QV from the Import yet

  if (COLDSTART) then

   if (FV_Atm(1)%flagstruct%grid_type == 4) then
       if ( FV_Atm(1)%flagstruct%make_hybrid_z ) then
         hybrid = .false.
       else
         hybrid = FV_Atm(1)%flagstruct%hybrid_z
       endif
       call init_double_periodic(FV_Atm(1)%u,FV_Atm(1)%v,FV_Atm(1)%w,FV_Atm(1)%pt,FV_Atm(1)%delp,FV_Atm(1)%q,FV_Atm(1)%phis, &
                                 FV_Atm(1)%ps,FV_Atm(1)%pe, FV_Atm(1)%peln, &
                                 FV_Atm(1)%pk,FV_Atm(1)%pkz, FV_Atm(1)%uc,FV_Atm(1)%vc, FV_Atm(1)%ua,FV_Atm(1)%va,        &
                                 FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%gridstruct,FV_Atm(1)%flagstruct, &
                                 FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ng, FV_Atm(1)%flagstruct%ncnst, FV_Atm(1)%flagstruct%nwat,  &
                                 FV_Atm(1)%flagstruct%ndims, FV_Atm(1)%flagstruct%ntiles, FV_Atm(1)%flagstruct%dry_mass, FV_Atm(1)%flagstruct%mountain, &
                                 FV_Atm(1)%flagstruct%moist_phys, FV_Atm(1)%flagstruct%hydrostatic, hybrid, FV_Atm(1)%delz, FV_Atm(1)%ze0, &
                                 FV_Atm(1)%ks, FV_Atm(1)%ptop, FV_Atm(1)%domain, tile_in, FV_Atm(1)%bd)
     ! Copy FV to internal State
       call FV_To_State ( STATE )
       if( gid==masterproc ) write(*,*) 'Doubly Periodic IC generated LAT:', FV_Atm(1)%flagstruct%deglat
   else
     ALLOCATE( UA(isc:iec  ,jsc:jec  ,1:FV_Atm(1)%npz) )
     ALLOCATE( VA(isc:iec  ,jsc:jec  ,1:FV_Atm(1)%npz) )
     ALLOCATE( UD(isc:iec  ,jsc:jec+1,1:FV_Atm(1)%npz) )
     ALLOCATE( VD(isc:iec+1,jsc:jec  ,1:FV_Atm(1)%npz) )
     UA(isc:iec,jsc:jec,:) = STATE%VARS%U(isc:iec,jsc:jec,:)
     VA(isc:iec,jsc:jec,:) = STATE%VARS%V(isc:iec,jsc:jec,:)
     call INTERP_AGRID_TO_DGRID( UA, VA, UD, VD )
     STATE%VARS%U(isc:iec,jsc:jec,:) = UD(isc:iec,jsc:jec,:)
     STATE%VARS%V(isc:iec,jsc:jec,:) = VD(isc:iec,jsc:jec,:)
     DEALLOCATE ( UA )
     DEALLOCATE ( VA )
     DEALLOCATE ( UD )
     DEALLOCATE ( VD )
     if ( .not. FV_HYDROSTATIC) then
        FV_Atm(1)%w = 0.0
        W   = 0.0
        PT = PT*PKZ
        call fv_getDELZ(DZ,PT,PE)
        PT = PT/PKZ
     endif
     call State_To_FV( STATE )
   endif ! doubly-periodic

  else ! COLDSTART

    !if ( (.not. FV_HYDROSTATIC) .and. (FV_Atm(1)%flagstruct%Make_NH) ) then
    !   FV_Atm(1)%w = 0.0
    !   W   = 0.0
    !   PT = PT*PKZ
    !   call fv_getDELZ(DZ,PT,PE)
    !   PT = PT/PKZ
    !endif
     call State_To_FV( STATE )

  endif

  if ( (gid==0)                              ) print*, ' '
  if ( (gid==0) .and. (COLDSTART)            ) print*, 'COLDSTARTING FV3'
  if ( (gid==0) .and. (ADIABATIC)            ) print*, 'FV3 being run Adiabatically'
  if ( (gid==0) .and. (.not. FV_HYDROSTATIC) ) print*, 'FV3 being run Non-Hydrostatic'
  if ( (gid==0) .and. (.not. FV_HYDROSTATIC) .and. (FV_Atm(1)%flagstruct%Make_NH) ) print*, 'FV3 Coldstarting Non-Hydrostatic W and DZ'
  if ( (gid==0) .and. (FV_HYDROSTATIC)       ) print*, 'FV3 being run Hydrostatic'
  if ( (gid==0) .and. (SW_DYNAMICS)          ) print*, 'FV3 being run as Shallow-Water Model: test_case=', test_case
  if ( (gid==0) .and. (FV_Atm(1)%flagstruct%grid_type == 4) ) print*, 'FV3 being run as Doubly-Periodic: test_case=', test_case
  if ( (gid==0)                              ) print*, ' '

  if (DEBUG) call debug_fv_state('DEBUG_RESTART',STATE)
!
! Write the vertical coordinate to STDOUT
!
  if( gid.eq.0 .and. .not. SW_DYNAMICS) then
        print *
        write(6,*) ' + denotes a layer within the "dz-filter" of the dynamics'
        write(6,*) ' * denotes a layer within the "sponge-layer" of the dynamics'
        write(6,100)
100     format(2x,' k ','      A(k)    ',2x,' B(k)   ',2x,'  Pref    ',2x,'  DelP',/, &
               1x,'----',3x,'----------',2x,'--------',2x,'----------',2x,'---------' )
          k=0
          if ( (FV_Atm(1)%flagstruct%fv_sg_adj > 0) .AND. (k<=FV_Atm(1)%flagstruct%n_zfilter) ) then
            if (k<=FV_Atm(1)%flagstruct%n_sponge) then
              write(6,101) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k)
            else
              write(6,102) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k)
            endif
          else
            if (k<=FV_Atm(1)%flagstruct%n_sponge) then
              write(6,105) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k)
            else
              write(6,106) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k)
            endif
          endif
          do k=1,ubound(ak,1)
            if ( (FV_Atm(1)%flagstruct%fv_sg_adj > 0) .AND. (k<=FV_Atm(1)%flagstruct%n_zfilter) ) then
              if (k<=FV_Atm(1)%flagstruct%n_sponge) then
                 write(6,103) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k), &
                              (ak(k)-ak(k-1))*0.01 + 1000.0*(bk(k)-bk(k-1))
              else
                 write(6,104) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k), &
                              (ak(k)-ak(k-1))*0.01 + 1000.0*(bk(k)-bk(k-1))
              endif
            else
              if (k<=FV_Atm(1)%flagstruct%n_sponge) then
                 write(6,107) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k), &
                              (ak(k)-ak(k-1))*0.01 + 1000.0*(bk(k)-bk(k-1))
              else
                 write(6,108) k+1,ak(k)*0.01, bk(k), ak(k)*0.01 + 1000.0*bk(k), &
                              (ak(k)-ak(k-1))*0.01 + 1000.0*(bk(k)-bk(k-1))
              endif
            endif
          enddo
        print *
101     format('*+',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
102     format(' +',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
103     format('*+',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
104     format(' +',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
105     format('* ',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
106     format('  ',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4)
107     format('* ',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
108     format('  ',2x,i3,2x,f10.6,2x,f8.4,2x,f10.4,3x,f8.4)
  endif

  call MAPL_MemUtilsWrite(VM, 'FV_StateMod: FV Initialize', RC=STATUS )
  VERIFY_(STATUS)

  RETURN_(ESMF_SUCCESS)

end subroutine FV_InitState

subroutine FV_Run (STATE, CLOCK, GC, RC)

  type (T_FVDYCORE_STATE),pointer              :: STATE

  type (ESMF_Clock), target,     intent(IN   ) :: CLOCK
  type (ESMF_GridComp)         , intent(INOUT) :: GC
  integer, optional            , intent(OUT  ) :: RC

! Local variables
  integer              :: status
  character(len=ESMF_MAXSTR)       :: IAm='FV:FV_Run'

  type (ESMF_Time) :: fv_time
  integer  :: days, seconds
  real(FVPRC) :: time_total, massD

  integer :: i,j,k,n,nn
  integer :: isc,iec,jsc,jec,ng
  integer :: isd,ied,jsd,jed
  integer :: npx,npy,npz

  real(FVPRC), allocatable :: QV(:,:,:)

  real(FVPRC), allocatable :: u_dt(:,:,:)
  real(FVPRC), allocatable :: v_dt(:,:,:)
  real(FVPRC), allocatable :: t_dt(:,:,:)
  real(FVPRC), allocatable :: q_dt(:,:,:,:)
  real(FVPRC), allocatable :: u_srf(:,:)
  real(FVPRC), allocatable :: v_srf(:,:)
  real(FVPRC), allocatable :: ts(:,:)

  real(FVPRC), allocatable :: mass(:,:), tqtot(:,:)
  real(REAL8), allocatable :: ratio(:)
  real(REAL8) :: dpd
  real(FVPRC) :: FQC

! Splitting for Pure Advection
  real(FVPRC) :: myDT, lnp, rdg, pek, ak1

! Convience variables
  integer :: nwat_tracers
  integer :: sphu = -1
  integer :: qliq = -1
  integer :: qice = -1
  integer :: qlls = -1
  integer :: qlcn = -1
  integer :: qils = -1
  integer :: qicn = -1
  integer :: clls = -1
  integer :: clcn = -1
  integer :: qcld = -1
  integer :: rain = -1
  integer :: snow = -1
  integer :: grpl = -1

  type(domain2D) :: domain

  type (MAPL_MetaComp),          pointer :: mapl  => NULL()

  logical :: SPHU_FILLED = .FALSE.
  logical :: QLIQ_FILLED = .FALSE.
  logical :: QICE_FILLED = .FALSE.
  logical :: RAIN_FILLED = .FALSE.
  logical :: SNOW_FILLED = .FALSE.
  logical :: GRPL_FILLED = .FALSE.
  logical :: QCLD_FILLED = .FALSE.
  logical :: QLLS_FILLED = .FALSE.
  logical :: QLCN_FILLED = .FALSE.
  logical :: QILS_FILLED = .FALSE.
  logical :: QICN_FILLED = .FALSE.
  logical :: CLLS_FILLED = .FALSE.
  logical :: CLCN_FILLED = .FALSE.

  logical :: NWAT_TEST

! Begin

! Retrieve the pointer to the state
! ---------------------------------

  call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
  VERIFY_(STATUS)

  call ESMF_ClockGet( CLOCK, currTime=fv_time, rc=STATUS ) 
  VERIFY_(STATUS)
  call ESMF_TimeGet( fv_time, dayOfYear=days, s=seconds, rc=STATUS )
  VERIFY_(STATUS)

  time_total = days*86400. + seconds

  isc = FV_Atm(1)%bd%isc
  iec = FV_Atm(1)%bd%iec
  jsc = FV_Atm(1)%bd%jsc
  jec = FV_Atm(1)%bd%jec
  isd = FV_Atm(1)%bd%isd
  ied = FV_Atm(1)%bd%ied
  jsd = FV_Atm(1)%bd%jsd
  jed = FV_Atm(1)%bd%jed
  npx = FV_Atm(1)%npx
  npy = FV_Atm(1)%npy
  npz = FV_Atm(1)%npz
  ng  = FV_Atm(1)%ng
  domain = FV_Atm(1)%domain

  ! Be sure we have the correct PHIS and number of tracers for this run
   if (fv_first_run) then
    ! Determine how many water species we have
     nwat_tracers = 0
     if (.not. ADIABATIC) then
       do n=1,STATE%GRID%NQ
         if (TRIM(state%vars%tracer(n)%tname) == 'Q'       ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QLCN'    ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QLLS'    ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QICN'    ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QILS'    ) nwat_tracers = nwat_tracers + 1
       enddo
      ! We must have these first 5 at a minimum
       ASSERT_(nwat_tracers == 5)
      ! Check for CLLS, CLCN, QRAIN, QSNOW, QGRAUPEL
       do n=1,STATE%GRID%NQ
         if (TRIM(state%vars%tracer(n)%tname) == 'CLLS'    ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'CLCN'    ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QRAIN'   ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QSNOW'   ) nwat_tracers = nwat_tracers + 1
         if (TRIM(state%vars%tracer(n)%tname) == 'QGRAUPEL') nwat_tracers = nwat_tracers + 1
       enddo
       if (FV_Atm(1)%flagstruct%nwat == 0) then
         if (nwat_tracers >=  5) FV_Atm(1)%flagstruct%nwat = 1 ! Tell FV3 about QV only
         if (.not. FV_Atm(1)%flagstruct%hydrostatic) then
           if (nwat_tracers >=  5) FV_Atm(1)%flagstruct%nwat = 3 ! Tell FV3 about QV, QLIQ, QICE
           if (nwat_tracers == 10) FV_Atm(1)%flagstruct%nwat = 6 ! Tell FV3 about QV, QLIQ, QICE, QRAIN, QSNOW, QGRAUPEL plus QCLD
         endif
       endif
       if (FV_Atm(1)%flagstruct%do_sat_adj) then
          ASSERT_(FV_Atm(1)%flagstruct%nwat == 6)
       endif
       STATE%VARS%nwat = FV_Atm(1)%flagstruct%nwat
     endif
    ! Set FV3 surface geopotential
     FV_Atm(1)%phis(isc:iec,jsc:jec) = real(phis,kind=FVPRC)
     call mpp_update_domains(FV_Atm(1)%phis, FV_Atm(1)%domain, complete=.true.)
    ! How many tracers do we really have?
     ! MAT GCC cannot handle multi-line asserts. For ease of reading,
     !     create a new variable to test
     NWAT_TEST = ( (FV_Atm(1)%flagstruct%nwat == 0) .OR. &
                   (FV_Atm(1)%flagstruct%nwat == 1) .OR. &
                   (FV_Atm(1)%flagstruct%nwat == 3) .OR. &
                   (FV_Atm(1)%flagstruct%nwat == 6) )
     ASSERT_( NWAT_TEST )
     select case ( FV_Atm(1)%flagstruct%nwat )
     case (6) 
          FV_Atm(1)%ncnst = STATE%GRID%NQ + 3 ! NQ + Combined QLIQ,QICE,QCLD
     case (3)
          FV_Atm(1)%ncnst = STATE%GRID%NQ + 2 ! NQ + Combined QLIQ,QICE
     case default
          FV_Atm(1)%ncnst = STATE%GRID%NQ
     end select
     deallocate( FV_Atm(1)%q )
     allocate  ( FV_Atm(1)%q(isd:ied  ,jsd:jed  ,npz, FV_Atm(1)%ncnst) )
    ! Echo FV3 setup
     call echo_fv3_setup()
   endif

   select case ( FV_Atm(1)%flagstruct%nwat )
  ! Assign Tracer Indices for FV3
   case (6)
    sphu = 1
    qliq = 2
    qice = 3
    rain = 4
    snow = 5
    grpl = 6
    qcld = 7
  ! Advect around split CN/LS species still for now...
    qlcn = 8 
    qlls = 9 
    qicn = 10
    qils = 11
    clcn = 12
    clls = 13
   case (3)
    sphu = 1
    qliq = 2
    qice = 3
  ! Advect around split CN/LS species still for now...
    qlcn = 4
    qlls = 5
    qicn = 6
    qils = 7
   case (1)
    sphu = 1
  ! Advect around split CN/LS species
    qlcn = 2
    qlls = 3
    qicn = 4
    qils = 5
   end select

 ! Pull Tracers
  nn = 0
  if (.not. ADIABATIC) then
    select case ( FV_Atm(1)%flagstruct%nwat )
    case (0)
       ASSERT_(FV_Atm(1)%ncnst == STATE%GRID%NQ)
    case (1)
       ASSERT_(FV_Atm(1)%ncnst >= 5)
       ASSERT_(FV_Atm(1)%ncnst == STATE%GRID%NQ)
    case (3)
       ASSERT_(FV_Atm(1)%ncnst >= 7)
       ASSERT_(FV_Atm(1)%ncnst == STATE%GRID%NQ + 2)
    case (6)
       ASSERT_(FV_Atm(1)%ncnst >= 13)
       ASSERT_(FV_Atm(1)%ncnst == STATE%GRID%NQ + 3)
    end select
    FV_Atm(1)%q(:,:,:,:) = 0.0
    if (FV_Atm(1)%flagstruct%nwat > 0) then
    do n=1,STATE%GRID%NQ
       if (TRIM(state%vars%tracer(n)%tname) == 'Q') then
          SPHU_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,sphu) = state%vars%tracer(n)%content_r4(:,:,:)
          else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,sphu) = state%vars%tracer(n)%content(:,:,:)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QLCN') then
         if (FV_Atm(1)%flagstruct%nwat >= 3) then ! QLIQ
           QLIQ_FILLED = .TRUE.
           nn = nn+1
           if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) + state%vars%tracer(n)%content_r4(:,:,:)
           else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) + state%vars%tracer(n)%content(:,:,:)
           endif
         endif
         QLCN_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlcn) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlcn) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QLLS') then
         if (FV_Atm(1)%flagstruct%nwat >= 3) then ! QLIQ
           QLIQ_FILLED = .TRUE.
          ! nn increment already handled in QLCN
           if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) + state%vars%tracer(n)%content_r4(:,:,:)
           else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qliq) + state%vars%tracer(n)%content(:,:,:)
           endif
         endif
         QLLS_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlls) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlls) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QICN') then
         if (FV_Atm(1)%flagstruct%nwat >= 3) then ! QICE
           QICE_FILLED = .TRUE.
           nn = nn+1
           if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) + state%vars%tracer(n)%content_r4(:,:,:)
           else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) + state%vars%tracer(n)%content(:,:,:)
           endif
         endif
         QICN_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qicn) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qicn) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QILS') then
         if (FV_Atm(1)%flagstruct%nwat >= 3) then ! QICE
           QICE_FILLED = .TRUE.
          ! nn increment already handled in QICN
           if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) + state%vars%tracer(n)%content_r4(:,:,:)
           else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qice) + state%vars%tracer(n)%content(:,:,:)
           endif
         endif
         QILS_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qils) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qils) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
     ! Extra species for 6-phase microphysics
       if (FV_Atm(1)%flagstruct%nwat == 6) then
       if (TRIM(state%vars%tracer(n)%tname) == 'QRAIN') then
         RAIN_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,rain) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,rain) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QSNOW') then
         SNOW_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,snow) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,snow) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QGRAUPEL') then
         GRPL_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,grpl) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,grpl) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'CLCN') then
         if (FV_Atm(1)%flagstruct%nwat > 0) then ! QCLD
           QCLD_FILLED = .TRUE.
           nn = nn+1
           if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) + state%vars%tracer(n)%content_r4(:,:,:)
           else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) + state%vars%tracer(n)%content(:,:,:)
           endif
         endif
         CLCN_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clcn) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clcn) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'CLLS') then
         if (FV_Atm(1)%flagstruct%nwat > 0) then ! QCLD
           QCLD_FILLED = .TRUE.
          ! nn increment already handled in CLCN
           if (state%vars%tracer(n)%is_r4) then
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) + state%vars%tracer(n)%content_r4(:,:,:)
           else
             FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qcld) + state%vars%tracer(n)%content(:,:,:)
           endif
         endif
         CLLS_FILLED = .TRUE.
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clls) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clls) = state%vars%tracer(n)%content(:,:,:)
         endif
       endif
       endif !nwat==6
    enddo
   ! Verify
    select case (FV_Atm(1)%flagstruct%nwat)
    case (6)
      ASSERT_(nn == 13) ! Q, QLCN, QLLS, QICN, QILS, CLLS, CLCN, QRAIN, QSNOW, QGRAUPEL, QLIQ, QICE, QCLD
      ASSERT_(SPHU_FILLED)
      ASSERT_(QLIQ_FILLED)
      ASSERT_(QICE_FILLED)
      ASSERT_(RAIN_FILLED)
      ASSERT_(SNOW_FILLED)
      ASSERT_(GRPL_FILLED)
      ASSERT_(QCLD_FILLED)
      ASSERT_(QLCN_FILLED)
      ASSERT_(QLLS_FILLED)
      ASSERT_(QICN_FILLED)
      ASSERT_(QILS_FILLED)
      ASSERT_(CLCN_FILLED)
      ASSERT_(CLLS_FILLED)
    case (3)
      ASSERT_(nn == 7) ! Q, QLCN, QLLS, QICN, QILS, QLIQ, QICE
      ASSERT_(SPHU_FILLED)
      ASSERT_(QLIQ_FILLED)
      ASSERT_(QICE_FILLED)
      ASSERT_(QLCN_FILLED)
      ASSERT_(QLLS_FILLED)
      ASSERT_(QICN_FILLED)
      ASSERT_(QILS_FILLED)
    case (1)
      ASSERT_(nn == 5) ! Q, QLCN, QLLS, QICN, QILS
      ASSERT_(SPHU_FILLED)
      ASSERT_(QLCN_FILLED)
      ASSERT_(QLLS_FILLED)
      ASSERT_(QICN_FILLED)
      ASSERT_(QILS_FILLED)
    end select
    endif !nwat > 0
      select case (FV_Atm(1)%flagstruct%nwat)
      case (0)
       do n=1,STATE%GRID%NQ
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content(:,:,:)
         endif
       enddo
      case (1)
       do n=1,STATE%GRID%NQ
         if ( (TRIM(state%vars%tracer(n)%tname) /= 'Q'       ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QLCN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QLLS'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QICN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QILS'    ) ) then
           nn=nn+1
           if (state%vars%tracer(n)%is_r4) then
              FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content_r4(:,:,:)
           else
              FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content(:,:,:)
           endif
         endif
       enddo
      case (3)
       do n=1,STATE%GRID%NQ
         if ( (TRIM(state%vars%tracer(n)%tname) /= 'Q'       ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QLCN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QLLS'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QICN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QILS'    ) ) then
           nn=nn+1
           if (state%vars%tracer(n)%is_r4) then
              FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content_r4(:,:,:)
           else
              FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content(:,:,:)
           endif
         endif
       enddo
      case (6)
       do n=1,STATE%GRID%NQ
         if ( (TRIM(state%vars%tracer(n)%tname) /= 'Q'       ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QLCN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QLLS'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QICN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QILS'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'CLCN'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'CLLS'    ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QRAIN'   ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QSNOW'   ) .and. &
              (TRIM(state%vars%tracer(n)%tname) /= 'QGRAUPEL') ) then
           nn=nn+1
           if (state%vars%tracer(n)%is_r4) then
              FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content_r4(:,:,:)
           else
              FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content(:,:,:)
           endif
         endif
       enddo
      end select
      ASSERT_(nn == FV_Atm(1)%ncnst)
  else
    if (mpp_pe()==0) print*, 'Running In Adiabatic Mode'
      do n=1,STATE%GRID%NQ
         nn = nn+1
         if (state%vars%tracer(n)%is_r4) then
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content_r4(:,:,:)
         else
            FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn) = state%vars%tracer(n)%content(:,:,:)
         endif
      enddo
      ASSERT_(nn == FV_Atm(1)%ncnst)
  endif

    myDT = state%dt

    elapsed_time = elapsed_time + myDT

    if (DEBUG) call debug_fv_state('Before Dynamics Execution',STATE)

! Update FV with Internal State
    call State_To_FV( STATE )

    ! Query for PSDRY from AGCM.rc and set to MAPL_PSDRY if not found
    call MAPL_GetResource( MAPL, massD0, 'PSDRY:', default=MAPL_PSDRY, RC=STATUS )
    VERIFY_(STATUS)
    FV_Atm(1)%flagstruct%dry_mass = massD0

    if (fv_first_run) then
     ! Make_NH
      if ( .not. FV_Atm(1)%flagstruct%hydrostatic ) then
        if (all(FV_Atm(1)%w(isc:iec,jsc:jec,:) == 0.0)) FV_Atm(1)%flagstruct%Make_NH = .true.
        if ( FV_Atm(1)%flagstruct%Make_NH ) then
          if (mpp_pe()==0) print*, 'fv_first_run: FV3 is making Non-Hydrostatic W and DZ'
          call p_var(FV_Atm(1)%npz,         isc,         iec,       jsc,     jec,  FV_Atm(1)%ptop,     ptop_min,  &
                     FV_Atm(1)%delp, FV_Atm(1)%delz, FV_Atm(1)%pt, FV_Atm(1)%ps, FV_Atm(1)%pe,  FV_Atm(1)%peln,   &
                     FV_Atm(1)%pk,   FV_Atm(1)%pkz, kappa, FV_Atm(1)%q, FV_Atm(1)%ng, &
                     FV_Atm(1)%ncnst, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%flagstruct%dry_mass,  &
                     FV_Atm(1)%flagstruct%adjust_dry_mass,  FV_Atm(1)%flagstruct%mountain, &
                     FV_Atm(1)%flagstruct%moist_phys,  FV_Atm(1)%flagstruct%hydrostatic, &
                     FV_Atm(1)%flagstruct%nwat, FV_Atm(1)%domain, FV_Atm(1)%flagstruct%make_nh)
          FV_Atm(1)%flagstruct%Make_NH=.false. 
        endif
      endif
     ! Mark FV setup complete
      fv_first_run = .false.
    endif

! Check Dry Mass (Apply fixer is option is enabled)
   if ( check_mass .OR. fix_mass ) then
      call MAPL_TimerOn(MAPL,"--MASS_FIX")

      if ( FV_Atm(1)%flagstruct%adjust_dry_mass .AND. &
            ((.not. FV_Atm(1)%flagstruct%hydrostatic) .OR. FV_Atm(1)%flagstruct%nwat==6)  ) then

         call p_var(FV_Atm(1)%npz,         isc,         iec,       jsc,     jec,  FV_Atm(1)%ptop,     ptop_min,  &
                    FV_Atm(1)%delp, FV_Atm(1)%delz, FV_Atm(1)%pt, FV_Atm(1)%ps, FV_Atm(1)%pe,  FV_Atm(1)%peln,   &
                    FV_Atm(1)%pk,   FV_Atm(1)%pkz, kappa, FV_Atm(1)%q, FV_Atm(1)%ng, &
                    FV_Atm(1)%ncnst, FV_Atm(1)%gridstruct%area_64, FV_Atm(1)%flagstruct%dry_mass,  &
                    FV_Atm(1)%flagstruct%adjust_dry_mass,  FV_Atm(1)%flagstruct%mountain, &
                    FV_Atm(1)%flagstruct%moist_phys,  FV_Atm(1)%flagstruct%hydrostatic, &
                    FV_Atm(1)%flagstruct%nwat, FV_Atm(1)%domain, FV_Atm(1)%flagstruct%make_nh)

      else

      allocate ( mass(isc:iec,jsc:jec))
      allocate (tqtot(isc:iec,jsc:jec))
      do j=jsc,jec
         do i=isc,iec
            mass(i,j) = FV_Atm(1)%pe(i,npz+1,j)
         enddo
      enddo
      tqtot = 0.0
      if ( (.not. ADIABATIC) .AND. (FV_Atm(1)%flagstruct%nwat /= 0) ) then
       if (FV_Atm(1)%flagstruct%nwat == 6) then
            tqtot(:,:) = tqtot(:,:) + ( &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,sphu) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,qliq) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,qice) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,rain) + & 
            FV_Atm(1)%q(isc:iec,jsc:jec,k,snow) + & 
            FV_Atm(1)%q(isc:iec,jsc:jec,k,grpl) ) * FV_Atm(1)%delp(isc:iec,jsc:jec,k)
       else
         do k=1,npz
            tqtot(:,:) = tqtot(:,:) + ( &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,sphu) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,qlcn) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,qlls) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,qicn) + &
            FV_Atm(1)%q(isc:iec,jsc:jec,k,qils) ) * FV_Atm(1)%delp(isc:iec,jsc:jec,k)
         enddo
       endif
      endif

      massD  = g_sum(FV_Atm(1)%domain, mass-tqtot, isc, iec, jsc, jec, state%grid%ng, fv_atm(1)%gridstruct%area_64, 1, reproduce = .true.)

      ! If PSDRY is negative, set to use the incoming drymass.
      ! NOTE: THIS WILL NOT TIME REGRESS
      if (massD0 < 0.0d0) then
         massD0 = massD
      end if

      if(ESMF_AlarmIsRinging(MASSALARM) .AND. check_mass) then
         if (ABS((massD-massD0)/massD0) >= epsilon(1.0_REAL4)) then
            if (mpp_pe()==mpp_root_pe()) then
               write(6,126) 
               write(6,127) massD0
               write(6,128) massD
               write(6,129) massD0/massD, (massD-massD0)/massD0
               126 format('Dry Mass Difference of > epsilon relative found!')
               127 format('Dry Mass Expected:'2x,g21.14)
               128 format('Dry Mass Calculated in FV3:'2x,g21.14)
               129 format('Dry Mass Scaling Factor and Relative Difference'2x,g21.14,2x,g21.14)
            end if
         end if
      end if

      ! Apply the fixer if asked for
      if (fix_mass) then

         if (ABS((massD-massD0)/massD0) >= epsilon(1.0_REAL4)) then
            if (mpp_pe()==mpp_root_pe()) then
               write(6,119) massD0, massD, massD0/massD, (massD-massD0)/massD0
               119 format('Dry Mass Violation (epsilon relative)!'2x,g21.14,2x,g21.14,2x,g21.14,2x,g21.14)
            end if
         end if

         ! Accumulate added dry mass from the fixer
         massADDED = massADDED + (massD-massD0)/massD0

         ! -----------------------------------------------
         ! Fix Dry Mass after increments have been applied
         ! -----------------------------------------------
         FV_Atm(1)%pe = FV_Atm(1)%pe*massD0/massD

         if(ESMF_AlarmIsRinging(MASSALARM) .AND. check_mass) then
            if (mpp_pe()==mpp_root_pe()) then
               write(6,109) massD0, massD, massD0/massD, (massD-massD0)/massD0, massADDED
               109    format('Dry Mass Fixer'2x,g21.14,2x,g21.14,2x,g21.14,2x,g21.14,2x,g21.14)
            end if
         end if
      else
         ! Check Dry Mass Conservation and write to log
         if(ESMF_AlarmIsRinging(MASSALARM) .AND. mpp_pe()==mpp_root_pe()) write(6,110) massD0, massD, massD0/massD, (massD-massD0)/massD0
         110    format('Dry Mass Check'2x,g21.14,2x,g21.14,2x,g21.14,2x,g21.14)
      endif

      deallocate (mass)
      deallocate (tqtot)

      endif

      call MAPL_TimerOff(MAPL,"--MASS_FIX")
   endif

          ! if (mpp_pe()==mpp_root_pe()) then
          !    write(6,*) 'Advecting tracers: ', FV_Atm(1)%ncnst, STATE%GRID%NQ
          ! endif

    call MAPL_TimerOn(MAPL,"--FV_DYNAMICS")
    if (.not. FV_OFF) then
    call set_domain(FV_Atm(1)%domain)  ! needed for diagnostic output done in fv_dynamics
    call fv_dynamics(FV_Atm(1)%npx, FV_Atm(1)%npy, FV_Atm(1)%npz, FV_Atm(1)%ncnst, FV_Atm(1)%ng,   &
                     myDT, FV_Atm(1)%flagstruct%consv_te, FV_Atm(1)%flagstruct%fill, FV_Atm(1)%flagstruct%reproduce_sum, kappa,   &
                     cp, zvir, FV_Atm(1)%ptop, FV_Atm(1)%ks, FV_Atm(1)%flagstruct%ncnst, FV_Atm(1)%flagstruct%n_split, FV_Atm(1)%flagstruct%q_split, &
                     FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%w, FV_Atm(1)%delz,       &
                     FV_Atm(1)%flagstruct%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%delp, FV_Atm(1)%q, FV_Atm(1)%ps,       &
                     FV_Atm(1)%pe, FV_Atm(1)%pk, FV_Atm(1)%peln, FV_Atm(1)%pkz,                         &
                     FV_Atm(1)%phis, FV_Atm(1)%q_con, FV_Atm(1)%omga, FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%uc, FV_Atm(1)%vc,  &
                     FV_Atm(1)%ak, FV_Atm(1)%bk, FV_Atm(1)%mfx, FV_Atm(1)%mfy, FV_Atm(1)%cx, FV_Atm(1)%cy,    &
                     FV_Atm(1)%ze0, FV_Atm(1)%flagstruct%hybrid_z, FV_Atm(1)%gridstruct, FV_Atm(1)%flagstruct, &
                     FV_Atm(1)%neststruct, FV_Atm(1)%idiag, FV_Atm(1)%bd, FV_Atm(1)%parent_grid, FV_Atm(1)%domain, FV_Atm(1)%diss_est, time_total)

    if ( FV_Atm(1)%flagstruct%fv_sg_adj > 0 ) then
         allocate ( u_dt(isd:ied,jsd:jed,npz) )
         allocate ( v_dt(isd:ied,jsd:jed,npz) )
         allocate ( t_dt(isc:iec,jsc:jec,npz) )
         u_dt(:,:,:) = 0.0
         v_dt(:,:,:) = 0.0
         t_dt(:,:,:) = 0.0
         call fv_subgrid_z(isd, ied, jsd, jed, isc, iec, jsc, jec, FV_Atm(1)%npz, &
                           FV_Atm(1)%ncnst, myDT, FV_Atm(1)%flagstruct%fv_sg_adj,      &
                           FV_Atm(1)%flagstruct%nwat, FV_Atm(1)%delp, FV_Atm(1)%pe,     &
                           FV_Atm(1)%peln, FV_Atm(1)%pkz, FV_Atm(1)%pt, FV_Atm(1)%q,       &
                           FV_Atm(1)%ua, FV_Atm(1)%va, FV_Atm(1)%flagstruct%hydrostatic,&
                           FV_Atm(1)%w, FV_Atm(1)%delz, u_dt, v_dt, t_dt, FV_Atm(1)%flagstruct%n_zfilter)
         deallocate ( u_dt )
         deallocate ( v_dt )
         deallocate ( t_dt )
    endif

    call nullify_domain()

    endif
    call MAPL_TimerOff(MAPL,"--FV_DYNAMICS")

! Copy FV to internal State
   call FV_To_State ( STATE )

  SPHU_FILLED = .FALSE.
  QLIQ_FILLED = .FALSE.
  QICE_FILLED = .FALSE.
  RAIN_FILLED = .FALSE.
  SNOW_FILLED = .FALSE.
  GRPL_FILLED = .FALSE.
  QCLD_FILLED = .FALSE.
  QLLS_FILLED = .FALSE.
  QLCN_FILLED = .FALSE.
  QILS_FILLED = .FALSE.
  QICN_FILLED = .FALSE.
  CLLS_FILLED = .FALSE.
  CLCN_FILLED = .FALSE.

 ! Push Tracers
  nn = 0
  if (.not. ADIABATIC) then

     ! Redistribute CN/LS liq, ice and cld condensate based on advected CN/LS species
     if (FV_Atm(1)%flagstruct%nwat >= 3) then
      do k=1,npz
         do j=jsc,jec
            do i=isc,iec
              ! LIQUID
               FQC = 0.0
               if ( FV_Atm(1)%q(i,j,k,qlcn)+FV_Atm(1)%q(i,j,k,qlls) > 0.0 ) then
                  FQC = MIN(1.0, MAX(0.0,FV_Atm(1)%q(i,j,k,qlcn)) / (FV_Atm(1)%q(i,j,k,qlcn)+FV_Atm(1)%q(i,j,k,qlls)))
               endif
               FV_Atm(1)%q(i,j,k,qlcn) = FV_Atm(1)%q(i,j,k,qliq)*(    FQC)
               FV_Atm(1)%q(i,j,k,qlls) = FV_Atm(1)%q(i,j,k,qliq)*(1.0-FQC)
              ! ICE
               FQC = 0.0
               if ( FV_Atm(1)%q(i,j,k,qicn)+FV_Atm(1)%q(i,j,k,qils) > 0.0 ) then
                  FQC = MIN(1.0, MAX(0.0,FV_Atm(1)%q(i,j,k,qicn)) / (FV_Atm(1)%q(i,j,k,qicn)+FV_Atm(1)%q(i,j,k,qils)))
               endif
               FV_Atm(1)%q(i,j,k,qicn) = FV_Atm(1)%q(i,j,k,qice)*(    FQC)
               FV_Atm(1)%q(i,j,k,qils) = FV_Atm(1)%q(i,j,k,qice)*(1.0-FQC)
              ! CLOUD
               if (FV_Atm(1)%flagstruct%nwat == 6) then
                  FQC = 0.0
                  if ( FV_Atm(1)%q(i,j,k,clcn)+FV_Atm(1)%q(i,j,k,clls) > 0.0 ) then
                     FQC = MIN(1.0, MAX(0.0,FV_Atm(1)%q(i,j,k,clcn)) / (FV_Atm(1)%q(i,j,k,clcn)+FV_Atm(1)%q(i,j,k,clls)))
                  endif
                  FV_Atm(1)%q(i,j,k,clcn) = FV_Atm(1)%q(i,j,k,qcld)*(    FQC)
                  FV_Atm(1)%q(i,j,k,clls) = FV_Atm(1)%q(i,j,k,qcld)*(1.0-FQC)
               endif
            enddo
         enddo
      enddo
      if (FV_Atm(1)%flagstruct%nwat == 3) then
        nn = nn+2
        QLIQ_FILLED = .TRUE.
        QICE_FILLED = .TRUE.
      endif
      if (FV_Atm(1)%flagstruct%nwat == 6) then
        nn = nn+3
        QLIQ_FILLED = .TRUE.
        QICE_FILLED = .TRUE.
        QCLD_FILLED = .TRUE.
      endif
     endif

     do n=1,STATE%GRID%NQ

       if (FV_Atm(1)%flagstruct%nwat >= 1) then
       if (TRIM(state%vars%tracer(n)%tname) == 'Q') then
          SPHU_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,sphu)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,sphu)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QLCN') then
          QLCN_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlcn)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlcn)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QLLS') then
          QLLS_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlls)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qlls)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QICN') then
          QICN_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qicn)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qicn)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QILS') then
          QILS_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qils)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,qils)
          endif
       endif
       endif ! nwat >= 1
       
       if (FV_Atm(1)%flagstruct%nwat == 6) then
       if (TRIM(state%vars%tracer(n)%tname) == 'QRAIN') then
          RAIN_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,rain)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,rain)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QSNOW') then
          SNOW_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,snow)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,snow)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'QGRAUPEL') then
          GRPL_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,grpl)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,grpl)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'CLCN') then
          CLCN_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clcn)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clcn)
          endif
       endif
       if (TRIM(state%vars%tracer(n)%tname) == 'CLLS') then
          CLLS_FILLED = .TRUE.
          nn = nn+1
          if (state%vars%tracer(n)%is_r4) then
             state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clls)
          else
                state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,clls)
          endif
       endif
       endif ! nwat == 6
     enddo

   ! Verify
    select case (FV_Atm(1)%flagstruct%nwat)
    case (6)
      ASSERT_(nn == 13) ! Q, QLCN, QLLS, QICN, QILS, CLLS, CLCN, QRAIN, QSNOW, QGRAUPEL, QLIQ, QICE, QCLD
      ASSERT_(SPHU_FILLED)
      ASSERT_(QLIQ_FILLED)
      ASSERT_(QICE_FILLED)
      ASSERT_(RAIN_FILLED)
      ASSERT_(SNOW_FILLED)
      ASSERT_(GRPL_FILLED)
      ASSERT_(QCLD_FILLED)
      ASSERT_(QLCN_FILLED)
      ASSERT_(QLLS_FILLED)
      ASSERT_(QICN_FILLED)
      ASSERT_(QILS_FILLED)
      ASSERT_(CLCN_FILLED)
      ASSERT_(CLLS_FILLED)
    case (3)
      ASSERT_(nn == 7) ! Q, QLCN, QLLS, QICN, QILS, QLIQ, QICE
      ASSERT_(SPHU_FILLED)
      ASSERT_(QLIQ_FILLED)
      ASSERT_(QICE_FILLED)
      ASSERT_(QLCN_FILLED)
      ASSERT_(QLLS_FILLED)
      ASSERT_(QICN_FILLED)
      ASSERT_(QILS_FILLED)
    case (1)
      ASSERT_(nn == 5) ! Q, QLCN, QLLS, QICN, QILS
      ASSERT_(SPHU_FILLED)
      ASSERT_(QLCN_FILLED)
      ASSERT_(QLLS_FILLED)
      ASSERT_(QICN_FILLED)
      ASSERT_(QILS_FILLED)
    end select

      select case(FV_Atm(1)%flagstruct%nwat)
      case (0)
       do n=1,STATE%GRID%NQ
         nn=nn+1
         if (state%vars%tracer(n)%is_r4) then
            state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         else
            state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         endif
       enddo
      case (1)
       do n=1,STATE%GRID%NQ
        if ((TRIM(state%vars%tracer(n)%tname) /= 'Q'       ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QLCN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QLLS'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QICN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QILS'    ) ) then
         nn=nn+1
         if (state%vars%tracer(n)%is_r4) then
            state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         else
            state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         endif
        endif
       enddo
      case (3)
       do n=1,STATE%GRID%NQ
        if ((TRIM(state%vars%tracer(n)%tname) /= 'Q'       ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QLCN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QLLS'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QICN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QILS'    ) ) then
         nn=nn+1
         if (state%vars%tracer(n)%is_r4) then
            state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         else
            state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         endif
        endif
       enddo
      case (6)
       do n=1,STATE%GRID%NQ
        if ((TRIM(state%vars%tracer(n)%tname) /= 'Q'       ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QLCN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QLLS'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QICN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QILS'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'CLCN'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'CLLS'    ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QRAIN'   ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QSNOW'   ) .and. &
            (TRIM(state%vars%tracer(n)%tname) /= 'QGRAUPEL') ) then
         nn=nn+1
         if (state%vars%tracer(n)%is_r4) then
            state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         else
            state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         endif
        endif
       enddo
      end select
      ASSERT_(nn == FV_Atm(1)%ncnst)
  else
      do n=1,STATE%GRID%NQ
         nn=nn+1
         if (state%vars%tracer(n)%is_r4) then
            state%vars%tracer(n)%content_r4(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         else
            state%vars%tracer(n)%content(:,:,:) = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,nn)
         endif
      enddo
      ASSERT_(nn == FV_Atm(1)%ncnst)
  endif

    if (DEBUG) call debug_fv_state('After Dynamics Execution',STATE)

    RETURN_(ESMF_SUCCESS)

end subroutine FV_Run

 subroutine FV_DA_Incs (LONS,LATS,IM,JM,KM,&
                        u_amb, v_amb, t_amb, dp_amb, q_amb, o3_amb, &
                        u_inc, v_inc, t_inc, dp_inc, q_inc, o3_inc)
 use fv_treat_da_inc_mod, only : geos_get_da_increments
! The ANA Lat-Lon Grid
 integer, intent(in) :: IM,JM,KM 
 real, intent(inout) :: LONS(IM), LATS(JM)
! ANA-BKG fields on the ANA Lat-Lon Grid
 real, dimension(IM,JM,KM), intent(inout) :: u_amb, v_amb, t_amb, dp_amb, q_amb, o3_amb
! increments are returned on Cubed A-Grid with Native Cubed-Sphere Wind Incs on the D-Grid
 real, dimension(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec, &
                 FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec, &
                 FV_Atm(1)%npz), intent(out) :: u_inc, v_inc, t_inc, dp_inc, q_inc, o3_inc

!call geos_get_da_increments(FV_Atm, FV_Atm(1)%domain, LONS,LATS,IM,JM,KM, &
!                     u_amb, v_amb, t_amb, dp_amb, q_amb, o3_amb,  &
!                     u_inc, v_inc, t_inc, dp_inc, q_inc, o3_inc)

 end subroutine FV_DA_Incs

 subroutine FV_Finalize (STATE)

  use fv_control_mod, only : fv_end

  type (T_FVDYCORE_STATE),pointer              :: STATE

  integer isc, iec, jsc, jec
  integer isd, ied, jsd, jed
  integer npz, ng

  isc = FV_Atm(1)%bd%isc
  iec = FV_Atm(1)%bd%iec
  jsc = FV_Atm(1)%bd%jsc
  jec = FV_Atm(1)%bd%jec
  isd = FV_Atm(1)%bd%isd
  ied = FV_Atm(1)%bd%ied
  jsd = FV_Atm(1)%bd%jsd
  jed = FV_Atm(1)%bd%jed
  npz = FV_Atm(1)%npz
  ng  = FV_Atm(1)%ng

    if (DEBUG) call debug_fv_state('FV_Finalize',STATE)

      call timing_off('TOTAL')
      call fv_end(FV_Atm, grids_on_this_pe ,.false.)

#if defined( MAPL_MODE )
!    call ESMF_GridDestroy  (STATE%GRID%GRID)
#endif

 end subroutine FV_Finalize

subroutine State_To_FV ( STATE )

! !INPUT PARAMETERS:

   type(T_FVDYCORE_STATE),      pointer   :: STATE

    integer               :: ISC,IEC, JSC,JEC
    integer               :: ISD,IED, JSD,JED 
    integer               :: KM, NG
    integer               :: I,J,K
    real(REAL8)              :: akap

    real(REAL8) :: uatemp(state%grid%isd:state%grid%ied, &
                          state%grid%jsd:state%grid%jed,state%grid%npz)
    real(REAL8) :: vatemp(state%grid%isd:state%grid%ied, &
                          state%grid%jsd:state%grid%jed,state%grid%npz)

    real(FVPRC) :: wbuffer(state%grid%js:state%grid%je,state%grid%npz)
    real(FVPRC) :: sbuffer(state%grid%is:state%grid%ie,state%grid%npz)
    real(FVPRC) :: ebuffer(state%grid%js:state%grid%je,state%grid%npz)
    real(FVPRC) :: nbuffer(state%grid%is:state%grid%ie,state%grid%npz)

    ISC = state%grid%is
    IEC = state%grid%ie
    JSC = state%grid%js
    JEC = state%grid%je
    ISD = state%grid%isd
    IED = state%grid%ied
    JSD = state%grid%jsd
    JED = state%grid%jed
    KM  = state%grid%npz
    NG  = state%grid%ng

    akap  = kappa
    if (SW_DYNAMICS) akap  = 1.

!------------
! Update Winds
!------------

! D-Grid
  FV_Atm(1)%u(:,:,:) = tiny_number
  FV_Atm(1)%v(:,:,:) = tiny_number
  if (FV_Atm(1)%flagstruct%grid_type>=4) then
  ! Doubly Periodic
    uatemp(isc:iec,jsc:jec,:) = STATE%VARS%U
    vatemp(isc:iec,jsc:jec,:) = STATE%VARS%V
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    FV_Atm(1)%u(isc:iec,jsc:jec+1,:) = uatemp(isc:iec,jsc:jec+1,:)
    FV_Atm(1)%v(isc:iec+1,jsc:jec,:) = vatemp(isc:iec+1,jsc:jec,:)
  else
    FV_Atm(1)%u(isc:iec,jsc:jec,:) = STATE%VARS%U
    FV_Atm(1)%v(isc:iec,jsc:jec,:) = STATE%VARS%V
    call mpp_get_boundary(FV_Atm(1)%u, FV_Atm(1)%v, FV_Atm(1)%domain, &
                          wbuffery=wbuffer, ebuffery=ebuffer, &
                          sbufferx=sbuffer, nbufferx=nbuffer, &
                          gridtype=DGRID_NE, complete=.true. )
    do k=1,km
       do i=isc,iec
          FV_Atm(1)%u(i,jec+1,k) = nbuffer(i,k)
       enddo
    enddo
    do k=1,km
       do j=jsc,jec
          FV_Atm(1)%v(iec+1,j,k) = ebuffer(j,k)
       enddo
    enddo
  endif

   if (.not. FV_Atm(1)%flagstruct%hydrostatic) FV_Atm(1)%w(isc:iec,jsc:jec,:) = STATE%VARS%W
 
!------------
! Update Pressures
!------------

   FV_Atm(1)%pe(:,:,:) = tiny_number
   if (SW_DYNAMICS) then
      do k=1,km+1
        do j=jsc,jec
          do i=isc,iec
            FV_Atm(1)%pe(i,k,j)   = STATE%VARS%PE(i,j,k)
          enddo
        enddo
      enddo
   else
      do k=1,km+1
        do j=jsc,jec
          do i=isc,iec
            FV_Atm(1)%pe(i,k,j)   = STATE%VARS%PE(i,j,k)
          enddo
        enddo
      enddo

      do k=1,km+1
        do j=jsc,jec
          do i=isc,iec
            FV_Atm(1)%peln(i,k,j) = log(FV_Atm(1)%pe(i,k,j))
          enddo
        enddo
      enddo

      do k=1,km+1
        do j=jsc,jec
          do i=isc,iec
            FV_Atm(1)%pk(i,j,k)   = exp( akap*FV_Atm(1)%peln(i,k,j) )
          enddo
        enddo
      enddo

      FV_Atm(1)%ps(isc:iec,jsc:jec) = FV_Atm(1)%pe(isc:iec,km+1,jsc:jec)

   endif

    FV_Atm(1)%delp(:,:,:) = tiny_number
    do k=1,km
      do j=jsc,jec
        do i=isc,iec
          FV_Atm(1)%delp(i,j,k) = FV_Atm(1)%pe(i,k+1,j) - FV_Atm(1)%pe(i,k,j) 
        enddo
      enddo
    enddo

    if (.not. SW_DYNAMICS) then

!-----------------------
! Copy PT and make Dry T 
!-----------------------
       FV_Atm(1)%pt(:,:,:) = tiny_number
       FV_Atm(1)%pt(isc:iec,jsc:jec,:) = STATE%VARS%PT*STATE%VARS%PKZ

!------------
! Get delz
!------------
       if (.not. FV_Atm(1)%flagstruct%hydrostatic) FV_Atm(1)%delz(isc:iec,jsc:jec,:) = STATE%VARS%DZ

!------------------------------------------------------------------------------
! Get pkz
!------------------------------------------------------------------------------
       FV_Atm(1)%pkz(isc:iec,jsc:jec,:) = STATE%VARS%PKZ

    endif


   return

end subroutine State_To_FV

subroutine FV_To_State ( STATE )

!
! !INPUT PARAMETERS:

   type(T_FVDYCORE_STATE),      pointer   :: STATE

    integer               :: ISC,IEC, JSC,JEC, KM
    integer               :: I,J

    ISC = state%grid%is
    IEC = state%grid%ie
    JSC = state%grid%js
    JEC = state%grid%je
    KM  = state%grid%npz

! Copy updated FV data to internal state
    STATE%VARS%U(:,:,:) = FV_Atm(1)%u(isc:iec,jsc:jec,:)
    STATE%VARS%V(:,:,:) = FV_Atm(1)%v(isc:iec,jsc:jec,:)
    if (.not. FV_Atm(1)%flagstruct%hydrostatic) STATE%VARS%W = FV_Atm(1)%w(isc:iec,jsc:jec,:)

    if (SW_DYNAMICS) then
       STATE%VARS%PE(:,:,1) = FV_Atm(1)%phis(isc:iec,jsc:jec)
       STATE%VARS%PE(:,:,2) = FV_Atm(1)%phis(isc:iec,jsc:jec) + FV_Atm(1)%delp(isc:iec,jsc:jec,1)
    else
       do j=jsc,jec
          do i=isc,iec
             STATE%VARS%PE(i,j,:) = FV_Atm(1)%pe(i,:,j)
          enddo
       enddo

!-----------------------------------
! Fill Dry Temperature to PT
!-----------------------------------
       STATE%VARS%PT  = FV_Atm(1)%pt(isc:iec,jsc:jec,:)

!------------------------------
! Get delz from FV3
!------------------------------
       if (.not. FV_Atm(1)%flagstruct%hydrostatic) STATE%VARS%DZ = FV_Atm(1)%delz(isc:iec,jsc:jec,:)
       
!--------------------------------
! Get pkz from FV3
!--------------------------------
       STATE%VARS%PKZ = FV_Atm(1)%pkz(isc:iec,jsc:jec,:)

!---------------------------------------------------------------------
! Convert to Dry Temperature to PT with hydrostatic pkz
!---------------------------------------------------------------------
       STATE%VARS%PT  = STATE%VARS%PT/STATE%VARS%PKZ
    endif

   return

end subroutine FV_To_State

subroutine fv_getDELZ(delz,temp,pe)
  real(REAL8), intent(OUT) :: delz(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent( IN) :: temp(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent( IN) ::   pe(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  real(REAL8) :: peln(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  real(REAL8) :: rdg
  peln = log(pe)
  rdg   = -rgas / grav
  delz = rdg*temp*(peln(:,:,2:)-peln(:,:,1:))
return 
end subroutine fv_getDELZ

subroutine fv_getQ(Q, qNAME)
  real(FVPRC), intent(OUT) :: Q(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  character(LEN=*), intent(IN) :: qNAME
  integer :: isc,iec, jsc,jec, npz
  isc = FV_Atm(1)%bd%isc
  iec = FV_Atm(1)%bd%iec
  jsc = FV_Atm(1)%bd%jsc
  jec = FV_Atm(1)%bd%jec
  npz = FV_Atm(1)%npz
  if (TRIM(qNAME) == 'Q'   ) Q = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,1)
  if (TRIM(qNAME) == 'QLCN') Q = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,2)
  if (TRIM(qNAME) == 'QLLS') Q = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,3)
  if (TRIM(qNAME) == 'QICN') Q = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,4)
  if (TRIM(qNAME) == 'QILS') Q = FV_Atm(1)%q(isc:iec,jsc:jec,1:npz,5)

return
end subroutine fv_getQ

subroutine fv_getPKZ(pkz,temp,qv,pe,delz,HYDROSTATIC)
  real(REAL8), intent(OUT) ::  pkz(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent( IN) :: temp(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(FVPRC), intent( IN) ::   qv(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent( IN) ::   pe(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  real(REAL8), intent( IN) :: delz(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  logical    , intent( IN) :: HYDROSTATIC
! Local
  real(REAL8) ::   pk(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  real(REAL8) :: peln(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  real(REAL8) :: delp(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8) :: rdg
  integer :: i,j,k, isc,iec, jsc,jec, npz

  isc = FV_Atm(1)%bd%isc
  iec = FV_Atm(1)%bd%iec
  jsc = FV_Atm(1)%bd%jsc
  jec = FV_Atm(1)%bd%jec
  npz = FV_Atm(1)%npz

  rdg  = -rgas / grav
  peln = log(pe)
  pk   = exp( kappa*peln )
  delp = pe(:,:,2:npz+1)-pe(:,:,1:npz)

!-------------------------------------------------------------------------
! Re-compute the full (nonhydrostatic) pressure due to temperature changes
!-------------------------------------------------------------------------
    if ( .not.hydrostatic ) then
!$omp parallel do default(shared)
      do k=1,npz
         do j=jsc,jec
            do i=isc,iec
! perfect gas law: p = density * rdgas * virtual_temperature
!              pkz(i,j,k) = ( rdg*delp(i,j,k)*pt(i,j,k)/delz(i,j,k) )**kappa
               pkz(i,j,k) = exp( kappa*log(rdg*delp(i,j,k)*temp(i,j,k)*    &
                                      (1.d0+zvir*qv(i,j,k))/delz(i,j,k)) )
            enddo
         enddo
      enddo
    else
!$omp parallel do default(shared)
      do k=1,npz
         do j=jsc,jec
            do i=isc,iec
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k)) /     &
                   (kappa*(peln(i,j,k+1)-peln(i,j,k)))
            enddo
         enddo
      enddo
    endif

return
end subroutine fv_getPKZ


subroutine a2d3d(ua, va, ud, vd)

! Move A-Grid winds/tendencies oriented on lat/lon to the D-grid cubed-sphere orientation

! !INPUT/OUTPUT PARAMETERS:
      real(REAL8)                :: ua(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec  ,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec  ,FV_Atm(1)%npz) ! U-Wind
      real(REAL8)                :: va(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec  ,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec  ,FV_Atm(1)%npz) ! V-Wind
      real(REAL8), intent(inout) :: ud(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec  ,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec+1,FV_Atm(1)%npz) ! U-Wind
      real(REAL8), intent(inout) :: vd(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec  ,FV_Atm(1)%npz) ! V-Wind
! !Local Variables
      integer :: is ,ie , js ,je 
      integer :: npx, npy, npz
      integer :: i,j,k, im2,jm2

      real(REAL8) :: uatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%npz)
      real(REAL8) :: vatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,FV_Atm(1)%npz)

      real(REAL8) :: v3(FV_Atm(1)%bd%isc-1:FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc-1:FV_Atm(1)%bd%jec+1,3)
      real(REAL8) :: ue(FV_Atm(1)%bd%isc-1:FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc  :FV_Atm(1)%bd%jec+1,3)    ! 3D winds at edges
      real(REAL8) :: ve(FV_Atm(1)%bd%isc  :FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc-1:FV_Atm(1)%bd%jec+1,3)    ! 3D winds at edges
      real(REAL8), dimension(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec):: ut1, ut2, ut3
      real(REAL8), dimension(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec):: vt1, vt2, vt3

      npx = FV_Atm(1)%npx
      npy = FV_Atm(1)%npy
      npz = FV_Atm(1)%npz
      is  = FV_Atm(1)%bd%isc
      ie  = FV_Atm(1)%bd%iec
      js  = FV_Atm(1)%bd%jsc
      je  = FV_Atm(1)%bd%jec

      im2 = (npx-1)/2
      jm2 = (npy-1)/2

    uatemp(:,:,:) = 0.0
    vatemp(:,:,:) = 0.0

    uatemp(is:ie,js:je,:) = ua
    vatemp(is:ie,js:je,:) = va

    if (FV_Atm(1)%flagstruct%grid_type<4) then
   ! Cubed-Sphere
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, complete=.true.)
    do k=1, npz
! Compute 3D wind tendency on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,1) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,1)
             v3(i,j,2) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,2) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,2)
             v3(i,j,3) = uatemp(i,j,k)*fv_atm(1)%gridstruct%vlon(i,j,3) + vatemp(i,j,k)*fv_atm(1)%gridstruct%vlat(i,j,3)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
             ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
             ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
             ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
             ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
          enddo
       enddo

! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
        else
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_e(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_e(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
             ut2(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
             ut3(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
             ut2(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
             ut3(i) = fv_atm(1)%gridstruct%edge_vect_s(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
             ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
             ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
             ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
             ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif

! Update:
       do j=js,je+1
          do i=is,ie
             ud(i,j,k) = 0.5*( ue(i,j,1)*fv_atm(1)%gridstruct%es(1,i,j,1) +  &
                               ue(i,j,2)*fv_atm(1)%gridstruct%es(2,i,j,1) +  &
                               ue(i,j,3)*fv_atm(1)%gridstruct%es(3,i,j,1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             vd(i,j,k) = 0.5*( ve(i,j,1)*fv_atm(1)%gridstruct%ew(1,i,j,2) +  &
                               ve(i,j,2)*fv_atm(1)%gridstruct%ew(2,i,j,2) +  &
                               ve(i,j,3)*fv_atm(1)%gridstruct%ew(3,i,j,2) )
          enddo
       enddo

    enddo         ! k-loop
   else
   ! Cartesian
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    do k=1,npz
       do j=js,je+1
          do i=is,ie
             ud(i,j,k) = 0.5*( uatemp(i,j,k) + uatemp(i,j-1,k) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             vd(i,j,k) = 0.5*( vatemp(i,j,k) + vatemp(i-1,j,k) )
          enddo
       enddo
    enddo         ! k-loop
   endif

end subroutine a2d3d

subroutine a2d2d(ua, va, ud, vd)

! Move A-Grid winds/tendencies oriented on lat/lon to the D-grid cubed-sphere orientation

! !INPUT/OUTPUT PARAMETERS:
      real(REAL8)                :: ua(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec  ,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec  ) ! U-Wind
      real(REAL8)                :: va(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec  ,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec  ) ! V-Wind
      real(REAL8), intent(inout) :: ud(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec  ,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec+1) ! U-Wind
      real(REAL8), intent(inout) :: vd(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec  ) ! V-Wind
! !Local Variables
      integer :: is ,ie , js ,je
      integer :: npx, npy
      integer :: i,j, im2,jm2

      real(REAL8) :: uatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
      real(REAL8) :: vatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)

      real(REAL8) :: v3(FV_Atm(1)%bd%isc-1:FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc-1:FV_Atm(1)%bd%jec+1,3)
      real(REAL8) :: ue(FV_Atm(1)%bd%isc-1:FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc  :FV_Atm(1)%bd%jec+1,3)    ! 3D winds at edges
      real(REAL8) :: ve(FV_Atm(1)%bd%isc  :FV_Atm(1)%bd%iec+1,FV_Atm(1)%bd%jsc-1:FV_Atm(1)%bd%jec+1,3)    ! 3D winds at edges
      real(REAL8), dimension(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec):: ut1, ut2, ut3
      real(REAL8), dimension(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec):: vt1, vt2, vt3

      npx = FV_Atm(1)%npx
      npy = FV_Atm(1)%npy
      is  = FV_Atm(1)%bd%isc
      ie  = FV_Atm(1)%bd%iec
      js  = FV_Atm(1)%bd%jsc
      je  = FV_Atm(1)%bd%jec

      im2 = (npx-1)/2
      jm2 = (npy-1)/2

    uatemp(:,:) = 0.0
    vatemp(:,:) = 0.0

    uatemp(is:ie,js:je) = ua
    vatemp(is:ie,js:je) = va

    call mpp_update_domains(uatemp, FV_Atm(1)%domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)

    if (FV_Atm(1)%flagstruct%grid_type<4) then
   ! Cubed-Sphere
! Compute 3D wind tendency on A grid
       do j=js-1,je+1
          do i=is-1,ie+1
             v3(i,j,1) = uatemp(i,j)*fv_atm(1)%gridstruct%vlon(1,i,j) + vatemp(i,j)*fv_atm(1)%gridstruct%vlat(1,i,j)
             v3(i,j,2) = uatemp(i,j)*fv_atm(1)%gridstruct%vlon(2,i,j) + vatemp(i,j)*fv_atm(1)%gridstruct%vlat(2,i,j)
             v3(i,j,3) = uatemp(i,j)*fv_atm(1)%gridstruct%vlon(3,i,j) + vatemp(i,j)*fv_atm(1)%gridstruct%vlat(3,i,j)
          enddo
       enddo

! A --> D
! Interpolate to cell edges
       do j=js,je+1
          do i=is-1,ie+1
             ue(i,j,1) = v3(i,j-1,1) + v3(i,j,1)
             ue(i,j,2) = v3(i,j-1,2) + v3(i,j,2)
             ue(i,j,3) = v3(i,j-1,3) + v3(i,j,3)
          enddo
       enddo

       do j=js-1,je+1
          do i=is,ie+1
             ve(i,j,1) = v3(i-1,j,1) + v3(i,j,1)
             ve(i,j,2) = v3(i-1,j,2) + v3(i,j,2)
             ve(i,j,3) = v3(i-1,j,3) + v3(i,j,3)
          enddo
       enddo

! --- E_W edges (for v-wind):
     if ( is==1 ) then
       i = 1
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
     if ( (ie+1)==npx ) then
       i = npx
       do j=js,je
        if ( j>jm2 ) then
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j-1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
        else
             vt1(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,1)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,1)
             vt2(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,2)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,2)
             vt3(j) = fv_atm(1)%gridstruct%edge_vect_w(j)*ve(i,j+1,3)+(1.-fv_atm(1)%gridstruct%edge_vect_w(j))*ve(i,j,3)
        endif
       enddo
       do j=js,je
          ve(i,j,1) = vt1(j)
          ve(i,j,2) = vt2(j)
          ve(i,j,3) = vt3(j)
       enddo
     endif
! N-S edges (for u-wind):
     if ( js==1 ) then
       j = 1
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = (i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
             ut2(i) = (i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
             ut3(i) = (i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
        else
             ut1(i) = (i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,1)
             ut2(i) = (i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,2)
             ut3(i) = (i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_s(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif
     if ( (je+1)==npy ) then
       j = npy
       do i=is,ie
        if ( i>im2 ) then
             ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
             ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
             ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i-1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
        else
             ut1(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,1)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,1)
             ut2(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,2)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,2)
             ut3(i) = fv_atm(1)%gridstruct%edge_vect_n(i)*ue(i+1,j,3)+(1.-fv_atm(1)%gridstruct%edge_vect_n(i))*ue(i,j,3)
        endif
       enddo
       do i=is,ie
          ue(i,j,1) = ut1(i)
          ue(i,j,2) = ut2(i)
          ue(i,j,3) = ut3(i)
       enddo
     endif

! Update:
       do j=js,je+1
          do i=is,ie
             ud(i,j) = 0.5*( ue(i,j,1)*fv_atm(1)%gridstruct%es(1,i,j,1) +  &
                             ue(i,j,2)*fv_atm(1)%gridstruct%es(2,i,j,1) +  &
                             ue(i,j,3)*fv_atm(1)%gridstruct%es(3,i,j,1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             vd(i,j) = 0.5*( ve(i,j,1)*fv_atm(1)%gridstruct%ew(1,i,j,2) +  &
                             ve(i,j,2)*fv_atm(1)%gridstruct%ew(2,i,j,2) +  &
                             ve(i,j,3)*fv_atm(1)%gridstruct%ew(3,i,j,2) )
          enddo
       enddo
   else
   ! Cartesian
       do j=js,je+1
          do i=is,ie
             ud(i,j) = 0.5*( uatemp(i,j) + uatemp(i,j-1) )
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             vd(i,j) = 0.5*( vatemp(i,j) + vatemp(i-1,j) )
          enddo
       enddo
   endif

end subroutine a2d2d

subroutine fv_computeMassFluxes_r4(ucI, vcI, ple, mfx, mfy, cx, cy, dt)
  use tp_core_mod,       only: fv_tp_2d
  use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed, mp_reduce_max
  real(REAL4), intent(IN   ) ::  ucI(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL4), intent(IN   ) ::  vcI(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL4), intent(IN   ) ::  ple(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz+1)
  real(REAL4), intent(INOUT) ::  mfx(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL4), intent(INOUT) ::  mfy(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL4), intent(INOUT) ::   cx(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL4), intent(INOUT) ::   cy(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC), intent(IN   ) :: dt
  integer i,j,k

! Local ghosted arrays
  real(FVPRC) ::  uc(isd:ied+1,jsd:jed  ,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC) ::  vc(isd:ied  ,jsd:jed+1,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC) ::  ut(isd:ied+1,jsd:jed  )
  real(FVPRC) ::  vt(isd:ied  ,jsd:jed+1)

  real(FVPRC) ::  crx(is :ie +1,jsd:jed  )
  real(FVPRC) ::  cry(isd:ied  ,js :je +1)
  real(FVPRC) ::  xfx(is :ie +1,jsd:jed  )
  real(FVPRC) ::  yfx(isd:ied  ,js :je +1)
  real(FVPRC) :: ra_x(is :ie   ,jsd:jed  )
  real(FVPRC) :: ra_y(isd:ied  ,js :je   )
  real(FVPRC) ::   fx(is :ie +1,js :je   )
  real(FVPRC) ::   fy(is :ie   ,js :je +1)

  real(FVPRC) ::  delp(isd:ied,jsd:jed)

  real(FVPRC) :: wbuffer(js:je,FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: sbuffer(is:ie,FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: ebuffer(js:je,FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: nbuffer(is:ie,FV_Atm(1)%flagstruct%npz)

  real(FVPRC) :: cmax, frac
  integer     :: it, nsplt

! Fill Ghosted arrays and update halos
  uc = 0.0
  vc = 0.0
  uc(is:ie,js:je,:) = ucI
  vc(is:ie,js:je,:) = vcI
  call mpp_get_boundary(uc, vc, FV_Atm(1)%domain, &
                        wbufferx=wbuffer, ebufferx=ebuffer, &
                        sbuffery=sbuffer, nbuffery=nbuffer, &
                        gridtype=CGRID_NE, complete=.true. )
  do k=1,FV_Atm(1)%flagstruct%npz
     do j=js,je
        uc(ie+1,j,k) = ebuffer(j,k)
     enddo
     do i=is,ie
        vc(i,je+1,k) = nbuffer(i,k)
     enddo
  enddo
  call mpp_update_domains( uc, vc, FV_Atm(1)%domain, gridtype=CGRID_NE, complete=.true.)

  do k=1,FV_Atm(1)%flagstruct%npz
    ! Prepare pressures for Mass Flux calculations
     delp(is:ie,js:je) = ple(:,:,k+1)-ple(:,:,k)

     call compute_utvt(uc(isd,jsd,k), vc(isd,jsd,k), ut(isd,jsd), vt(isd,jsd), dt)
     do j=jsd,jed
        do i=is,ie+1
           xfx(i,j) = dt*ut(i,j)
        enddo
     enddo
     do j=js,je+1
        do i=isd,ied
           yfx(i,j) = dt*vt(i,j)
        enddo
     enddo
     do j=jsd,jed
        do i=is,ie+1
           if ( xfx(i,j) > 0. ) then
              crx(i,j) = xfx(i,j) * fv_atm(1)%gridstruct%rdxa(i-1,j)
              xfx(i,j) = fv_atm(1)%gridstruct%dy(i,j)*xfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i-1,j,3)
           else
              crx(i,j) = xfx(i,j) * fv_atm(1)%gridstruct%rdxa(i,j)
              xfx(i,j) = fv_atm(1)%gridstruct%dy(i,j)*xfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j,1)
           endif
        enddo
     enddo
     do j=js,je+1
        do i=isd,ied
           if ( yfx(i,j) > 0. ) then
              cry(i,j) = yfx(i,j) * fv_atm(1)%gridstruct%rdya(i,j-1)
              yfx(i,j) = fv_atm(1)%gridstruct%dx(i,j)*yfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j-1,4)
           else
              cry(i,j) = yfx(i,j) * fv_atm(1)%gridstruct%rdya(i,j)
              yfx(i,j) = fv_atm(1)%gridstruct%dx(i,j)*yfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j,2)
           endif
        enddo
     enddo

!#define SIMPLE_FV_MASS_FLUXES
#ifndef SIMPLE_FV_MASS_FLUXES
     if ( FV_Atm(1)%flagstruct%n_split==0 ) then
! Determine nsplt for tracer advection
        cmax = 0.
        do j=js,je
           do i=is,ie
              cmax = max(abs(crx(i,j))+(1.-fv_atm(1)%gridstruct%sina_u(i,j)),     &
                         abs(cry(i,j))+(1.-fv_atm(1)%gridstruct%sina_v(i,j)), cmax)
           enddo
        enddo
        call mp_reduce_max(cmax)
        nsplt = int(1.01 + cmax)
        if ( mpp_pe() == 0 )  write(6,*) k, 'Tracer_2d_split=', nsplt, cmax
     else
        nsplt = FV_Atm(1)%flagstruct%n_split
     endif
     frac  = 1. / real(nsplt)
     crx = crx*frac
     cry = cry*frac
     xfx = xfx*frac
     yfx = yfx*frac
     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = fv_atm(1)%gridstruct%area(i,j) + xfx(i,j) - xfx(i+1,j)
        enddo
     enddo
     do j=js,je
        do i=isd,ied
           ra_y(i,j) = fv_atm(1)%gridstruct%area(i,j) + yfx(i,j) - yfx(i,j+1)
        enddo
     enddo
! Zero out accumulated mass fluxes and courant numbers
      cx(:,:,k) = 0.0
     mfx(:,:,k) = 0.0
      cy(:,:,k) = 0.0
     mfy(:,:,k) = 0.0
! Compute mass fluxes using FV3 advection
     do it=1,nsplt
        call mpp_update_domains( delp, FV_Atm(1)%domain, complete=.true. )
        call fv_tp_2d(delp, crx, cry, FV_Atm(1)%flagstruct%npx, FV_Atm(1)%flagstruct%npy, &
                      FV_Atm(1)%flagstruct%hord_dp, fx, fy, xfx, yfx, &
                      FV_Atm(1)%gridstruct,FV_Atm(1)%bd, ra_x, ra_y, FV_Atm(1)%flagstruct%lim_fac)
! Update delp
        do j=js,je
           do i=is,ie
              delp(i,j) = delp(i,j) + (fx(i,j) - fx(i+1,j) +  &
                                       fy(i,j) - fy(i,j+1)) * fv_atm(1)%gridstruct%rarea(i,j)
           enddo
        enddo
! Accumulate Mass Fluxes and Courant Number outputs
         cx(:,:,k) =  cx(:,:,k) + crx(is:ie,js:je)
        mfx(:,:,k) = mfx(:,:,k) +  fx(is:ie,js:je)
         cy(:,:,k) =  cy(:,:,k) + cry(is:ie,js:je)
        mfy(:,:,k) = mfy(:,:,k) +  fy(is:ie,js:je)
     enddo
#else
! Fill DELP ghost zones
  call mpp_update_domains(delp, FV_Atm(1)%domain, complete=.true.)
! Compute Mass Fluxes and Fill Courant number outputs
     do j=js,je
        do i=is,ie
       ! X-Dir 
         if (crx(i,j) > 0.) then
            fx(i,j) = xfx(i,j) * delp(i-1,j)
         else
            fx(i,j) = xfx(i,j) * delp(i,j)
         endif
       ! Y-Dir 
         if (cry(i,j) > 0.) then
            fy(i,j) = yfx(i,j) * delp(i,j-1)
         else
            fy(i,j) = yfx(i,j) * delp(i,j)
         endif
      enddo
    enddo
! Fill Mass Fluxes and Fill Courant number outputs
      cx(:,:,k) = crx(is:ie,js:je)
     mfx(:,:,k) =  fx(is:ie,js:je)
      cy(:,:,k) = cry(is:ie,js:je)
     mfy(:,:,k) =  fy(is:ie,js:je)
#endif
  enddo

return
end subroutine fv_computeMassFluxes_r4


subroutine fv_computeMassFluxes_r8(ucI, vcI, ple, mfx, mfy, cx, cy, dt)
  use tp_core_mod,       only: fv_tp_2d
  use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed, mp_reduce_max
  real(REAL8), intent(IN   ) ::  ucI(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(IN   ) ::  vcI(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(IN   ) ::  ple(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz+1)
  real(REAL8), intent(INOUT) ::  mfx(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(INOUT) ::  mfy(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(INOUT) ::   cx(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(INOUT) ::   cy(is:ie,js:je,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC), intent(IN   ) :: dt
  integer i,j,k

! Local ghosted arrays
  real(FVPRC) ::  uc(isd:ied+1,jsd:jed  ,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC) ::  vc(isd:ied  ,jsd:jed+1,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC) ::  ut(isd:ied+1,jsd:jed  )
  real(FVPRC) ::  vt(isd:ied  ,jsd:jed+1)

  real(FVPRC) ::  crx(is :ie +1,jsd:jed  )
  real(FVPRC) ::  cry(isd:ied  ,js :je +1)
  real(FVPRC) ::  xfx(is :ie +1,jsd:jed  )
  real(FVPRC) ::  yfx(isd:ied  ,js :je +1)
  real(FVPRC) :: ra_x(is :ie   ,jsd:jed  )
  real(FVPRC) :: ra_y(isd:ied  ,js :je   )
  real(FVPRC) ::   fx(is :ie +1,js :je   )
  real(FVPRC) ::   fy(is :ie   ,js :je +1)

  real(FVPRC) ::  delp(isd:ied,jsd:jed)

  real(FVPRC) :: wbuffer(js:je,FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: sbuffer(is:ie,FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: ebuffer(js:je,FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: nbuffer(is:ie,FV_Atm(1)%flagstruct%npz)

  real(FVPRC) :: cmax, frac
  integer     :: it, nsplt

! Fill Ghosted arrays and update halos
  uc(is:ie,js:je,:) = ucI
  vc(is:ie,js:je,:) = vcI
  call mpp_get_boundary(uc, vc, FV_Atm(1)%domain, &
                        wbufferx=wbuffer, ebufferx=ebuffer, &
                        sbuffery=sbuffer, nbuffery=nbuffer, &
                        gridtype=CGRID_NE, complete=.true. )
  do k=1,FV_Atm(1)%flagstruct%npz
     do j=js,je
        uc(ie+1,j,k) = ebuffer(j,k)
     enddo
     do i=is,ie
        vc(i,je+1,k) = nbuffer(i,k)
     enddo
  enddo
  call mpp_update_domains( uc, vc, FV_Atm(1)%domain, gridtype=CGRID_NE, complete=.true.)

  do k=1,FV_Atm(1)%flagstruct%npz
    ! Prepare pressures for Mass Flux calculations
     delp(is:ie,js:je) = ple(:,:,k+1)-ple(:,:,k)

     call compute_utvt(uc(isd,jsd,k), vc(isd,jsd,k), ut(isd,jsd), vt(isd,jsd), dt)
     do j=jsd,jed
        do i=is,ie+1
           xfx(i,j) = dt*ut(i,j)
        enddo
     enddo
     do j=js,je+1
        do i=isd,ied
           yfx(i,j) = dt*vt(i,j)
        enddo
     enddo
     do j=jsd,jed
        do i=is,ie+1
           if ( xfx(i,j) > 0. ) then
              crx(i,j) = xfx(i,j) * fv_atm(1)%gridstruct%rdxa(i-1,j)
              xfx(i,j) = fv_atm(1)%gridstruct%dy(i,j)*xfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i-1,j,3)
           else
              crx(i,j) = xfx(i,j) * fv_atm(1)%gridstruct%rdxa(i,j)
              xfx(i,j) = fv_atm(1)%gridstruct%dy(i,j)*xfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j,1)
           endif
        enddo
     enddo
     do j=js,je+1
        do i=isd,ied
           if ( yfx(i,j) > 0. ) then
              cry(i,j) = yfx(i,j) * fv_atm(1)%gridstruct%rdya(i,j-1)
              yfx(i,j) = fv_atm(1)%gridstruct%dx(i,j)*yfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j-1,4)
           else
              cry(i,j) = yfx(i,j) * fv_atm(1)%gridstruct%rdya(i,j)
              yfx(i,j) = fv_atm(1)%gridstruct%dx(i,j)*yfx(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j,2)
           endif
        enddo
     enddo
!#define SIMPLE_FV_MASS_FLUXES
#ifndef SIMPLE_FV_MASS_FLUXES
     if ( FV_Atm(1)%flagstruct%n_split==0 ) then
! Determine nsplt for tracer advection
        cmax = 0.
        do j=js,je
           do i=is,ie
              cmax = max(abs(crx(i,j))+(1.-fv_atm(1)%gridstruct%sina_u(i,j)),     &
                         abs(cry(i,j))+(1.-fv_atm(1)%gridstruct%sina_v(i,j)), cmax)
           enddo
        enddo
        call mp_reduce_max(cmax)
        nsplt = int(1.01 + cmax)
        if ( mpp_pe() == 0 )  write(6,*) k, 'Tracer_2d_split=', nsplt, cmax
     else
        nsplt = FV_Atm(1)%flagstruct%n_split
     endif
     frac  = 1. / real(nsplt)
     crx = crx*frac
     cry = cry*frac
     xfx = xfx*frac
     yfx = yfx*frac
     do j=jsd,jed
        do i=is,ie
           ra_x(i,j) = fv_atm(1)%gridstruct%area(i,j) + xfx(i,j) - xfx(i+1,j)
        enddo
     enddo
     do j=js,je
        do i=isd,ied
           ra_y(i,j) = fv_atm(1)%gridstruct%area(i,j) + yfx(i,j) - yfx(i,j+1)
        enddo
     enddo
! Zero out accumulated mass fluxes and courant numbers
      cx(:,:,k) = 0.0
     mfx(:,:,k) = 0.0
      cy(:,:,k) = 0.0
     mfy(:,:,k) = 0.0
! Compute mass fluxes using FV3 advection
     do it=1,nsplt
        call mpp_update_domains( delp, FV_Atm(1)%domain, complete=.true. )
        call fv_tp_2d(delp, crx, cry, FV_Atm(1)%flagstruct%npx, FV_Atm(1)%flagstruct%npy, &
                      FV_Atm(1)%flagstruct%hord_dp, fx, fy, xfx, yfx, &
                      FV_Atm(1)%gridstruct,FV_Atm(1)%bd, ra_x, ra_y, FV_Atm(1)%flagstruct%lim_fac)
! Update delp
        do j=js,je
           do i=is,ie
              delp(i,j) = delp(i,j) + (fx(i,j) - fx(i+1,j) +  &
                                       fy(i,j) - fy(i,j+1)) * fv_atm(1)%gridstruct%rarea(i,j)
           enddo
        enddo
! Accumulate Mass Fluxes and Courant Number outputs
         cx(:,:,k) =  cx(:,:,k) + crx(is:ie,js:je)
        mfx(:,:,k) = mfx(:,:,k) +  fx(is:ie,js:je)
         cy(:,:,k) =  cy(:,:,k) + cry(is:ie,js:je)
        mfy(:,:,k) = mfy(:,:,k) +  fy(is:ie,js:je)
     enddo
#else
! Fill DELP ghost zones
  call mpp_update_domains(delp, FV_Atm(1)%domain, complete=.true.)
! Compute Mass Fluxes and Fill Courant number outputs
     do j=js,je
        do i=is,ie
       ! X-Dir 
         if (crx(i,j) > 0.) then
            fx(i,j) = xfx(i,j) * delp(i-1,j)
         else
            fx(i,j) = xfx(i,j) * delp(i,j)
         endif
       ! Y-Dir 
         if (cry(i,j) > 0.) then
            fy(i,j) = yfx(i,j) * delp(i,j-1)
         else
            fy(i,j) = yfx(i,j) * delp(i,j)
         endif
      enddo
    enddo
! Fill Mass Fluxes and Fill Courant number outputs
      cx(:,:,k) = crx(is:ie,js:je)
     mfx(:,:,k) =  fx(is:ie,js:je)
      cy(:,:,k) = cry(is:ie,js:je)
     mfy(:,:,k) =  fy(is:ie,js:je)
#endif
  enddo

return
end subroutine fv_computeMassFluxes_r8

subroutine fv_fillMassFluxes(mfx, mfy, cx, cy)
  real(REAL8), intent(OUT) :: mfx(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(OUT) :: mfy(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(OUT) ::  cx(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(OUT) ::  cy(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  integer isc,iec,jsc,jec

  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec

  mfx(:,:,:) = FV_Atm(1)%mfx(isc:iec,jsc:jec,:)
  mfy(:,:,:) = FV_Atm(1)%mfy(isc:iec,jsc:jec,:)
   cx(:,:,:) = FV_Atm(1)%cx(isc:iec,jsc:jec,:) 
   cy(:,:,:) = FV_Atm(1)%cy(isc:iec,jsc:jec,:) 

return
end subroutine fv_fillMassFluxes

subroutine fv_getVerticalMassFlux(mfx, mfy, mfz, dt)
  real(REAL8), intent(IN   ) :: mfx(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(IN   ) :: mfy(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(  OUT) :: mfz(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  real(FVPRC), intent(IN)  :: dt

  real(REAL8) :: conv(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8) :: pit(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)

  real(REAL8) :: fac

  real(REAL8) :: wbuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(REAL8) :: sbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)
  real(REAL8) :: ebuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(REAL8) :: nbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)
  real(REAL8) :: xfx(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)
  real(REAL8) :: yfx(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1,1:FV_Atm(1)%npz)
  real(REAL8) :: xfxtemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%npz)
  real(REAL8) :: yfxtemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%npz)

  integer isc,iec,jsc,jec,npz,i,j,k
  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  npz=FV_Atm(1)%npz
! Fill Ghosted arrays and update halos
  xfx=0.0d0
  yfx=0.0d0
  xfxtemp=0.0d0
  yfxtemp=0.0d0
  if (FV_Atm(1)%flagstruct%grid_type>=4) then
    xfxtemp(isc:iec,jsc:jec,:) = mfx
    yfxtemp(isc:iec,jsc:jec,:) = mfy

   ! Doubly Periodic
    call mpp_update_domains(xfxtemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(yfxtemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    xfx(isc:iec+1,jsc:jec,:) = xfxtemp(isc:iec+1,jsc:jec,:)
    yfx(isc:iec,jsc:jec+1,:) = yfxtemp(isc:iec,jsc:jec+1,:)
  else
     xfx(isc:iec,jsc:jec,:) = mfx
     yfx(isc:iec,jsc:jec,:) = mfy
     call mpp_get_boundary(xfx, yfx, FV_Atm(1)%domain, &
                           wbufferx=wbuffer, ebufferx=ebuffer, &
                           sbuffery=sbuffer, nbuffery=nbuffer, &
                           gridtype=CGRID_NE, complete=.true. )
     do k=1,npz
        do j=jsc,jec
           xfx(iec+1,j,k) = ebuffer(j,k)
        enddo
        do i=isc,iec
           yfx(i,jec+1,k) = nbuffer(i,k)
        enddo
     enddo
  end if

  fac = 1.0/(dt*MAPL_GRAV)
!
! Compute the vertical mass flux
!
!   Compute Convergence of the horizontal Mass flux
    do k=1,npz
       do j=jsc,jec
          do i=isc,iec
             conv(i,j,k) = ( xfx(i,j,k) - xfx(i+1,j,k) +  &
                             yfx(i,j,k) - yfx(i,j+1,k) ) * fac
          enddo
       enddo
    enddo
!   Surface pressure tendency
    pit(:,:) = 0.0
    do k=1,npz
       do j=jsc,jec
          do i=isc,iec
             pit(i,j) = pit(i,j) + conv(i,j,k)
          enddo
       enddo
    enddo
!   Sum over levels
    do k=2,npz
       do j=jsc,jec
          do i=isc,iec
             conv(i,j,k) = conv(i,j,k) + conv(i,j,k-1)
          enddo
       enddo
    enddo
    mfz(:,:,:) = 0.0
    do k=2,npz
       do j=jsc,jec
          do i=isc,iec
             mfz(i,j,k) = ( conv(i,j,k-1)  - FV_Atm(1)%bk(k)*pit(i,j) )/(MAPL_GRAV*fv_atm(1)%gridstruct%area(i,j))  ! Kg/m^2/s
          enddo
       enddo
    enddo

return
end subroutine fv_getVerticalMassFlux

subroutine compute_utvt(uc, vc, ut, vt, dt)
 use fv_mp_mod,         only: is,js,ie,je, isd,jsd,ied,jed
  real(FVPRC), intent(IN   ) ::  uc(isd:ied+1,jsd:jed  )
  real(FVPRC), intent(IN   ) ::  vc(isd:ied  ,jsd:jed+1)
  real(FVPRC), intent(INOUT) ::  ut(isd:ied+1,jsd:jed  )
  real(FVPRC), intent(INOUT) ::  vt(isd:ied  ,jsd:jed+1)
  real(FVPRC), intent(IN   ) :: dt
! Local vars
  real(FVPRC) :: damp
  integer i,j,npx,npy
  
  npx = FV_Atm(1)%flagstruct%npx
  npy = FV_Atm(1)%flagstruct%npy

  if ( FV_Atm(1)%flagstruct%grid_type < 3 ) then

! Center of Cube Faces
        do j=jsd,jed
           if(j/=0 .and. j/=1 .and. j/=(npy-1) .and. j/=npy) then
             do i=is-1,ie+2
                ut(i,j) = ( uc(i,j) - 0.25 * fv_atm(1)%gridstruct%cosa_u(i,j) *     &
                    (vc(i-1,j)+vc(i,j)+vc(i-1,j+1)+vc(i,j+1)))*fv_atm(1)%gridstruct%rsin_u(i,j)
             enddo
           endif
        enddo
        do j=js-1,je+2
           if( j/=1 .and. j/=npy ) then
              do i=isd,ied
                 vt(i,j) = ( vc(i,j) - 0.25 * fv_atm(1)%gridstruct%cosa_v(i,j) *     &
                    (uc(i,j-1)+uc(i+1,j-1)+uc(i,j)+uc(i+1,j)))*fv_atm(1)%gridstruct%rsin_v(i,j)
              enddo
           endif
        enddo

! West edge:
       if ( is==1 ) then
          do j=jsd,jed
             if ( uc(1,j)*dt > 0. ) then
                ut(1,j) = uc(1,j) / fv_atm(1)%gridstruct%sin_sg(0,j,3)
             else
                ut(1,j) = uc(1,j) / fv_atm(1)%gridstruct%sin_sg(1,j,1)
             endif
          enddo
          do j=max(3,js), min(npy-2,je+1)
             vt(0,j) = vc(0,j) - 0.25*fv_atm(1)%gridstruct%cosa_v(0,j)*   &
                  (ut(0,j-1)+ut(1,j-1)+ut(0,j)+ut(1,j))
             vt(1,j) = vc(1,j) - 0.25*fv_atm(1)%gridstruct%cosa_v(1,j)*   &
                  (ut(1,j-1)+ut(2,j-1)+ut(1,j)+ut(2,j))
          enddo
       endif   ! West face

! East edge:
       if ( (ie+1)==npx ) then
          do j=jsd,jed
             if ( uc(npx,j)*dt > 0. ) then
                ut(npx,j) = uc(npx,j) / fv_atm(1)%gridstruct%sin_sg(npx-1,j,3)
             else
                ut(npx,j) = uc(npx,j) / fv_atm(1)%gridstruct%sin_sg(npx,j,1)
             endif
          enddo

           do j=max(3,js), min(npy-2,je+1)
              vt(npx-1,j) = vc(npx-1,j) - 0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,j)*   &
                           (ut(npx-1,j-1)+ut(npx,j-1)+ut(npx-1,j)+ut(npx,j))
              vt(npx,j) = vc(npx,j) - 0.25*fv_atm(1)%gridstruct%cosa_v(npx,j)*   &
                         (ut(npx,j-1)+ut(npx+1,j-1)+ut(npx,j)+ut(npx+1,j))
           enddo
       endif

! South (Bottom) edge:
       if ( js==1 ) then

           do i=isd,ied
              if ( vc(i,1)*dt > 0. ) then
                   vt(i,1) = vc(i,1) / fv_atm(1)%gridstruct%sin_sg(i,0,4)
              else
                   vt(i,1) = vc(i,1) / fv_atm(1)%gridstruct%sin_sg(i,1,2)
              endif
           enddo

           do i=max(3,is),min(npx-2,ie+1)
              ut(i,0) = uc(i,0) - 0.25*fv_atm(1)%gridstruct%cosa_u(i,0)*   &
                       (vt(i-1,0)+vt(i,0)+vt(i-1,1)+vt(i,1))
              ut(i,1) = uc(i,1) - 0.25*fv_atm(1)%gridstruct%cosa_u(i,1)*   &
                       (vt(i-1,1)+vt(i,1)+vt(i-1,2)+vt(i,2))
           enddo
       endif

! North edge:
       if ( (je+1)==npy ) then
           do i=isd,ied
              if ( vc(i,npy)*dt > 0. ) then
                   vt(i,npy) = vc(i,npy) / fv_atm(1)%gridstruct%sin_sg(i,npy-1,4)
              else
                   vt(i,npy) = vc(i,npy) / fv_atm(1)%gridstruct%sin_sg(i,npy,2)
              endif
           enddo
           do i=max(3,is),min(npx-2,ie+1)
              ut(i,npy-1) = uc(i,npy-1) - 0.25*fv_atm(1)%gridstruct%cosa_u(i,npy-1)*   &
                           (vt(i-1,npy-1)+vt(i,npy-1)+vt(i-1,npy)+vt(i,npy))
              ut(i,npy) = uc(i,npy) - 0.25*fv_atm(1)%gridstruct%cosa_u(i,npy)*   &
                         (vt(i-1,npy)+vt(i,npy)+vt(i-1,npy+1)+vt(i,npy+1))
           enddo
       endif

!The following code solves a 2x2 system to get the interior parallel-to-edge 
! uc,vc values near the corners (ex: for the sw corner ut(2,1) and vt(1,2) are solved for simultaneously).
! It then computes the halo uc, vc values so as to be consistent with the computations on the facing panel.

       !The system solved is:
       !  ut(2,1) = uc(2,1) - avg(vt)*fv_atm(1)%gridstruct%cosa_u(2,1)
       !  vt(1,2) = vc(1,2) - avg(ut)*fv_atm(1)%gridstruct%cosa_v(1,2)
       ! in which avg(vt) includes vt(1,2) and avg(ut) includes ut(2,1)

        if( fv_atm(1)%gridstruct%sw_corner ) then
            damp = 1. / (1.-0.0625*fv_atm(1)%gridstruct%cosa_u(2,0)*fv_atm(1)%gridstruct%cosa_v(1,0))
            ut(2,0) = (uc(2,0)-0.25*fv_atm(1)%gridstruct%cosa_u(2,0)*(vt(1,1)+vt(2,1)+vt(2,0) +vc(1,0) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_v(1,0)*(ut(1,0)+ut(1,-1)+ut(2,-1))) ) * damp
            damp = 1. / (1.-0.0625*fv_atm(1)%gridstruct%cosa_u(0,1)*fv_atm(1)%gridstruct%cosa_v(0,2))
            vt(0,2) = (vc(0,2)-0.25*fv_atm(1)%gridstruct%cosa_v(0,2)*(ut(1,1)+ut(1,2)+ut(0,2)+uc(0,1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_u(0,1)*(vt(0,1)+vt(-1,1)+vt(-1,2))) ) * damp

            damp = 1. / (1.-0.0625*fv_atm(1)%gridstruct%cosa_u(2,1)*fv_atm(1)%gridstruct%cosa_v(1,2))
            ut(2,1) = (uc(2,1)-0.25*fv_atm(1)%gridstruct%cosa_u(2,1)*(vt(1,1)+vt(2,1)+vt(2,2)+vc(1,2) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_v(1,2)*(ut(1,1)+ut(1,2)+ut(2,2))) ) * damp

            vt(1,2) = (vc(1,2)-0.25*fv_atm(1)%gridstruct%cosa_v(1,2)*(ut(1,1)+ut(1,2)+ut(2,2)+uc(2,1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_u(2,1)*(vt(1,1)+vt(2,1)+vt(2,2))) ) * damp
        endif

        if( fv_atm(1)%gridstruct%se_corner ) then
            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(npx-1,0)*fv_atm(1)%gridstruct%cosa_v(npx-1,0))
            ut(npx-1,0) = ( uc(npx-1,0)-0.25*fv_atm(1)%gridstruct%cosa_u(npx-1,0)*(   &
                            vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,0)+vc(npx-1,0) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,0)*(ut(npx,0)+ut(npx,-1)+ut(npx-1,-1))) ) * damp
            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(npx+1,1)*fv_atm(1)%gridstruct%cosa_v(npx,2))
            vt(npx,  2) = ( vc(npx,2)-0.25*fv_atm(1)%gridstruct%cosa_v(npx,2)*(  &
                            ut(npx,1)+ut(npx,2)+ut(npx+1,2)+uc(npx+1,1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_u(npx+1,1)*(vt(npx,1)+vt(npx+1,1)+vt(npx+1,2))) ) * damp

            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(npx-1,1)*fv_atm(1)%gridstruct%cosa_v(npx-1,2))
            ut(npx-1,1) = ( uc(npx-1,1)-0.25*fv_atm(1)%gridstruct%cosa_u(npx-1,1)*(  &
                            vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,2)+vc(npx-1,2) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,2)*(ut(npx,1)+ut(npx,2)+ut(npx-1,2))) ) * damp
            vt(npx-1,2) = ( vc(npx-1,2)-0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,2)*(  &
                            ut(npx,1)+ut(npx,2)+ut(npx-1,2)+uc(npx-1,1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_u(npx-1,1)*(vt(npx-1,1)+vt(npx-2,1)+vt(npx-2,2))) ) * damp
        endif

        if( fv_atm(1)%gridstruct%ne_corner ) then
            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(npx-1,npy)*fv_atm(1)%gridstruct%cosa_v(npx-1,npy+1))
            ut(npx-1,npy) = ( uc(npx-1,npy)-0.25*fv_atm(1)%gridstruct%cosa_u(npx-1,npy)*(   &
                              vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy+1)+vc(npx-1,npy+1) -   &
                0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,npy+1)*(ut(npx,npy)+ut(npx,npy+1)+ut(npx-1,npy+1))) ) * damp
            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(npx+1,npy-1)*fv_atm(1)%gridstruct%cosa_v(npx,npy-1))
            vt(npx,  npy-1) = ( vc(npx,npy-1)-0.25*fv_atm(1)%gridstruct%cosa_v(npx,npy-1)*(   &
                                ut(npx,npy-1)+ut(npx,npy-2)+ut(npx+1,npy-2)+uc(npx+1,npy-1) -   &
                0.25*fv_atm(1)%gridstruct%cosa_u(npx+1,npy-1)*(vt(npx,npy)+vt(npx+1,npy)+vt(npx+1,npy-1))) ) * damp

            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(npx-1,npy-1)*fv_atm(1)%gridstruct%cosa_v(npx-1,npy-1))
            ut(npx-1,npy-1) = ( uc(npx-1,npy-1)-0.25*fv_atm(1)%gridstruct%cosa_u(npx-1,npy-1)*(  &
                                vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy-1)+vc(npx-1,npy-1) -  &
                0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,npy-1)*(ut(npx,npy-1)+ut(npx,npy-2)+ut(npx-1,npy-2))) ) * damp
            vt(npx-1,npy-1) = ( vc(npx-1,npy-1)-0.25*fv_atm(1)%gridstruct%cosa_v(npx-1,npy-1)*(  &
                                ut(npx,npy-1)+ut(npx,npy-2)+ut(npx-1,npy-2)+uc(npx-1,npy-1) -  &
                0.25*fv_atm(1)%gridstruct%cosa_u(npx-1,npy-1)*(vt(npx-1,npy)+vt(npx-2,npy)+vt(npx-2,npy-1))) ) * damp
        endif

        if( fv_atm(1)%gridstruct%nw_corner ) then
            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(2,npy)*fv_atm(1)%gridstruct%cosa_v(1,npy+1))
            ut(2,npy) = ( uc(2,npy)-0.25*fv_atm(1)%gridstruct%cosa_u(2,npy)*(   &
                          vt(1,npy)+vt(2,npy)+vt(2,npy+1)+vc(1,npy+1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_v(1,npy+1)*(ut(1,npy)+ut(1,npy+1)+ut(2,npy+1))) ) * damp
            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(0,npy-1)*fv_atm(1)%gridstruct%cosa_v(0,npy-1))
            vt(0,npy-1) = ( vc(0,npy-1)-0.25*fv_atm(1)%gridstruct%cosa_v(0,npy-1)*(  &
                            ut(1,npy-1)+ut(1,npy-2)+ut(0,npy-2)+uc(0,npy-1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_u(0,npy-1)*(vt(0,npy)+vt(-1,npy)+vt(-1,npy-1))) ) * damp

            damp = 1. / (1. - 0.0625*fv_atm(1)%gridstruct%cosa_u(2,npy-1)*fv_atm(1)%gridstruct%cosa_v(1,npy-1))
            ut(2,npy-1) = ( uc(2,npy-1)-0.25*fv_atm(1)%gridstruct%cosa_u(2,npy-1)*(  &
                            vt(1,npy)+vt(2,npy)+vt(2,npy-1)+vc(1,npy-1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_v(1,npy-1)*(ut(1,npy-1)+ut(1,npy-2)+ut(2,npy-2))) ) * damp

            vt(1,npy-1) = ( vc(1,npy-1)-0.25*fv_atm(1)%gridstruct%cosa_v(1,npy-1)*(  &
                            ut(1,npy-1)+ut(1,npy-2)+ut(2,npy-2)+uc(2,npy-1) -   &
                      0.25*fv_atm(1)%gridstruct%cosa_u(2,npy-1)*(vt(1,npy)+vt(2,npy)+vt(2,npy-1))) ) * damp
        endif
           
 else
! grid_type >= 3

        do j=jsd,jed
           do i=is,ie+1
              ut(i,j) =  uc(i,j)
           enddo
        enddo

        do j=js,je+1
           do i=isd,ied
              vt(i,j) = vc(i,j)
           enddo
        enddo
 endif

end subroutine compute_utvt

subroutine fv_getOmega(omga)
  real(REAL8), intent(OUT) :: omga(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  integer isc,iec,jsc,jec

  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec

  omga(:,:,:) = FV_Atm(1)%omga(isc:iec,jsc:jec,:)

return
end subroutine fv_getOmega

subroutine fv_getPK(pkxyz)
  real(REAL8), intent(OUT) :: pkxyz(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz+1)
  integer isc,iec,jsc,jec
  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  pkxyz(:,:,:) = FV_Atm(1)%pk(isc:iec,jsc:jec,:)
  return
end subroutine fv_getPK

subroutine fv_getVorticity(u, v, vort)
  real(REAL8), intent(IN)  ::     u(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(IN)  ::     v(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(FVPRC), intent(OUT) ::  vort(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)

  real(REAL8) ::  utemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1,1:FV_Atm(1)%npz)
  real(REAL8) ::  vtemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)

  real(REAL8) :: uatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%npz)
  real(REAL8) :: vatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%npz)

  real(REAL8) :: wbuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(REAL8) :: sbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)
  real(REAL8) :: ebuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(REAL8) :: nbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)

  integer isc,iec,jsc,jec
  integer npz
  integer i,j,k

  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  npz = FV_Atm(1)%npz

  utemp  = 0d0
  vtemp  = 0d0
  if (FV_Atm(1)%flagstruct%grid_type>=4) then
  ! Doubly Periodic
    uatemp(isc:iec,jsc:jec,:) = u
    vatemp(isc:iec,jsc:jec,:) = v
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    utemp(isc:iec,jsc:jec+1,:) = uatemp(isc:iec,jsc:jec+1,:)
    vtemp(isc:iec+1,jsc:jec,:) = vatemp(isc:iec+1,jsc:jec,:)
  else
    utemp(isc:iec,jsc:jec,:) = u
    vtemp(isc:iec,jsc:jec,:) = v
    call mpp_get_boundary(utemp, vtemp, FV_Atm(1)%domain, &
                          wbuffery=wbuffer, ebuffery=ebuffer, &
                          sbufferx=sbuffer, nbufferx=nbuffer, &
                          gridtype=DGRID_NE, complete=.true. )
    do k=1,npz
       do i=isc,iec
          utemp(i,jec+1,k) = nbuffer(i,k)
       enddo
       do j=jsc,jec
          vtemp(iec+1,j,k) = ebuffer(j,k)
       enddo
    enddo
  endif
! Calc Vorticity
    do k=1,npz
       do j=jsc,jec
          do i=isc,iec
             vort(i,j,k) = fv_atm(1)%gridstruct%rarea(i,j)*(utemp(i,j,k)*fv_atm(1)%gridstruct%dx(i,j)-utemp(i,j+1,k)*fv_atm(1)%gridstruct%dx(i,j+1) - &
                                       vtemp(i,j,k)*fv_atm(1)%gridstruct%dy(i,j)+vtemp(i+1,j,k)*fv_atm(1)%gridstruct%dy(i+1,j))
         enddo
       enddo
    enddo
end subroutine fv_getVorticity

subroutine fv_getDivergence(uc, vc, divg)
  real(REAL8), intent(IN)  ::    uc(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(IN)  ::    vc(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(FVPRC), intent(OUT) ::  divg(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)

  real(FVPRC) :: uctemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)
  real(FVPRC) :: vctemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1,1:FV_Atm(1)%npz)

  real(FVPRC) :: uatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%npz)
  real(FVPRC) :: vatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%npz)

  real(FVPRC) :: ut(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  )
  real(FVPRC) :: vt(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1)

  real(FVPRC) :: wbuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(FVPRC) :: sbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)
  real(FVPRC) :: ebuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(FVPRC) :: nbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)

  integer isd,ied,jsd,jed
  integer isc,iec,jsc,jec
  integer npz
  integer i,j,k
  real(FVPRC) :: dt

  isd=FV_Atm(1)%bd%isd ; ied=FV_Atm(1)%bd%ied
  jsd=FV_Atm(1)%bd%jsd ; jed=FV_Atm(1)%bd%jed
  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  npz = FV_Atm(1)%npz

  uctemp  = 0d0
  vctemp  = 0d0
  if (FV_Atm(1)%flagstruct%grid_type>=4) then
  ! Doubly Periodic
    uatemp(isc:iec,jsc:jec,:) = uc
    vatemp(isc:iec,jsc:jec,:) = vc
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    uctemp(isc:iec+1,jsc:jec,:) = uatemp(isc:iec+1,jsc:jec,:)
    vctemp(isc:iec,jsc:jec+1,:) = vatemp(isc:iec,jsc:jec+1,:)
  else
    uctemp(isc:iec,jsc:jec,:) = uc
    vctemp(isc:iec,jsc:jec,:) = vc
    call mpp_get_boundary(uctemp, vctemp, FV_Atm(1)%domain, &
                          wbufferx=wbuffer, ebufferx=ebuffer, &
                          sbuffery=sbuffer, nbuffery=nbuffer, &
                          gridtype=CGRID_NE, complete=.true.)
    do k=1,npz
       do j=jsc,jec
          uctemp(iec+1,j,k) = ebuffer(j,k)
       enddo
       do i=isc,iec
          vctemp(i,jec+1,k) = nbuffer(i,k)
       enddo
    enddo
  endif
  call mpp_update_domains( uctemp, vctemp, FV_Atm(1)%domain, gridtype=CGRID_NE, complete=.true.)
! Calc Divergence
    dt = 1.0
    do k=1,npz
        call compute_utvt(uctemp(isd,jsd,k), vctemp(isd,jsd,k), ut(isd,jsd), vt(isd,jsd), dt)
        do j=jsc,jec
           do i=isc,iec+1
              if ( ut(i,j) > 0. ) then
                   ut(i,j) = fv_atm(1)%gridstruct%dy(i,j)*ut(i,j)*fv_atm(1)%gridstruct%sin_sg(i-1,j,3)
              else
                   ut(i,j) = fv_atm(1)%gridstruct%dy(i,j)*ut(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j,1)
             endif
           enddo
        enddo
        do j=jsc,jec+1
           do i=isc,iec
              if ( vt(i,j) > 0. ) then
                   vt(i,j) = fv_atm(1)%gridstruct%dx(i,j)*vt(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j-1,4)
              else
                   vt(i,j) = fv_atm(1)%gridstruct%dx(i,j)*vt(i,j)*fv_atm(1)%gridstruct%sin_sg(i,j,2)
              endif
           enddo
        enddo
        do j=jsc,jec
           do i=isc,iec
              divg(i,j,k) = fv_atm(1)%gridstruct%rarea(i,j)*( ut(i+1,j)-ut(i,j) + &
                                         vt(i,j+1)-vt(i,j) )
           enddo
        enddo
    enddo
end subroutine fv_getDivergence

subroutine fv_getUpdraftHelicity(uh25)
   use constants_mod, only: fms_grav=>grav
   use fv_diagnostics_mod, only: get_vorticity, updraft_helicity
   real(FVPRC), intent(OUT) :: uh25(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
   integer :: sphum=1
   real(FVPRC) :: vort(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
   call get_vorticity(FV_Atm(1)%bd%isc, FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc, FV_Atm(1)%bd%jec, &
                      FV_Atm(1)%bd%isd, FV_Atm(1)%bd%ied, FV_Atm(1)%bd%jsd, FV_Atm(1)%bd%jed, &
                      FV_Atm(1)%npz, FV_Atm(1)%u, FV_Atm(1)%v, vort, &
                      FV_Atm(1)%gridstruct%dx, FV_Atm(1)%gridstruct%dy, FV_Atm(1)%gridstruct%rarea)
   call updraft_helicity(FV_Atm(1)%bd%isc, FV_Atm(1)%bd%iec, FV_Atm(1)%bd%jsc, FV_Atm(1)%bd%jec, FV_Atm(1)%ng, FV_Atm(1)%npz, &
                     zvir, sphum, uh25, &
                     FV_Atm(1)%w, vort, FV_Atm(1)%delz, FV_Atm(1)%q,   &
                     FV_Atm(1)%flagstruct%hydrostatic, FV_Atm(1)%pt, FV_Atm(1)%peln, FV_Atm(1)%phis, fms_grav, 2.e3, 5.e3)
end subroutine fv_getUpdraftHelicity

subroutine fv_getEPV(pt, vort, ua, va, epv)
  real(REAL8), intent(IN)  ::    pt(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC), intent(IN)  ::  vort(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(IN)  ::    ua(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(IN)  ::    va(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz)
  real(REAL8), intent(OUT) ::   epv(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz)

  real(FVPRC) :: dz_g(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%flagstruct%npz)

  real(FVPRC) :: pt_g(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed,1:FV_Atm(1)%flagstruct%npz)
  real(FVPRC) :: pt_e(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz+1)

  real(FVPRC) :: ua_e(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz+1)
  real(FVPRC) :: va_e(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%flagstruct%npz+1)

  real(FVPRC) :: vtcorr
  real(FVPRC) :: pt_l, pt_r, pt_b, pt_t, dptdx, dptdy, dptdp
  real(FVPRC) :: w_l, w_r, w_b, w_t, dwdx, dwdy, dudp, dvdp

  real(FVPRC) :: dudz, dvdz, dptdz, area_dxdz, area_dydz, area_dxdy
  real(FVPRC) :: w_im1, w_jm1, w_ip1, w_jp1
  real(FVPRC) :: pt_im1, pt_jm1, pt_ip1, pt_jp1
  real(FVPRC) :: dz_im1, dz_jm1, dz_ip1, dz_jp1
  real(FVPRC) :: vort_i, vort_j, vort_k

  integer isc,iec,jsc,jec
  integer npz
  integer i,j,k
  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  npz = FV_Atm(1)%npz

   pt_g(isc:iec,jsc:jec,:) = pt(isc:iec,jsc:jec,:)
   call mpp_update_domains(pt_g, FV_Atm(1)%domain, complete=.true.)
! Get PT/UA/VA at layer edges
   do j=jsc,jec
      call ppme(pt(isc:iec,j,:),pt_e(isc:iec,j,:),FV_Atm(1)%delp(isc:iec,j,:),iec-isc+1,npz)
      call ppme(ua(isc:iec,j,:),ua_e(isc:iec,j,:),FV_Atm(1)%delp(isc:iec,j,:),iec-isc+1,npz)
      call ppme(va(isc:iec,j,:),va_e(isc:iec,j,:),FV_Atm(1)%delp(isc:iec,j,:),iec-isc+1,npz)
   enddo

   if (.not. FV_HYDROSTATIC) then
      dz_g(isc:iec,jsc:jec,:) = FV_Atm(1)%delz(isc:iec,jsc:jec,:)
      call mpp_update_domains(dz_g, FV_Atm(1)%domain, complete=.false.)
      call mpp_update_domains(FV_Atm(1)%w,  FV_Atm(1)%domain, complete=.true.)
      do k=1,npz
        do j=jsc,jec
          do i=isc,iec
! get all dz and areas needed
             dz_im1 = 0.5*(dz_g(i,j,k) + dz_g(i-1,j,k))
             dz_ip1 = 0.5*(dz_g(i,j,k) + dz_g(i+1,j,k))
             dz_jm1 = 0.5*(dz_g(i,j,k) + dz_g(i,j-1,k))
             dz_jp1 = 0.5*(dz_g(i,j,k) + dz_g(i,j+1,k))
             area_dxdz = fv_atm(1)%gridstruct%dxa(i,j) * 0.5 * ( dz_im1 + dz_ip1 )
             area_dydz = fv_atm(1)%gridstruct%dya(i,j) * 0.5 * ( dz_jm1 + dz_jp1 )
             area_dxdy = fv_atm(1)%gridstruct%area(i,j)
! Get W on center of horizontal-cell edges
             w_im1 = (fv_atm(1)%gridstruct%dxa(i-1,j)*FV_Atm(1)%w(i-1,j,k) + fv_atm(1)%gridstruct%dxa(i  ,j)*FV_Atm(1)%w(i  ,j,k))/(fv_atm(1)%gridstruct%dxa(i-1,j)+fv_atm(1)%gridstruct%dxa(i  ,j))
             w_ip1 = (fv_atm(1)%gridstruct%dxa(i  ,j)*FV_Atm(1)%w(i  ,j,k) + fv_atm(1)%gridstruct%dxa(i+1,j)*FV_Atm(1)%w(i+1,j,k))/(fv_atm(1)%gridstruct%dxa(i  ,j)+fv_atm(1)%gridstruct%dxa(i+1,j))
             w_jm1 = (fv_atm(1)%gridstruct%dya(i,j-1)*FV_Atm(1)%w(i,j-1,k) + fv_atm(1)%gridstruct%dya(i,j  )*FV_Atm(1)%w(i,j  ,k))/(fv_atm(1)%gridstruct%dya(i,j-1)+fv_atm(1)%gridstruct%dya(i,j  ))
             w_jp1 = (fv_atm(1)%gridstruct%dya(i,j  )*FV_Atm(1)%w(i,j  ,k) + fv_atm(1)%gridstruct%dya(i,j+1)*FV_Atm(1)%w(i,j+1,k))/(fv_atm(1)%gridstruct%dya(i,j  )+fv_atm(1)%gridstruct%dya(i,j+1))
! Get PT on center of horizontal-cell edges
             pt_im1 = (fv_atm(1)%gridstruct%dxa(i-1,j)*pt_g(i-1,j,k) + fv_atm(1)%gridstruct%dxa(i  ,j)*pt_g(i  ,j,k))/(fv_atm(1)%gridstruct%dxa(i-1,j)+fv_atm(1)%gridstruct%dxa(i  ,j))
             pt_ip1 = (fv_atm(1)%gridstruct%dxa(i  ,j)*pt_g(i  ,j,k) + fv_atm(1)%gridstruct%dxa(i+1,j)*pt_g(i+1,j,k))/(fv_atm(1)%gridstruct%dxa(i  ,j)+fv_atm(1)%gridstruct%dxa(i+1,j))
             pt_jm1 = (fv_atm(1)%gridstruct%dya(i,j-1)*pt_g(i,j-1,k) + fv_atm(1)%gridstruct%dya(i,j  )*pt_g(i,j  ,k))/(fv_atm(1)%gridstruct%dya(i,j-1)+fv_atm(1)%gridstruct%dya(i,j  ))
             pt_jp1 = (fv_atm(1)%gridstruct%dya(i,j  )*pt_g(i,j  ,k) + fv_atm(1)%gridstruct%dya(i,j+1)*pt_g(i,j+1,k))/(fv_atm(1)%gridstruct%dya(i,j  )+fv_atm(1)%gridstruct%dya(i,j+1))
! d(pt)/dx
             dptdx = (pt_ip1      - pt_im1       )/fv_atm(1)%gridstruct%dxa(i,j)
! d(pt)/dy
             dptdy = (pt_jp1      - pt_jm1       )/fv_atm(1)%gridstruct%dya(i,j)
! d(pt)/dz
             dptdz = (pt_e(i,j,k) - pt_e(i,j,k+1))/dz_g(i,j,k)
! i-component of EPV
             vort_i = (1./area_dydz) * ( va_e(i,j,k+1)*fv_atm(1)%gridstruct%dya(i,j) - va_e(i,j,k)*fv_atm(1)%gridstruct%dya(i,j) - &
                                                 w_jm1*dz_jm1   +       w_jp1*dz_jp1 )
! j-component of EPV 
             vort_j = (1./area_dxdz) * ( ua_e(i,j,k)*fv_atm(1)%gridstruct%dxa(i,j)   - ua_e(i,j,k+1)*fv_atm(1)%gridstruct%dxa(i,j) - &
                                                 w_ip1*dz_ip1   +       w_im1*dz_im1 ) + &
                                         2.*MAPL_OMEGA*COS(FV_Atm(1)%gridstruct%agrid(i,j,2))
! k-component of EPV (pre-computed and passed in to this subroutine)
             vort_k = (vort(i,j,k)+FV_Atm(1)%gridstruct%f0(i,j))
           ! vort_k = (1./area_dxdy) * (FV_Atm(1)%u(i,j,k)*dx(i,j)-FV_Atm(1)%u(i,j+1,k)*dx(i,j+1)  - &
           !                            FV_Atm(1)%v(i,j,k)*dy(i,j)+FV_Atm(1)%v(i+1,j,k)*dy(i+1,j)) + &
           !                            f0(i,j)
! Complete full non-hydrostatic EPV calculation
             epv(i,j,k) = MAPL_GRAV * ( FV_Atm(1)%delz(i,j,k)/FV_Atm(1)%delp(i,j,k) ) * &
                                      ( vort_i*dptdx + vort_j*dptdy + vort_k*dptdz )
          enddo
        enddo
      enddo
   else
      do k=1,npz
        do j=jsc,jec
          do i=isc,iec
! d(u)/dp
             dudp = (ua_e(i,j,k)-ua_e(i,j,k+1))/FV_Atm(1)%delp(i,j,k)
! d(v)/dp
             dvdp = (va_e(i,j,k)-va_e(i,j,k+1))/FV_Atm(1)%delp(i,j,k)
! d(pt)/dx
             pt_l = (fv_atm(1)%gridstruct%dxa(i-1,j)*pt_g(i-1,j,k) + fv_atm(1)%gridstruct%dxa(i  ,j)*pt_g(i  ,j,k))/(fv_atm(1)%gridstruct%dxa(i-1,j)+fv_atm(1)%gridstruct%dxa(i  ,j))
             pt_r = (fv_atm(1)%gridstruct%dxa(i  ,j)*pt_g(i  ,j,k) + fv_atm(1)%gridstruct%dxa(i+1,j)*pt_g(i+1,j,k))/(fv_atm(1)%gridstruct%dxa(i  ,j)+fv_atm(1)%gridstruct%dxa(i+1,j))
             dptdx = (pt_r-pt_l)/fv_atm(1)%gridstruct%dxa(i,j)
! d(pt)/dy
             pt_b = (fv_atm(1)%gridstruct%dya(i,j-1)*pt_g(i,j-1,k) + fv_atm(1)%gridstruct%dya(i,j  )*pt_g(i,j  ,k))/(fv_atm(1)%gridstruct%dya(i,j-1)+fv_atm(1)%gridstruct%dya(i,j  ))
             pt_t = (fv_atm(1)%gridstruct%dya(i,j  )*pt_g(i,j  ,k) + fv_atm(1)%gridstruct%dya(i,j+1)*pt_g(i,j+1,k))/(fv_atm(1)%gridstruct%dya(i,j  )+fv_atm(1)%gridstruct%dya(i,j+1))
             dptdy = (pt_t-pt_b)/fv_atm(1)%gridstruct%dya(i,j)
! d(pt)/dp
             dptdp = (pt_e(i,j,k)-pt_e(i,j,k+1))/FV_Atm(1)%delp(i,j,k)
! Vorticity Correction to account for changes in theta along isobaric surfaces
             vtcorr = -dvdp*dptdx + dudp*dptdy + (vort(i,j,k)+FV_Atm(1)%gridstruct%f0(i,j))*dptdp
             epv(i,j,k) = MAPL_GRAV*vtcorr
          enddo
        enddo
      enddo
   endif

end subroutine fv_getEPV

!------------------------------------------------------------------------------
!BOP         
!
! !IROUTINE: fv_getAgridWinds_3D
!
! !INTERFACE:
!    
subroutine fv_getAgridWinds_3D(u, v, ua, va, uc, vc, rotate)

! !INPUT PARAMETERS:
  real(REAL8), intent(IN)  ::  u(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), intent(IN)  ::  v(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  logical, optional, intent(IN) :: rotate
! 
! !OUTPUT PARAMETERS:
  real(REAL8),           intent(OUT) :: ua(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8),           intent(OUT) :: va(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), optional, intent(OUT) :: uc(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
  real(REAL8), optional, intent(OUT) :: vc(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,1:FV_Atm(1)%npz)
!
! !DESCRIPTION:
! 
! !LOCAL VARIABLES:
  real(FVPRC) :: wbuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(FVPRC) :: sbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)
  real(FVPRC) :: ebuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec,FV_Atm(1)%npz)
  real(FVPRC) :: nbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%npz)

  integer isc,iec,jsc,jec
  integer isd,ied,jsd,jed
  integer npz
  integer i,j,k
  
  real(FVPRC) :: ut(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied, FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
  real(FVPRC) :: vt(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied, FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)

  real(FVPRC) ::  utemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1,1:FV_Atm(1)%npz)
  real(FVPRC) ::  vtemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)
  real(FVPRC) :: uatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)
  real(FVPRC) :: vatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)
  real(FVPRC) :: uctemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed  ,1:FV_Atm(1)%npz)
  real(FVPRC) :: vctemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied  ,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1,1:FV_Atm(1)%npz)
!EOP
!------------------------------------------------------------------------------
!BOC
  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  isd=FV_Atm(1)%bd%isd ; ied=FV_Atm(1)%bd%ied
  jsd=FV_Atm(1)%bd%jsd ; jed=FV_Atm(1)%bd%jed
  npz = FV_Atm(1)%npz
  
  utemp  = 0
  vtemp  = 0
  uatemp = 0
  vatemp = 0
  uctemp = 0
  vctemp = 0

  if (FV_Atm(1)%flagstruct%grid_type>=4) then
  ! Doubly Periodic
    uatemp(isc:iec,jsc:jec,:) = u
    vatemp(isc:iec,jsc:jec,:) = v
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    utemp(isc:iec,jsc:jec+1,:) = uatemp(isc:iec,jsc:jec+1,:)
    vtemp(isc:iec+1,jsc:jec,:) = vatemp(isc:iec+1,jsc:jec,:)
  else
    utemp(isc:iec,jsc:jec,:) = u
    vtemp(isc:iec,jsc:jec,:) = v
  ! update shared edges
    call mpp_get_boundary(utemp, vtemp, FV_Atm(1)%domain, &
                          wbuffery=wbuffer, ebuffery=ebuffer, &
                          sbufferx=sbuffer, nbufferx=nbuffer, &
                          gridtype=DGRID_NE, complete=.true. )
    do k=1,npz
       do i=isc,iec
          utemp(i,jec+1,k) = nbuffer(i,k)
       enddo  
       do j=jsc,jec
          vtemp(iec+1,j,k) = ebuffer(j,k)
       enddo  
    enddo   
  endif

  call mpp_update_domains(utemp, vtemp, FV_Atm(1)%domain, gridtype=DGRID_NE, complete=.true.)
  do k=1,npz
   call d2a2c_vect(utemp(:,:,k),  vtemp(:,:,k), &
                   uatemp(:,:,k), vatemp(:,:,k), &
                   uctemp(:,:,k), vctemp(:,:,k), ut, vt, .true., &
                   FV_Atm(1)%gridstruct,FV_Atm(1)%bd, FV_Atm(1)%flagstruct%npx, FV_Atm(1)%flagstruct%npy, &
                   FV_Atm(1)%gridstruct%nested, FV_Atm(1)%gridstruct%grid_type)
  enddo
  if (FV_Atm(1)%flagstruct%grid_type<4 .AND. present(rotate)) then 
   if (rotate) call cubed_to_latlon(utemp  , vtemp  , &
                                    uatemp , vatemp , &
                                    FV_Atm(1)%gridstruct, &
                                    FV_Atm(1)%flagstruct%npx, FV_Atm(1)%flagstruct%npy, FV_Atm(1)%flagstruct%npz, -1, &
                                    FV_Atm(1)%gridstruct%grid_type, &
                                    FV_Atm(1)%domain,FV_Atm(1)%gridstruct%nested,FV_Atm(1)%flagstruct%c2l_ord,FV_Atm(1)%bd)
  endif

  ua(:,:,:) = uatemp(isc:iec,jsc:jec,:)
  va(:,:,:) = vatemp(isc:iec,jsc:jec,:)
  if (present(uc)) uc(:,:,:) = uctemp(isc:iec,jsc:jec,:)
  if (present(vc)) vc(:,:,:) = vctemp(isc:iec,jsc:jec,:)

  return
end subroutine fv_getAgridWinds_3D
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fv_getAgridWinds_2D
!
! !INTERFACE:
!
subroutine fv_getAgridWinds_2D(u, v, ua, va, rotate)

!
! !INPUT PARAMETERS:
  real(REAL8), intent(IN)  ::  u(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
  real(REAL8), intent(IN)  ::  v(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
  logical, optional, intent(IN) :: rotate
!
! !OUTPUT PARAMETERS:
  real(REAL8), intent(OUT) :: ua(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
  real(REAL8), intent(OUT) :: va(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec,FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
!
! !DESCRIPTION:
!
! !LOCAL VARIABLES:
  real(FVPRC) :: wbuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
  real(FVPRC) :: sbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec)
  real(FVPRC) :: ebuffer(FV_Atm(1)%bd%jsc:FV_Atm(1)%bd%jec)
  real(FVPRC) :: nbuffer(FV_Atm(1)%bd%isc:FV_Atm(1)%bd%iec)

  integer isc,iec,jsc,jec
  integer i,j, npz

  real(FVPRC) :: ut(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied, FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
  real(FVPRC) :: vt(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied, FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)

  real(FVPRC) ::  utemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1)
  real(FVPRC) ::  vtemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
  real(FVPRC) :: uatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
  real(FVPRC) :: vatemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
  real(FVPRC) :: uctemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied+1,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed)
  real(FVPRC) :: vctemp(FV_Atm(1)%bd%isd:FV_Atm(1)%bd%ied,FV_Atm(1)%bd%jsd:FV_Atm(1)%bd%jed+1)
!EOP
!------------------------------------------------------------------------------
!BOC
  isc=FV_Atm(1)%bd%isc ; iec=FV_Atm(1)%bd%iec
  jsc=FV_Atm(1)%bd%jsc ; jec=FV_Atm(1)%bd%jec
  npz = 1

  utemp  = 0d0
  vtemp  = 0d0
  uatemp = 0d0
  vatemp = 0d0
  uctemp = 0d0
  vctemp = 0d0

  if (FV_Atm(1)%flagstruct%grid_type>=4) then
  ! Doubly Periodic
    uatemp(isc:iec,jsc:jec) = u
    vatemp(isc:iec,jsc:jec) = v
    call mpp_update_domains(uatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.false.)
    call mpp_update_domains(vatemp, FV_Atm(1)%domain, &
                            whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
    utemp(isc:iec,jsc:jec+1) = uatemp(isc:iec,jsc:jec+1)
    vtemp(isc:iec+1,jsc:jec) = vatemp(isc:iec+1,jsc:jec)
  else
    utemp(isc:iec,jsc:jec) = u
    vtemp(isc:iec,jsc:jec) = v
  ! update shared edges
    call mpp_get_boundary(utemp, vtemp, FV_Atm(1)%domain, &
                          wbuffery=wbuffer, ebuffery=ebuffer, &
                          sbufferx=sbuffer, nbufferx=nbuffer, &
                          gridtype=DGRID_NE, complete=.true. )
    do i=isc,iec
       utemp(i,jec+1) = nbuffer(i)
    enddo
    do j=jsc,jec
       vtemp(iec+1,j) = ebuffer(j)
    enddo
  endif

  call mpp_update_domains(utemp, vtemp, FV_Atm(1)%domain, gridtype=DGRID_NE, complete=.true.)
  call d2a2c_vect( utemp(:,:),  vtemp(:,:), &
                  uatemp(:,:), vatemp(:,:), &
                  uctemp(:,:), vctemp(:,:), ut, vt, .true., &
                  FV_Atm(1)%gridstruct,FV_Atm(1)%bd, FV_Atm(1)%flagstruct%npx, FV_Atm(1)%flagstruct%npy, &
                  FV_Atm(1)%gridstruct%nested, FV_Atm(1)%gridstruct%grid_type)

  if (FV_Atm(1)%flagstruct%grid_type<4 .AND. present(rotate)) then 
   if (rotate) call cubed_to_latlon(utemp  , vtemp  , &
                                    uatemp , vatemp , &
                                    FV_Atm(1)%gridstruct, &
                                    FV_Atm(1)%flagstruct%npx, FV_Atm(1)%flagstruct%npy, 1, -1, &
                                    FV_Atm(1)%gridstruct%grid_type, &
                                    FV_Atm(1)%domain,FV_Atm(1)%gridstruct%nested,FV_Atm(1)%flagstruct%c2l_ord,FV_Atm(1)%bd)
  endif

  ua(:,:) = uatemp(isc:iec,jsc:jec)
  va(:,:) = vatemp(isc:iec,jsc:jec)

  return
end subroutine fv_getAgridWinds_2D
!EOC
!------------------------------------------------------------------------------

  subroutine CREATE_VARS (I1, IN, J1, JN, K1, KN, KP, &
       U, V, PT, PE, PKZ, DZ, W, VARS )

    integer, intent(IN   ) :: I1, IN, J1, JN, K1, KN, KP
    real(REAL8), target ::   U(I1:IN,J1:JN,K1:KN  )
    real(REAL8), target ::   V(I1:IN,J1:JN,K1:KN  )
    real(REAL8), target ::  PT(I1:IN,J1:JN,K1:KN  )
    real(REAL8), target ::  PE(I1:IN,J1:JN,K1:KP  )
    real(REAL8), target :: PKZ(I1:IN,J1:JN,K1:KN  )
    real(REAL8), target ::  DZ(I1:IN,J1:JN,K1:KN  )
    real(REAL8), target ::   W(I1:IN,J1:JN,K1:KN  )

    type (T_FVDYCORE_VARS), intent(INOUT) :: VARS

    VARS%U => U
    VARS%V => V
    VARS%PT => PT
    VARS%PE => PE
    VARS%PKZ => PKZ
    VARS%DZ => DZ
    VARS%W => W

    return
  end subroutine CREATE_VARS

  subroutine debug_fv_state(debug_txt,STATE)
  character(LEN=*), intent(IN) :: debug_txt
  type (T_FVDYCORE_STATE),pointer :: STATE
  integer   :: isc,iec, jsc,jec, npz   !  Local dims
  real(FVPRC) :: fac1    = 1.0
  real(FVPRC) :: fac1em2 = 1.e-2
  real(FVPRC), allocatable           :: DEBUG_ARRAY(:,:,:)
  ISC    = FV_Atm(1)%bd%isc
  IEC    = FV_Atm(1)%bd%iec
  JSC    = FV_Atm(1)%bd%jsc
  JEC    = FV_Atm(1)%bd%jec
  NPZ    = FV_Atm(1)%npz
  prt_minmax     = FV_Atm(1)%flagstruct%fv_debug
  allocate( DEBUG_ARRAY(ISC:IEC,JSC:JEC,NPZ+1) )
  if (mpp_pe()==0) print*,''
  if (mpp_pe()==0) print*,'--------------', TRIM(debug_txt), '--------------'
  DEBUG_ARRAY(:,:,1) = FV_Atm(1)%phis(isc:iec,jsc:jec)
  call prt_maxmin('PHIS', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0,   1, fac1   )
  DEBUG_ARRAY(:,:,1) = STATE%VARS%PE(:,:,NPZ+1)
  call prt_maxmin('PS  ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0,   1, fac1em2)
  DEBUG_ARRAY(:,:,1:npz) = STATE%VARS%U
  call prt_maxmin('U   ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
  DEBUG_ARRAY(:,:,1:npz) = STATE%VARS%V
  call prt_maxmin('V   ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
  DEBUG_ARRAY(:,:,1:npz) = FV_Atm(1)%ua(isc:iec,jsc:jec,1:npz)
  call prt_maxmin('UA  ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
  DEBUG_ARRAY(:,:,1:npz) = FV_Atm(1)%va(isc:iec,jsc:jec,1:npz)
  call prt_maxmin('VA  ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
  if (.not. SW_DYNAMICS) then
  DEBUG_ARRAY(:,:,1:npz) = STATE%VARS%PT*STATE%VARS%PKZ
  call prt_maxmin('TA  ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
 !DEBUG_ARRAY(:,:,1:npz) = FV_Atm(1)%q(isc:iec,jsc:jec,:,1)
 !call prt_maxmin('Q1  ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
  DEBUG_ARRAY(:,:,1:npz) = STATE%VARS%PKZ
  call prt_maxmin('PZ  ', DEBUG_ARRAY  , isc, iec  , jsc, jec  , 0, npz, fac1   )
  if ( .not. FV_Atm(1)%flagstruct%hydrostatic) then
      DEBUG_ARRAY(:,:,1:npz) = STATE%VARS%W
      call prt_maxmin('W   ', DEBUG_ARRAY  , isc, iec, jsc, jec, 0, npz, fac1   )
      DEBUG_ARRAY(:,:,1:npz) = STATE%VARS%DZ
      call prt_maxmin('DZ  ', DEBUG_ARRAY  , isc, iec, jsc, jec, 0, npz, fac1   )
  endif
  endif
  if (mpp_pe()==0) print*,'--------------', TRIM(debug_txt), '--------------'
  if (mpp_pe()==0) print*,''
  prt_minmax     = .false.
  deallocate( DEBUG_ARRAY)
  end subroutine debug_fv_state

      subroutine get_latlon_vector (pp, elon, elat)
      real(REAL8), intent(IN)  :: pp(2)
      real(REAL8), intent(OUT) :: elon(3), elat(3)

         elon(1) = -SIN(pp(1))
         elon(2) =  COS(pp(1))
         elon(3) =  0.0
         elat(1) = -SIN(pp(2))*COS(pp(1))
         elat(2) = -SIN(pp(2))*SIN(pp(1))
#ifdef RIGHT_HAND
         elat(3) =  COS(pp(2))
#else
! Left-hand system needed to be consistent with rest of the codes
         elat(3) = -COS(pp(2))
#endif

      end subroutine get_latlon_vector

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  ppme --- PPM scheme at vertical edges
!
! !INTERFACE:
      subroutine ppme(p,qe,delp,im,km)
! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer,  intent(in)     ::  im, km
      real(REAL8)         , intent(in)     ::  p(im,km)
      real(FVPRC), intent(in)     ::  delp(im,km)

! !INPUT/OUTPUT PARAMETERS:
      real(FVPRC), intent(out)    ::  qe(im,km+1)

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!    05.06.13   Sawyer    Inserted file ppme.F90 here, added ProTeX
!
!EOP
!-----------------------------------------------------------------------
!BOC

      integer  km1
      integer  i, k
! local arrays.
      real(REAL8) dc(im,km),delq(im,km), a6(im,km)
      real(REAL8) c1, c2, c3, tmp, qmax, qmin
      real(REAL8) a1, a2, s1, s2, s3, s4, ss3, s32, s34, s42
      real(REAL8) a3, b2, sc, dm, d1, d2, f1, f2, f3, f4
      real(REAL8) qm, dq

      real(REAL8), parameter ::  D1EM14                  =  1.0e-14
      real(REAL8), parameter ::  D3_0                    =  3.0
      real(REAL8), parameter ::  D5_0                    =  5.0
      real(REAL8), parameter ::  D8_0                    =  8.0

      km1 = km - 1

      do 500 k=2,km
      do 500 i=1,im
500   a6(i,k) = delp(i,k-1) + delp(i,k)

      do 1000 k=1,km1
      do 1000 i=1,im
      delq(i,k) = p(i,k+1) - p(i,k)
1000  continue

      do 1220 k=2,km1
      do 1220 i=1,im
      c1 = (delp(i,k-1)+D0_5*delp(i,k))/a6(i,k+1)
      c2 = (delp(i,k+1)+D0_5*delp(i,k))/a6(i,k)
      tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /    &
                                    (a6(i,k)+delp(i,k+1))
      qmax = max(p(i,k-1),p(i,k),p(i,k+1)) - p(i,k)
      qmin = p(i,k) - min(p(i,k-1),p(i,k),p(i,k+1))
      dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp)
1220  continue

!****6***0*********0*********0*********0*********0*********0**********72
! 4th order interpolation of the provisional cell edge value
!****6***0*********0*********0*********0*********0*********0**********72

      do 12 k=3,km1
      do 12 i=1,im
      c1 = delq(i,k-1)*delp(i,k-1) / a6(i,k)
      a1 = a6(i,k-1) / (a6(i,k) + delp(i,k-1))
      a2 = a6(i,k+1) / (a6(i,k) + delp(i,k))
      qe(i,k) = p(i,k-1) + c1 + D2_0/(a6(i,k-1)+a6(i,k+1)) *        &
                ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -         &
                                delp(i,k-1)*a1*dc(i,k  ) )
12    continue

! three-cell parabolic subgrid distribution at model top

      do 10 i=1,im
! three-cell PP-distribution
! Compute a,b, and c of q = aP**2 + bP + c using cell averages and delp
! a3 = a / 3
! b2 = b / 2
      s1 = delp(i,1)
      s2 = delp(i,2) + s1
!
      s3 = delp(i,2) + delp(i,3)
      s4 = s3 + delp(i,4)
      ss3 =  s3 + s1
      s32 = s3*s3
      s42 = s4*s4
      s34 = s3*s4
! model top
      a3 = (delq(i,2) - delq(i,1)*s3/s2) / (s3*ss3)
!
      if(abs(a3) .gt. D1EM14) then
         b2 =  delq(i,1)/s2 - a3*(s1+s2)
         sc = -b2/(D3_0*a3)
         if(sc .lt. D0_0 .or. sc .gt. s1) then
             qe(i,1) = p(i,1) - s1*(a3*s1 + b2)
         else
             qe(i,1) = p(i,1) - delq(i,1)*s1/s2
         endif
      else
! Linear
         qe(i,1) = p(i,1) - delq(i,1)*s1/s2
      endif
      dc(i,1) = p(i,1) - qe(i,1)
! compute coef. for the off-centered area preserving cubic poly.
      dm = delp(i,1) / (s34*ss3*(delp(i,2)+s3)*(s4+delp(i,1)))
      f1 = delp(i,2)*s34 / ( s2*ss3*(s4+delp(i,1)) )
      f2 = (delp(i,2)+s3) * (ss3*(delp(i,2)*s3+s34+delp(i,2)*s4)   &
            + s42*(delp(i,2)+s3+s32/s2))
      f3 = -delp(i,2)*( ss3*(s32*(s3+s4)/(s4-delp(i,2))            &
            + (delp(i,2)*s3+s34+delp(i,2)*s4))                     &
            + s42*(delp(i,2)+s3) )
      f4 = ss3*delp(i,2)*s32*(delp(i,2)+s3) / (s4-delp(i,2))
      qe(i,2) = f1*p(i,1)+(f2*p(i,2)+f3*p(i,3)+f4*p(i,4))*dm
10    continue

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do 15 i=1,im
      d1 = delp(i,km)
      d2 = delp(i,km1)
      qm = (d2*p(i,km)+d1*p(i,km1)) / (d1+d2)
      dq = D2_0*(p(i,km1)-p(i,km)) / (d1+d2)
      c1 = (qe(i,km1)-qm-d2*dq) / (d2*(D2_0*d2*d2+d1*(d2+D3_0*d1)))
      c3 = dq - D2_0*c1*(d2*(D5_0*d1+d2)-D3_0*d1**2)
      qe(i,km  ) = qm - c1*d1*d2*(d2+D3_0*d1)
      qe(i,km+1) = d1*(D8_0*c1*d1**2-c3) + qe(i,km)
15    continue
      return
!EOC
      end subroutine ppme
!----------------------------------------------------------------------- 

subroutine echo_fv3_setup()

   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%grid_name ,format='("FV3 grid_name: ",(A))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%grid_file ,format='("FV3 grid_file: ",(A))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%grid_type ,format='("FV3 grid_type: ",(I2))' )
   call WRITE_PARALLEL ( 'FV3 Momentum (or KE) options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%hord_mt ,format='("FV3 hord_mt: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%kord_mt ,format='("FV3 kord_mt: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%kord_wz ,format='("FV3 kord_wz: ",(I3))' )
   call WRITE_PARALLEL ( 'FV3 Vorticity & w transport options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%hord_vt ,format='("FV3 hord_vt: ",(I3))' )
   call WRITE_PARALLEL ( 'FV3 Heat & air mass (delp) transport options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%hord_tm ,format='("FV3 hord_tm: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%hord_dp ,format='("FV3 hord_dp: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%kord_tm ,format='("FV3 kord_tm: ",(I3))' )
   call WRITE_PARALLEL ( 'FV3 Tracer transport options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%ncnst ,format='("FV3 ncnst: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%hord_tr ,format='("FV3 hord_tr: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%kord_tr ,format='("FV3 kord_tr: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%lim_fac ,format='("FV3 lim_fac: ",(F7.5))' )
!  real(FVPRC)    :: scale_z = 0.   ! diff_z = scale_z**2 * 0.25
!  real(FVPRC)    :: w_max = 75.    ! max w (m/s) threshold for hydostatiic adjustment 
!  real(FVPRC)    :: z_min = 0.05   ! min ratio of dz_nonhydrostatic/dz_hydrostatic
   call WRITE_PARALLEL ( 'FV3 Damping options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%nord ,format='("FV3 nord: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%dddmp ,format='("FV3 dddmp: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d2_bg ,format='("FV3 d2_bg: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d4_bg ,format='("FV3 d4_bg: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%vtdm4 ,format='("FV3 vtdm4: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d2_bg_k1 ,format='("FV3 d2_bg_k1: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d2_bg_k2 ,format='("FV3 d2_bg_k2: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d2_divg_max_k1 ,format='("FV3 d2_divg_max_k1: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d2_divg_max_k2 ,format='("FV3 d2_divg_max_k2: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%damp_k_k1 ,format='("FV3 damp_k_k1: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%damp_k_k2 ,format='("FV3 damp_k_k2: ",(F7.5))' )
!! Additional (after the fact) terrain filter (to further smooth the terrain after cold start)
!   integer ::    n_zs_filter=0      !  number of application of the terrain filter
!   integer :: nord_zs_filter=4      !  use del-2 (2) OR del-4 (4)
!! Additional FV3 options
!   logical :: consv_am  = .false.   ! Apply Angular Momentum Correction (to zonal wind component)
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%do_sat_adj ,format='("FV3 do_sat_adj: ",(A))' )
!   logical :: do_f3d    = .false.   ! 
!   logical :: no_dycore = .false.   ! skip the dycore
!   logical :: convert_ke = .false. 
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%do_vort_damp ,format='("FV3 do_vort_damp: ",(A))' )
!   logical :: use_old_omega = .true. 
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%beta ,format='("FV3 beta: ",(F7.4))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%n_zfilter ,format='("FV3 n_zfilter: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%n_sponge ,format='("FV3 n_sponge: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d_ext ,format='("FV3 d_ext: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%nwat ,format='("FV3 nwat: ",(I3))' )
!  logical :: warm_start = .false. 
!  logical :: inline_q = .true.
!  logical :: adiabatic = .true.     ! Run without physics (full or idealized).
! Grid options:
!  real(FVPRC) :: shift_fac   =  18.   ! shift west by 180/shift_fac = 10 degrees
!  logical :: do_schmidt = .false. 
!  real(kind=R_GRID) :: stretch_fac =   1.   ! No stretching
!  real(kind=R_GRID) :: target_lat  = -90.   ! -90: no grid rotation 
!  real(kind=R_GRID) :: target_lon  =   0.   ! 
! NH core and splitting options
!   logical :: reset_eta = .false. 
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%p_fac ,format='("FV3 p_fac: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%a_imp ,format='("FV3 a_imp: ",(F7.5))' )
   call WRITE_PARALLEL ( 'FV3 Splitting options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%n_split ,format='("FV3 n_split: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%m_split ,format='("FV3 m_split: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%k_split ,format='("FV3 k_split: ",(I3))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%q_split ,format='("FV3 q_split: ",(I3))' )
!   logical :: use_logp = .false.
!   integer :: print_freq = 0 ! Print max/min of selected fields
   call WRITE_PARALLEL ( 'FV3 Grid options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%npx ,format='("FV3 npx: ",(I4))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%npy ,format='("FV3 npy: ",(I4))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%npz ,format='("FV3 npz: ",(I4))' )
!   integer :: npz_rst = 0             ! Original Vertical Levels (in the restart)
!   integer :: pnats = 0               ! Number of non-advected consituents
!   integer :: dnats = 0               ! Number of non-advected consituents (as seen by dynamics)
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%ntiles ,format='("FV3 ntiles: ",(I4))' )
!   integer :: ndims = 2     ! Lat-Lon Dims for Grid in Radians
   call WRITE_PARALLEL ( 'FV3 Additional options:' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%nf_omega ,format='("FV3 nf_omega: ",(I7))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%fv_sg_adj ,format='("FV3 fv_sg_adj: ",(I7))' )
!   integer :: na_init = 0             ! Perform adiabatic initialization
!   real(FVPRC)    :: p_ref = 1.E5
!   real(FVPRC)    :: dry_mass = 98290.
!   integer :: nt_prog = 0
!   integer :: nt_phys = 0
!   real(FVPRC)    :: tau_h2o = 0.            ! Time scale (days) for ch4_chem
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%d_con ,format='("FV3 d_con: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%consv_te ,format='("FV3 consv_te: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%delt_max ,format='("FV3 delt_max: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%ke_bg ,format='("FV3 ke_bg: ",(F7.5))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%tau ,format='("FV3 tau: ",(F7.4))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%rf_cutoff ,format='("FV3 rf_cutoff: ",(F7.1))' )
   call WRITE_PARALLEL ( 'FV3 Logical options:' )
!   logical :: filter_phys = .false.
!   logical :: dwind_2d = .false.
!   logical :: breed_vortex_inline = .false.
!   logical :: range_warn = .false.
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%fill ,format='("FV3 fill: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%fill_dp ,format='("FV3 fill_dp: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%fill_wz ,format='("FV3 fill_wz: ",(A))' )
!   logical :: check_negative = .false.
!   logical :: non_ortho = .true.
!   logical :: moist_phys = .true.     ! Run with moist physics
!   logical :: do_Held_Suarez = .false.
!   logical :: do_reed_physics = .false.
!   logical :: reed_cond_only = .false.
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%reproduce_sum ,format='("FV3 reproduce_sum: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%adjust_dry_mass ,format='("FV3 adjust_dry_mass: ",(A))' )
!   logical :: fv_debug  = .false.
!   logical :: srf_init  = .false.
!   logical :: mountain  = .true.
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%remap_option ,format='("FV3 remap_option: ",(I4))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%z_tracer ,format='("FV3 z_tracer: ",(A))' )
!   logical :: old_divg_damp = .false. ! parameter to revert damping parameters back to values
!   logical :: fv_land = .false.       ! To cold starting the model with USGS terrain
!   logical :: nudge = .false.         ! Perform nudging
!   logical :: nudge_ic = .false.      ! Perform nudging on IC
!   logical :: ncep_ic = .false.       ! use NCEP ICs 
!   logical :: fv_diag_ic = .false.    ! reconstruct IC from fv_diagnostics on lat-lon grid
!   logical :: external_ic = .false.   ! use ICs from external sources; e.g. lat-lon FV core
!   character(len=128) :: res_latlon_dynamics = 'INPUT/fv_rst.res.nc'
!   character(len=128) :: res_latlon_tracers  = 'INPUT/atmos_tracers.res.nc'
! Parameters related to non-hydrostatic dynamics:
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%hydrostatic ,format='("FV3 hydrostatic: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%phys_hydrostatic ,format='("FV3 phys_hydrostatic: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%hybrid_z ,format='("FV3 hybrid_z: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%Make_NH ,format='("FV3 Make_NH: ",(A))' )
   call WRITE_PARALLEL_L ( FV_Atm(1)%flagstruct%make_hybrid_z ,format='("FV3 make_hybrid_z: ",(A))' )
!   real(FVPRC)    :: add_noise = -1.            !Amplitude of random noise added upon model startup; <=0 means no noise added
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%a2b_ord ,format='("FV3 a2b_ord: ",(I4))' )
   call WRITE_PARALLEL ( FV_Atm(1)%flagstruct%c2l_ord ,format='("FV3 c2l_ord: ",(I4))' )
!! Doubly preiodic options:
!  real(kind=R_GRID) :: dx_const = 1000.    ! spatial resolution for double periodic boundary configuration [m]
!  real(kind=R_GRID) :: dy_const = 1000.
!  real(kind=R_GRID) :: deglat=15.
!  real(kind=R_GRID) :: deglon_start = -30., deglon_stop = 30., &  ! boundaries of latlon patch
!                       deglat_start = -30., deglat_stop = 30.
!! Convenience pointers
!  integer, pointer :: grid_number

end subroutine echo_fv3_setup

subroutine WRITE_PARALLEL_L ( field, format )
  logical, intent(in) :: field
  character(len=*), intent(in ), optional :: format

  if (field) then
   call WRITE_PARALLEL ( 'T' ,format=format )
  else
   call WRITE_PARALLEL ( 'F' ,format=format )
  endif
end subroutine WRITE_PARALLEL_L

end module FV_StateMod

