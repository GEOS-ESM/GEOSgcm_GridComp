! $Id$

! VERIFY_ and RETURN_ macros for error handling.

!#define UWDIAG 1

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GFS_Moist -- A Module to compute moist processes from CCPP

! !INTERFACE:

module GFS_MoistGridCompMod

  ! !USES:
  use ESMF
  use MAPL, r8 => MAPL_R8
  use GEOS_UtilsMod

   use gfdl_cloud_microphys, only: gfdl_cloud_microphys_init
   use get_phi_fv3, only: get_phi_fv3_run
   use GFS_suite_interstitial_1, only: GFS_suite_interstitial_1_run
   use GFS_suite_interstitial_3, only: GFS_suite_interstitial_3_run
   use GFS_DCNV_generic_pre, only: GFS_DCNV_generic_pre_run
   use samfdeepcnv, only: samfdeepcnv_run
   use GFS_DCNV_generic_post, only: GFS_DCNV_generic_post_run
   use GFS_SCNV_generic_pre, only: GFS_SCNV_generic_pre_run
   use samfshalcnv, only: samfshalcnv_run
   use GFS_SCNV_generic_post, only: GFS_SCNV_generic_post_run
   use GFS_suite_interstitial_4, only: GFS_suite_interstitial_4_run
   use cnvc90, only: cnvc90_run
   use GFS_MP_generic_pre, only: GFS_MP_generic_pre_run
   use gfdl_cloud_microphys, only: gfdl_cloud_microphys_run
   use GFS_MP_generic_post, only: GFS_MP_generic_post_run
   use maximum_hourly_diagnostics, only: maximum_hourly_diagnostics_run
   use phys_tend, only: phys_tend_run
   use machine, only: kind_phys
   use gfdl_cloud_microphys, only: gfdl_cloud_microphys_finalize

  implicit none

  private

  ! !PUBLIC MEMBER FUNCTIONS:

  character(len=ESMF_MAXSTR) :: HYDROSTATIC ! TRUE or FALSE
  logical                    :: LHYDROSTATIC

  character(len=ESMF_MAXSTR) :: PHYS_HYDROSTATIC ! TRUE or FALSE
  logical                    :: LPHYS_HYDROSTATIC

  !integer :: DOSHLW

#include "namelist_phy.inc" ! Physics namelist read from input.nml

  public SetServices

  ! !DESCRIPTION:
  ! 
  !   {\tt GEOS\_MoistGridCompMod} implements moist processes in GEOS-5. These
  !   include all processes that involve phase changes in the atmosphere, such
  !   as large-scale condensation, convective clouds, and all rain and cloud
  !   formation. Its state consists of water vapor, various types of condensate,
  !   and fractions of various cloud types. 
  !   two moment cloud microphysics (Barahona et al., GMD, 2014.) can be run by setting CLDMACRO==2MOMENT. 
  !   When using 2-moment microphysics the number concentration of ice crystals and cloud droplets 
  !   are part of the state.
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

    ! !DESCRIPTION:  {\tt GEOS\_MoistGridCompMod} uses the default Initialize and Finalize 
    !                services, but registers its own Run method.

    !EOP

    !=============================================================================
    !
    ! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    ! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: STATE 
    type (ESMF_Config      )            :: CF

    integer      :: RFRSHINT
    integer      :: AVRGNINT
    !integer      :: IQVAINC
    real         :: DT
    
    character(len=ESMF_MAXSTR) :: FRIENDLIES_NCPL , FRIENDLIES_NCPI , &
                                  FRIENDLIES_NRAIN, FRIENDLIES_NSNOW, FRIENDLIES_NGRAUPEL
    character(len=ESMF_MAXSTR) :: FRIENDLIES_QRAIN, FRIENDLIES_QSNOW, FRIENDLIES_QGRAUPEL

    !=============================================================================

    ! Begin...

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


    ! Set the Init entry point
    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
    VERIFY_(STATUS)

    ! Get the configuration from the component
    !-----------------------------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

    ! Set the Run entry point
    ! -----------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,  &
         RC=STATUS)
    VERIFY_(STATUS)

    ! Set the state variable specs.
    ! -----------------------------

    call ESMF_ConfigGetAttribute ( CF, DT, Label="RUN_DT:",                                   RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, RFRSHINT, Label="REFRESH_INTERVAL:",  default=nint(DT), RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, AVRGNINT, Label='AVERAGING_INTERVAL:',default=RFRSHINT, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_ConfigGetAttribute( CF, HYDROSTATIC, Label="HYDROSTATIC:",  default="TRUE", RC=STATUS)
    VERIFY_(STATUS)
    if (adjustl(HYDROSTATIC)=="TRUE" ) LHYDROSTATIC=.true.
    if (adjustl(HYDROSTATIC)=="FALSE") LHYDROSTATIC=.false.

    call ESMF_ConfigGetAttribute( CF, PHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default="TRUE", RC=STATUS)
    VERIFY_(STATUS)
    if (adjustl(PHYS_HYDROSTATIC)=="TRUE" ) LPHYS_HYDROSTATIC=.true.
    if (adjustl(PHYS_HYDROSTATIC)=="FALSE") LPHYS_HYDROSTATIC=.false.

    FRIENDLIES_NCPI     = trim(COMP_NAME)
    FRIENDLIES_NCPL     = trim(COMP_NAME)
    FRIENDLIES_NRAIN    = trim(COMP_NAME)    
    FRIENDLIES_NSNOW    = trim(COMP_NAME)
    FRIENDLIES_NGRAUPEL = trim(COMP_NAME)
    FRIENDLIES_QRAIN = 'DYNAMICS:TURBULENCE'
    FRIENDLIES_QSNOW = 'DYNAMICS:TURBULENCE'
    FRIENDLIES_QGRAUPEL = 'DYNAMICS:TURBULENCE'

#include "AddSpec.inc"  ! All MAPL_Add[Export/Import/Internal]Spec

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,  Finalize,  &
         RC=STATUS)
    VERIFY_(STATUS)


    ! Set generic init and final methods
    ! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!-!-!-!!!!!!!!!!!!!!!

  !EOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOP

  ! !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component

  ! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

    use, intrinsic :: iso_fortran_env, only : output_unit
    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    ! !DESCRIPTION: The Initialize method of the Moist Physics Gridded Component first
    !   calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
    !   create the ESMF GRID, which will then be used to set the GRID associated with the
    !   SuperDyn Composite Component itself.  It should be noted that the
    !   SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
    !   topography related quantities.

    !EOP

    ! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    type(ESMF_VM) :: vm
    integer :: me, nlunit
    character(len=ESMF_MAXSTR) :: input_nml_file(1)=""
    integer, parameter :: imp_physics_gfdl = 11
    character(len=ESMF_MAXSTR) :: fn_nml="input_GFS_v16beta.nml"
    character(len=ESMF_MAXSTR) :: errmsg

    Iam = 'Initialize'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Call Generic Initialize
!------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGetCurrent(vm, rc=status)
    VERIFY_(STATUS)
    call ESMF_VMGet(vm, localPet=me, rc=status)
    VERIFY_(STATUS)
    nlunit = get_file_unit(99)
! Read GFS physics namelist. 
    open(unit=nlunit, file=trim(fn_nml))
    rewind(nlunit)
    read (nlunit, nml=gfs_physics_nml)
    close (nlunit)

    call gfdl_cloud_microphys_init(me=me,master=MAPL_Root,nlunit=nlunit, &
                  input_nml_file=input_nml_file(:),logunit=output_unit, &
                  fn_nml=fn_nml,imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl, &
                  do_shoc=do_shoc,errmsg=errmsg,errflg=status)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
    contains
      integer function get_file_unit (lu_max)
!
!   get_file_unit returns a unit number that is not in use
      integer :: lu_max,  lu, m, iostat
      logical ::  opened
!
      m = lu_max  ;  if (m < 1) m = 97
      do lu = m,1,-1
         inquire (unit=lu, opened=opened, iostat=iostat)
         if (iostat.ne.0) cycle
         if (.not.opened) exit
      end do
!
      get_file_unit = lu
      return
      end function get_file_unit
   end subroutine Initialize


  !===================================================================================

  !BOP

  ! !IROUTINE: RUN -- Run method for the CONVECT component

  ! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

    use physcons ! constants (TODO: check overlap with MAPL constants)
    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
    !                the Initialize and Finalize services, as well as allocating

    !EOP


    ! ErrLog Variables

    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: STATE
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM

    ! Local variables

    type(ESMF_VM) :: vm
    integer :: me
    integer                         :: IM,JM,LM
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS

    ! Begin...

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
    VERIFY_(STATUS)

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( STATE, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         CF       = CF,                &
         LONS     = LONS,              &
         LATS     = LATS,              &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

    ! retrieve proc id
    call ESMF_VMGetCurrent(vm, rc=status)
    VERIFY_(STATUS)
    call ESMF_VMGet(vm, localPet=me, rc=status)
    VERIFY_(STATUS)

    ! If its time, calculate convective tendencies
    ! --------------------------------------------

    if ( ESMF_AlarmIsRinging( ALARM, RC=status) ) then
       VERIFY_(STATUS)
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)
       call MOIST_DRIVER(IM,JM,LM, RC=STATUS)
       VERIFY_(STATUS)
    endif

    RETURN_(ESMF_SUCCESS)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MOIST_DRIVER(IM,JM,LM, RC)
      integer,           intent(IN ) :: IM, JM, LM
      integer, optional, intent(OUT) :: RC

      !locals
      real    :: HEARTBEAT, DT_MOIST
      real(ESMF_KIND_R8)  :: DT_R8
      type (ESMF_TimeInterval)  :: TINT
      real, pointer, dimension(:,:  ) :: FRLANDICE, FRLAND
      real, pointer, dimension(:,:  ) :: AREA, KPBLIN
      real, pointer, dimension(:,:,:) :: PLE, T, U, V, TH, OMEGA
      real, pointer, dimension(:,:,:) :: Q, QRAIN, QSNOW, QGRAUPEL, QLLS, QLCN, CLLS, CLCN, QILS, QICN
      logical, save :: compute_firsttime = .true.
      integer :: i, j, k, ij
#ifdef IMJM
#undef IMJM
#endif
#define IMJM IM*JM      
#include "phy_var_declarations.inc"
      real, pointer, dimension(:,:) :: slmsk => Null()
      real, dimension(IM,JM,LM,ntrac) :: gq0
      real, dimension(IM,JM,LM)  ::  FQAL, FQAI, FQA
      real, dimension(IM,JM,0:LM):: PKE
      real, dimension(IM,JM,0:LM):: ZLE  !geopotential at model layer interfaces

      call ESMF_ConfigGetAttribute (CF, HEARTBEAT, Label="RUN_DT:", RC=STATUS); VERIFY_(STATUS)

      ! Get the time step from the alarm
      !---------------------------------

      call ESMF_AlarmGet( ALARM, RingInterval=TINT, RC=STATUS)
      VERIFY_(STATUS)
      call ESMF_TimeIntervalGet(TINT, S_R8=DT_R8, RC=STATUS)
      VERIFY_(STATUS)

      DT_MOIST = DT_R8

      ! Pointers to internals
      !----------------------

      call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)

      ! Pointers to imports
      !--------------------
      call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, AREA,    'AREA'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FRLANDICE,  'FRLANDICE'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, FRLAND,  'FRLAND'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, T,       'T'       ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, KPBLIN,  'KPBL'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, TH,      'TH'      ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   ,RC=STATUS); VERIFY_(STATUS)

      if (compute_firsttime) then
         allocate(slmsk(IM,JM), stat=status)
         VERIFY_(STATUS)
         slmsk = zero
         where(FRLAND >= 0.5) slmsk=1.0
         where(FRLANDICE >= 0.5) slmsk=2.0

         tem     = con_rerth*con_rerth*(con_pi+con_pi)*con_pi
         dxmax  = log(tem/(max_lon*max_lat))
         dxmin  = log(tem/(min_lon*min_lat))
         dxinv  = 1.0d0 / (dxmax-dxmin)
         compute_firsttime = .false.
      end if

      call GFS_suite_interstitial_1_run(                       &
                  !ntent(in)
                  im=IMJM,levs=LM,ntrac=ntrac, &
                  dtf=HEARTBEAT,dtp=DT_MOIST,slmsk=slmsk, &
                  area=AREA,dxmin=dxmin,dxinv=dxinv, &
                  pgr=PLE(:,:,LM),                           &
                  !ntent(out)
                  islmsk=islmsk,work1=work1, &
                  work2=work2,psurf=psurf,dudt=dudt, &
                  dvdt=dvdt,dtdt=dtdt,dtdtc=dtdtc, &
                  dqdt=dqdt,errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)
    
      gq0(:,:,:,ntqv)    = Q
      gq0(:,:,:,ntcw)    = QLLS+QLCN
      gq0(:,:,:,ntiw)    = QILS+QICN
      gq0(:,:,:,ntrw)    = QRAIN
      gq0(:,:,:,ntsw)    = QSNOW
      gq0(:,:,:,ntgl)    = QGRAUPEL
      gq0(:,:,:,ntclamt) = MIN(CLCN+CLLS,1.0)
      FQA =  MIN(1.0,MAX(CLCN/MAX(gq0(:,:,:,ntclamt),1.e-5),0.0))
      FQAL =  MIN(1.0,MAX(QLCN/MAX(gq0(:,:,:,ntcw),1.E-8),0.0))
      FQAI =  MIN(1.0,MAX(QICN/MAX(gq0(:,:,:,ntiw),1.E-8),0.0))

      do j = 1, JM
         do i = 1, IM
            ij =  i + (j-1)*IM
            prsl(ij,1:LM) = 0.5*(PLE(i,j,0:LM-1) + PLE(i,j,1:LM))
            delp(ij,1:LM) = PLE(i,j,0:LM-1) - PLE(i,j,1:LM)
            prslk(ij,1:LM) = (prsl(ij,1:LM)/MAPL_P00)**MAPL_KAPPA
            kpbl(ij) = int(KPBLIN(i,j))
         enddo
      enddo

      kinver = 0   ! not used
      clw = zero
      clw(:,:,2) = -999.9
      rhc = zero
      save_q = zero
      save_tcp = zero
      call GFS_suite_interstitial_3_run(   &
                  !ntent(in)
                  im=IMJM,levs=LM,nn=nn, &
                  cscnv=cscnv,satmedmf=satmedmf,trans_trac=trans_trac, &
                  do_shoc=do_shoc,ltaerosol=ltaerosol,ntrac=ntrac, &
                  ntcw=ntcw,ntiw=ntiw,ntclamt=ntclamt, &
                  ntrw=ntrw,ntsw=ntsw,ntrnc=ntrnc, &
                  ntsnc=ntsnc,ntgl=ntgl,ntgnc=ntgnc, &
                  xlon=LONS,xlat=LATS,gt0=T, &
                  gq0=gq0,imp_physics=imp_physics,imp_physics_mg=imp_physics_mg, &
                  imp_physics_zhao_carr=imp_physics_zhao_carr,imp_physics_zhao_carr_pdf=imp_physics_zhao_carr_pdf, &
                  imp_physics_gfdl=imp_physics_gfdl,imp_physics_thompson=imp_physics_thompson, &
                  imp_physics_wsm6=imp_physics_wsm6,imp_physics_fer_hires=imp_physics_fer_hires, &
                  prsi=PLE,prsl=prsl,prslk=prslk, &
                  rhcbot=crtrh(1),rhcpbl=crtrh(2),rhctop=crtrh(3), &
                  rhcmax=rhcmax,islmsk=islmsk,work1=work1, &
                  work2=work2,kpbl=kpbl,kinver=kinver, &
                  ras=ras,me=me,  &
                  !ntent(inout)
                  clw=clw, &
                  rhc=rhc,save_qc=save_q(:,:,ntcw), &
                  save_qi=save_q(:,:,ntiw),save_tcp=save_tcp, &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      call GFS_DCNV_generic_pre_run(   &
                  !ntent(in)
                  im=IMJM,levs=LM,ldiag3d=ldiag3d, &
                  qdiag3d=qdiag3d,do_cnvgwd=do_cnvgwd,cplchm=cplchm, &
                  gu0=U,gv0=V,gt0=T, &
                  !ntent(inout)
                  gq0_water_vapor=gq0(:,:,:,ntqv),save_u=save_u, &
                  save_v=save_v,save_t=save_t,save_qv=save_q(:,:,ntqv), &
                  dqdti=dqdti,   &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      !AOO Adapted from GEOS_TurbulenceGridComp.F90
      ! Compute the edge heights using Arakawa-Suarez hydrostatic equation
      !---------------------------------------------------------------------------

      PKE = (PLE/MAPL_P00)**MAPL_KAPPA
      ZLE(:,:,LM) = zero
      do k = LM, 1, -1
         ZLE(:,:,k-1) = ZLE(:,:,k) + MAPL_CP*TH(:,:,k)*(PKE(:,:,k)-PKE(:,:,k-1)) !m2 s-2
      end do

      do j = 1, JM
         do i = 1, IM
            ij =  i + (j-1)*IM
            phil(ij,1:LM) = 0.5*(ZLE(i,j,1:LM-1) + ZLE(i,j,1:LM))
         enddo
      enddo
      kcnv = 0
      qlcn_mg=0.; qicn_mg=0.; w_upi=0.; cf_upi=0.; cnv_mfd=0.   ! not used for GFDL MP
      cnv_dqldt=0.; clcn_mg=0.; cnv_fice=0.; cnv_ndrop=0.; cnv_nice=0.
      ca_deep=0.   ! no cellular automata
      call samfdeepcnv_run(           &
                  !ntent(in)
                  im=IMJM,km=LM,itc=itc, &
                  ntc=ntchm,cliq=con_cliq,cp=con_cp,cvap=con_cvap,eps=con_eps, &
                  epsm1=con_epsm1,fv=con_fvirt,grav=con_g,hvap=con_hvap,rd=con_rd,rv=con_rv, &
                  t0c=con_t0c,delt=DT_MOIST,ntk=ntk,ntr=nsamftrac, &
                  delp=delp,prslp=prsl,psp=PLE(:,:,LM), &
                  phil=phil,             &
                  !ntent(inout)
                  qtr=clw,q1=gq0(:,:,:,ntqv), &
                  t1=T,u1=U,v1=U, &
                  !ntent(in)
                  fscav=fscav,hwrf_samfdeep=hwrf_samfdeep, &
                  !ntent(out)
                  cldwrk=cld1d, &
                  rn=raincd,kbot=kbot,ktop=ktop, &
                  !ntent(inout)
                  kcnv=kcnv, &
                  !ntent(in)
                  islimsk=islmsk,garea=AREA, &
                  dot=OMEGA,ncloud=ncld, &
                  !ntent(out)
                  ud_mf=ud_mf, &
                  dd_mf=dd_mf,dt_mf=dt_mf, &
                  !ntent(inout)
                  cnvw=cnvw, &
                  cnvc=cnvc,qlcn=qlcn_mg,qicn=qicn_mg, &
                  w_upi=w_upi,cf_upi=cf_upi,cnv_mfd=cnv_mfd, &
                  cnv_dqldt=cnv_dqldt,clcn=clcn_mg, &
                  cnv_fice=cnv_fice,cnv_ndrop=cnv_ndrop, &
                  cnv_nice=cnv_nice,mp_phys=imp_physics, &
                  !ntent(in)
                  mp_phys_mg=imp_physics_mg,clam=clam_deep,c0s=c0s_deep, &
                  c1=c1_deep,betal=betal_deep,betas=betas_deep, &
                  evfact=evfact_deep,evfactl=evfactl_deep,pgcon=pgcon_deep, &
                  asolfac=asolfac_deep,do_ca=do_ca,ca_closure=ca_closure, &
                  ca_entr=ca_entr,ca_trigger=ca_trigger,nthresh=nthresh, &
                  ca_deep=ca_deep, &
                  !ntent(out)
                  rainevap=condition,errmsg=errmsg, &
                  errflg=status)
      VERIFY_(STATUS)

      frain  = HEARTBEAT/DT_MOIST
      cldwrk = zero
      dt3dt  = zero
      dq3dt  = zero
      du3dt  = zero
      dv3dt  = zero
      call GFS_DCNV_generic_post_run(    &
                  !ntent(in)
                  im=IMJM,levs=LM,lssav=lssav, &
                  ldiag3d=ldiag3d,qdiag3d=qdiag3d,ras=ras, &
                  cscnv=cscnv,frain=frain,rain1=raincd, &
                  dtf=HEARTBEAT,cld1d=cld1d,save_u=save_u, &
                  save_v=save_v,save_t=save_t,save_qv=save_q(:,:,ntqv), &
                  gu0=U,gv0=V,gt0=T, &
                  gq0_water_vapor=gq0(:,:,:,ntqv),ud_mf=ud_mf, &
                  dd_mf=dd_mf,dt_mf=dt_mf,con_g=con_g, &
                  npdf3d=npdf3d,num_p3d=num_p3d,ncnvcld3d=ncnvcld3d, &
                  !ntent(inout)
                  rainc=rainc,cldwrk=cldwrk,dt3dt=dt3dt(:,:,4), &
                  dq3dt=dq3dt(:,:,2),du3dt=du3dt(:,:,3),dv3dt=dv3dt(:,:,3), &
                  upd_mf=upd_mf,dwn_mf=dwn_mf,det_mf=det_mf, &
                  cnvw=cnvw,cnvc=cnvc,cnvw_phy_f3d=phy_f3d(:,:,ncnvw), &
                  cnvc_phy_f3d=phy_f3d(:,:,ncnvc), &
                  !ntent(in)
                  flag_for_dcnv_generic_tend=flag_for_dcnv_generic_tend, &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)

      call GFS_SCNV_generic_pre_run(       &
                  !ntent(in)
                  im=IMJM,levs=LM,ldiag3d=ldiag3d, &
                  qdiag3d=qdiag3d,gu0=U,gv0=V, &
                  gt0=T,gq0_water_vapor=gq0(:,:,:,ntqv), &
                  !ntent(out)
                  save_u=save_u,save_v=save_v,save_t=save_t, &
                  save_qv=save_q(:,:,ntqv),   &
                  !ntent(in)
                  flag_for_scnv_generic_tend=flag_for_scnv_generic_tend, &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      call samfshalcnv_run(     &
                  !ntent(in)
                  im=IMJM,km=LM,itc=itc, &
                  ntc=ntchm,cliq=con_cliq,cp=con_cp,cvap=con_cvap,eps=con_eps, &
                  epsm1=con_epsm1,fv=con_fvirt,grav=con_g,hvap=con_hvap,rd=con_rd,rv=con_rv, &
                  t0c=con_t0c,delt=DT_MOIST,ntk=ntk,ntr=nsamftrac, &
                  delp=delp,prslp=prsl,psp=PLE(:,:,LM), &
                  phil=phil,       &
                  !ntent(inout)
                  qtr=clw,q1=gq0(:,:,:,ntqv), &
                  t1=T,u1=U,v1=V, &
                  !ntent(in)
                  fscav=fscav,     &
                  !ntent(out)
                  rn=raincs,kbot=kbot,  ktop=ktop,  &
                  !ntent(inout)
                  kcnv=kcnv,   &
                  !ntent(in)
                  islimsk=islmsk, &
                  garea=AREA,dot=OMEGA,ncloud=ncld, hpbl=hpbl,  &
                  !ntent(out)
                  ud_mf=ud_mf,dt_mf=dt_mf, &
                  cnvw=cnvw,cnvc=cnvc,     &
                  !ntent(in)
                  clam=clam_shal, &
                  c0s=c0s_shal,c1=c1_shal,pgcon=pgcon_shal, &
                  asolfac=asolfac_shal,hwrf_samfshal=hwrf_samfshal, &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      call GFS_SCNV_generic_post_run(                  &
                  !ntent(in)
                  im=IMJM,levs=LM,nn=nn, &
                  lssav=lssav,qdiag3d=qdiag3d,ldiag3d=ldiag3d, &
                  cplchm=cplchm,frain=frain,gu0=U, &
                  gv0=V,gt0=T,gq0_water_vapor=gq0(:,:,:,ntqv), &
                  save_u=save_u,save_v=save_v,save_t=save_t, &
                  save_qv=save_q(:,:,ntqv),  &
                  !ntent(inout)
                  dqdti=dqdti, &
                  du3dt=du3dt(:,:,6),dv3dt=dv3dt(:,:,6),dt3dt=dt3dt(:,:,5), &
                  dq3dt=dq3dt(:,:,3),clw=clw, &
                  !ntent(in)
                  shcnvcw=shcnvcw, &
                  rain1=raincs,npdf3d=npdf3d,num_p3d=num_p3d, &
                  ncnvcld3d=ncnvcld3d,cnvc=cnvc,cnvw=cnvw, &
                  !ntent(inout)
                  rainc=rainc,cnvprcp=cnvprcp,cnvprcpb=cnvprcpb, &
                  cnvw_phy_f3d=phy_f3d(:,:,ncnvw),cnvc_phy_f3d=phy_f3d(:,:,ncnvc), &
                  !ntent(in)
                  flag_for_scnv_generic_tend=flag_for_scnv_generic_tend,imfshalcnv=imfshalcnv, &
                  imfshalcnv_sas=imfshalcnv_sas,imfshalcnv_samf=imfshalcnv_samf, &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      call GFS_suite_interstitial_4_run(     &
                  !ntent(in)
                  im=IMJM,levs=LM,ltaerosol=ltaerosol, &
                  cplchm=cplchm,tracers_total=tracers_total, &
                  ntrac=ntrac,ntcw=ntcw,ntiw=ntiw, &
                  ntclamt=ntclamt,ntrw=ntrw,ntsw=ntsw, &
                  ntrnc=ntrnc,ntsnc=ntsnc,ntgl=ntgl, &
                  ntgnc=ntgnc,ntlnc=ntlnc,ntinc=ntinc, &
                  nn=nn,imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl, &
                  imp_physics_thompson=imp_physics_thompson,imp_physics_zhao_carr=imp_physics_zhao_carr, &
                  imp_physics_zhao_carr_pdf=imp_physics_zhao_carr_pdf,dtf=HEARTBEAT, &
                  save_qc=save_q(:,:,ntcw),save_qi=save_q(:,:,ntiw), &
                  con_pi=con_pi, &
                  !ntent(inout)
                  gq0=gq0,clw=clw,  &
                  !ntent(in)
                  prsl=prsl, &
                  save_tcp=save_tcp,con_rd=con_rd,nwfa=qgrs(:,:,ntwa), &
                  !ntent(inout)
                  spechum=gq0(:,:,:,ntqv),dqdti=dqdti, &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      call GFS_MP_generic_pre_run(           &
                  !ntent(in)
                  im=IMJM,levs=LM,ldiag3d=ldiag3d, &
                  qdiag3d=qdiag3d,do_aw=do_aw,ntcw=ntcw, &
                  nncl=nncl,ntrac=ntrac,gt0=T, &
                  gq0=gq0,     &
                  !ntent(inout)
                  save_t=save_t,     &
                  !ntent(out)
                  save_qv=save_q(:,:,ntqv), &
                  !ntent(inout)
                  save_q=save_q,                     &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      print *, IMJM, LHYDROSTATIC, LPHYS_HYDROSTATIC, lradar, reset, effr_in
      call gfdl_cloud_microphys_run(     &
                  !ntent(in)
                  levs=LM,im=IMJM,con_g=con_g, &
                  con_fvirt=con_fvirt,con_rd=con_rd,frland=FRLAND,garea=AREA, &
                  islmsk=islmsk,     &
                  !ntent(inout)
                  gq0=gq0(:,:,:,ntqv), &
                  gq0_ntcw=gq0(:,:,:,ntcw),gq0_ntrw=gq0(:,:,:,ntrw), &
                  gq0_ntiw=gq0(:,:,:,ntiw),gq0_ntsw=gq0(:,:,:,ntsw), &
                  gq0_ntgl=gq0(:,:,:,ntgl),gq0_ntclamt=gq0(:,:,:,ntclamt), &
                  gt0=T,gu0=U,gv0=V, &
                  !ntent(in)
                  vvl=OMEGA,prsl=prsl,phii=ZLE, &
                  del=delp,        &
                  !ntent(out)
                  rain0=rainmp,ice0=icemp, &
                  snow0=snowmp,graupel0=graupelmp, &
                  prcp0=prcpmp,sr=sr,  &
                  !ntent(in)
                  dtp=DT_MOIST, &
                  hydrostatic=LHYDROSTATIC,phys_hydrostatic=LPHYS_HYDROSTATIC, &
                  lradar=lradar,    &
                  !ntent(inout)
                  refl_10cm=refl_10cm,  &
                  !ntent(in)
                  reset=reset, &
                  effr_in=effr_in,   &
                  !ntent(inout)
                  rew=phy_f3d(:,:,nleffr), &
                  rei=phy_f3d(:,:,nieffr),rer=phy_f3d(:,:,nreffr), &
                  res=phy_f3d(:,:,nseffr),reg=phy_f3d(:,:,ngeffr), &
                  !ntent(out)
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      call ESMF_Finalize()
      stop
      if (me == 0) print *, "phy_f3d(:,:,nleffr) ... ", nleffr
      if (me == 0) print *, phy_f3d(:,:,nleffr)
      if (me == 0) print *, "rainmp ... "
      if (me == 0) print *, rainmp

      call GFS_MP_generic_post_run(im=IMJM,levs=LM,kdt=kdt, &
                  nrcm=nrcm,ncld=ncld,nncl=nncl, &
                  ntcw=ntcw,ntrac=ntrac,imp_physics=imp_physics, &
                  imp_physics_gfdl=imp_physics_gfdl,imp_physics_thompson=imp_physics_thompson, &
                  imp_physics_mg=imp_physics_mg,imp_physics_fer_hires=imp_physics_fer_hires, &
                  cal_pre=cal_pre,lssav=lssav,ldiag3d=ldiag3d, &
                  qdiag3d=qdiag3d,cplflx=cplflx,cplchm=cplchm, &
                  con_g=con_g,dtf=HEARTBEAT,frain=frain,rainc=rainc, &
                  rain1=prcpmp,rann=rann,xlat=LATS, &
                  xlon=LONS,gt0=T,gq0=gq0, &
                  prsl=prsl,prsi=prsi,phii=phii, &
                  tsfc=tsfc,ice=ice,snow=snow,graupel=graupel, &
                  save_t=save_t,save_qv=save_q(:,:,ntqv), &
                  rain0=rainmp,ice0=icemp,snow0=snowmp, &
                  graupel0=graupelmp,del=del,rain=rain, &
                  domr_diag=tdomr,domzr_diag=tdomzr,domip_diag=tdomip, &
                  doms_diag=tdoms,tprcp=tprcp,srflag=srflag, &
                  sr=sr,cnvprcp=cnvprcp,totprcp=totprcp, &
                  totice=totice,totsnw=totsnw,totgrp=totgrp, &
                  cnvprcpb=cnvprcpb,totprcpb=totprcpb,toticeb=toticeb, &
                  totsnwb=totsnwb,totgrpb=totgrpb,dt3dt=dt3dt(:,:,6), &
                  dq3dt=dq3dt(:,:,4),rain_cpl=rain_cpl,rainc_cpl=rainc_cpl, &
                  snow_cpl=snow_cpl,pwat=pwat,do_sppt=do_sppt, &
                  ca_global=ca_global,dtdtr=dtdtr,dtdtc=dtdtc, &
                  drain_cpl=drain_cpl,dsnow_cpl=dsnow_cpl,lsm=lsm, &
                  lsm_ruc=lsm_ruc,lsm_noahmp=lsm_noahmp,raincprv=raincprv, &
                  rainncprv=rainncprv,iceprv=iceprv,snowprv=snowprv, &
                  graupelprv=graupelprv,draincprv=draincprv, &
                  drainncprv=drainncprv,diceprv=diceprv,dsnowprv=dsnowprv, &
                  dgraupelprv=dgraupelprv,dtp=DT_MOIST,errmsg=errmsg, &
                  errflg=status)
      VERIFY_(STATUS)



      call phys_tend_run(ldiag3d=ldiag3d,qdiag3d=qdiag3d,du3dt_pbl=du3dt(:,:,1), &
                  du3dt_orogwd=du3dt(:,:,2),du3dt_deepcnv=du3dt(:,:,3), &
                  du3dt_congwd=du3dt(:,:,4),du3dt_rdamp=du3dt(:,:,5), &
                  du3dt_shalcnv=du3dt(:,:,6),du3dt_phys=du3dt(:,:,7), &
                  dv3dt_pbl=dv3dt(:,:,1),dv3dt_orogwd=dv3dt(:,:,2), &
                  dv3dt_deepcnv=dv3dt(:,:,3),dv3dt_congwd=dv3dt(:,:,4), &
                  dv3dt_rdamp=dv3dt(:,:,5),dv3dt_shalcnv=dv3dt(:,:,6), &
                  dv3dt_phys=dv3dt(:,:,7),dt3dt_lw=dt3dt(:,:,1), &
                  dt3dt_sw=dt3dt(:,:,2),dt3dt_pbl=dt3dt(:,:,3), &
                  dt3dt_deepcnv=dt3dt(:,:,4),dt3dt_shalcnv=dt3dt(:,:,5), &
                  dt3dt_mp=dt3dt(:,:,6),dt3dt_orogwd=dt3dt(:,:,7), &
                  dt3dt_rdamp=dt3dt(:,:,8),dt3dt_congwd=dt3dt(:,:,9), &
                  dt3dt_phys=dt3dt(:,:,10),dq3dt_pbl=dq3dt(:,:,1), &
                  dq3dt_deepcnv=dq3dt(:,:,2),dq3dt_shalcnv=dq3dt(:,:,3), &
                  dq3dt_mp=dq3dt(:,:,4),dq3dt_o3pbl=dq3dt(:,:,5), &
                  dq3dt_o3prodloss=dq3dt(:,:,6),dq3dt_o3mix=dq3dt(:,:,7), &
                  dq3dt_o3tmp=dq3dt(:,:,8),dq3dt_o3column=dq3dt(:,:,9), &
                  dq3dt_phys=dq3dt(:,:,10),dq3dt_o3phys=dq3dt(:,:,11), &
                  errmsg=errmsg,errflg=status)
      VERIFY_(STATUS)

      ! Redistribute cloud/liquid/ice species...
      CLCN     = gq0(:,:,:,ntclamt)*     FQA
      CLLS     = gq0(:,:,:,ntclamt)*(1.0-FQA)
      QLCN     = gq0(:,:,:,ntcw)*     FQAL
      QLLS     = gq0(:,:,:,ntcw)*(1.0-FQAL)
      QICN     = gq0(:,:,:,ntiw)*     FQAI
      QILS     = gq0(:,:,:,ntiw)*(1.0-FQAI)
      Q        = gq0(:,:,:,ntqv)
      QRAIN    = gq0(:,:,:,ntrw)
      QSNOW    = gq0(:,:,:,ntsw)
      QGRAUPEL = gq0(:,:,:,ntgl)

      RETURN_(ESMF_SUCCESS)

    end subroutine MOIST_DRIVER

!!!!!!!!-!-!-!!!!!!!!
   end subroutine Run


  !===================================================================================

  !BOP

  ! !IROUTINE: RUN -- Run method for the CONVECT component

  ! !INTERFACE:

  subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
    !                the Initialize and Finalize services, as well as allocating

    !EOP


    ! ErrLog Variables

    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME
    character(len=ESMF_MAXSTR)      :: errmsg

    call gfdl_cloud_microphys_finalize(errmsg=errmsg,errflg=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

   end subroutine Finalize

end module GFS_MoistGridCompMod
