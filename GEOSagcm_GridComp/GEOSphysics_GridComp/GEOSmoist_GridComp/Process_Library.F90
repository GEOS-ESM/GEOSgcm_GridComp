! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

module GEOSmoist_Process_Library

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use Aer_Actv_Single_Moment
  use aer_cloud

  implicit none
  private

  character(len=ESMF_MAXSTR)              :: IAm="GEOSmoist_Process_Library"
  integer                                 :: STATUS

  interface MELTFRZ
    module procedure MELTFRZ_3D
    module procedure MELTFRZ_2D
    module procedure MELTFRZ_1D
    module procedure MELTFRZ_SC
  end interface MELTFRZ
  interface ICE_FRACTION
    module procedure ICE_FRACTION_3D
    module procedure ICE_FRACTION_2D
    module procedure ICE_FRACTION_1D
    module procedure ICE_FRACTION_SC
  end interface ICE_FRACTION
  ! ICE_FRACTION constants
   ! In anvil/convective clouds
   real, parameter :: aT_ICE_ALL = 252.16
   real, parameter :: aT_ICE_MAX = 268.16
   real, parameter :: aICEFRPWR  = 2.0
   ! Over snow SRF_TYPE = 2 and over ice SRF_TYPE = 3
   real, parameter :: iT_ICE_ALL = 236.16
   real, parameter :: iT_ICE_MAX = 261.16
   real, parameter :: iICEFRPWR  = 5.0
   ! Over Land     SRF_TYPE = 1
   real, parameter :: lT_ICE_ALL = 239.16
   real, parameter :: lT_ICE_MAX = 261.16
   real, parameter :: lICEFRPWR  = 2.0
   ! Over Oceans   SRF_TYPE = 0
   real, parameter :: oT_ICE_ALL = 238.16
   real, parameter :: oT_ICE_MAX = 263.16
   real, parameter :: oICEFRPWR  = 4.0
   ! Jason
   ! In anvil/convective clouds
   real, parameter :: JaT_ICE_ALL = 245.16
   real, parameter :: JaT_ICE_MAX = 261.16
   real, parameter :: JaICEFRPWR  = 2.0
        ! Over snow/ice
   real, parameter :: JiT_ICE_ALL = MAPL_TICE-40.0
   real, parameter :: JiT_ICE_MAX = MAPL_TICE
   real, parameter :: JiICEFRPWR  = 4.0

 ! parameters
  real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
  real, parameter :: K_COND  =  2.4e-2    ! J m**-1 s**-1 K**-1
  real, parameter :: DIFFU   =  2.2e-5    ! m**2 s**-1
  real, parameter :: taufrz  =  450.0
  real, parameter :: dQCmax  =  1.e-4
  ! LDRADIUS4
  ! Jason
  real, parameter :: abeta = 0.07
  real, parameter :: r13bbeta = 1./3. - 0.14
  real, parameter :: bx = 100.* (3./(4.*MAPL_PI))**(1./3.)
  ! Liquid  based on DOI 10.1088/1748-9326/3/4/045021
  real, parameter :: RHO_W   = 1000.0  ! Density of liquid water in kg/m^3
  real, parameter :: Ldiss   = 0.07    ! tunable dispersion effect
  real, parameter :: Lk      = 0.75    ! tunable shape effect (0.5:1)
  real, parameter :: Lbe     = 1./3. - 0.14
  real, parameter :: Lbx     = Ldiss*1.e3*(3./(4.*MAPL_PI*Lk*RHO_W*1.e-3))**(1./3.)
                             ! LDRADIUS eqs are in cgs units
  ! Ice
  real, parameter :: RHO_I   =  916.8  ! Density of ice crystal in kg/m^3

  ! combined constants
  real, parameter :: cpbgrav = MAPL_CP/MAPL_GRAV
  real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP
  real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
  real, parameter :: alhfbcp = MAPL_ALHF/MAPL_CP
  real, parameter :: alhsbcp = MAPL_ALHS/MAPL_CP

  ! base grid length for sigma calculation
  real :: SIGMA_DX  = 500.0
  real :: SIGMA_EXP = 1.0

  ! control for order of plumes
  logical :: SH_MD_DP = .FALSE.

  ! Radar parameter
  integer :: DBZ_VAR_INTERCP=1 ! use variable intercept parameters 
  integer :: DBZ_LIQUID_SKIN=1 ! use liquid skin on snow/ice in warm environments

  ! option for cloud liq/ice radii
  integer :: LIQ_RADII_PARAM = 1
  integer :: ICE_RADII_PARAM = 1

  ! defined to determine CNV_FRACTION
  real    :: CNV_FRACTION_MIN
  real    :: CNV_FRACTION_MAX
  real    :: CNV_FRACTION_EXP

  ! Storage of aerosol properties for activation
  type(AerPropsNew) :: AeroPropsNew(nsmx_par)
  type(AerProps), allocatable, dimension (:,:,:) :: AeroProps

  ! Tracer Bundle things for convection
  type CNV_Tracer_Type
      real, pointer              :: Q(:,:,:) => null()
      real                       :: fscav = 0.0
      real                       :: Vect_Hcts(4)
      real                       :: KcScal(3)
      real                       :: convfaci2g
      real                       :: retfactor
      real                       :: liq_and_gas
      real                       :: online_cldliq
      real                       :: online_vud
      real                       :: ftemp_threshold
      logical                    :: use_gcc_washout
      logical                    :: use_gocart
      logical                    :: is_wetdep
      character(len=ESMF_MAXSTR) :: QNAME ! Tracer Name
      character(len=ESMF_MAXSTR) :: CNAME ! Component Name
  end type CNV_Tracer_Type
  type(CNV_Tracer_Type), allocatable :: CNV_Tracers(:)

  public :: AeroProps
  public :: AeroPropsNew
  public :: CNV_Tracer_Type, CNV_Tracers, CNV_Tracers_Init
  public :: ICE_FRACTION, EVAP3, SUBL3, LDRADIUS4, BUOYANCY, BUOYANCY2
  public :: REDISTRIBUTE_CLOUDS, RADCOUPLE, FIX_UP_CLOUDS
  public :: hystpdf, fix_up_clouds_2M
  public :: FILLQ2ZERO, FILLQ2ZERO1
  public :: MELTFRZ
  public :: DIAGNOSE_PRECIP_TYPE
  public :: VertInterp, cs_interpolator
  public :: find_l, FIND_EIS, FIND_KLCL
  public :: find_cldtop, find_cldbase, gw_prof
  public :: make_IceNumber, make_DropletNumber, make_RainNumber
  public :: dissipative_ke_heating
  public :: pdffrac, pdfcondensate, partition_dblgss
  public :: SIGMA_DX, SIGMA_EXP
  public :: CNV_FRACTION_MIN, CNV_FRACTION_MAX, CNV_FRACTION_EXP
  public :: SH_MD_DP, DBZ_VAR_INTERCP, DBZ_LIQUID_SKIN, LIQ_RADII_PARAM, ICE_RADII_PARAM
  public :: update_cld, meltfrz_inst2M
  public :: FIX_NEGATIVE_PRECIP
  public :: FIND_KLID
  public :: sigma
  public :: pdf_alpha

  contains

  subroutine CNV_Tracers_Init(TR, RC)
    type (ESMF_FieldBundle), intent(inout) :: TR
    integer,       optional, intent(inout) :: RC
   ! Local
    type (ESMF_Field) :: FIELD
    integer :: TotalTracers, FriendlyTracers
    logical :: isPresent, isFriendly
    integer :: ind, N, F
    real    :: rtmp
    character(len=ESMF_MAXSTR), pointer, dimension(:) :: QNAMES
    character(len=ESMF_MAXSTR) :: QNAME, STR_CNV_TRACER

    call ESMF_FieldBundleGet(TR, FieldCount=TotalTracers, RC=STATUS); VERIFY_(STATUS)
    allocate(QNAMES(TotalTracers), stat=STATUS); VERIFY_(STATUS)
    call ESMF_FieldBundleGet(TR, fieldNameList=QNAMES, RC=STATUS); VERIFY_(STATUS)
    FriendlyTracers = 0
    do N=1,TotalTracers
       QNAME = trim(QNAMES(N))
       call ESMF_FieldBundleGet(TR, fieldName=trim(QNAME), Field=FIELD, RC=STATUS); VERIFY_(STATUS)
       call ESMF_AttributeGet  (FIELD, "FriendlyToMOIST",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
       if(isPresent) then
          call ESMF_AttributeGet(FIELD, "FriendlyToMOIST", isFriendly, RC=STATUS); VERIFY_(STATUS)
          if (isFriendly) FriendlyTracers = FriendlyTracers + 1
       end if
    enddo

    ! see if we need to allocate
    if (allocated(CNV_Tracers)) then
      ASSERT_( size(CNV_Tracers) == FriendlyTracers )
    else
      call WRITE_PARALLEL ("List of species friendly to MoistGridComp:")
      ! fill CNV_Tracers
      allocate( CNV_Tracers(FriendlyTracers), stat=STATUS); VERIFY_(STATUS)
      F = 0
      do N=1,TotalTracers
         QNAME = trim(QNAMES(N))
         call ESMF_FieldBundleGet(TR, fieldName=trim(QNAME), Field=FIELD, RC=STATUS); VERIFY_(STATUS)
         call ESMF_AttributeGet  (FIELD, "FriendlyToMOIST",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
         if(isPresent) then
            call ESMF_AttributeGet(FIELD, "FriendlyToMOIST", isFriendly, RC=STATUS); VERIFY_(STATUS)
            if (isFriendly) then
               ! Iterate the friendly index
               !-------------------------------
               F = F + 1
               ! Get items scavenging fraction
               !-------------------------------
               CNV_Tracers(F)%fscav = 0.0
               call ESMF_AttributeGet(FIELD, "ScavengingFractionPerKm", isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
               if(isPresent) then
                  call ESMF_AttributeGet(FIELD, "ScavengingFractionPerKm", CNV_Tracers(F)%fscav, RC=STATUS); VERIFY_(STATUS)
               end if
               ! Get component and tracer names
               !-------------------------------------------------------------------------------------
               ind= index(QNAME, '::')
               if (ind > 0) then
                  CNV_Tracers(F)%CNAME = trim(QNAME(1:ind-1))  ! Component name (e.g., GOCART, CARMA)
                  CNV_Tracers(F)%QNAME = trim(QNAME(ind+2:))
               end if
               ! Get items for the wet removal parameterization for gases based on the Henry's Law
               !-------------------------------------------------------------------------------------
               CNV_Tracers(F)%Vect_Hcts(:)=-99.
               call ESMF_AttributeGet(FIELD, "SetofHenryLawCts", isPresent=isPresent,  RC=STATUS); VERIFY_(STATUS)
               if (isPresent) then
                  call ESMF_AttributeGet(FIELD, "SetofHenryLawCts", CNV_Tracers(F)%Vect_Hcts,  RC=STATUS); VERIFY_(STATUS)
               end if
               ! Additional items, needed for GEOS-Chem washout parameterization
               !-------------------------------------------------------------------------------------
               ! Defaults
               CNV_Tracers(F)%is_wetdep       = .FALSE.
               CNV_Tracers(F)%use_gcc_washout = .FALSE.
               CNV_Tracers(F)%KcScal(:)       = 1.0
               CNV_Tracers(F)%retfactor       = 1.0
               CNV_Tracers(F)%liq_and_gas     = 0.0
               CNV_Tracers(F)%convfaci2g      = 0.0
               CNV_Tracers(F)%online_cldliq   = 0.0
               CNV_Tracers(F)%online_vud      = 1.0
               CNV_Tracers(F)%use_gocart      = .FALSE.
               CNV_Tracers(F)%ftemp_threshold = -999.0
               ! Check if GEOS-Chem washout should be used. Assume this is the case if Kc scale factors are
               ! present
               call ESMF_AttributeGet(FIELD, "SetofKcScalFactors", isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
               CNV_Tracers(F)%use_gcc_washout = isPresent
               ! If using GEOS-Chem parameterization, retrieve all necessary parameter
               if ( CNV_Tracers(F)%use_gcc_washout ) then
                  ! KC scale factors
                  call ESMF_AttributeGet(FIELD, "SetofKcScalFactors", isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet(FIELD, "SetofKcScalFactors", CNV_Tracers(F)%KcScal, RC=STATUS); VERIFY_(STATUS)
                  end if
                  ! is this a wetdep species?
                  call ESMF_AttributeGet(FIELD, "IsWetDep", isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet(FIELD, "IsWetDep", rtmp, RC=STATUS); VERIFY_(STATUS)
                     CNV_Tracers(F)%is_wetdep = ( rtmp == 1.0 )
                  end if
                  ! Gas-phase washout parameter for GEOS-Chem
                  call ESMF_AttributeGet (FIELD, "RetentionFactor",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD, "RetentionFactor", CNV_Tracers(F)%retfactor, RC=STATUS); VERIFY_(STATUS)
                  endif
                  call ESMF_AttributeGet (FIELD, "LiqAndGas",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD, "LiqAndGas", CNV_Tracers(F)%liq_and_gas, RC=STATUS); VERIFY_(STATUS)
                  endif
                  call ESMF_AttributeGet (FIELD, "ConvFacI2G",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD, "ConvFacI2G", CNV_Tracers(F)%convfaci2g, RC=STATUS); VERIFY_(STATUS)
                  endif
                  call ESMF_AttributeGet (FIELD, "OnlineCLDLIQ",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD, "OnlineCLDLIQ", CNV_Tracers(F)%online_cldliq, RC=STATUS); VERIFY_(STATUS)
                  endif
                  call ESMF_AttributeGet (FIELD,"OnlineVUD",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD,"OnlineVUD", CNV_Tracers(F)%online_vud, RC=STATUS); VERIFY_(STATUS)
                  endif
                  call ESMF_AttributeGet (FIELD,"UseGOCART",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD,"UseGOCART", rtmp, RC=STATUS); VERIFY_(STATUS)
                     CNV_Tracers(F)%use_gocart = ( rtmp == 1.0 )
                  endif
                  call ESMF_AttributeGet (FIELD,"GOCARTfTempThreshold",isPresent=isPresent, RC=STATUS); VERIFY_(STATUS)
                  if (isPresent) then
                     call ESMF_AttributeGet (FIELD,"GOCARTfTempThreshold", CNV_Tracers(F)%ftemp_threshold, RC=STATUS); VERIFY_(STATUS)
                  endif
               end if ! use_gcc_washout
               ! Get pointer to friendly tracers
               !-----------------------------------------
               call ESMFL_BundleGetPointerToData(TR, trim(QNAME), CNV_Tracers(F)%Q, RC=STATUS); VERIFY_(STATUS)
               ! Report tracer status
               !-----------------------------------------
               if (CNV_Tracers(F)%fscav > 1.e-6) then
                   WRITE(STR_CNV_TRACER,101) TRIM(QNAME), CNV_Tracers(F)%fscav
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
               elseif (CNV_Tracers(F)%Vect_Hcts(1)>1.e-6) then
                   WRITE(STR_CNV_TRACER,102) TRIM(QNAME), CNV_Tracers(F)%Vect_Hcts
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
               else
                   WRITE(STR_CNV_TRACER,103) TRIM(QNAME)
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
               endif
101            FORMAT(a,' ScavengingFractionPerKm:',1(1x,f3.1))
102            FORMAT(a,' SetofHenryLawCts:',4(1x,es9.2))
103            FORMAT(a,' is transported by Moist')
               ! Additional information for GEOS-Chem washout species
               !-----------------------------------------------------
               if (CNV_Tracers(F)%use_gcc_washout .and. CNV_Tracers(F)%is_wetdep) then
                   STR_CNV_TRACER = TRIM(QNAME)//": will use GEOS-Chem washout formulation"
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   WRITE(STR_CNV_TRACER,104) TRIM(QNAME), CNV_Tracers(F)%KcScal
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   WRITE(STR_CNV_TRACER,105) TRIM(QNAME), CNV_Tracers(F)%retfactor
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   WRITE(STR_CNV_TRACER,106) TRIM(QNAME), CNV_Tracers(F)%liq_and_gas
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   WRITE(STR_CNV_TRACER,107) TRIM(QNAME), CNV_Tracers(F)%convfaci2g
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   WRITE(STR_CNV_TRACER,108) TRIM(QNAME), CNV_Tracers(F)%online_cldliq
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   WRITE(STR_CNV_TRACER,109) TRIM(QNAME), CNV_Tracers(F)%online_vud
                   call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   if (CNV_Tracers(F)%use_gocart)then
                       STR_CNV_TRACER = TRIM(QNAME)//": will treat like GOCART aerosol"
                       call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                       WRITE(STR_CNV_TRACER,110) TRIM(QNAME), CNV_Tracers(F)%ftemp_threshold
                       call WRITE_PARALLEL (trim(STR_CNV_TRACER))
                   endif
               endif
104            FORMAT(a,' KcScaleFactors:',3(1x,es9.2))
105            FORMAT(a,' RetentionFactor:',1(1x,es9.2))
106            FORMAT(a,' Liq_and_gas:',1(1x,es9.2))
107            FORMAT(a,' ConvFacI2G:',1(1x,es9.2))
108            FORMAT(a,' online_cldliq:',1(1x,es9.2))
109            FORMAT(a,' online_vud:',1(1x,es9.2))
110            FORMAT(a,' ftemp_threshold:',1(1x,es9.2))
            end if
         end if
      enddo
    end if

    deallocate(QNAMES)

  end subroutine CNV_Tracers_Init

  real function sigma (dx, BASE_DX, BASE_EXP)
      real, intent(in) :: dx
      real, optional , intent(in) :: BASE_DX, BASE_EXP
      real                :: tmp_exp
                             tmp_exp = SIGMA_EXP
      if (present(BASE_EXP)) tmp_exp = BASE_EXP
     ! Arakawa 2011 based sigma function
      if (present(BASE_DX)) then
        sigma = (1.0-0.9839*exp(-0.09835*(dx/ BASE_DX)))**tmp_exp
      else
        sigma = (1.0-0.9839*exp(-0.09835*(dx/SIGMA_DX)))**tmp_exp
      endif
  end function sigma

  function ICE_FRACTION_3D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
      real, intent(in) :: TEMP(:,:,:),CNV_FRACTION(:,:),SRF_TYPE(:,:)
      real :: ICEFRCT(size(TEMP,1),size(TEMP,2),size(TEMP,3))
      integer :: i,j,l
      do l=1,size(TEMP,3)
      do j=1,size(TEMP,2)
      do i=1,size(TEMP,1)
        ICEFRCT(i,j,l) = ICE_FRACTION_SC(TEMP(i,j,l),CNV_FRACTION(i,j),SRF_TYPE(i,j))
      enddo
      enddo
      enddo
  end function ICE_FRACTION_3D

  function ICE_FRACTION_2D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
      real, intent(in) :: TEMP(:,:),CNV_FRACTION(:,:),SRF_TYPE(:,:)
      real :: ICEFRCT(size(TEMP,1),size(TEMP,2))
      integer :: i,j
      do j=1,size(TEMP,2)
      do i=1,size(TEMP,1)
        ICEFRCT(i,j) = ICE_FRACTION_SC(TEMP(i,j),CNV_FRACTION(i,j),SRF_TYPE(i,j))
      enddo
      enddo
  end function ICE_FRACTION_2D

  function ICE_FRACTION_1D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
      real, intent(in) :: TEMP(:),CNV_FRACTION(:),SRF_TYPE(:)
      real :: ICEFRCT(size(TEMP))
      integer :: i
      do i=1,size(TEMP)
        ICEFRCT(i) = ICE_FRACTION_SC(TEMP(i),CNV_FRACTION(i),SRF_TYPE(i))
      enddo
  end function ICE_FRACTION_1D

  function ICE_FRACTION_SC (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
      real, intent(in) :: TEMP,CNV_FRACTION,SRF_TYPE
      real             :: ICEFRCT
      real             :: tc, ptc
      real             :: ICEFRCT_C, ICEFRCT_M

#ifdef USE_MODIS_ICE_POLY
     ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
      tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
      ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
      ICEFRCT = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
#else
     ! Use sigmoidal functions based on surface type from Hu et al, DOI: (10.1029/2009JD012384)
     ! Anvil clouds
     ! Anvil-Convective sigmoidal function like figure 6(right)
     ! Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
      if (ICE_RADII_PARAM == 1) then
        ! Jason formula
        ICEFRCT_C  = 0.00
        if ( TEMP <= JaT_ICE_ALL ) then
           ICEFRCT_C = 1.000
        else if ( (TEMP > JaT_ICE_ALL) .AND. (TEMP <= JaT_ICE_MAX) ) then
           ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - JaT_ICE_ALL ) / ( JaT_ICE_MAX - JaT_ICE_ALL ) ) )
        end if
      else
        ICEFRCT_C  = 0.00
        if ( TEMP <= aT_ICE_ALL ) then
           ICEFRCT_C = 1.000
        else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
           ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
        end if
      end if
      ICEFRCT_C = MIN(ICEFRCT_C,1.00)
      ICEFRCT_C = MAX(ICEFRCT_C,0.00)
      ICEFRCT_C = ICEFRCT_C**aICEFRPWR
     ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
      if (SRF_TYPE >= 2.0) then
        ! Over snow (SRF_TYPE == 2.0) and ice (SRF_TYPE == 3.0)
        if (ICE_RADII_PARAM == 1) then
          ! Jason formula
          ICEFRCT_M  = 0.00
          if ( TEMP <= JiT_ICE_ALL ) then
             ICEFRCT_M = 1.000
          else if ( (TEMP > JiT_ICE_ALL) .AND. (TEMP <= JiT_ICE_MAX) ) then
             ICEFRCT_M = 1.00 -  ( TEMP - JiT_ICE_ALL ) / ( JiT_ICE_MAX - JiT_ICE_ALL )
          end if
        else
          ICEFRCT_M  = 0.00
          if ( TEMP <= iT_ICE_ALL ) then
             ICEFRCT_M = 1.000
          else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
             ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL ) ) )
          end if
        end if
        ICEFRCT_M = MIN(ICEFRCT_M,1.00)
        ICEFRCT_M = MAX(ICEFRCT_M,0.00)
        ICEFRCT_M = ICEFRCT_M**iICEFRPWR
      else if (SRF_TYPE > 1.0) then
        ! Over Land
        ICEFRCT_M  = 0.00
        if ( TEMP <= lT_ICE_ALL ) then
           ICEFRCT_M = 1.000
        else if ( (TEMP > lT_ICE_ALL) .AND. (TEMP <= lT_ICE_MAX) ) then
           ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - lT_ICE_ALL ) / ( lT_ICE_MAX - lT_ICE_ALL ) ) )
        end if
        ICEFRCT_M = MIN(ICEFRCT_M,1.00)
        ICEFRCT_M = MAX(ICEFRCT_M,0.00)
        ICEFRCT_M = ICEFRCT_M**lICEFRPWR
      else
        ! Over Oceans
        ICEFRCT_M  = 0.00
        if ( TEMP <= oT_ICE_ALL ) then
           ICEFRCT_M = 1.000
        else if ( (TEMP > oT_ICE_ALL) .AND. (TEMP <= oT_ICE_MAX) ) then
           ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - oT_ICE_ALL ) / ( oT_ICE_MAX - oT_ICE_ALL ) ) )
        end if
        ICEFRCT_M = MIN(ICEFRCT_M,1.00)
        ICEFRCT_M = MAX(ICEFRCT_M,0.00)
        ICEFRCT_M = ICEFRCT_M**oICEFRPWR
      endif
      ! Combine the Convective and MODIS functions
        ICEFRCT  = ICEFRCT_M*(1.0-CNV_FRACTION) + ICEFRCT_C*(CNV_FRACTION)
#endif

  end function ICE_FRACTION_SC

   subroutine EVAP3(&
         DT      , &
         A_EFF   , &
         RHCR    , &
         PL      , &
         TE      , &
         QV      , &
         QL      , &
         QI      , &
         F       , &
         NL      , &
         NI      , &
         QS        )

      real, intent(in   ) :: DT
      real, intent(in   ) :: A_EFF
      real, intent(in   ) :: RHCR
      real, intent(in   ) :: PL
      real, intent(inout) :: TE
      real, intent(inout) :: QV
      real, intent(inout) :: QL,QI
      real, intent(inout) :: F
      real, intent(in   ) :: NL,NI
      real, intent(in   ) :: QS

      real :: ES,RADIUS,K1,K2,QCm,EVAP,RHx,QC

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         EVAPORATION OF CLOUD WATER.             !!
      !!                                                 !!
      !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
      !!  formulation  (Eq.s 15-17)                      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )

      ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100's <-^ convert from mbar to Pa)

      RHx = MIN( QV/QS , 1.00 )

      K1 = (MAPL_ALHL**2) * RHO_W / ( K_COND*MAPL_RVAP*(TE**2))

      K2 = MAPL_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Here DIFFU is given for 1000 mb  !!
      !! so 1000./PR accounts for inc-    !!
      !! reased diffusivity at lower      !!
      !! pressure.                        !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( ( F > 0.) .and. ( QL > 0. ) ) then
         QCm=QL/F
      else
         QCm=0.
      end if

      RADIUS = LDRADIUS4(PL,TE,QCm,NL,NI,1)

      if ( (RHx < RHCR ) .and. (RADIUS > 0.0) ) then
         EVAP = A_EFF*QL*DT*(RHCR - RHx) / ((K1+K2)*RADIUS**2)
         EVAP = MAX(0.0, MIN( EVAP , QL  ))
      else
         EVAP = 0.0
      end if

      QC=QL+QI
      if (QC > 0.) then
         F = F * ( QC - EVAP ) / QC
      end if

      QV = QV + EVAP
      QL = QL - EVAP
      TE = TE - alhlbcp*EVAP

   end subroutine EVAP3

   subroutine SUBL3( &
         DT        , &
         A_EFF     , &
         RHCR      , &
         PL        , &
         TE        , &
         QV        , &
         QL        , &
         QI        , &
         F         , &
         NL        , &
         NI        , &
         QS        )

      real, intent(in   ) :: DT
      real, intent(in   ) :: A_EFF
      real, intent(in   ) :: RHCR
      real, intent(in   ) :: PL
      real, intent(inout) :: TE
      real, intent(inout) :: QV
      real, intent(inout) :: QL,QI
      real, intent(inout) :: F
      real, intent(in   ) :: NL,NI
      real, intent(in   ) :: QS

      real :: ES,RADIUS,K1,K2,TEFF,QCm,SUBL,RHx,QC

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         SUBLORATION OF CLOUD WATER.             !!
      !!                                                 !!
      !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
      !!  formulation  (Eq.s 15-17)                      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )

      ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100s <-^ convert from mbar to Pa)

      RHx = MIN( QV/QS , 1.00 )

      K1 = (MAPL_ALHL**2) * RHO_I / ( K_COND*MAPL_RVAP*(TE**2))

      K2 = MAPL_RVAP * TE * RHO_I / ( DIFFU * (1000./PL) * ES )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Here DIFFU is given for 1000 mb  !!
      !! so 1000./PR accounts for inc-    !!
      !! reased diffusivity at lower      !!
      !! pressure.                        !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( ( F > 0.) .and. ( QI > 0. ) ) then
         QCm=QI/F
      else
         QCm=0.
      end if

      RADIUS = LDRADIUS4(PL,TE,QCm,NL,NI,2)

      if ( (RHx < RHCR) .and.(RADIUS > 0.0) ) then
         SUBL = A_EFF*QI*DT*(RHCR - RHx) / ((K1+K2)*RADIUS**2)
         SUBL = MAX(0.0, MIN( SUBL , QI  ))
      else
         SUBL = 0.0
      end if

      QC=QL+QI
      if (QC > 0.) then
         F = F * ( QC - SUBL ) / QC
      end if

      QV = QV + SUBL
      QI = QI - SUBL
      TE = TE - alhsbcp*SUBL

   end subroutine SUBL3

   function LDRADIUS4(PL,TE,QC,NNL,NNI,ITYPE) RESULT(RADIUS)

       REAL   , INTENT(IN) :: TE,PL,QC
       REAL   , INTENT(IN) :: NNL,NNI ! #/m^3
       INTEGER, INTENT(IN) :: ITYPE
       REAL  :: RADIUS
       INTEGER, PARAMETER  :: LIQUID=1, ICE=2
       REAL :: NNX,RHO,BB,WC
       REAL :: TC,AA

       !- air density (kg/m^3)
       RHO = (100.*PL) / (MAPL_RGAS*TE )
       IF(ITYPE == LIQUID) THEN

       !- liquid cloud effective radius -----
          !- liquid water content
          WC = 1.e3*RHO*QC ! air density [g/m3] * liquid cloud mixing ratio [kg/kg]
          !- cloud drop number concentration
          !- from the aerosol model + ....
          NNX = max(NNL*1.e-6, 10.0)
          !- radius in meters
          if (LIQ_RADII_PARAM == 1) then
            !- Jason Version
            RADIUS= MIN(60.e-6,MAX(2.5e-6, 1.e-6*bx*(WC/NNX)**r13bbeta*abeta*6.92))
          else
            !- [liu&daum, 2000 and 2005. liu et al 2008]
            RADIUS = MIN(60.e-6,MAX(2.5e-6, 1.e-6*Lbx*(WC/NNX)**Lbe))
          endif

       ELSEIF(ITYPE == ICE) THEN

       !- ice cloud effective radius -----
          !- ice water content
          WC = 1.e3*RHO*QC ! air density [g/m3] * ice cloud mixing ratio [kg/kg]
          !- radius in meters
          if (ICE_RADII_PARAM == 1) then
            !------ice cloud effective radius ----- [klaus wyser, 1998]
            if(TE>MAPL_TICE .or. QC < 1.e-9) then
              BB = -2.
            else
              BB = -2. + log10(WC/50.)*(1.e-3*(MAPL_TICE-TE)**1.5)
            endif
            BB     = MIN((MAX(BB,-6.)),-2.)
            RADIUS = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
            RADIUS = MIN(150.e-6,MAX(5.e-6, 1.e-6*RADIUS))
          else
            !------ice cloud effective radius ----- [Sun, 2001]
            ! https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GL102521
            TC = TE - MAPL_TICE
            AA = 45.8966 * (WC**0.2214)
            BB = 0.79570 * (WC**0.2535) * (TE - 83.15)
            RADIUS = MIN(155.0  ,MAX(30.0  , (1.2351 + 0.0105*TC) * (AA + BB)))
            RADIUS = MIN(150.e-6,MAX( 5.e-6, 1.e-6*0.64952*RADIUS))
          endif

      ELSE
        STOP "WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2"
      ENDIF

   end function LDRADIUS4


  subroutine BUOYANCY2( IM, JM, LM, T, Q, QS, DQS, DZ, ZLO, PLO, PS, SBCAPE, MLCAPE, MUCAPE, &
                        SBCIN, MLCIN, MUCIN, BYNCY, LFC, LNB )

    ! Computes surface-based (SB), mixed-layer (ML) and most unstable (MU) versions
    ! of CAPE and CIN.

    integer,                intent(in)  :: IM, JM, LM
    real, dimension(:,:,:), intent(in)  :: T, Q, QS, DQS, DZ, ZLO, PLO
    real, dimension(:,:,:), intent(out) :: BYNCY
    real, pointer, dimension(:,:)       :: MLCAPE, MUCAPE, MLCIN, MUCIN
    real, dimension(:,:)                :: SBCAPE, SBCIN, LFC, LNB
    real, dimension(:,:),   intent(in)  :: PS

    real, dimension(IM,JM,LM)           :: Tve
    real, dimension(IM,JM)              :: MSEp, Qp, tmp1, tmp2
    integer, dimension(IM,JM)           :: Lev0

    integer :: I, J, L

    Tve   = T*(1.+MAPL_VIREPS*Q)
    BYNCY = MAPL_UNDEF
    MSEp  = 0.
    Qp    = 0.

    ! Mixed-layer calculation. Parcel properties averaged over lowest 90 hPa
    if ( associated(MLCAPE) .and. associated(MLCIN) ) then
       BYNCY = MAPL_UNDEF
       tmp1 = 0.
       Lev0 = LM
       do L = LM,1,-1
         where (PS-PLO(:,:,L).lt.90.)
            MSEp = MSEp + (T(:,:,L) + gravbcp*ZLO(:,:,L) + alhlbcp*Q(:,:,L))*DZ(:,:,L)
            Qp   = Qp   + Q(:,:,L)*DZ(:,:,L)
            tmp1 = tmp1 + DZ(:,:,L)
            Lev0 = L
         end where
         if (all(PS-PLO(:,:,L).gt.90.)) exit
       end do
       where (tmp1.gt.0.)   ! average
          MSEp = MSEp / tmp1
          Qp = Qp / tmp1
       end where
       do I = 1,IM
          do J = 1,JM
             call RETURN_CAPE_CIN( ZLO(I,J,1:Lev0(I,J)), PLO(I,J,1:Lev0(I,J)), DZ(I,J,1:Lev0(I,J)),      &
                                   MSEp(I,J), Qp(I,J), Tve(I,J,1:Lev0(I,J)), QS(I,J,1:Lev0(I,J)), DQS(I,J,1:Lev0(I,J)),       &
                                   MLCAPE(I,J), MLCIN(I,J), BYNCY(I,J,1:Lev0(I,J)), LFC(I,J), LNB(I,J) )
          end do
       end do
       where (MLCAPE.le.0.)
          MLCAPE = MAPL_UNDEF
          MLCIN  = MAPL_UNDEF
       end where
    end if

    ! Most unstable calculation. Parcel in lowest 255 hPa with largest CAPE
    if ( associated(MUCAPE) .and. associated(MUCIN) ) then
       MUCAPE = 0.
       MUCIN  = 0.
       BYNCY = MAPL_UNDEF
       LFC = MAPL_UNDEF
       LNB = MAPL_UNDEF
       do I = 1,IM
          do J = 1,JM
             do L = LM,1,-1
                if (PS(I,J)-PLO(I,J,L).gt.255.) exit
                MSEp(I,J) = T(I,J,L) + gravbcp*ZLO(I,J,L) + alhlbcp*Q(I,J,L)
                Qp(I,J)   = Q(I,J,L)
                call RETURN_CAPE_CIN( ZLO(I,J,1:L), PLO(I,J,1:L), DZ(I,J,1:L),      &
                                      MSEp(I,J), Qp(I,J), Tve(I,J,1:L), QS(I,J,1:L), DQS(I,J,1:L),       &
                                      tmp1(I,J), tmp2(I,J), BYNCY(I,J,1:L), LFC(I,J), LNB(I,J) )
                if (tmp1(I,J) .gt. MUCAPE(I,J)) then
                   MUCAPE(I,J) = tmp1(I,J)
                   MUCIN(I,J)  = tmp2(I,J)
                end if
             end do
          end do
       end do

       where (MUCAPE.le.0.)
          MUCAPE = MAPL_UNDEF
          MUCIN  = MAPL_UNDEF
       end where
    end if

    ! Surface-based calculation
    MSEp = T(:,:,LM) + gravbcp*ZLO(:,:,LM) + alhlbcp*Q(:,:,LM)  ! parcel moist static energy
    Qp   = Q(:,:,LM)                                            ! parcel specific humidity
    do I = 1,IM
       do J = 1,JM
          call RETURN_CAPE_CIN( ZLO(I,J,:), PLO(I,J,:), DZ(I,J,:),      &
                                MSEp(I,J), Qp(I,J), Tve(I,J,:), QS(I,J,:), DQS(I,J,:),       &
                                SBCAPE(I,J), SBCIN(I,J), BYNCY(I,J,:), LFC(I,J), LNB(I,J) )
       end do
    end do
    where (SBCAPE.le.0.)
       SBCAPE = MAPL_UNDEF
       SBCIN  = MAPL_UNDEF
    end where

  end subroutine BUOYANCY2


  subroutine RETURN_CAPE_CIN( ZLO, PLO, DZ, MSEp, Qp, Tve, Qsate, DQS, CAPE, CIN, BYNCY, LFC, LNB )

    real,               intent(in)  :: MSEp, Qp
    real, dimension(:), intent(in)  :: ZLO, PLO, DZ, Tve, Qsate, DQS
    real,               intent(out) :: CAPE, CIN, LFC, LNB
    real, dimension(:), intent(out) :: BYNCY

    integer :: I, L, LM, KLNB, KLFC
    real    :: Qpnew, Tp, Tvp, Tlcl, Buoy, dq
    logical :: aboveLNB, aboveLFC, aboveLCL

    LM = size(ZLO,1)

    aboveLNB = .false.
    aboveLFC = .false.

    Qpnew = Qp

    CAPE = 0.
    CIN  = 0.
    BYNCY = 0.
    LFC = MAPL_UNDEF
    LNB = MAPL_UNDEF

    Tp = MSEp - gravbcp*ZLO(LM) - alhlbcp*Qp  ! initial parcel temp at source level LM
    Tlcl = find_tlcl( Tp, 100.*Qp/QSATE(LM) )
    aboveLCL = (Tp.lt.Tlcl)

    do L = LM-1,1,-1   ! start at level above source air

      ! determine parcel Qp, Tp
      if ( .not. aboveLCL ) then
         Tp = Tp - gravbcp*(ZLO(L)-ZLO(L+1))                ! new parcel temperature w/o condensation
         if (Tp.lt.Tlcl) then
            Tp = Tp + gravbcp*(ZLO(L)-ZLO(L+1))             ! if cross LCL, revert Tp and go to aboveLCL below
            aboveLCL = .true.
         end if
      end if
      if ( aboveLCL .and. Qpnew*alhlbcp.gt.0.01 ) then
         Tp = Tp - gravbcp*( ZLO(L)-ZLO(L+1) ) / ( 1.+alhlbcp*DQS(L) )     ! initial guess including condensation
         DO I = 1,10                                                       ! iterate until Qp=qsat(Tp)
            dq = Qpnew - GEOS_QSAT( Tp, PLO(L) )
            if (abs(dq*alhlbcp)<0.01) then
               exit
            end if
            Tp = Tp + dq*alhlbcp/(1.+alhlbcp*DQS(L))
            Qpnew = Qpnew - dq/(1.+alhlbcp*DQS(L))
         END DO
      end if
      Tp = MSEp - gravbcp*ZLO(L) - alhlbcp*Qpnew
      !  Qc = qp - qpnew.   ! condensate (not used for pseudoadiabatic ascent)

      Tvp = Tp*(1.+MAPL_VIREPS*Qpnew)              ! parcel virtual temp
    !  Tvp = Tp*(1.+0.61*Qpnew - Qc) ! condensate loading

      BYNCY(L) = MAPL_GRAV*(Tvp-Tve(L))/Tve(L)         ! parcel buoyancy

    end do

    ! if surface parcel immediately buoyant, scan upward to find first elevated
    ! B>0 level above a B<0 level, label it LFC.  If no such level, set LFC at surface.
    KLFC = LM
    KLNB = LM
    aboveLFC = .false.
    if (BYNCY(LM-1).gt.0.) then
       do L = LM-2,1,-1   ! scan up to find elevated LFC
          if (BYNCY(L).gt.0. .and. BYNCY(L+1).le.0.) then
             KLFC = L
             aboveLFC = .true.
          end if
          if (aboveLFC .and. BYNCY(L).lt.0. ) then
             KLNB = L
             exit
          end if
       end do
    else   ! if surface parcel not immediately buoyant, LFC is first B>0 level
       do L = LM-1,1,-1
          if (BYNCY(L).gt.0. .and. .not.aboveLFC) then
             KLFC = L
             aboveLFC = .true.
          end if
          if (aboveLFC .and. BYNCY(L).lt.0.) then
             KLNB = L
             exit
          end if
       end do
    end if
    LFC = ZLO(KLFC)
    LNB = ZLO(KLNB)

    CIN = SUM( min(0.,BYNCY(KLFC:)*DZ(KLFC:)) )        ! define CIN as negative
!    CAPE = SUM( max(0.,BYNCY(KLNB:KLFC)*DZ(KLNB:KLFC)) )
    CAPE = SUM( max(0.,BYNCY(:)*DZ(:)) )

  end subroutine RETURN_CAPE_CIN


  subroutine BUOYANCY( T, Q, QS, DQS, DZ, ZLO, BUOY, CAPE, INHB)


    ! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
    !  for a parcel raised from the surface. $T_c$ is the virtual temperature of
    !  the parcel and $T_e$ is the virtual temperature of the environment.

    real, dimension(:,:,:),   intent(in)  :: T, Q, QS, DQS, DZ, ZLO
    real, dimension(:,:,:),   intent(out) :: BUOY
    real, dimension(:,:),     intent(out) :: CAPE, INHB

    integer :: L, LM

    LM = size(T,3)

    BUOY(:,:,LM) =  T(:,:,LM) + gravbcp*ZLO(:,:,LM) + alhlbcp*Q(:,:,LM)

    do L=LM-1,1,-1
       BUOY(:,:,L) = BUOY(:,:,LM) - (T(:,:,L) + gravbcp*ZLO(:,:,L) + alhlbcp*QS(:,:,L))
       BUOY(:,:,L) = MAPL_GRAV*BUOY(:,:,L) / ( (1.+ alhlbcp*DQS(:,:,L))*T(:,:,L) )
    enddo

    BUOY(:,:,LM) = 0.0

    CAPE = 0.
    INHB = 0.

    do L=1,LM-1
       where(BUOY(:,:,L)>0.)
          CAPE = CAPE + BUOY(:,:,L)*DZ(:,:,L)
       end where
       where(BUOY(:,:,L)<0.)
          INHB = INHB - BUOY(:,:,L)*DZ(:,:,L)
       end where
    end do

    where(CAPE <= 0.0)
       CAPE=MAPL_UNDEF
       INHB=MAPL_UNDEF
    end where

  end subroutine BUOYANCY

   subroutine RADCOUPLE(  &
         TE,              &
         PL,              &
         CF,              &
         AF,              &
         QV,              &
         QClLS,           &
         QCiLS,           &
         QClAN,           &
         QCiAN,           &
         QRN_ALL,         &
         QSN_ALL,         &
         QGR_ALL,         &
         NL,              &
         NI,              &
         RAD_QV,          &
         RAD_QL,          &
         RAD_QI,          &
         RAD_QR,          &
         RAD_QS,          &
         RAD_QG,          &
         RAD_CF,          &
         RAD_RL,          &
         RAD_RI,          &
         FAC_RL, MIN_RL, MAX_RL, &
         FAC_RI, MIN_RI, MAX_RI)

      real, intent(in ) :: TE
      real, intent(in ) :: PL
      real, intent(in ) :: AF,CF, QV, QClAN, QCiAN, QClLS, QCiLS
      real, intent(in ) :: QRN_ALL, QSN_ALL, QGR_ALL
      real, intent(in ) :: NL,NI
      real, intent(out) :: RAD_QV,RAD_QL,RAD_QI,RAD_QR,RAD_QS,RAD_QG,RAD_CF,RAD_RL,RAD_RI
      real, intent(in ) :: FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI

      ! Limits on Radii needed to ensure
      ! correct behavior of cloud optical
      ! properties currently calculated in
      ! sorad and irrad (1e-6 m = micron)

      ! water vapor
      RAD_QV = QV

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Total cloud fraction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RAD_CF = MAX(MIN(CF+AF,1.0),0.0)
      if ( RAD_CF >= 1.e-5 ) then
        ! Total In-cloud liquid
        if ( (QClLS + QClAN) >= 1.e-8 ) then
           RAD_QL = ( QClLS + QClAN ) / RAD_CF
        else
           RAD_QL = 0.0
        end if
        ! Total In-cloud ice
        if ( (QCiLS + QCiAN) >= 1.e-8 ) then
           RAD_QI = ( QCiLS + QCiAN ) / RAD_CF
        else
           RAD_QI = 0.0
        end if
        ! Total In-cloud precipitation
        if (QRN_ALL >= 1.e-8 ) then
           RAD_QR = ( QRN_ALL ) / RAD_CF
        else
           RAD_QR = 0.0
        end if
        if (QSN_ALL >= 1.e-8 ) then
           RAD_QS = ( QSN_ALL ) / RAD_CF
        else
           RAD_QS = 0.0
        end if
        if (QGR_ALL >= 1.e-8 ) then
           RAD_QG = ( QGR_ALL ) / RAD_CF
        else
           RAD_QG = 0.0
        end if
      else
        RAD_CF = 0.0
        RAD_QL = 0.0
        RAD_QI = 0.0
        RAD_QR = 0.0
        RAD_QS = 0.0
        RAD_QG = 0.0
      end if
     ! cap the high end of condensates
      RAD_QL = MIN( RAD_QL, 0.01 )
      RAD_QI = MIN( RAD_QI, 0.01 )
      RAD_QR = MIN( RAD_QR, 0.01 )
      RAD_QS = MIN( RAD_QS, 0.01 )
      RAD_QG = MIN( RAD_QG, 0.01 )

     ! LIQUID RADII
      if (RAD_QL > 1.e-8) then
        !-BRAMS formulation
        RAD_RL = LDRADIUS4(PL,TE,RAD_QL,NL,NI,1)
        ! apply limits
        RAD_RL = MAX( MIN_RL, MIN(RAD_RL*FAC_RL, MAX_RL) )
      else
        RAD_RL = MAPL_UNDEF
      end if

    ! ICE RADII
      if (RAD_QI > 1.e-8) then
        !-BRAMS formulation
        RAD_RI = LDRADIUS4(PL,TE,RAD_QI,NL,NI,2)
        ! apply limits
        RAD_RI = MAX( MIN_RI, MIN(RAD_RI*FAC_RI, MAX_RI) )
      else
        RAD_RI = MAPL_UNDEF
      end if

   end subroutine RADCOUPLE

   subroutine  FIX_UP_CLOUDS( &
         QV, &
         TE, &
         QLC,&
         QIC,&
         CF, &
         QLA,&
         QIA,&
         AF ,&
         REMOVE_CLOUDS )

      real, intent(inout) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA
      logical, optional, intent(IN) :: REMOVE_CLOUDS
      real :: FCLD
      logical :: RM_CLDS

                                  RM_CLDS = .false.
      if (present(REMOVE_CLOUDS)) RM_CLDS = REMOVE_CLOUDS

      if (RM_CLDS) then

      ! Remove ALL cloud quants above the klid
         QV = QV + QLA + QIA + QLC + QIC
         TE = TE - (alhlbcp)*(QLA+QLC) - (alhsbcp)*(QIA+QIC)
         AF  = 0.
         QLA = 0.
         QIA = 0.
         CF  = 0.
         QLC = 0.
         QIC = 0.

      else

      ! Ensure total cloud fraction <= 1.0
      FCLD = CF + AF
      if (FCLD > 1.0) then
         CF = CF*(1.0/FCLD)
         AF = AF*(1.0/FCLD)
      end if

      ! Fix if Anvil cloud fraction too small
      if (AF < 1.E-5) then
         QV  = QV + QLA + QIA
         TE  = TE - (alhlbcp)*QLA - (alhsbcp)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if

      ! Fix if LS cloud fraction too small
      if ( CF < 1.E-5 ) then
         QV = QV + QLC + QIC
         TE = TE - (alhlbcp)*QLC - (alhsbcp)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if

      ! LS LIQUID too small
      if ( QLC  < 1.E-8 ) then
         QV = QV + QLC
         TE = TE - (alhlbcp)*QLC
         QLC = 0.
      end if
      ! LS ICE too small
      if ( QIC  < 1.E-8 ) then
         QV = QV + QIC
         TE = TE - (alhsbcp)*QIC
         QIC = 0.
      end if

      ! Anvil LIQUID too small
      if ( QLA  < 1.E-8 ) then
         QV = QV + QLA
         TE = TE - (alhlbcp)*QLA
         QLA = 0.
      end if
      ! Anvil ICE too small
      if ( QIA  < 1.E-8 ) then
         QV = QV + QIA
         TE = TE - (alhsbcp)*QIA
         QIA = 0.
      end if

      ! Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
      if ( ( QLA + QIA ) < 1.E-8 ) then
         QV = QV + QLA + QIA
         TE = TE - (alhlbcp)*QLA - (alhsbcp)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if
      ! Ditto if LS cloud LIQUID+ICE too small
      if ( ( QLC + QIC ) < 1.E-8 ) then
         QV = QV + QLC + QIC
         TE = TE - (alhlbcp)*QLC - (alhsbcp)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if

      end if

   end subroutine FIX_UP_CLOUDS

   subroutine fix_up_clouds_2M( &
         QV, &
         TE, &
         QLC,&
         QIC,&
         CF, &
         QLA,&
         QIA,&
         AF, &
         NL, &
         NI, &
         QR, &
         QS, &
         QG, &
         NR, &
         NS, &
         NG, &
         MASS, & 
         TMP2D)

      real, intent(inout), dimension(:,:,:) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA, QR, QS, QG
      real, intent(inout), dimension(:,:,:) :: NI, NL, NS, NR, NG
      real, dimension(:,:,:),   intent(in)     :: MASS
      real, dimension(:,:),     intent(  out)  :: TMP2D
      integer :: IM, JM, LM

      real, parameter  :: qmin  = 1.0e-12
      real, parameter :: cfmin  = 1.0e-4
      real, parameter :: nmin  = 100.0




      ! Fix if Anvil cloud fraction too small
      where (AF < cfmin)
         QV  = QV + QLA + QIA
         TE  = TE - (MAPL_ALHL/MAPL_CP)*QLA - (MAPL_ALHS/MAPL_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end where

      ! Fix if LS cloud fraction too small
      where ( CF < cfmin)
         QV = QV + QLC + QIC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC - (MAPL_ALHS/MAPL_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end where

      ! LS LIQUID too small
      where ( QLC  < qmin )
         QV = QV + QLC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC
         QLC = 0.
      end where
      ! LS ICE too small
      where ( QIC  < qmin)
         QV = QV + QIC
         TE = TE - (MAPL_ALHS/MAPL_CP)*QIC
         QIC = 0.
      end where

      ! Anvil LIQUID too small
      where ( QLA  < qmin )
         QV = QV + QLA
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLA
         QLA = 0.
      end where
      ! Anvil ICE too small
      where ( QIA  < qmin)
         QV = QV + QIA
         TE = TE - (MAPL_ALHS/MAPL_CP)*QIA
         QIA = 0.
      end where

      ! Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
      where ( ( QLA + QIA ) < qmin)
         QV = QV + QLA + QIA
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLA - (MAPL_ALHS/MAPL_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end where
      ! Ditto if LS cloud LIQUID+ICE too small
      where ( ( QLC + QIC ) < qmin )
         QV = QV + QLC + QIC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC - (MAPL_ALHS/MAPL_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end where

        IM = SIZE( QV, 1 )
    	JM = SIZE( QV, 2 )
    	LM = SIZE( QV, 3 )


      !make sure QI , NI stay within T limits 
         call meltfrz_inst2M  ( IM, JM, LM,    &
              TE              , &
              QLC          , &
              QLA         , &
              QIC           , &
              QIA          , &               
              NL         , &
              NI          )
              
      !make sure no negative number concentrations are passed
      !and that N goes to minimum defaults in the microphysics when mass is too small

      NL =  max(NL, 0.)
      NI =  max(NI, 0.)
      NR =  max(NR, 0.)
      NS =  max(NS, 0.)
      NG =  max(NG, 0.)

      where ((QLA+QLC) .le. qmin) NL = 0.0

      where ((QIA+QIC) .le. qmin) NI = 0.0

      where (QR .le. qmin) NR = 0.

      where (QS .le. qmin) NS = 0.

      where (QG .le. qmin) NG = 0.

      ! need to clean up small negative values. MG does can't handle them
          call FILLQ2ZERO( QV, MASS, TMP2D) 
          call FILLQ2ZERO( QG, MASS, TMP2D) 
          call FILLQ2ZERO( QR, MASS, TMP2D) 
          call FILLQ2ZERO( QS, MASS, TMP2D) 
          call FILLQ2ZERO( QLC, MASS, TMP2D)
          call FILLQ2ZERO( QLA, MASS, TMP2D)  
          call FILLQ2ZERO( QIC, MASS, TMP2D)
          call FILLQ2ZERO( QIA, MASS, TMP2D)
          call FILLQ2ZERO( CF, MASS, TMP2D)
          call FILLQ2ZERO( AF, MASS, TMP2D)

   end subroutine fix_up_clouds_2M

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
      subroutine pdffrac (flag,qtmean,sigmaqt1,sigmaqt2,qstar,clfrac)
      implicit none

      integer flag            ! flag to indicate shape of pdf
                              ! 1 for tophat, 2 for triangular, 3 for Gaussian
      real qtmean             ! Grid box value of q total
      real sigmaqt1           ! width of distribution (sigma)
      real sigmaqt2           ! width of distribution (sigma)
      real qstar              ! saturation q at grid box avg T
      real clfrac             ! cloud fraction (area under pdf from qs)

      real :: qtmode, qtmin, qtmax

      if(flag.eq.1) then
       if((qtmean+sigmaqt1).lt.qstar) then
        clfrac = 0.
       else
        if(sigmaqt1.gt.0.) then
        clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
        else
        clfrac = 1.
        endif
       endif
      elseif(flag.eq.2) then
       qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.
       qtmin = max(qtmode-sigmaqt1,0.)
       qtmax = qtmode + sigmaqt2
       if(qtmax.le.qstar) then
        clfrac = 0.
       elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
        clfrac = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
       elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
        clfrac = 1. - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
       elseif ( qstar.le.qtmin ) then
        clfrac = 1.
       endif
      endif

      return
      end subroutine pdffrac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif

      subroutine pdfcondensate (flag,qtmean4,sigmaqt14,sigmaqt24,qstar4,condensate4)
      implicit none

      integer flag            ! flag to indicate shape of pdf
                              ! 1 for tophat, 2 for triangular
      real qtmean4            ! Grid box value of q total
      real sigmaqt14          ! width of distribution (to left)
      real sigmaqt24          ! width of distribution (to right)
      real qstar4             ! saturation q at grid box avg T
      real condensate4        ! condensate (area under (q*-qt)*pdf from qs)

      real *8 :: qtmode, qtmin, qtmax, constA, constB, cloudf
      real *8 :: term1, term2, term3
      real *8 :: qtmean, sigmaqt1, sigmaqt2, qstar, condensate

      qtmean = dble(qtmean4)
      sigmaqt1 = dble(sigmaqt14)
      sigmaqt2 = dble(sigmaqt24)
      qstar = dble(qstar4)

      if(flag.eq.1) then
       if(qtmean+sigmaqt1.lt.qstar) then
        condensate = 0.d0
       elseif(qstar.gt.qtmean-sigmaqt1)then
        if(sigmaqt1.gt.0.d0) then
         condensate = (min(qtmean + sigmaqt1 - qstar,2.d0*sigmaqt1)**2)/ (4.d0*sigmaqt1)
        else
         condensate = qtmean-qstar
        endif
       else
        condensate = qtmean-qstar
       endif
      elseif(flag.eq.2) then
       qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.d0
       qtmin = max(qtmode-sigmaqt1,0.d0)
       qtmax = qtmode + sigmaqt2
       if ( qtmax.le.qstar ) then
        condensate = 0.d0
       elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
        constB = 2.d0 / ( (qtmax - qtmin)*(qtmax-qtmode) )
        cloudf = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
        term1 = (qstar*qstar*qstar)/3.d0
        term2 = (qtmax*qstar*qstar)/2.d0
        term3 = (qtmax*qtmax*qtmax)/6.d0
        condensate = constB * (term1-term2+term3) - qstar*cloudf
       elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
        constA = 2.d0 / ( (qtmax - qtmin)*(qtmode-qtmin) )
        cloudf = 1.d0 - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
        term1 = (qstar*qstar*qstar)/3.d0
        term2 = (qtmin*qstar*qstar)/2.d0
        term3 = (qtmin*qtmin*qtmin)/6.d0
        condensate = qtmean - ( constA * (term1-term2+term3) ) - qstar*cloudf
       elseif ( qstar.le.qtmin ) then
        condensate = qtmean-qstar
       endif
      endif
      condensate4 = real(condensate)

      return
      end subroutine pdfcondensate

!------------------------------------------------------------------------!
!                   Subroutine partition_dblgss                          !
!------------------------------------------------------------------------!
! Compute SGS cloud fraction and SGS condensation                        !
! using assumed analytic DOUBLE-GAUSSIAN PDF for SGS vertical velocity,  !
! moisture, and  liquid/ice water static energy, based on the            !
! general approach of  Larson et al 2002, JAS, 59, 3519-3539,            !
! and Golaz et al 2002, JAS, 59, 3540-3551                               !
! References in the comments in this code are given to                   !
! the Appendix A of Pete Bogenschutz's dissertation.                     !
!------------------------------------------------------------------------!
 subroutine partition_dblgss( fQi,          &  ! IN
                              tabs,         &
                              qwv,          &
                              qc,           &  ! OUT
!                              qi,           &
                              omega,        &  ! IN
                              zl,           &
                              pval,         &
                              total_water,  &
                              thl_first,    &
                              wthlsec,      &
                              wqwsec,       &
!                              wqtfac,       & ! inter-gaussian qt flux
!                              whlfac,       & ! inter-gaussian hl flux
                              thlsec,       &
                              qwsec,        &
                              qwthlsec,     &
                              w3var,        &
                              w_sec,        &
                              qt3,          &
                              hl3,          &
!                              mffrc,        &
                              PDF_A,        &  ! INOUT
#ifdef PDFDIAG
                              PDF_SIGW1,    &  ! OUT - diagnostic only
                              PDF_SIGW2,    &
                              PDF_W1,       &
                              PDF_W2,       &
                              PDF_SIGTH1,   &
                              PDF_SIGTH2,   &
                              PDF_TH1,      &
                              PDF_TH2,      &
                              PDF_SIGQT1,   &
                              PDF_SIGQT2,   &
                              PDF_QT1,      &
                              PDF_QT2,      &
                              PDF_RQTTH,    &
                              PDF_RWTH,     &
                              PDF_RWQT,     &
#endif
                              wthv_sec,     &  ! OUT - needed elsewhere
                              wqls,         &
                              cld_sgs)

 use MAPL_ConstantsMod, only: ggr    => MAPL_GRAV,   &
                              cp     => MAPL_CP,     &
                              rgas   => MAPL_RGAS,   &
                              rv     => MAPL_RVAP,   &
                              lcond  => MAPL_ALHL,   &
                              lfus   => MAPL_ALHS,   &
                              pi     => MAPL_PI,     &
                              MAPL_H2OMW, MAPL_AIRMW
 use MAPL_SatVaporMod,  only: MAPL_EQsat

   real, intent(in   )  :: fQi         ! ice fraction
!   real, intent(in   )  :: DT          ! timestep [s]
   real, intent(in   )  :: tabs        ! absolute temperature [K]
   real, intent(in   )  :: qwv         ! specific humidity [kg kg-1]
   real, intent(  out)  :: qc          ! liquid+ice condensate [kg kg-1]
   real, intent(in   )  :: omega       ! resolved pressure velocity
   real, intent(in   )  :: zl          ! layer heights [m]
   real, intent(in   )  :: pval        ! layer pressure [Pa]
   real, intent(in   )  :: total_water ! total water [kg kg-1]
   real, intent(in   )  :: thl_first   ! liquid water potential temperature [K]
   real, intent(in   )  :: wthlsec     ! thl flux [K m s-1]
   real, intent(in   )  :: wqwsec      ! total water flux [kg kg-1 m s-1]
!   real, intent(in   )  :: wqtfac      !
!   real, intent(in   )  :: whlfac      !
   real, intent(in   )  :: thlsec
   real, intent(in   )  :: qwsec
   real, intent(in   )  :: qwthlsec
   real, intent(in   )  :: w3var       ! 3rd moment vertical velocity [m3 s-3]
   real, intent(in   )  :: qt3         ! 3rd moment qt from mass flux
   real, intent(in   )  :: hl3         ! 3rd moment hl from mass flux
   real, intent(in   )  :: w_sec       ! 2nd moment vertical velocity [m2 s-2]
!   real, intent(in   )  :: mffrc       ! total EDMF updraft fraction
!   real, intent(inout)  :: qi         ! ice condensate [kg kg-1]
   real, intent(  out)  :: cld_sgs     ! cloud fraction
   real, intent(in   )  ::    PDF_A           ! fractional area of 1st gaussian
#ifdef PDFDIAG
   real, intent(  out)  ::    PDF_SIGW1,    & ! std dev w of 1st gaussian [m s-1]
                              PDF_SIGW2,    & ! std dev w of 2nd gaussian
                              PDF_W1,       & ! mean vertical velocity of 1st gaussian [m s-1]
                              PDF_W2,       & ! mean vertical velocity of 2nd gaussian [m s-1]
                              PDF_SIGTH1,   & ! std dev pot temp of 1st gaussian [K]
                              PDF_SIGTH2,   & ! std dev pot temp of 2nd gaussian [K]
                              PDF_TH1,      & ! mean pot temp of 1st gaussian [K]
                              PDF_TH2,      & ! mean pot temp of 2nd gaussian [K]
                              PDF_SIGQT1,   & ! std dev total water of 1st gaussian [kg kg-1]
                              PDF_SIGQT2,   & ! std dev total water of 2nd gaussian [kg kg-1]
                              PDF_QT1,      & ! mean total water of 1st gaussian [kg kg-1]
                              PDF_QT2,      & ! mean total water of 2nd gaussian [kg kg-1]
                              PDF_RQTTH,    & ! QT-TH correlation coeff
                              PDF_RWTH,     & ! W-TH correlation
                              PDF_RWQT        ! W-QT correlation
#endif
   real, intent(  out)  :: wthv_sec
   real, intent(  out)  :: wqls


! Local variables

   integer i,j,k,ku,kd
   real wrk, wrk1, wrk2, wrk3, wrk4, bastoeps
   real gamaz, thv, rwqt, rwthl, wql1, wql2
   real pkap, diag_qn, diag_frac, diag_ql, diag_qi,w_first,                     &
        sqrtw2, sqrtthl, sqrtqt, w1_1, w1_2, w2_1, w2_2, thl1_1, thl1_2,        &
        thl2_1, thl2_2, qw1_1, qw1_2, qw2_1, qw2_2, aterm, onema, sm,           &
        km1, skew_w, skew_qw, skew_thl, cond_w, sqrtw2t,                                 &
        sqrtthl2_1, sqrtthl2_2, sqrtqw2_1, sqrtqw2_2, corrtest1, corrtest2,     &
        tsign, testvar, r_qwthl_1, Tl1_1, Tl1_2, esval1_1, esval1_2, esval2_1,  &
        esval2_2, om1, om2, lstarn1, lstarn2, qs1, qs2, beta1, beta2, cqt1,     &
        cqt2, s1, s2, cthl1, cthl2, std_s1, std_s2, qn1, qn2, C1, C2, ql1, ql2, &
        qi1, qi2, wqis, wqtntrgs, whlntrgs


! Set constants and parameters
   real, parameter :: sqrt2 = sqrt(2.0)
   real, parameter :: sqrtpii = 1.0/sqrt(pi+pi)
   real, parameter :: tbgmin = 233.16
   real, parameter :: tbgmax = MAPL_TICE
   real, parameter :: a_bg   = 1.0/(tbgmax-tbgmin)
   real, parameter :: thl_tol = 1.e-2
   real, parameter :: w_thresh = 0.001
   real, parameter :: rt_tol = 1.e-4
   real, parameter :: w_tol_sqd = 4.0e-04   ! Min vlaue of second moment of w
   real, parameter :: onebrvcp = 1.0/(rv*cp)
   real, parameter :: skew_facw = 1.2
   real, parameter :: skew_fact = 0.5
   real, parameter :: lsub = lcond+lfus
   real, parameter :: fac_cond = lcond/cp
   real, parameter :: fac_sub = lsub/cp
   real, parameter :: fac_fus = lfus/cp
   real, parameter :: gocp = ggr/cp
   real, parameter :: rog = rgas / ggr
   real, parameter :: kapa = rgas / cp
   real, parameter :: epsv=MAPL_H2OMW/MAPL_AIRMW

   real, parameter :: use_aterm_memory = 1.
   real, parameter :: tauskew = 2400.

! define conserved variables
   gamaz = gocp * zl
   thv   = tabs * (1.0+epsv*qwv)
   thv   = thv*(100000.0/pval) ** kapa

   w_first = - rog * omega * thv / pval

! Initialize cloud variables to zero
   diag_qn   = 0.0
   diag_frac = 0.0
   diag_ql   = 0.0
   diag_qi   = 0.0

   pkap = (pval/100000.0) ** kapa


! Compute square roots of some variables so we don't have to do it again
          if (w_sec > w_thresh*w_thresh) then
            sqrtw2   = sqrt(w_sec)
            Skew_w   = w3var / (sqrtw2*sqrtw2*sqrtw2)
          else
            sqrtw2   = w_thresh
            Skew_w   = 0.
          endif
          if (thlsec > 1e-6) then
            sqrtthl  = sqrt(thlsec)
            skew_thl = hl3 / sqrtthl**3
          else
            sqrtthl  = 1e-3
            skew_thl = 0.
          endif
          if (qwsec > 1e-8*total_water*total_water) then
            sqrtqt   = sqrt(qwsec)
            skew_qw =  qt3/sqrtqt**3
          else
            sqrtqt   = 1e-4*total_water
            skew_qw  = 0.
          endif

! Find parameters of the double Gaussian PDF of vertical velocity

!          aterm = pdf_a

!         if (use_aterm_memory/=0) then   ! use memory in aterm and qt skewness
!          aterm = pdf_a

!           pdf_a = min(0.5,max(0.,(pdf_a+mffrc)/(1.+DT/tauskew)))
!          if (mffrc>=1e-3) then                ! if active updraft this timestep
!            if (pdf_a>1e-3) then                ! if distribution is skewed (recent updrafts)
!              pdf_a = (pdf_a+mffrc)/(1.+DT/tauskew) !max(mffrc,aterm*max(1.-DT/tauskew,0.0))
!            else                               ! if distribution unskewed
!              pdf_a = mffrc
!            end if
!          else                                 ! if no active updraft
!            if (pdf_a.gt.1e-3) then            ! but there is residual skewness
!              pdf_a = pdf_a*max(1.-DT/tauskew,0.0)
!            else
!              pdf_a = 0.0
!            end if
!          end if
!         else  ! don't use memory in aterm and qt skewness
!           pdf_a = max(0.,min(0.5,mffrc))
!         end if

         if (pdf_a>1e-3 .and. pdf_a<0.5) then
           aterm = pdf_a
         else
           aterm = 0.5
         end if
         onema = 1.0 - aterm


! If variance of w is too small or no skewness then
!          IF (w_sec <= w_tol_sqd .or. mffrc.lt.0.01) THEN ! If variance of w is too small then
          IF (w_sec <= w_tol_sqd) THEN ! If variance of w is too small then
            Skew_w = 0.
            w1_1   = 0.
            w1_2   = 0.
            w2_1   = w_sec
            w2_2   = w_sec
!            aterm  = 0.5
!            onema  = 0.5
          ELSE

! Proportionality coefficients between widths of each vertical velocity
! gaussian and the sqrt of the second moment of w
 !           w2_1 = 0.4
 !           w2_2 = 0.4

! analytic double gaussian 2, variable sigma_w

            wrk2 = 0.667*abs(Skew_w)**0.333    ! m below A.24
! not used     wrk = (1+wrk2*wrk2)**3/((3.+wrk2*wrk2)*wrk2)**2  ! M in A.24

            w2_1 = (onema/(aterm*(1.+wrk2**2)))**0.5
            w2_2 = (aterm/(onema*(1.+wrk2**2)))**0.5

            w1_1 = wrk2*w2_1             ! w1_tilde in A.23
            w1_2 = -wrk2*w2_2

! Compute realtive weight of the first PDF "plume"
! See Eq A4 in Pete's dissertaion -  Ensure 0.01 < a < 0.99

!            wrk = 1.0 - w2_1    ! 1-sigw2tilde = 1-0.4
!            aterm = max(0.01,min(0.5*(1.-Skew_w*sqrt(1./(4.*wrk*wrk*wrk+Skew_w*Skew_w))),0.99))

!            sqrtw2t = sqrt(wrk)

! Eq. A.5-A.6
!            wrk  =   sqrt(onema/aterm)
!            w1_1 =   sqrtw2t * wrk  ! w1tilde (A.5)
!            w1_2 = - sqrtw2t / wrk  ! w2tilde (A.6)

!            w2_1 = w2_1 * w_sec  ! sigma_w1 **2
!            w2_2 = w2_2 * w_sec  ! sigma_w2 **2

          ENDIF


!  Find parameters of the PDF of liquid/ice static energy

          ! inter-gaussian flux limited to 2x total flux
!          whlntrgs = max(min(whlfac,2.*abs(wthlsec)),-2.*abs(wthlsec))

          IF (thlsec <= thl_tol*thl_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
            thl1_1     = thl_first
            thl1_2     = thl_first
            thl2_1     = thlsec
            thl2_2     = thlsec
            sqrtthl2_1 = sqrt(thlsec)
            sqrtthl2_2 = sqrtthl2_1

          ELSE

!            corrtest1 = max(-1.0,min(1.0,whlntrgs/(sqrtw2*sqrtthl)))
            corrtest1 = max(-1.0,min(1.0,wthlsec/(sqrtw2*sqrtthl)))

            thl1_1 = -corrtest1 / w1_2       ! A.7
            thl1_2 = -corrtest1 / w1_1       ! A.8

!            thl1_1 = -whlntrgs / (w1_2*sqrtthl)   !   normalized
!            thl1_2 = -whlntrgs / (w1_1*sqrtthl)

            wrk1   = thl1_1 * thl1_1
            wrk2   = thl1_2 * thl1_2
            wrk3   = 1.0 - aterm*wrk1 - onema*wrk2
            wrk4   = skew_thl - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2
            wrk    = 3. * (thl1_2-thl1_1)
            if (wrk /= 0.0) then
              thl2_1 = thlsec * min(100.,max(0.,( 3.*thl1_2*wrk3-wrk4)/(aterm*wrk))) ! A.9
              thl2_2 = thlsec * min(100.,max(0.,(-3.*thl1_1*wrk3+wrk4)/(onema*wrk))) ! A.10
            else
!              thl2_1 = 0.0
!              thl2_2 = 0.0
              thl2_1 = thlsec
              thl2_2 = thlsec
            endif

            thl1_1 = thl1_1*sqrtthl + thl_first    ! convert to physical units
            thl1_2 = thl1_2*sqrtthl + thl_first

            sqrtthl2_1 = sqrt(thl2_1)
            sqrtthl2_2 = sqrt(thl2_2)

          ENDIF

          ! implied correlation coefficient
#ifdef PDFDIAG
          PDF_RWTH = max(-1.,min(1.,( wthlsec/sqrtw2-aterm*(thl1_1-thl_first)*(w1_1-w_first) &
                     -onema*(thl1_2-thl_first)*(w1_2-w_first) )               &
                     / (aterm*sqrt(thl2_1*w2_1)+onema*sqrt(thl2_2*w2_2)) ))
#endif

!  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

          ! inter-gaussian flux, limited to 2x total flux
!          wqtntrgs = max(min(wqtfac,2.*abs(wqwsec)),-2.*abs(wqwsec))

          IF (qwsec <= rt_tol*rt_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN ! if no active updrafts

            if (aterm .lt. 1e-3 .or. aterm.gt.0.499 .or. abs(Skew_qw).lt.1e-8) then ! if no residual skewness
              qw1_1     = total_water
              qw1_2     = total_water
              qw2_1     = qwsec
              qw2_2     = qwsec
              sqrtqw2_1 = sqrt(qw2_1)
              sqrtqw2_2 = sqrt(qw2_2)

            else
!              qw1_1     = total_water
!              qw1_2     = total_water
!              qw2_1     = qwsec
!              qw2_2     = qwsec
!              wrk1 = max(min(10.,skew_qw*sqrtqt**3)   ! third moment qt
              wrk1 = qt3
              qw1_1 = total_water + (wrk1/(2.*aterm-aterm**3/onema**2))**(1./3.)
              qw1_2 = (total_water -aterm*qw1_1)/onema

              qw2_1 = qwsec - min(0.5*qwsec,max(0.,(aterm/onema)*(qw1_1-total_water)**2))
              qw2_2 = qw2_1
              sqrtqw2_1 = sqrt(qw2_1)
              sqrtqw2_2 = sqrt(qw2_2)
            end if

          ELSE  ! active updrafts

!            corrtest2 = max(-1.0,min(1.0,wqtntrgs/(sqrtw2*sqrtqt)))
            corrtest2 = max(-1.0,min(1.0,0.5*wqwsec/(sqrtw2*sqrtqt)))

            qw1_1 = - corrtest2 / w1_2            ! A.7
            qw1_2 = - corrtest2 / w1_1            ! A.8

            tsign = abs(qw1_2-qw1_1)

            wrk1  = qw1_1 * qw1_1
            wrk2  = qw1_2 * qw1_2
            wrk3  = 1.      - aterm*wrk1       - onema*wrk2
            wrk4  = Skew_qw - aterm*wrk1*qw1_1 - onema*wrk2*qw1_2
            wrk   = 3. * (qw1_2-qw1_1)

            if (wrk /= 0.0) then
              qw2_1 = qwsec * min(100.,max(0.,( 3.*qw1_2*wrk3-wrk4)/(aterm*wrk))) ! A.10
              qw2_2 = qwsec * min(100.,max(0.,(-3.*qw1_1*wrk3+wrk4)/(onema*wrk))) ! A.11
            else
!              qw2_1 = 0.0
!              qw2_2 = 0.0
              qw2_1 = qwsec
              qw2_2 = qwsec
            endif

            qw1_1 = qw1_1*sqrtqt + total_water
            qw1_2 = qw1_2*sqrtqt + total_water

            sqrtqw2_1 = sqrt(qw2_1)
            sqrtqw2_2 = sqrt(qw2_2)

          ENDIF   ! if qwsec small

          ! implied correlation coefficient
#ifdef PDFDIAG
          PDF_RWQT = max(-1.,min(1.,( wqwsec/sqrtw2-aterm*(qw1_1-total_water)*(w1_1-w_first) &
                     -onema*(qw1_2-total_water)*(w1_2-w_first) )              &
                     / (aterm*sqrt(qw2_1*w2_1)+onema*sqrt(qw2_2*w2_2)) ))
#endif

!  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

          w1_1 = w1_1*sqrtw2 + w_first    ! using A.5 and A.6
          w1_2 = w1_2*sqrtw2 + w_first    ! note: this is already done for w2_x


!=== Assign PDF diagnostics ===!

!          pdf_a = aterm

#ifdef PDFDIAG
          pdf_th1 = thl1_1
          pdf_th2 = thl1_2
          pdf_sigth1 = sqrtthl2_1
          pdf_sigth2 = sqrtthl2_2

          pdf_qt1 = qw1_1
          pdf_qt2 = qw1_2
          pdf_sigqt1 = sqrtqw2_1
          pdf_sigqt2 = sqrtqw2_2

          pdf_w1 = w1_1
          pdf_w2 = w1_2
          if (w2_1.ne.0.) then
            pdf_sigw1 = w2_1*sqrtw2
            pdf_sigw2 = w2_2*sqrtw2
          else
            pdf_sigw1 = 0.0
            pdf_sigw2 = 0.0
          end if
#endif

!==============================!


!  FIND WITHIN-PLUME CORRELATIONS

          testvar = aterm*sqrtqw2_1*sqrtthl2_1 + onema*sqrtqw2_2*sqrtthl2_2

          IF (testvar == 0) THEN
            r_qwthl_1 = 0.
          ELSE
            r_qwthl_1 = max(-1.0,min(1.0,(qwthlsec-aterm*(qw1_1-total_water)*(thl1_1-thl_first)-onema*(qw1_2-total_water)*(thl1_2-thl_first))/testvar)) ! A.12
          ENDIF

#ifdef PDFDIAG
          pdf_rqtth = r_qwthl_1
#endif


!  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS
! This section follows Bogenschutz thesis Appendix A, based on
! Sommeria and Deardorff (1977) and Lewellen and Yoh (1993).

          Tl1_1 = thl1_1 - gamaz
          Tl1_2 = thl1_2 - gamaz

! Now compute qs

! I assume that temperature difference is small between plumes, and
! om1=om2=ice_fraction() from hystpdf
          esval1_1 = 0.
          esval1_2 = 0.
          esval2_1 = 0.
          esval2_2 = 0.
          om1      = 1.
          om2      = 1.

! Partition based on temperature for the first plume

!          IF (Tl1_1 >= tbgmax) THEN
!            esval1_1 = MAPL_EQsat(Tl1_1)
!            lstarn1  = lcond
!          ELSE IF (Tl1_1 < tbgmin) THEN
!            esval1_1 = MAPL_EQsat(Tl1_1,OverIce=.TRUE.)
!            lstarn1  = lsub
!          ELSE
            esval1_1 = MAPL_EQsat(Tl1_1)
            esval2_1 = MAPL_EQsat(Tl1_1,OverIce=.TRUE.)
            om1      = max(0.,min(1.,1.-fQi)) !max(0.,min(1.,a_bg*(Tl1_1-tbgmin)))  ! may be inconsistent with hystpdf ice fraction
            lstarn1  = lcond + (1.-om1)*lfus
!          ENDIF

          ! this is qs evaluated at Tl
          qs1   =     om1  * (0.622*esval1_1/max(esval1_1,pval-0.378*esval1_1))      &
                + (1.-om1) * (0.622*esval2_1/max(esval2_1,pval-0.378*esval2_1))

          beta1 = (lstarn1*lstarn1*onebrvcp) / (Tl1_1*Tl1_1)

! Are the two plumes equal?  If so then set qs and beta
! in each column to each other to save computation
          IF (Tl1_1 == Tl1_2) THEN
            qs2   = qs1
            beta2 = beta1
          ELSE

              esval1_2 = MAPL_EQsat(Tl1_2)
              esval2_2 = MAPL_EQsat(Tl1_2,OverIce=.TRUE.)
              om2      = max(0.,min(1.,1.-fQi)) !max(0.,min(1.,a_bg*(Tl1_2-tbgmin)))
              lstarn2  = lcond + (1.-om2)*lfus

            qs2   =     om2  * (0.622*esval1_2/max(esval1_2,pval-0.378*esval1_2))    &
                  + (1.-om2) * (0.622*esval2_2/max(esval2_2,pval-0.378*esval2_2))
!            qs2 = GEOS_QSAT( Tl1_2, pval*0.01 )

            beta2 = (lstarn2*lstarn2*onebrvcp) / (Tl1_2*Tl1_2)              ! A.18

          ENDIF


!  Now compute cloud stuff -  compute s term

          cqt1  = 1.0 / (1.0+beta1*qs1)                                     ! A.19
          wrk   = (1.0+beta1*qw1_1) * cqt1

          s1    = qw1_1 - qs1* wrk                                          ! A.17
          cthl1 = cqt1*wrk*(cp/lcond)*beta1*qs1*pkap                        ! A.20

          wrk1   = cthl1 * cthl1
          wrk2   = cqt1  * cqt1
          std_s1 = sqrt(max(0.,wrk1*thl2_1+wrk2*qw2_1-2.*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))

          qn1 = 0.
          C1  = 0.

          IF (std_s1 /= 0) THEN
            wrk = s1 / (std_s1*sqrt2)
            C1 = 0.5*(1.+erf(wrk))                                         ! A.15
            IF (C1 /= 0) qn1 = s1*C1 + (std_s1*sqrtpii)*exp(-wrk*wrk)      ! A.16
          ELSEIF (s1 > 0) THEN
            C1  = 1.0
            qn1 = s1
          ENDIF

! now compute non-precipitating cloud condensate

! If two plumes exactly equal, then just set many of these
! variables to themselves to save on computation.
          IF (qw1_1 == qw1_2 .and. thl2_1 == thl2_2 .and. qs1 == qs2) THEN
            s2     = s1
            cthl2  = cthl1
            cqt2   = cqt1
            std_s2 = std_s1
            C2     = C1
            qn2    = qn1
          ELSE

            cqt2   = 1.0 / (1.0+beta2*qs2)
            wrk    = (1.0+beta2*qw1_2) * cqt2
            s2     = qw1_2 - qs2*wrk
            cthl2  = wrk*cqt2*(cp/lcond)*beta2*qs2*pkap
            wrk1   = cthl2 * cthl2
            wrk2   = cqt2  * cqt2
            std_s2 = sqrt(max(0.,wrk1*thl2_2+wrk2*qw2_2-2.*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

            qn2 = 0.
            C2  = 0.

            IF (std_s2 /= 0) THEN
              wrk = s2 / (std_s2*sqrt2)
              C2  = 0.5*(1.+erf(wrk))
              IF (C2 /= 0) qn2 = s2*C2 + (std_s2*sqrtpii)*exp(-wrk*wrk)
            ELSEIF (s2 > 0) THEN
              C2  = 1.0
              qn2 = s2
            ENDIF

          ENDIF


! finally, compute the SGS cloud fraction
          cld_sgs = aterm*C1 + onema*C2

!          om1 = max(0.,min(1.,(Tl1_1-tbgmin)*a_bg))
!          om2 = max(0.,min(1.,(Tl1_2-tbgmin)*a_bg))

          qn1 = min(qn1,qw1_1)
          qn2 = min(qn2,qw1_2)

          ql1 = qn1*om1
          ql2 = qn2*om2

          qi1 = qn1 - ql1
          qi2 = qn2 - ql2

          qc = min(max(0.0, aterm*qn1 + onema*qn2), total_water)


! Update temperature variable based on diagnosed cloud properties
!          om1         = 1.-fQi !max(0.,min(1.,(tabs-tbgmin)*a_bg))
!          lstarn1     = lcond + (1.-om1)*lfus
!          tabs = thl_first - gamaz + fac_cond*(diag_ql) &
!                            + fac_sub *(diag_qi) !&
                    !  + tkesbdiss(i,j,k) * (dtn/cp)      ! tke dissipative heating
! Update moisture fields

         if (sqrtqt>1e-4*total_water .AND. sqrtw2>w_thresh) then
            rwqt = 0.5*wqwsec/(sqrtqt*sqrtw2)
!            rwqt = max(-1.,min(1.,pdf_rwqt))
         else
            rwqt = 0.0
         end if
         if (sqrtthl>1e-3 .AND. sqrtw2>w_thresh) then
            rwthl = wthlsec/(sqrtthl*sqrtw2)
!            rwthl = max(-1.,min(1.,pdf_rwth))
         else
            rwthl = 0.0
         end if

         wql1 = C1*(cqt1*sqrt(w2_1)*sqrt(qw2_1)*rwqt-cthl1*sqrt(w2_1)*sqrt(thl2_1)*rwthl)
         wql2 = C2*(cqt2*sqrt(w2_2)*sqrt(qw2_2)*rwqt-cthl2*sqrt(w2_2)*sqrt(thl2_2)*rwthl)

! Compute the liquid water flux
          wqls = aterm * ((w1_1-w_first)*ql1+wql1) + onema * ((w1_2-w_first)*ql2+wql2)
          wqis = aterm * ((w1_1-w_first)*qi1) + onema * ((w1_2-w_first)*qi2)

! diagnostic buoyancy flux.  Includes effects from liquid water, ice
! condensate, liquid & ice precipitation
          wrk = epsv * thv

          bastoeps = (rv/rgas) * thv   ! thetav / epsilon

          wthv_sec = wthlsec + wrk*wqwsec                                     &
                   + (fac_cond-bastoeps)*wqls                                 &
                   + (fac_sub-bastoeps) *wqis
!                   + ((lstarn1/cp)-thv(i,j,k))*0.5*(wqp_sec(i,j,kd)+wqp_sec(i,j,ku))

  end subroutine partition_dblgss

   subroutine hystpdf( &
         DT          , &
         ALPHA       , &
         PDFSHAPE    , &
         CNVFRC      , &
         SRF_TYPE    , &
         PL          , &
         ZL          , &
         QV          , &
         QLLS        , &
         QLCN        , &
         QILS        , &
         QICN        , &
         TE          , &
         CLLS        , &
         CLCN        , &
         NL          , &
         NI          , &
         WHL         , &
         WQT         , &
         HL2         , &
         QT2         , &
         HLQT        , &
         W3          , &
         W2          , &
         MFQT3       , &
         MFHL3       , &
         PDF_A       , &  ! can remove these after development
         PDFITERS    , &
#ifdef PDFDIAG
         PDF_SIGW1,  &
         PDF_SIGW2,  &
         PDF_W1,     &
         PDF_W2,     &
         PDF_SIGHL1, &
         PDF_SIGHL2, &
         PDF_HL1,    &
         PDF_HL2,    &
         PDF_SIGQT1, &
         PDF_SIGQT2, &
         PDF_QT1,    &
         PDF_QT2,    &
         PDF_RHLQT,  &
         PDF_RWHL,   &
         PDF_RWQT,   &
#endif
         WTHV2,      &
         WQL,        &
         needs_preexisting, &
         USE_BERGERON, &
         SC_ICE )

      real, intent(in)    :: DT,ALPHA,PL,ZL
      integer, intent(in) :: PDFSHAPE
      real, intent(inout) :: TE,QV,QLLS,QILS,CLLS,QLCN,QICN,CLCN,PDF_A
      real, intent(in)    :: NL,NI,CNVFRC,SRF_TYPE
      real, intent(in)    :: WHL,WQT,HL2,QT2,HLQT,W3,W2,MFQT3,MFHL3
#ifdef PDFDIAG
      real, intent(out)   :: PDF_SIGW1, PDF_SIGW2, PDF_W1, PDF_W2, &
                             PDF_SIGHL1, PDF_SIGHL2, PDF_HL1, PDF_HL2, &
                             PDF_SIGQT1, PDF_SIGQT2, PDF_QT1, PDF_QT2, &
                             PDF_RHLQT,  PDF_RWHL, PDF_RWQT
#endif
      real, intent(out)   :: WTHV2, WQL
      real, intent(out)   :: PDFITERS
      logical, intent(in) :: needs_preexisting, USE_BERGERON
      real, optional , intent(in) :: SC_ICE

      ! internal arrays
      real :: TAU,HL
      real :: QT, sigmaqt1, sigmaqt2, scice

      real :: QSx,DQsx,QS,DQs

      real :: TEp, QSp, CFp, QVp, QCp
      real :: TEn, QSn, CFn, QVn, QCn

      real :: QCx, QC, fQi, QCi, qsnx
      real :: dQICN, dQLCN, dQILS, dQLLS, Nfac, NLv, NIv

      real :: tmpARR
      real :: alhxbcp, DQCALL
      ! internal scalars
      integer :: N, nmax

      character*(10) :: Iam='Process_Library:hystpdf'

      scice =  1.0

                      tmpARR = 0.0
      if (CLCN < 1.0) tmpARR = 1.0/(1.0-CLCN)

      CFn = (CLLS       )*tmpARR
      QCn = (QLLS + QILS)*tmpARR
      QCi = (QILS)*tmpARR
      TEn = TE

      DQS = GEOS_DQSAT( TEn, PL, QSAT=QSx )
      QVn = ( QV - QSx*CLCN )*tmpARR

      QT = QCn + QVn  !Total LS water after microphysics

      nmax = 20
      do n=1,nmax

         QVp = QVn
         QCp = QCn
         CFp = CFn
         TEp = TEn
         DQS = GEOS_DQSAT( TEn, PL, QSAT=QSn )

         if(present(SC_ICE)) then
            scice = min(max(SC_ICE, 1.0), 1.7)
            qsnx= Qsn*scice !
            if ((QCi .ge. 0.0) .and. (Qsn .gt. Qt))  QSn=Qsnx !this way we do not evaporate preexisting ice but maintain supersat
         end if

         if(PDFSHAPE.lt.2) then  ! top-hat
            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn
         elseif(PDFSHAPE.eq.2) then  ! triangular
            ! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*QSn (alpha is half width)
            ! for triangular, skewed r : sigmaqt1 < sigmaqt2
            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn
         elseif(PDFSHAPE .eq. 3) then ! single gaussian
            ! missing
         elseif(PDFSHAPE .eq. 4) then !lognormal (sigma is dimmensionless)
            sigmaqt1 =  max(ALPHA/sqrt(3.0), 0.001)
         endif

         if (PDFSHAPE.lt.5) then
           call pdffrac(PDFSHAPE,QT,sigmaqt1,sigmaqt2,QSn,CFn)
           call pdfcondensate(PDFSHAPE,QT,sigmaqt1,sigmaqt2,QSn,QCn)
         elseif (PDFSHAPE.eq.5) then

           ! Update the liquid water static energy
           fQi = ice_fraction( TEn, CNVFRC,SRF_TYPE )
           alhxbcp = (1.0-fQi)*alhlbcp + fQi*alhsbcp
           HL = TEn + gravbcp*ZL - alhxbcp*QCn

           call partition_dblgss(fQi,          &
                                 TEn,          &
                                 QVn,          &
                                 QCn,          &
                                 0.0,          & ! assume OMEGA=0
                                 ZL,           &
                                 PL*100.,      &
                                 QT,           &
                                 HL,           &
                                 WHL,          &
                                 WQT,          &
                                 HL2,          &
                                 QT2,          &
                                 HLQT,         &
                                 W3,           &
                                 W2,           &
                                 MFQT3,        &
                                 MFHL3,        &
                                 PDF_A,        &
#ifdef PDFDIAG
                                 PDF_SIGW1,    &
                                 PDF_SIGW2,    &
                                 PDF_W1,       &
                                 PDF_W2,       &
                                 PDF_SIGHL1,   &
                                 PDF_SIGHL2,   &
                                 PDF_HL1,      &
                                 PDF_HL2,      &
                                 PDF_SIGQT1,   &
                                 PDF_SIGQT2,   &
                                 PDF_QT1,      &
                                 PDF_QT2,      &
                                 PDF_RHLQT,    &
                                 PDF_RWHL,     &
                                 PDF_RWQT,     &
#endif
                                 WTHV2,        &
                                 WQL,          &
                                 CFn)
         endif

         IF(USE_BERGERON) THEN
           DQCALL = QCn - QCp
           CLLS = CFn * (1.-CLCN)
           Nfac = 100.*PL*R_AIR/TEn !density times conversion factor
           NLv = NL/Nfac
           NIv = NI/Nfac
           call Bergeron_Partition( &         !Microphysically-based partitions the new condensate
                 DT               , &
                 PL               , &
                 TEn              , &
                 QT               , &
                 QILS             , &
                 QICN             , &
                 QLLS             , &
                 QLCN             , &
                 CLLS             , &
                 CLCN             , &
                 NLv              , &
                 NIv              , &
                 DQCALL           , &
                 fQi              , &
                 CNVFRC,SRF_TYPE  , &
                 needs_preexisting)
         ELSE
           fQi = ice_fraction( TEn, CNVFRC,SRF_TYPE )
         ENDIF

         alhxbcp = (1.0-fQi)*alhlbcp + fQi*alhsbcp
         if(PDFSHAPE.eq.1) then
            QCn = QCp +     (QCn-QCp)/(1.-(CFn*(ALPHA-1.)-(QCn/QSn))*DQS*alhxbcp)
         elseif(PDFSHAPE.eq.2) then
            ! This next line needs correcting - need proper d(del qc)/dT derivative for triangular
            ! for now, just use relaxation of 1/2 of top-hat.
            QCn = QCp + 0.5*(QCn-QCp)/(1.-(CFn*(ALPHA-1.)-(QCn/QSn))*DQS*alhxbcp)
         elseif(PDFSHAPE.eq.5) then
            QCn = QCp + 0.5*(QCn-QCp)
         endif

         QVn = QVp - (QCn - QCp)
         TEn = TEp + (1.0-fQi)*(alhlbcp)*(QCn - QCp)*(1.-CLCN)  &
                   +      fQi *(alhsbcp)*(QCn - QCp)*(1.-CLCN)

         PDFITERS = n
         if (abs(TEn - TEp) .lt. 0.00001) exit

      enddo ! qsat iteration

      ! Now take {\em New} condensate and partition into ice and liquid

      ! large-scale
      CLLS = CFn * (1.-CLCN)
      QCn  = QCn * (1.-CLCN)
      QCx  = QCn - (QLLS+QILS)
      if (QCx .lt. 0.0) then  !net evaporation
         dQLLS = max(QCx        , -QLLS) ! Water evaporates first
         dQILS = max(QCx - dQLLS, -QILS) ! Then sublimation
      else
         dQLLS = (1.0-fQi)*QCx
         dQILS =      fQi *QCx
      end if

      ! Clean-up cloud if fractions are too small
      if ( CLLS < 1.e-5 ) then
         dQILS = -QILS
         dQLLS = -QLLS
      end if

      QILS   = QILS + dQILS
      QLLS   = QLLS + dQLLS
      QV     = QV -         (dQILS+dQLLS)
      TE     = TE + alhlbcp*(dQILS+dQLLS) + alhfbcp*(dQILS)

   end subroutine hystpdf

!==========Estimate RHcrit========================
!==============================
 subroutine pdf_alpha(PP,P_LM, ALPHA, FRLAND, MINRHCRIT, TURNRHCRIT, TURNRHCRIT_UPPER, EIS, RHC_OPTION)

      real,    intent(in)  :: PP, P_LM !mbar
      real,    intent(out) :: ALPHA
      real,    intent(in)  :: FRLAND
      real,    intent(in)  :: MINRHCRIT, TURNRHCRIT, EIS, TURNRHCRIT_UPPER
      integer, intent(in)  :: RHC_OPTION !0-Slingo(1985), 1-QUAAS (2012)
      real :: dw_land = 0.20 !< base value for subgrid deviation / variability over land
      real :: dw_ocean = 0.10 !< base value for ocean
      real :: sloperhcrit =20.
      !real :: TURNRHCRIT_UPPER = 300.
      real ::  aux1, aux2, maxalpha

      IF (RHC_OPTION .lt. 1) then

          !  Use Slingo-Ritter (1985) formulation for critical relative humidity
          !Reformulated by Donifan Barahona

          maxalpha=1.0-MINRHCRIT
          aux1 = min(max((pp- TURNRHCRIT)/sloperhcrit, -20.0), 20.0)
          aux2 = min(max((TURNRHCRIT_UPPER - pp)/sloperhcrit, -20.0), 20.0)

          if (FRLAND > 0.05)  then
             aux1=1.0
          else
             aux1 = 1.0/(1.0+exp(aux1)) !this function reproduces the old Sligo function.
          end if

           if (TURNRHCRIT_UPPER .gt. 0.0) then 
          	 aux2= 1.0/(1.0+exp(aux2)) !this function reverses the profile P< TURNRHCRIT_UPPER
           else
           aux2=1.0
           end if 
           
           ALPHA  = min(maxalpha*aux1*aux2, 0.4)

       ELSE
           ! based on Quass 2012 https://doi.org/10.1029/2012JD017495
             if (EIS > 5.0) then ! Stable
                ALPHA = 1.0 - ((1.0-dw_land ) + (0.99 - (1.0-dw_land ))*exp(1.0-(P_LM/PP)**2))
             else ! Unstable
                ALPHA = 1.0 - ((1.0-dw_ocean) + (0.99 - (1.0-dw_ocean))*exp(1.0-(P_LM/PP)**4))
             endif
       END IF

   end subroutine pdf_alpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Parititions DQ into ice and liquid. Follows Barahona et al. GMD. 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine Bergeron_Partition (           &
         DTIME            , &
         PL               , &
         TE               , &
         QV               , &
         QILS             , &
         QICN             , &
         QLLS             , &
         QLCN             , &
         CF               , &
         AF               , &
         NL               , &
         NI               , &
         DQALL            , &
         FQI              , &
         CNVFRC, SRF_TYPE , &
         needs_preexisting )

      real ,  intent(in   )    :: DTIME, PL, TE       !, RHCR
      real ,  intent(inout   )    ::  DQALL
      real ,  intent(in)    :: QV, QLLS, QLCN, QICN, QILS
      real ,  intent(in)    :: CF, AF, NL, NI
      real, intent (out) :: FQI
      real, intent(in) :: CNVFRC, SRF_TYPE
      logical, intent (in)  :: needs_preexisting

      real  :: DC, TEFF,QCm,DEP, &
            QC, QS, RHCR, DQSL, DQSI, QI, TC, &
            DIFF, DENAIR, DENICE, AUX, &
            DCF, QTOT, LHCORR,  QL, DQI, DQL, &
            QVINC, QSLIQ, CFALL,  new_QI, new_QL, &
            QSICE, fQI_0, QS_0, DQS_0, FQA, NIX

      DIFF = 0.0
      DEP  = 0.0
      QI   = QILS + QICN !neccesary because NI is for convective and large scale
      QL   = QLLS + QLCN
      QTOT = QI+QL
      FQA  = 0.0
      if (QTOT .gt. 0.0) FQA = (QICN+QILS)/QTOT
      NIX  = (1.0-FQA)*NI

      DQALL = DQALL/DTIME
      CFALL = min(CF+AF, 1.0)
      TC    = TE-MAPL_TICE
      fQI_0 = fQI

      !Completely glaciated cloud:
      if (TE .ge. iT_ICE_MAX) then   !liquid cloud
         FQI   = 0.0
      elseif(TE .le. iT_ICE_ALL) then !ice cloud
         FQI   = 1.0
      else !mixed phase cloud
         FQI   = 0.0
         if (QILS .le. 0.0) then
              if (needs_preexisting) then
              ! new 0518 this line ensures that only preexisting ice can grow by deposition.
              ! Only works if explicit ice nucleation is available (2 moment muphysics and up)
              else
                  fQi = ice_fraction( TE, CNVFRC, SRF_TYPE )
              end if
              return
         end if

         QVINC =  QV
         QSLIQ = GEOS_QsatLQU( TE, PL*100.0 , DQ=DQSL )
         QSICE = GEOS_QsatICE( TE, PL*100.0 , DQ=DQSI )
         QVINC = MIN(QVINC, QSLIQ) !limit to below water saturation

         ! Calculate deposition onto preexisting ice

         DIFF=(0.211*1013.25/(PL+0.1))*(((TE+0.1)/MAPL_TICE)**1.94)*1e-4  !From Seinfeld and Pandis 2006
         DENAIR=PL*100.0/MAPL_RGAS/TE
         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
         LHcorr = ( 1.0 + DQSI*MAPL_ALHS/MAPL_CP) !must be ice deposition

         if  ((NIX .gt. 1.0) .and. (QILS .gt. 1.0e-10)) then
            DC=max((QILS/(NIX*DENICE*MAPL_PI))**(0.333), 20.0e-6) !Assumme monodisperse size dsitribution
         else
            DC=20.0e-6
         end if

         TEFF= NIX*DENAIR*2.0*MAPL_PI*DIFF*DC/LHcorr ! 1/Dep time scale

         DEP=0.0
         if ((TEFF .gt. 0.0) .and. (QILS .gt. 1.0e-14)) then
            AUX=max(min(DTIME*TEFF, 20.0), 0.0)
            DEP=(QVINC-QSICE)*(1.0-EXP(-AUX))/DTIME
         end if
         DEP=MAX(DEP, -QILS/DTIME) !only existing ice can be sublimated

         DQI = 0.0
         DQL = 0.0
         FQI = 0.0
         !Partition DQALL accounting for Bergeron-Findensen process
         if  (DQALL .ge. 0.0) then !net condensation. Note: do not allow bergeron with QLCN
            if (DEP .gt. 0.0) then
               DQI = min(DEP, DQALL + QLLS/DTIME)
               DQL = DQALL - DQI
            else
               DQL = DQALL ! could happen because the PDF allows condensation in subsaturated conditions
               DQI = 0.0
            end if
         end if
         if  (DQALL .lt. 0.0) then  !net evaporation. Water evaporates first regaardless of DEP
            DQL = max(DQALL, -QLLS/DTIME)
            DQI = max(DQALL - DQL, -QILS/DTIME)
         end if
         if (DQALL .ne. 0.0)  FQI=max(min(DQI/DQALL, 1.0), 0.0)

      end if !=====

   end subroutine Bergeron_Partition

   subroutine MELTFRZ_3D ( DT, CNVFRC, SRFTYPE, TE, QL, QI )
      real, intent(in   ) :: DT, CNVFRC(:,:),SRFTYPE(:,:)
      real, intent(inout) :: TE(:,:,:), QL(:,:,:), QI(:,:,:)
      integer :: i,j,l
      do l=1,size(TE,3)
      do j=1,size(TE,2)
      do i=1,size(TE,1)
        call MELTFRZ_SC(DT, CNVFRC(i,j),SRFTYPE(i,j), TE(i,j,l), QL(i,j,l), QI(i,j,l))
      enddo
      enddo
      enddo
   end subroutine MELTFRZ_3D

   subroutine MELTFRZ_2D ( DT, CNVFRC, SRFTYPE, TE, QL, QI )
      real, intent(in   ) :: DT, CNVFRC(:,:),SRFTYPE(:,:)
      real, intent(inout) :: TE(:,:), QL(:,:), QI(:,:)
      integer :: i,j
      do j=1,size(TE,2)
      do i=1,size(TE,1)
        call MELTFRZ_SC(DT, CNVFRC(i,j),SRFTYPE(i,j), TE(i,j), QL(i,j), QI(i,j))
      enddo
      enddo
   end subroutine MELTFRZ_2D

   subroutine MELTFRZ_1D ( DT, CNVFRC, SRFTYPE, TE, QL, QI )
      real, intent(in   ) :: DT, CNVFRC(:),SRFTYPE(:)
      real, intent(inout) :: TE(:), QL(:), QI(:)
      integer :: i
      do i=1,size(TE)
        call MELTFRZ_SC(DT, CNVFRC(i),SRFTYPE(i), TE(i), QL(i), QI(i))
      enddo
   end subroutine MELTFRZ_1D

   subroutine MELTFRZ_SC( DT, CNVFRC, SRFTYPE, TE, QL, QI )
      real, intent(in   ) :: DT, CNVFRC, SRFTYPE
      real, intent(inout) :: TE,QL,QI
      real  :: fQi,dQil
      integer :: K
      ! freeze liquid first
      if ( TE <= MAPL_TICE ) then
         fQi  = ice_fraction( TE, CNVFRC, SRFTYPE )
         dQil = Ql *(1.0 - EXP( -Dt * fQi / taufrz ) )
         dQil = max(  0., dQil )
         Qi   = Qi + dQil
         Ql   = Ql - dQil
         TE   = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if
      ! melt ice instantly above 0^C
      if ( TE > MAPL_TICE ) then
         dQil = -Qi
         dQil = min(  0., dQil )
         Qi   = Qi + dQil
         Ql   = Ql - dQil
         TE   = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if
   end subroutine MELTFRZ_SC

  subroutine FILLQ2ZERO( Q, MASS, FILLQ  )

    ! New algorithm to fill the negative q values in a mass conserving way.
    ! Conservation of TPW was checked. Donifan Barahona
    ! Updated from FILLQ2ZERO, avoids the usage of scalars

    real, dimension(:,:,:),   intent(inout)  :: Q
    real, dimension(:,:,:),   intent(in)     :: MASS
    real, dimension(:,:),     intent(  out)  :: FILLQ
    real, dimension(:,:), allocatable        :: TPW1, TPW2, TPWC
    integer                                  :: IM,JM,LM, l

    IM = SIZE( Q, 1 )
    JM = SIZE( Q, 2 )
    LM = SIZE( Q, 3 )

    ALLOCATE(TPW1(IM, JM))
    ALLOCATE(TPW2(IM, JM))
    ALLOCATE(TPWC(IM, JM))

    TPW2 =0.0
    TPWC= 0.0
    TPW1 = SUM( Q*MASS, 3 )

    WHERE (Q < 0.0)
       Q=0.0
    END WHERE

    TPW2 = SUM( Q*MASS, 3 )

    WHERE (TPW2 > 0.0)
       TPWC=(TPW2-TPW1)/TPW2
    END WHERE

    do l=1,LM
       Q(:, :, l)= Q(:, :, l)*(1.0-TPWC) !reduce Q proportionally to the increase in TPW
    end do

    FILLQ = TPW2-TPW1

    DEALLOCATE(TPW1)
    DEALLOCATE(TPW2)
    DEALLOCATE(TPWC)
  end subroutine FILLQ2ZERO

  subroutine FILLQ2ZERO1( Q, MASS, FILLQ  )
    real, dimension(:,:,:),   intent(inout)  :: Q
    real, dimension(:,:,:),   intent(in)     :: MASS
    real, dimension(:,:),     intent(  out)  :: FILLQ
    integer                                  :: IM,JM,LM
    integer                                  :: I,J,K,L
    real                                     :: TPW, NEGTPW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Fills in negative q values in a mass conserving way.
    ! Conservation of TPW was checked.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IM = SIZE( Q, 1 )
    JM = SIZE( Q, 2 )
    LM = SIZE( Q, 3 )
    do j=1,JM
       do i=1,IM
          TPW = SUM( Q(i,j,:)*MASS(i,j,:) )
          NEGTPW = 0.
          do l=1,LM
             if ( Q(i,j,l) < 0.0 ) then
                NEGTPW   = NEGTPW + ( Q(i,j,l)*MASS( i,j,l ) )
                Q(i,j,l) = 0.0
             endif
          enddo
          do l=1,LM
             if ( Q(i,j,l) >= 0.0 ) then
                Q(i,j,l) = Q(i,j,l)*( 1.0+NEGTPW/(TPW-NEGTPW) )
             endif
          enddo
          FILLQ(i,j) = -NEGTPW
       end do
    end do
  end subroutine FILLQ2ZERO1

  subroutine DIAGNOSE_PRECIP_TYPE(IM, JM, LM, TPREC, RAIN_LS, RAIN_CU, RAIN, SNOW, ICE, FRZR, PTYPE, PLE, TH, PK, PKE, ZL0, LUPDATE_PRECIP_TYPE)
    integer,                        intent(in   )  :: IM, JM, LM
    real,    dimension(IM,JM),      intent(inout)  :: TPREC, RAIN_LS, RAIN_CU, RAIN, SNOW, ICE, FRZR
    real,    dimension(IM,JM),      intent(  out)  :: PTYPE
    real,    dimension(IM,JM,LM),   intent(in   )  :: TH, PK, ZL0
    real,    dimension(IM,JM,0:LM), intent(in   )  :: PLE, PKE
    logical,                        intent(in   )  :: LUPDATE_PRECIP_TYPE

    integer :: I,J,L,KTOP
    real    :: TL, NA, PA, PA2, TH_TOP, TH_BOT, TL_MEAN, Z_LAYER, ZTHICK

         PTYPE(:,:) = 0 ! default PTYPE to rain
         ! Surface Precip Type diagnostic
         !
         !   PTYPE = 0  => Rain
         !   PTYPE = 1  => Freezing Rain
         !   PTYPE = 2  => Ice Pellets (sleet)
         !   PTYPE = 3  => Rain mixed with Snow
         !   PTYPE = 4  => Snow
         !
         ! Based on Pierre Bourgouin, 1999, "A Method to Determine Precipitation Types"
         !       in WEATHER AND FORECASTING Vol 15 pp 583-592
         do J=1,JM
         do I=1,IM
          if (SNOW(I,J)+ICE(I,J)+FRZR(I,J) > 0.0) then ! Only diagnose where model has frozen precip
           PTYPE(I,J) = 4 ! Start as snow
           PA2 = -999
           ! Sweep down the column from ~300mb looking for freezing/melting layers
           KTOP = max(1,count(PLE(I,J,:) < 30000.))
           do while (KTOP < LM)
            NA = 0.0
            PA = 0.0
            TH_TOP  = TH(I,J,KTOP) ! Potential Temp at top of layer
            TL_MEAN = 0.0 ! Layer mean absolute temperature
            Z_LAYER = 0.0 ! Layer thickness
            if (TH(I,J,KTOP)*PK(I,J,KTOP) > MAPL_TICE) then
               do L=KTOP,LM-1
                  KTOP = L
                  TL = TH(I,J,L)*PK(I,J,L)
                  if (TL > MAPL_TICE) then
                      ZTHICK = TH(I,J,L) * (PKE(I,J,L) - PKE(I,J,L-1)) * cpbgrav
                      TL_MEAN = TL_MEAN + TL*ZTHICK
                      Z_LAYER = Z_LAYER + ZTHICK
                   else
                      if (Z_LAYER > 0) then
                         TL_MEAN = TL_MEAN/Z_LAYER
                         TH_BOT = TH(I,J,L)
                        ! Determine depth of the warm layer [Positive Area (PA)]
                         PA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                      endif
                      EXIT
                   endif
               enddo
            else
               do L=KTOP,LM-1
                  KTOP = L
                  TL = TH(I,J,L)*PK(I,J,L)
                  if (TL <= MAPL_TICE) then
                      ZTHICK = TH(I,J,L) * (PKE(I,J,L) - PKE(I,J,L-1)) * cpbgrav
                      TL_MEAN = TL_MEAN + TL*ZTHICK
                      Z_LAYER = Z_LAYER + ZTHICK
                  else
                     if (Z_LAYER > 0) then
                        TL_MEAN = TL_MEAN/Z_LAYER
                        TH_BOT = TH(I,J,L)
                       ! Determine depth of the freezing layer [Negative Area (NA)]
                        NA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                     endif
                     EXIT
                  endif
               enddo
            endif
            if (KTOP == LM-1) then
               TL     = TH(I,J,LM)*PK(I,J,LM)
               TH_BOT = (TL + gravbcp*ZL0(I,J,LM)) / PKE(I,J,LM)
               ZTHICK = TH(I,J,LM) * (PKE(I,J,LM) - PKE(I,J,LM-1)) * cpbgrav
               TL_MEAN = TL_MEAN + TL*ZTHICK
               Z_LAYER = Z_LAYER + ZTHICK
               if (TL > MAPL_TICE) then
                  if (Z_LAYER > 0) then
                     TL_MEAN = TL_MEAN/Z_LAYER
                     ! Determine depth of the warm layer [Positive Area (PA)]
                     PA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                  endif
               else
                  if (Z_LAYER > 0) then
                     TL_MEAN = TL_MEAN/Z_LAYER
                     ! Determine depth of the freezing layer [Negative Area (NA)]
                     NA = MAPL_CP*TL_MEAN*log( TH_TOP/TH_BOT )
                  endif
               endif
               if (PTYPE(I,J) == 4) then ! No Warm layer found above the surface yet
                 if (PA <  5.6) PTYPE(I,J) = 4 ! SNOW
                 if (PA >= 5.6) PTYPE(I,J) = 3 ! Mix of Snow Ice and Rain
                 if (PA > 13.2) PTYPE(I,J) = 0 ! RAIN
               else
                 if ( NA <  (46.0 + 0.66*PA2) ) PTYPE(I,J) = 1   ! Freezing Rain
                 if ( NA >= (46.0 + 0.66*PA2) ) PTYPE(I,J) = 1.5 ! Freezing Rain & Ice Pellets (sleet)
                 if ( NA >  (66.0 + 0.66*PA2) ) PTYPE(I,J) = 2   ! Ice Pellets (sleet)
               endif
               KTOP = LM
            else
               if (PA > 2.0) then ! Found a warm layer...
                  if (PA <  5.6) PTYPE(I,J) = 4 ! SNOW
                  if (PA >= 5.6) PTYPE(I,J) = 3 ! Mix of Snow and Rain
                  if (PA > 13.2) PTYPE(I,J) = 0 ! RAIN
                  PA2 = PA
               else ! Found a freezing layer
                  PA2 = 0
                  if ( (PTYPE(I,J) <= 3) ) then
                     if (NA > 50.0 ) PTYPE(I,J) = 1 ! Freezing Rain
                     if (NA > 200.0) PTYPE(I,J) = 2 ! Ice Pellets (sleet)
                  endif
               endif
            endif
           enddo
          endif
         enddo
         enddo

      UPDATE_PTYPE: if (LUPDATE_PRECIP_TYPE) then
         SNOW = 0.0
         WHERE (PTYPE > 2)
            SNOW = TPREC
         END WHERE
         ICE = 0.0
         WHERE ( (PTYPE == 2) .OR. (PTYPE == 1.5) )
            ICE = TPREC
         END WHERE
         FRZR = 0.0
         WHERE (PTYPE == 1)
            FRZR = TPREC
         END WHERE
         RAIN = 0.0
         WHERE ( PTYPE < 1 )
            RAIN = TPREC
            RAIN_LS = MIN(RAIN_LS,TPREC)
            RAIN_CU = MAX(TPREC-RAIN_LS,0.0)
         ELSEWHERE
            RAIN_LS = 0.0
            RAIN_CU = 0.0
         END WHERE
      endif UPDATE_PTYPE

  end subroutine DIAGNOSE_PRECIP_TYPE

  subroutine VertInterp(v2,v3,ple,pp,rc)

    real    , intent(OUT) :: v2(:,:)
    real    , intent(IN ) :: v3(:,:,:)
    real    , intent(IN ) :: ple(:,:,:)
    real    , intent(IN ) :: pp
    integer, optional, intent(OUT) :: rc

    real, dimension(size(v2,1),size(v2,2)) :: al,PT,PB
    integer k,km
    logical edge

    character*(10) :: Iam='VertInterp'

    km   = size(ple,3)-1
    edge = size(v3,3)==km+1

    _ASSERT(edge .or. size(v3,3)==km,'needs informative message')

    v2   = MAPL_UNDEF

    if(EDGE) then
       pb   = ple(:,:,km+1)
       do k=km,1,-1
          pt = ple(:,:,k)
          if(all(pb<pp)) exit
          where(pp>pt .and. pp<=pb)
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k)*al + v3(:,:,k+1)*(1.0-al)
          end where
          pb = pt
       end do
    else
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       do k=km,2,-1
          pt = 0.5*(ple(:,:,k-1)+ple(:,:,k))
          if(all(pb<pp)) exit
          where( (pp>pt.and.pp<=pb) )
             al = (pb-pp)/(pb-pt)
             v2 = v3(:,:,k-1)*al + v3(:,:,k)*(1.0-al)
          end where
          pb = pt
       end do
       pt = 0.5*(ple(:,:,km)+ple(:,:,km-1))
       pb = 0.5*(ple(:,:,km)+ple(:,:,km+1))
       where( (pp>pb.and.pp<=ple(:,:,km+1)) )
          v2 = v3(:,:,km)
       end where
      ! final protection to avoid undef
       where( v2 .eq. MAPL_UNDEF)
          v2 = v3(:,:,km)
       end where
    end if

    RETURN_(ESMF_SUCCESS)
  end subroutine VertInterp

!finds the level closets to X2=criteria
  subroutine find_l(kp, X, crit, im, jm, lm, kmin, kmax)

      real, intent(in):: crit, X(im, jm, lm)
      integer, intent (in) :: im, jm, lm, kmin, kmax
      integer, intent (out) :: kp(im, jm)
      integer :: i, j , k

        DO j = 1, jm
             DO i = 1, im
                   DO k = lm-1, kmin, -1
                     if ((X(i, j, k) .lt. crit) .and.  (X(i, j, k+1) .gt. crit))then
                       kp(i, j) =  max(min(k  + 1, kmax), kmin)
                       exit
                     end if
                    end do
               end do
       end do
  end subroutine find_l

  subroutine FIND_EIS(TH, QSAT, T, ZL0, PLE, KLCL, IM, JM, LM, LTS, EIS)
   ! !DESCRIPTION:  Returns Estimated Inversion Strength (K) according to Wood and Betherton, J.Climate, 2006
   ! Written by Donifan Barahona

    integer, intent(in) :: IM,JM,LM
    real   , dimension(IM,JM,LM)  , intent(in)  :: TH, QSAT, T, ZL0
    real   , dimension(IM,JM,0:LM), intent(in)  :: PLE
    integer, dimension(IM,JM)     , intent(in)  :: KLCL
    real   , dimension(IM,JM)     , intent(out) :: EIS, LTS

    real   , dimension(IM,JM) :: TH700, T700, Z700
    real    ::  ZLCL, QS850, T850, GAMMA850
    integer :: I, J, K

    call VertInterp(TH700, TH,log(100.0*PLE),log(70000.))
    call VertInterp( T700,  T,log(100.0*PLE),log(70000.))
    call VertInterp( Z700,ZL0,log(100.0*PLE),log(70000.))

    do J = 1, JM
       do I = 1, IM
          LTS(I,J) = TH700(I,J)-TH(I,J,LM)
          ZLCL = ZL0(I,J,KLCL(I,J)-1)

          ! Simplified single adiabat eq4 of https://doi.org/10.1175/JCLI3988.1
           T850 = 0.5*(T(I,J,LM)+T700(I,J))
          QS850 = GEOS_QSAT(T850, 850.0)

          GAMMA850 =  (1.0+(          MAPL_ALHL*QS850/(        MAPL_RGAS*T850     )))/ &
                      (1.0+(MAPL_ALHL*MAPL_ALHL*QS850/(MAPL_CP*MAPL_RVAP*T850*T850)))
          GAMMA850 =  gravbcp*(1.0-GAMMA850)

          EIS(I,J) =  LTS(I,J) - GAMMA850*(Z700(I,J) - ZLCL)

       end do
    end do

  end subroutine FIND_EIS

  FUNCTION FIND_TLCL ( tk, rh ) result( tlcl )
  ! Description:
  !    This function calculates the temperature of a parcel of air would have
  !    if lifed dry adiabatically to it's lifting condensation level (lcl).
  ! References:
  !    Bolton (1980), Monthly Weather Review, pg. 1048, Eq. 22
    IMPLICIT NONE
    REAL, INTENT ( IN ) :: tK   !~ Temperature ( K )
    REAL, INTENT ( IN ) :: rh   !~ Relative Humidity ( % )
    REAL                :: tlcl
    REAL :: denom, term1, term2
    term1 = 1.0 / ( tK - 55.0 )
    term2 = ( LOG (max(0.1,rh)/100.0)  / 2840.0 )
    denom = term1 - term2
    tlcl = ( 1.0 / denom ) + 55.0
  END FUNCTION FIND_TLCL

  function FIND_KLCL( T, Q, PL, IM, JM, LM ) result( KLCL )
    ! !DESCRIPTION:
    integer,                      intent(in) :: IM, JM, LM
    real,    dimension(IM,JM,LM), intent(in) :: T, Q, PL ! T in K, PL in mb
    integer, dimension(IM,JM)             :: KLCL

    real    :: RHSFC, TLCL, Rm, Cpm, PLCL
    integer :: I, J, L, KOFFSET

    do J=1,JM
       do I=1,IM
          RHSFC = 100.0*Q(I,J,LM)/GEOS_QSAT( T(I,J,LM), PL(I,J,LM) ) ! surface RH %
          TLCL  = FIND_TLCL(T(I,J,LM),RHSFC) ! T at LCL
          Rm    = (1.0-Q(I,J,LM))*MAPL_RGAS  + Q(I,J,LM)*MAPL_RVAP
          Cpm   = (1.0-Q(I,J,LM))*MAPL_CPDRY + Q(I,J,LM)*MAPL_CPVAP
          PLCL  = PL(I,J,LM) * ( (TLCL/T(I,J,LM))**(Cpm/Rm) ) ! P at LCL
          do L=LM,1,-1
             KLCL(I,J) = L
             if (PL(I,J,L) <= PLCL) exit
          end do
       end do
    end do

  end function FIND_KLCL

  !Find cloud top based on cloud fraction

  subroutine find_cldtop(ncol, pver, cf, kcldtop)

    integer, intent(in)  :: pver , ncol ! number of vertical layers
    real,    intent(in)  :: cf(ncol,pver)     ! midpoint potential temperatures
    integer, intent(out) ::  kcldtop
    integer              ::  kuppest, ibot, k
    real                 ::  stab,  cfcrit, cf00, cfp1


    ibot = pver-1
    kcldtop  = ibot+1
    kuppest = 20
    cfcrit = 1.0e-2


    do k =  kuppest , ibot
       cfp1 = cf(ncol, k+1)  ! qc one level down

       if ( ( cfp1  .ge. cfcrit ) ) then
          kcldtop  = k +1
          exit
       end if
    end do

    if (kcldtop .ge. ibot) then
       kcldtop = pver
       return
    endif


  end subroutine find_cldtop


!Find cloud base  for a given cloud fraction

  subroutine find_cldbase(ncol, pver, cf, kcldbase)

    integer, intent(in)  :: pver , ncol ! number of vertical layers
    real,    intent(in)  :: cf(ncol,pver)     ! midpoint potential temperatures
    integer, intent(out) ::  kcldbase
    integer              ::  kuppest, ibot, k
    real                 ::  stab,  cfcrit, cf00, cfp1


    ibot = pver-1
    kcldbase  = 20
    kuppest = 20
    cfcrit = 1.0e-3


    do k =  ibot, kuppest, -1
       cfp1 = cf(ncol, k)  !

       if ( ( cfp1  .ge. cfcrit ) ) then
          kcldbase  = k
          exit
       end if
    end do

    if (kcldbase .le. kuppest) then
       kcldbase = 1
       return
    endif


  end subroutine find_cldbase

  !DONIF Calculate the Brunt_Vaisala frequency

  !===============================================================================
  subroutine gw_prof (pcols, pver, ncol, t, pm, pi, rhoi, ni, ti, nm)
    !-----------------------------------------------------------------------
    ! Compute profiles of background state quantities for the multiple
    ! gravity wave drag parameterization.
    !
    ! The parameterization is assumed to operate only where water vapor
    ! concentrations are negligible in determining the density.
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    integer, intent(in)  :: ncol               ! number of atmospheric columns
    integer, intent(in)  :: pcols              ! number of atmospheric columns
    integer, intent(in)  :: pver               ! number of vertical layers

    !real,    intent(in)  :: u(pcols,pver)      ! midpoint zonal wind
    !real,    intent(in)  :: v(pcols,pver)      ! midpoint meridional wind
    real,    intent(in)  :: t(pcols,pver)      ! midpoint temperatures
    real,    intent(in)  :: pm(pcols,pver)     ! midpoint pressures
    real,    intent(in)  :: pi(pcols,0:pver)   ! interface pressures

    real,    intent(out) :: rhoi(pcols,0:pver) ! interface density
    real,    intent(out) :: ni(pcols,0:pver)   ! interface Brunt-Vaisalla frequency
    real,    intent(out) :: ti(pcols,0:pver)   ! interface temperature
    real,    intent(out) :: nm(pcols,pver)     ! midpoint Brunt-Vaisalla frequency

    !---------------------------Local storage-------------------------------
    integer :: ix,kx                            ! loop indexes

    real    :: dtdp
    real    :: n2, cpair, r,g                              ! Brunt-Vaisalla frequency squared
    real :: n2min   = 1.e-8
    r=MAPL_RGAS
    cpair=MAPL_CP
    g=MAPL_GRAV

    !-----------------------------------------------------------------------------
    ! Determine the interface densities and Brunt-Vaisala frequencies.
    !-----------------------------------------------------------------------------

    ! The top interface values are calculated assuming an isothermal atmosphere
    ! above the top level.
    kx = 0
    do ix = 1, ncol
       ti(ix,kx) = t(ix,kx+1)
       rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx))
       ni(ix,kx) = sqrt (g*g / (cpair*ti(ix,kx)))
    end do

    ! Interior points use centered differences
    do kx = 1, pver-1
       do ix = 1, ncol
          ti(ix,kx) = 0.5 * (t(ix,kx) + t(ix,kx+1))
          rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx))
          dtdp = (t(ix,kx+1)-t(ix,kx)) / (pm(ix,kx+1)-pm(ix,kx))
          n2 = g*g/ti(ix,kx) * (1./cpair - rhoi(ix,kx)*dtdp)
          ni(ix,kx) = sqrt (max (n2min, n2))
       end do
    end do

    ! Bottom interface uses bottom level temperature, density; next interface
    ! B-V frequency.
    kx = pver
    do ix = 1, ncol
       ti(ix,kx) = t(ix,kx)
       rhoi(ix,kx) = pi(ix,kx) / (r*ti(ix,kx))
       ni(ix,kx) = ni(ix,kx-1)
    end do

    !-----------------------------------------------------------------------------
    ! Determine the midpoint Brunt-Vaisala frequencies.
    !-----------------------------------------------------------------------------
    do kx=1,pver
       do ix=1,ncol
          nm(ix,kx) = 0.5 * (ni(ix,kx-1) + ni(ix,kx))
       end do
    end do

    return
  end subroutine gw_prof


!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
!----- module_mp_thompson_make_number_concentrations
!- Developed by H. Barnes @ NOAA/OAR/ESRL/GSL Earth Prediction Advancement Division
!-----------------------------------------------------------------------
!      Q_ice              is cloud ice mixing ratio, units of kg/m3
!      Q_cloud            is cloud water mixing ratio, units of kg/m3
!      Q_rain             is rain mixing ratio, units of kg/m3
!      temp               is air temperature in Kelvin
!      make_IceNumber     is cloud droplet number mixing ratio, units of number per m3
!      make_DropletNumber is rain number mixing ratio, units of number per kg of m3
!      make_RainNumber    is rain number mixing ratio, units of number per kg of m3
!      qnwfa              is number of water-friendly aerosols in number per kg

!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
      elemental real function make_IceNumber (Q_ice, temp)

      implicit none
      real, parameter:: ice_density = 890.0
      real, intent(in):: q_ice, temp
      integer idx_rei
      real corr, reice, deice, mui, k,  TC, lambdai

      double precision lambda

!+---+-----------------------------------------------------------------+
!..Table of lookup values of radiative effective radius of ice crystals
!.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
!.. radiation code where it is attributed to Jon Egill Kristjansson
!.. and coauthors.
!+---+-----------------------------------------------------------------+

      real, dimension(95), parameter:: retab = (/                       &
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)

      if (Q_ice == 0) then
         make_IceNumber = 0
         return
      end if

!+---+-----------------------------------------------------------------+
!..From the model 3D temperature field, subtract 179K for which
!.. index value of retab as a start.  Value of corr is for
!.. interpolating between neighboring values in the table.
!+---+-----------------------------------------------------------------+

      idx_rei = int(temp-179.)
      idx_rei = min(max(idx_rei,1),94)
      corr = temp - FLOOR(temp)
      reice = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
      deice = 2.*reice * 1.E-6

!+---+-----------------------------------------------------------------+
!..Now we have the final radiative effective size of ice (as function
!.. of temperature only).  This size represents 3rd moment divided by
!.. second moment of the ice size distribution, so we can compute a
!.. number concentration from the mean size and mass mixing ratio.
!.. The mean (radiative effective) diameter is 3./Slope for an inverse
!.. exponential size distribution.  So, starting with slope, work
!.. backwords to get number concentration.
!+---+-----------------------------------------------------------------+

      lambda = 3.0 / deice


      ! value of the dispersion parameter according to Heymsfield et al 2002, Table3.
            TC=temp-273.15

            TC=MIN(MAX(TC, -70.0), -15.0)

            if (TC .gt. -27.0) then
               lambdai=6.8*exp(-0.096*TC)
            else
               lambdai=24.8*exp(-0.049*TC)
            end if

            mui=(0.13*(lambdai**0.64))-2.




      k =  (mui+3)*(mui*3)/(mui+2)/(mui+1)
      make_IceNumber = k*Q_ice * lambda*lambda*lambda / (MAPL_PI*Ice_density)

!+---+-----------------------------------------------------------------+
!..Example1: Common ice size coming from Thompson scheme is about 30 microns.
!.. An example ice mixing ratio could be 0.001 g/kg for a temperature of -50C.
!.. Remember to convert both into MKS units.  This gives N_ice=357652 per kg.
!..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
!.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
!.. which is 28 crystals per liter of air if the air density is 1.0.
!+---+-----------------------------------------------------------------+

      return
      end function make_IceNumber

!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+

      elemental real function make_DropletNumber (Q_cloud, qnwfa)

      implicit none
      real, intent(in):: q_cloud, qnwfa
      real, parameter:: am_r = MAPL_PI*1000./6.
      real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
                      504,720,990,1320,1716,2184,2730,3360,4080,4896/)
      double precision:: lambda, qnc
      real:: q_nwfa, x1, xDc
      integer:: nu_c

      if (Q_cloud .le. 0.) then
         make_DropletNumber = 0.0
         return
      end if

!+---+

      q_nwfa = MAX(99.E6, MIN(qnwfa,5.E10))
      nu_c = MAX(2, MIN(NINT(2.5E10/q_nwfa), 15))

      x1 = MAX(1., MIN(q_nwfa*1.E-9, 10.)) - 1.
      xDc = (30. - x1*20./9.) * 1.E-6

      lambda = (4.0D0 + nu_c) / xDc
      qnc = Q_cloud / g_ratio(nu_c) * lambda*lambda*lambda / am_r
      make_DropletNumber = SNGL(qnc)


      return
      end function make_DropletNumber

!+---+-----------------------------------------------------------------+

      elemental real function make_RainNumber (Q_rain, temp)

      IMPLICIT NONE

      real, intent(in):: Q_rain, temp
      double precision:: lambda, N0, qnr
      real, parameter:: am_r = MAPL_PI*1000./6.

      if (Q_rain == 0) then
         make_RainNumber = 0
         return
      end if

      !+---+-----------------------------------------------------------------+
      !.. Not thrilled with it, but set Y-intercept parameter to Marshal-Palmer value
      !.. that basically assumes melting snow becomes typical rain. However, for
      !.. -2C < T < 0C, make linear increase in exponent to attempt to keep
      !.. supercooled collision-coalescence (warm-rain) similar to drizzle rather
      !.. than bigger rain drops.  While this could also exist at T>0C, it is
      !.. more difficult to assume it directly from having mass and not number.
      !+---+-----------------------------------------------------------------+

      N0 = 8.E6

      if (temp .le. 271.15) then
         N0 = 8.E8
      elseif (temp .gt. 271.15 .and. temp.lt.273.15) then
         N0 = 8. * 10**(279.15-temp)
      endif

      lambda = SQRT(SQRT(N0*am_r*6.0/Q_rain))
      qnr = Q_rain / 6.0 * lambda*lambda*lambda / am_r
      make_RainNumber = SNGL(qnr)

      return
      end function make_RainNumber

!+---+-----------------------------------------------------------------+

   SUBROUTINE  dissipative_ke_heating(IM,JM,LM &
                             ,mass,us,vs,du,dv,ttend)

   implicit none
   integer                       ,intent (in ) :: IM,JM,LM
   real   , dimension (IM,JM,LM) ,intent (in ) :: mass  ! delp/gravity
   real   , dimension (IM,JM,LM) ,intent (in ) :: us,vs ! winds prior to cumulus process
   real   , dimension (IM,JM,LM) ,intent (in ) :: du,dv ! wind tendency due to cumulus process
   real   , dimension (IM,JM,LM) ,intent (out) :: ttend ! K/s

   real :: dts,fp,dp,fpi,ke(LM)
   integer ::i,j,l

! since kinetic energy is being dissipated, add heating accordingly (from ECMWF)
!
   ttend = 0.0
   do j=1,JM
     do i=1,IM
          dts=0.
          fpi=0.
          do l=1,LM
             !total KE dissiptaion estimate
             dts = dts - (du(i,j,l)*us(i,j,l)+dv(i,j,l)*vs(i,j,l))*mass(i,j,l)
             !
             ! fpi needed for calcualtion of conversion to pot. energyintegrated
             ke(l) = sqrt(du(i,j,l)*du(i,j,l) + dv(i,j,l)*dv(i,j,l))
             fpi = fpi + ke(l)*mass(i,j,l)
          enddo
          if(fpi.gt.0.)then
             do l=1,LM
                ttend(i,j,l) = (ke(l)/fpi)*dts*(1.0/MAPL_CP)
             enddo
          endif
     enddo
   enddo

  end SUBROUTINE  dissipative_ke_heating





!+---+-----------------------------------------------------------------+

subroutine update_cld( &
         DT          , &
         ALPHA       , &
         PDFFLAG     , &
         CNVFRC      , &
         SRFTYPE     , &
         PL          , &
         QV          , &
         QCl         , &
         QAl         , &
         QCi         , &
         QAi         , &
         TE          , &
         CF          , &
         AF          , &
         SCICE       , &
         NI          , &
         NL          , &
         RHcmicro)

      real, intent(in)    :: DT,ALPHA,PL,CNVFRC,SRFTYPE
      integer, intent(in) :: pdfflag
      real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF,SCICE,NI,NL,RHCmicro

      ! internal arrays
      real :: CFO
      real :: QT

      real :: QSx,DQsx

      real :: QCx, QC, QA

      real :: QX, QSLIQ, QSICE, CFALL, DQx, FQA, DELQ

      real :: SHOM
      real :: maxalpha =  0.4


    !  maxalpha=1.0-minrhcrit

      QC = QCl + QCi
      QA = QAl + QAi
      QT  =  QC  + QA + QV  !Total water after microphysics
      CFALL  = AF+CF
      FQA = 0.0
      if (QA+QC .gt. tiny(1.0))  FQA=QA/(QA+QC)

      SHOM=2.349-(TE/259.0) !hom threeshold Si according to Ren & McKenzie, 2005

      !================================================
      ! First find the cloud fraction that would correspond to the current condensate
      QSLIQ  = GEOS_QsatLQU( TE, PL*100.0 , DQ=DQx )
      QSICE  = GEOS_QsatICE( TE, PL*100.0 , DQ=DQX )

      
      IF (QCl + QAl .gt. 0.) then 
        QSx =  QSLIQ
      ELSEIF (QCi + QAi.gt. 0.) then 
        QSx =  QSICE        
      ELSE
		 DQSx  = GEOS_DQSAT( TE, PL, QSAT=QSx )
      end if

      if (TE .gt. 240.0)   SCICE = 1.0
      QCx=QC+QA
      QX=QT-QSx*SCICE
      CFo=0.

      !====== recalculate QX if too low and SCICE<SHOM
        if ((QX .gt. QCx) .and. (QCx .gt. 0.0)) then
           QX=QT-QSx*SHOM
       end if

      !=======================

     DELQ=max(min(2.0*maxalpha*QSx, 0.5*QT), 1.0e-12)

      if  ((QX .le. QCx)  .and. (QCx .gt. tiny(1.0)))  then
         CFo =  (1.0+SQRT(1.0-(QX/QCx)))
         if (CFo .gt. 1.e-6) then
            CFo = min(1.0/CFo, 1.0)
            DELQ=  2.0*QCx/(CFo*CFo)
         else
            CFo = 0.0
         end if
      elseif (QCx .gt. tiny(1.0)) then
         !   CFo = 1.0  !Outside of distribution but still with condensate
        DELQ=max(min(2.0*maxalpha*QSx, 0.5*QT), 1.0e-12)
        CFo = SQRT(2.0*QCx/DELQ)
      else
        CFo = 0.0
      end if

      if  (QSx .gt. tiny(1.0)) then
         RHCmicro = SCICE - 0.5*DELQ/Qsx
      else
         RHCmicro = 1.0-ALPHA
      end if

      RHCmicro =  max(min(RHCmicro, 0.99), 0.6)

      CFALL   = max(CFo, 0.0)
      CFALL   = min(CFo, 1.0)

      CF=CFALL*(1.0-FQA)
      AF=CFALL*FQA


   end subroutine update_cld




   subroutine meltfrz_inst2M  ( IM, JM, LM,    &
         TE       , &
         QCL       , &
         QAL       , &
         QCI        , &
         QAI        , &
         NL       , &
         NI            )

      real ,   intent(inout), dimension(:,:,:)   :: TE,QCL,QCI, QAL, QAI, NI, NL
      integer, intent(in) :: IM, JM, LM

      real ,   dimension(im,jm,lm)              :: dQil, DQmax, QLTOT, QITOT, dNil, FQA
      real :: T_ICE_ALL =  240.
      real :: T_ICE_MAX =  273.

      QITOT= QCI+QAI
      QLTOT=QCL + QAL
      FQA = 0.0

      where (QITOT+QLTOT .gt. tiny(0.0))
         FQA= (QAI+QAL)/(QITOT+QLTOT)
      end where


      dQil = 0.0
      dNil =0.0
      DQmax  = 0.0

      ! freeze liquid instantaneosly below -40 C
      where( TE <= T_ICE_ALL )
         DQmax = (T_ICE_ALL - TE)*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil = min(QLTOT , DQmax)
      end where

      where ((dQil .le. DQmax) .and. (dQil .gt. 0.0))
         dNil = NL
      end where

      where ((dQil .gt. DQmax) .and. (dQil .gt. tiny(0.0)))
         dNil  =  NL*DQmax/dQil
      end where

      dQil = max(  0., dQil )
      QITOT = max(QITOT + dQil, 0.0)
      QLTOT= max(QLTOT -  dQil, 0.0)
      NL  = NL - dNil
      NI   = NI  + dNil
      TE   = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP

      dQil = 0.0
      dNil =0.0
      DQmax  = 0.0

      ! melt ice instantly above 0^C
      where( TE > T_ICE_MAX )
         DQmax =  (TE-T_ICE_MAX) *MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil = min(QITOT, DQmax)
      endwhere

      where ((dQil .le. DQmax) .and. (dQil .gt. 0.0))
         dNil = NI
      end where
      where ((dQil .gt. DQmax) .and. (dQil .gt. tiny(0.0)))
         dNil  =  NI*DQmax/dQil
      end where
      dQil = max(  0., dQil )
      QLTOT =  max(QLTOT+ dQil, 0.)
      QITOT = max(QITOT - dQil, 0.)
      NL  = NL + dNil
      NI   = NI  - dNil

      TE   = TE - (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP

      QCI   = QITOT*(1.0-FQA)
      QAI = QITOT*FQA
      QCL   = QLTOT*(1.0-FQA)
      QAL = QLTOT*FQA

   end subroutine meltfrz_inst2M

      subroutine FIX_NEGATIVE_PRECIP(QRAIN, QSNOW, QGRAUPEL)
          real, dimension(:,:,:), intent(inout) :: QRAIN, QSNOW, QGRAUPEL

          WHERE (QRAIN < 1.e-8)
            QRAIN = 0.0
          END WHERE

          WHERE (QSNOW < 1.e-8)
            QSNOW = 0.0
          END WHERE

          WHERE (QGRAUPEL < 1.e-8)
            QGRAUPEL = 0.0
          END WHERE

    end subroutine FIX_NEGATIVE_PRECIP

   subroutine REDISTRIBUTE_CLOUDS(CF, QL, QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, QV, TE)
      real, dimension(:,:,:), intent(inout) :: CF, QL, QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, QV, TE

      ! Liquid
      QLLS = QLLS + (QL - (QLCN+QLLS))
      WHERE (QLLS < 0.0)
        QLCN = QLCN + QLLS
        QLLS = 0.0
      END WHERE            
      WHERE (QLCN < 1.E-8)
       ! QLCN is negative so the signs here -/+ are reversed
        QV = QV - QLCN
        TE = TE + (alhlbcp)*QLCN
        QLCN = 0.0
      END WHERE

      ! Ice
      QILS = QILS + (QI - (QICN+QILS))
      WHERE (QILS < 0.0)
        QICN = QICN + QILS
        QILS = 0.0
      END WHERE
      WHERE (QICN < 1.E-8)
       ! QLCN is negative so the signs here -/+ are reversed
        QV = QV - QICN
        TE = TE + (alhsbcp)*QICN
        QICN = 0.0
      END WHERE

      ! Cloud
      CLLS = min(1.0,CLLS + (CF - (CLCN+CLLS)))
      WHERE (CLLS < 0.0)
        CLCN = min(1.0,CLCN + CLLS)
        CLLS = 0.0
      END WHERE
      WHERE (CLCN < 1.E-8)
        CLCN = 0.
      END WHERE

      ! Evaporate liquid/ice where clouds are gone
      WHERE (CLLS < 1.E-8)
        QV = QV + QLLS + QILS
        TE = TE - (alhlbcp)*QLLS - (alhsbcp)*QILS
        CLLS = 0.
        QLLS = 0.
        QILS = 0.
      END WHERE
      WHERE (CLCN < 1.E-8)
        QV = QV + QLCN + QICN
        TE = TE - (alhlbcp)*QLCN - (alhsbcp)*QICN
        CLCN = 0.
        QLCN = 0.
        QICN = 0.
      END WHERE

   end subroutine REDISTRIBUTE_CLOUDS

 subroutine cs_interpolator(is, ie, js, je, km, qin, zout, wz, qout, qmin)
 integer,  intent(in):: is, ie, js, je, km
 real, intent(in):: zout, qmin
 real, intent(in):: qin(is:ie,js:je,km)
 real, intent(in):: wz(is:ie,js:je,km+1)
 real, intent(out):: qout(is:ie,js:je)
! local:
 real:: qe(is:ie,km+1)
 real, dimension(is:ie,km):: q2, dz
 real:: s0, a6
 integer:: i,j,k

 do j=js,je

   do i=is,ie
      do k=1,km
         dz(i,k) = wz(i,j,k) - wz(i,j,k+1)
         q2(i,k) = qin(i,j,k)
      enddo
   enddo

   call cs_prof(q2, dz, qe, km, is, ie, 1)

   do i=is,ie
      if( zout >= wz(i,j,1) ) then
! Higher than the top:
          qout(i,j) = qe(i,1)
      elseif ( zout <= wz(i,j,km+1) ) then
          qout(i,j) = qe(i,km+1)
      else
          do k=1,km
             if ( zout<=wz(i,j,k) .and. zout >= wz(i,j,k+1) ) then
! PPM distribution: f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s <= 1 )
                  a6 = 3.*(2.*q2(i,k) - (qe(i,k)+qe(i,k+1)))
                  s0 = (wz(i,j,k)-zout) / dz(i,k)
                  qout(i,j) = qe(i,k) + s0*(qe(i,k+1)-qe(i,k)+a6*(1.-s0))
                  go to 500
             endif
          enddo
      endif
500   qout(i,j) = max(qmin, qout(i,j))
   enddo
 enddo

! Send_data here

 end subroutine cs_interpolator

 subroutine cs_prof(q2, delp, q, km, i1, i2, iv)
! Latest: Dec 2015 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2, km
 integer, intent(in):: iv
 real, intent(in)   :: q2(i1:i2,km)
 real, intent(in)   :: delp(i1:i2,km)     ! layer pressure thickness
 real, intent(out):: q(i1:i2,km+1)
!-----------------------------------------------------------------------
 real  gam(i1:i2,km)
 real   d4(i1:i2)
 real   bet, a_bot, grat
 integer i, k

  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*q2(i,1) + q2(i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(q2(i,k-1)+d4(i)*q2(i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo

  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*q2(i,km)+q2(i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(q2(i,1), q2(i,2)) )
     q(i,2) = max( q(i,2), min(q2(i,1), q2(i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = q2(i,k) - q2(i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(q2(i,k-1),q2(i,k)) )
             q(i,k) = max( q(i,k), min(q2(i,k-1),q2(i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(q2(i,k-1),q2(i,k)))
          else
! There exists a local min
               q(i,k) = min(q(i,k), max(q2(i,k-1),q2(i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(q2(i,km-1), q2(i,km)) )
     q(i,km) = max( q(i,km), min(q2(i,km-1), q2(i,km)) )
  enddo

 end subroutine cs_prof

   integer function FIND_KLID (plid, ple, rc)  RESULT(klid)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)     :: plid ! pressure lid [hPa]
   real, dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]

! !OUTPUT PARAMETERS:
   integer, intent(out) :: rc ! return code; 0 - all is good
!                                            1 - bad

! !DESCRIPTION: Finds corresponding vertical index for defined pressure lid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k, j, i
   real :: plid_, diff, refDiff
   real, allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]

!EOP
!----------------------------------------------------------------------------------
!  Begin...
   klid = 1
   rc = 0

!  convert from hPa to Pa
   plid_ = plid*100.0

   allocate(pres(ubound(ple,3)))

!  find pressure at each model level
   do k = 1, ubound(ple,3)
      pres(k) = ple(1,1,k)
   end do

!  find smallest absolute difference between plid and average pressure at each model level
   refDiff = 150000.0
   do k = 1, ubound(ple,3)
      diff = abs(pres(k) - plid_)
      if (diff < refDiff) then
         klid = k
         refDiff = diff
      end if
   end do

!  Check to make sure that all pressures at (i,j) were the same
   do j = 1, ubound(ple,2)
      do i = 1, ubound(ple,1)
         if (pres(klid) /= ple(i,j,klid)) then
            rc = 1
            return
         end if
      end do
   end do

   end function FIND_KLID


#ifdef THOM_REF
!+---+-----------------------------------------------------------------+
!>\ingroup aathompson
!! Compute radar reflectivity assuming 10 cm wavelength radar and using
!! Rayleigh approximation.  Only complication is melted snow/graupel
!! which we treat as water-coated ice spheres and use Uli Blahak's
!! library of routines.  The meltwater fraction is simply the amount
!! of frozen species remaining from what initially existed at the
!! melting level interface.

      subroutine calc_refl10cm (qv1d, qr1d, nr1d, qs1d, qg1d, &
               t1d, p1d, dBZ, rand1, kts, kte, ii, jj, melti,       &
               vt_dBZ, first_time_step, ktopin, kbotin)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, INTENT(IN):: rand1
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                          qv1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ
      REAL, DIMENSION(kts:kte), OPTIONAL, INTENT(INOUT):: vt_dBZ
      LOGICAL, OPTIONAL, INTENT(IN) :: first_time_step
      INTEGER, OPTIONAL, INTENT(IN) :: ktopin, kbotin

!..Local variables
      LOGICAL :: do_vt_dBZ
      LOGICAL :: allow_wet_graupel
      LOGICAL :: allow_wet_snow
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho, rhof
      REAL, DIMENSION(kts:kte):: rr, nr, rs, rg

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
      REAL, DIMENSION(kts:kte):: mvd_r
      REAL, DIMENSION(kts:kte):: smob, smo2, smoc, smoz
      REAL:: oM3, M0, Mrat, slam1, slam2, xDs
      REAL:: ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
      REAL:: vtr_dbz_wt, vts_dbz_wt, vtg_dbz_wt

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

      DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamr, lamg
      REAL:: a_, b_, loga_, tc0, SR
      DOUBLE PRECISION:: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0, ktop, kbot, kdwn, n
      LOGICAL, INTENT(IN):: melti
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL:: xslw1, ygra1, zans1

      REAL, PARAMETER, PRIVATE :: R1 = 1.E-12
      REAL, PARAMETER, PRIVATE :: R2 = 1.E-6
!..Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
      REAL, PARAMETER, PRIVATE:: rho_not = 101325.0/(287.05*298.0)

!+---+
      if (present(ktopin) .and. present(kbotin)) then
         ktop=ktopin
         kbot=kbotin
         if (ktop < kbot) then
           kdwn= 1
         else
           kdwn=-1
         endif
      else
         ktop=kte
         kbot=kts
         kdwn=-1
      endif

      if (present(vt_dBZ) .and. present(first_time_step)) then
         do_vt_dBZ = .true.
         if (first_time_step) then
!           no bright banding, to be consistent with hydrometeor retrieval in GSI
            allow_wet_snow = .false.
         else
            allow_wet_snow = .true.
         endif
         allow_wet_graupel = .false.
      else
         do_vt_dBZ = .false.
         allow_wet_snow = .true.
         allow_wet_graupel = .false.
      endif

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rhof(k) = SQRT(RHO_NOT/rho(k))
         if (qr1d(k) .gt. R1) then
            rr(k) = qr1d(k)*rho(k)
            nr(k) = MAX(R2, nr1d(k)*rho(k))
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            ilamr(k) = 1./lamr
            N0_r(k) = nr(k)*org2*lamr**cre(2)
            mvd_r(k) = (3.0 + mu_r + 0.672) * ilamr(k)
            L_qr(k) = .true.
         else
            rr(k) = R1
            nr(k) = R1
            mvd_r(k) = 50.E-6
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R2) then
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R2) then
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         smo2(k) = 0.
         smob(k) = 0.
         smoc(k) = 0.
         smoz(k) = 0.
      enddo
      if (ANY(L_qs .eqv. .true.)) then
      do k = kts, kte
         if (.not. L_qs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
     &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
     &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
     &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
     &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
     &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
     &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
     &         + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif

!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
     &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
     &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
     &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_

!..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
     &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
     &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(3)*cse(3)*cse(3)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
     &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
         smoz(k) = a_ * smo2(k)**b_
      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

      if (ANY(L_qg .eqv. .true.)) then
      do k = ktop, kbot, kdwn
         ygra1 = alog10(max(1.E-9, rg(k)))
         zans1 = 3.4 + 2./7.*(ygra1+8.) + rand1
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      k_0 = kbot
      if ( melti ) then
        K_LOOP:do k = ktop+kdwn, kbot, kdwn
          if ((temp(k).gt.273.15) .and. L_qr(k)                         &
     &                            .and. (L_qs(k-kdwn).or.L_qg(k-kdwn)) ) then
             if (kdwn < 0) then
                k_0 = MAX(k-kdwn, k_0)
             else
                k_0 = MIN(k-kdwn, k_0)
             endif
             EXIT K_LOOP
          endif
        enddo K_LOOP
      endif
!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
     &                           * (am_s/900.0)*(am_s/900.0)*smoz(k)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
     &                              * (am_g/900.0)*(am_g/900.0)         &
     &                              * N0_g(k)*cgg(4)*ilamg(k)**cge(4)
      enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (.not. iiwarm .and. melti .and. k_0.ge.2) then
       do k = k_0+kdwn, kbot, kdwn

!..Reflectivity contributed by melting snow
          if (allow_wet_snow .and. L_qs(k) .and. L_qs(k_0) ) then
           SR = MAX(0.01, MIN(1.0 - rs(k)/(rs(k) + rr(k)), 0.99))
           fmelt_s = DBLE(SR*SR)
           eta = 0.d0
           oM3 = 1./smoc(k)
           M0 = (smob(k)*oM3)
           Mrat = smob(k)*M0*M0*M0
           slam1 = M0 * Lam0
           slam2 = M0 * Lam1
           do n = 1, nrbins
              x = am_s * xxDs(n)**bm_s
              call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
     &              fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
     &              CBACK, mixingrulestring_s, matrixstring_s,          &
     &              inclusionstring_s, hoststring_s,                    &
     &              hostmatrixstring_s, hostinclusionstring_s)
              f_d = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
     &              + Kap1*(M0*xxDs(n))**mu_s * DEXP(-slam2*xxDs(n)))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

!..Reflectivity contributed by melting graupel
          if (allow_wet_graupel .and. L_qg(k) .and. L_qg(k_0) ) then
           SR = MAX(0.01, MIN(1.0 - rg(k)/(rg(k) + rr(k)), 0.99))
           fmelt_g = DBLE(SR*SR)
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = am_g * xxDg(n)**bm_g
              call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
     &              fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
     &              CBACK, mixingrulestring_g, matrixstring_g,          &
     &              inclusionstring_g, hoststring_g,                    &
     &              hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**mu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = ktop, kbot, kdwn
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo

!..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).
      if (do_vt_dBZ) then
         do k = ktop, kbot, kdwn
            vt_dBZ(k) = 1.E-3
            if (rs(k).gt.R2) then
             Mrat = smob(k) / smoc(k)
             ils1 = 1./(Mrat*Lam0 + fv_s)
             ils2 = 1./(Mrat*Lam1 + fv_s)
             t1_vts = Kap0*csg(5)*ils1**cse(5)
             t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
             ils1 = 1./(Mrat*Lam0)
             ils2 = 1./(Mrat*Lam1)
             t3_vts = Kap0*csg(6)*ils1**cse(6)
             t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
             vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
             if (temp(k).ge.273.15 .and. temp(k).lt.275.15) then
                vts_dbz_wt = vts_dbz_wt*1.5
             elseif (temp(k).ge.275.15) then
                vts_dbz_wt = vts_dbz_wt*2.0
             endif
            else
             vts_dbz_wt = 1.E-3
            endif

            if (rr(k).gt.R1) then
             lamr = 1./ilamr(k)
             vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
                        / (crg(4)*lamr**(-cre(4)))
            else
             vtr_dbz_wt = 1.E-3
            endif

            if (rg(k).gt.R2) then
             lamg = 1./ilamg(k)
             vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
                        / (cgg(4)*lamg**(-cge(4)))
            else
             vtg_dbz_wt = 1.E-3
            endif

            vt_dBZ(k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
                         + vtg_dbz_wt*ze_graupel(k))                        &
                         / (ze_rain(k)+ze_snow(k)+ze_graupel(k))
         enddo
      endif

      end subroutine calc_refl10cm
#endif

end module GEOSmoist_Process_Library
