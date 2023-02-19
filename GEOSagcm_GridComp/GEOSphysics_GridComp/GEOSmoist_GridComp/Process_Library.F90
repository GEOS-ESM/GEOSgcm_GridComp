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
  ! In anvil/convective clouds
  real, parameter :: aT_ICE_ALL = 252.16
  real, parameter :: aT_ICE_MAX = 268.16
  real, parameter :: aICEFRPWR  = 2.0
  ! Over snow/ice SRF_TYPE = 2
  real, parameter :: iT_ICE_ALL = 236.16
  real, parameter :: iT_ICE_MAX = 261.16
  real, parameter :: iICEFRPWR  = 6.0
  ! Over Land     SRF_TYPE = 1
  real, parameter :: lT_ICE_ALL = 239.16
  real, parameter :: lT_ICE_MAX = 261.16
  real, parameter :: lICEFRPWR  = 2.0
  ! Over Oceans   SRF_TYPE = 0
  real, parameter :: oT_ICE_ALL = 238.16
  real, parameter :: oT_ICE_MAX = 263.16
  real, parameter :: oICEFRPWR  = 4.0

 ! parameters
  real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
  real, parameter :: K_COND  =  2.4e-2    ! J m**-1 s**-1 K**-1
  real, parameter :: DIFFU   =  2.2e-5    ! m**2 s**-1
 ! LDRADIUS4
  ! Liquid  based on DOI 10.1088/1748-9326/3/4/045021
  real, parameter :: RHO_W   = 1000.0  ! Density of liquid water in kg/m^3
  real, parameter :: Lbe     = 1./3. - 0.14
  real, parameter :: Lbx     = 0.07*1.e3*(3./(4.*MAPL_PI*RHO_W*1.e-3))**(1./3.)
                             ! 0.07 is a dispersion effect and eqs are in cgs units
  ! Ice
  real, parameter :: RHO_I   =  916.8  ! Density of ice crystal in kg/m^3

  ! combined constantc
  real, parameter :: cpbgrav = MAPL_CP/MAPL_GRAV
  real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP
  real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
  real, parameter :: alhfbcp = MAPL_ALHF/MAPL_CP
  real, parameter :: alhsbcp = MAPL_ALHS/MAPL_CP

  ! defined to determine CNV_FRACTION
  real    :: CNV_FRACTION_MIN
  real    :: CNV_FRACTION_MAX
  real    :: CNV_FRACTION_EXP

 ! Storage of aerosol properties for activation
  type(AerProps), allocatable, dimension (:,:,:) :: AeroProps

  ! Tracer Bundle things for convection
  type CNV_Tracer_Type
      real, pointer              :: Q(:,:,:) => null()
      real                       :: fscav = 0.0
      real                       :: Vect_Hcts(4)
      character(len=ESMF_MAXSTR) :: QNAME ! Tracer Name
      character(len=ESMF_MAXSTR) :: CNAME ! Component Name
  end type CNV_Tracer_Type
  type(CNV_Tracer_Type), allocatable :: CNV_Tracers(:)

  public :: AeroProps
  public :: CNV_Tracer_Type, CNV_Tracers, CNV_Tracers_Init
  public :: ICE_FRACTION, EVAP3, SUBL3, LDRADIUS4, BUOYANCY, RADCOUPLE, FIX_UP_CLOUDS
  public :: hystpdf, fix_up_clouds_2M
  public :: FILLQ2ZERO, FILLQ2ZERO1
  public :: MELTFRZ
  public :: DIAGNOSE_PRECIP_TYPE
  public :: VertInterp
  public :: find_l, FIND_EIS, FIND_KLCL
  public :: find_cldtop, find_cldbase, gw_prof
  public :: make_IceNumber, make_DropletNumber
  public :: dissipative_ke_heating
  public :: pdffrac, pdfcondensate, partition_dblgss
  public :: CNV_FRACTION_MIN, CNV_FRACTION_MAX, CNV_FRACTION_EXP
 
  contains

  subroutine CNV_Tracers_Init(TR, RC)
    type (ESMF_FieldBundle), intent(inout) :: TR
    integer,       optional, intent(inout) :: RC
   ! Local
    type (ESMF_Field) :: FIELD
    integer :: TotalTracers, FriendlyTracers
    logical :: isPresent, isFriendly
    integer :: ind, N, F
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
               ! Get items for the wet removal parameterization for gases based on the Henry's Law
               !-------------------------------------------------------------------------------------
               CNV_Tracers(F)%Vect_Hcts(:)=-99.
               call ESMF_AttributeGet(FIELD, "SetofHenryLawCts", isPresent=isPresent,  RC=STATUS); VERIFY_(STATUS)
               if (isPresent) then
                  call ESMF_AttributeGet(FIELD, "SetofHenryLawCts", CNV_Tracers(F)%Vect_Hcts,  RC=STATUS); VERIFY_(STATUS)
               end if
               ! Get component and tracer names
               !-------------------------------------------------------------------------------------
               ind= index(QNAME, '::')
               if (ind > 0) then
                  CNV_Tracers(F)%CNAME = trim(QNAME(1:ind-1))  ! Component name (e.g., GOCART, CARMA)
                  CNV_Tracers(F)%QNAME = trim(QNAME(ind+2:))
               end if
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
            end if
         end if
      enddo
    end if

    deallocate(QNAMES)

  end subroutine CNV_Tracers_Init

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

      ! Anvil clouds
      ! Anvil-Convective sigmoidal function like figure 6(right)
      ! Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
        ICEFRCT_C  = 0.00
        if ( TEMP <= aT_ICE_ALL ) then
           ICEFRCT_C = 1.000
        else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
           ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
        end if
        ICEFRCT_C = MIN(ICEFRCT_C,1.00)
        ICEFRCT_C = MAX(ICEFRCT_C,0.00)
        ICEFRCT_C = ICEFRCT_C**aICEFRPWR
#ifdef MODIS_ICE_POLY
      ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384) 
        tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
        ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
        ICEFRCT_M = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
#else
  ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
      if (SRF_TYPE == 2.0) then
        ! Over snow/ice
        ICEFRCT_M  = 0.00
        if ( TEMP <= iT_ICE_ALL ) then
           ICEFRCT_M = 1.000
        else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
           ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL ) ) )
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
#endif
      ! Combine the Convective and MODIS functions
        ICEFRCT  = ICEFRCT_M*(1.0-CNV_FRACTION) + ICEFRCT_C*(CNV_FRACTION)

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
         EVAP = MIN( EVAP , QL  )
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
         SUBL = MIN( SUBL , QI  )
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

       !- air density (kg/m^3)
       RHO = (100.*PL) / (MAPL_RGAS*TE )
       IF(ITYPE == LIQUID) THEN

       !- liquid cloud effective radius ----- 
          !- liquid water content
          WC = 1.e3*RHO*QC ! air density [g/m3] * liquid cloud mixing ratio [kg/kg]
          !- cloud drop number concentration
          !- from the aerosol model + ....
          NNX = NNL*1.e-6 ! #/cm3
          !- radius in meters
          !- [liu&daum, 2000 and 2005. liu et al 2008]
          RADIUS = MIN(60.e-6,MAX(2.5e-6, 1.e-6*Lbx*(WC/NNX)**Lbe))

       ELSEIF(ITYPE == ICE) THEN

       !- ice cloud effective radius ----- 
          !- ice water content
          WC = 1.e3*RHO*QC ! air density [g/m3] * ice cloud mixing ratio [kg/kg]
          !- cloud ice number concentration #/m3
          !- from the aerosol model + ....
          NNX = NNI*1.e-6 ! #/m3
          !- radius in meters
          !------ice cloud effective radius ----- [klaus wyser, 1998]
          if(TE>MAPL_TICE .or. QC <=0.) then
            BB = -2.
          else
            BB = -2. + log10(WC/50.)*(1.e-3*(MAPL_TICE-TE)**1.5)
          endif
          BB     = MIN((MAX(BB,-6.)),-2.)
          RADIUS = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
          RADIUS = MIN(150.e-6,MAX(5.e-6, 1.e-6*RADIUS))

      ELSE
        STOP "WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2"
      ENDIF

   end function LDRADIUS4

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
      !-BRAMS formulation     
      RAD_RL = LDRADIUS4(PL,TE,RAD_QL,NL,NI,1)
     ! apply limits
      RAD_RL = MAX( MIN_RL, MIN(RAD_RL*FAC_RL, MAX_RL) )

    ! ICE RADII
     !-BRAMS formulation  
      RAD_RI = LDRADIUS4(PL,TE,RAD_QI,NL,NI,2)
    ! apply limits
      RAD_RI = MAX( MIN_RI, MIN(RAD_RI*FAC_RI, MAX_RI) )

   end subroutine RADCOUPLE

   subroutine  FIX_UP_CLOUDS( &
         QV, &
         TE, &
         QLC,&
         QIC,&
         CF, &
         QLA,&
         QIA,&
         AF  )

      real, intent(inout) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA

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
         NG)

      real, intent(inout), dimension(:,:,:) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA, QR, QS, QG 
      real, intent(inout), dimension(:,:,:) :: NI, NL, NS, NR, NG
    
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
 subroutine partition_dblgss( dt,           &  ! IN
                              tabs,         &  ! INOUT
                              qwv,          &
                              qc,           &
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
                              mffrc,        &
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

   real, intent(in   )  :: DT          ! timestep [s]
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
   real, intent(in   )  :: mffrc       ! total EDMF updraft fraction
!   real, intent(inout)  :: qi         ! ice condensate [kg kg-1]
   real, intent(  out)  :: cld_sgs     ! cloud fraction
   real, intent(inout)  ::    PDF_A           ! fractional area of 1st gaussian
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
   real, parameter :: tauskew = 3600. 

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
          if (w_sec > 0.0) then
            sqrtw2   = sqrt(w_sec)
            Skew_w   = w3var / (sqrtw2*sqrtw2*sqrtw2)
          else
            sqrtw2   = w_thresh
            Skew_w   = 0.
          endif
          if (thlsec > 0.0) then
            sqrtthl  = sqrt(thlsec)
            skew_thl = hl3 / sqrtthl**3
          else
            sqrtthl  = 0.0
            skew_thl = 0.
          endif
          if (qwsec > 0.0) then
            sqrtqt   = sqrt(qwsec)
            skew_qw =  qt3/sqrtqt**3
          else
            sqrtqt   = 1e-4*total_water
            skew_qw  = 0.
          endif

! Find parameters of the double Gaussian PDF of vertical velocity

!          aterm = pdf_a

         if (use_aterm_memory/=0) then   ! use memory in aterm and qt skewness
          aterm = pdf_a

          if (mffrc>=1e-3) then                ! if active updraft this timestep
            if (aterm<0.5) then                ! if distribution is skewed (recent updrafts)
              aterm = max(mffrc,aterm*max(1.-DT/tauskew,0.0))
            else                               ! if distribution unskewed
              aterm = mffrc
            end if
          else                                 ! if no active updraft
            if (aterm.lt.0.5 .and. aterm.gt.1e-3) then  ! but there is residual skewness
              aterm = aterm*max(1.-DT/tauskew,0.0)
            else
              aterm = 0.5
            end if
          end if

         else  ! don't use memory in aterm and qt skewness

           aterm = mffrc
           aterm = max(1e-3,min(0.99,aterm))
           if (mffrc.le.1e-3) aterm = 0.5
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

            if (aterm .lt. 1e-3 .or. aterm.gt.0.499 .or. Skew_qw.eq.0.) then ! if no residual skewness
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
              wrk1 = min(10.,skew_qw*sqrtqt**3)   ! third moment qt
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

          pdf_a = aterm

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

          esval1_1 = 0.
          esval1_2 = 0.
          esval2_1 = 0.
          esval2_2 = 0.
          om1      = 1.
          om2      = 1.

! Partition based on temperature for the first plume

          IF (Tl1_1 >= tbgmax) THEN
            esval1_1 = MAPL_EQsat(Tl1_1)
            lstarn1  = lcond
          ELSE IF (Tl1_1 < tbgmin) THEN
            esval1_1 = MAPL_EQsat(Tl1_1,OverIce=.TRUE.)
            lstarn1  = lsub
          ELSE
            esval1_1 = MAPL_EQsat(Tl1_1)
            esval2_1 = MAPL_EQsat(Tl1_1,OverIce=.TRUE.)
            om1      = max(0.,min(1.,a_bg*(Tl1_1-tbgmin)))
            lstarn1  = lcond + (1.-om1)*lfus
          ENDIF

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

            IF (Tl1_2 < tbgmin) THEN
              esval1_2 = MAPL_EQsat(Tl1_2,OverIce=.TRUE.)
              lstarn2  = lsub
            ELSE IF (Tl1_2 >= tbgmax) THEN
              esval1_2 = MAPL_EQsat(Tl1_2)
              lstarn2  = lcond
            ELSE
              esval1_2 = MAPL_EQsat(Tl1_2)
              esval2_2 = MAPL_EQsat(Tl1_2,OverIce=.TRUE.)
              om2      = max(0.,min(1.,a_bg*(Tl1_2-tbgmin)))
              lstarn2  = lcond + (1.-om2)*lfus
            ENDIF

            qs2   =     om2  * (0.622*esval1_2/max(esval1_2,pval-0.378*esval1_2))    &
                  + (1.-om2) * (0.622*esval2_2/max(esval2_2,pval-0.378*esval2_2))

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
          diag_frac = aterm*C1 + onema*C2

          om1 = max(0.,min(1.,(Tl1_1-tbgmin)*a_bg))
          om2 = max(0.,min(1.,(Tl1_2-tbgmin)*a_bg))

          qn1 = min(qn1,qw1_1)
          qn2 = min(qn2,qw1_2)

          ql1 = qn1*om1
          ql2 = qn2*om2

          qi1 = qn1 - ql1
          qi2 = qn2 - ql2

          diag_qn = min(max(0.0, aterm*qn1 + onema*qn2), total_water)
          diag_ql = min(max(0.0, aterm*ql1 + onema*ql2), diag_qn)
          diag_qi = diag_qn - diag_ql

!!! temporary
!          if (abs(qc-diag_qn)>0.001) print *,'SHOC: t=',tabs,' s1=',s1,' qn1=',qn1,' qs1=',qs1,' qt1=',qw1_1


! Update temperature variable based on diagnosed cloud properties
          om1         = max(0.,min(1.,(tabs-tbgmin)*a_bg))
          lstarn1     = lcond + (1.-om1)*lfus
!          tabs = thl_first - gamaz + fac_cond*(diag_ql) &
!                            + fac_sub *(diag_qi) !&
                    !  + tkesbdiss(i,j,k) * (dtn/cp)      ! tke dissipative heating
! Update moisture fields



         qc      = diag_ql + diag_qi
!         qi      = diag_qi
!         qwv     = total_water - diag_qn
         cld_sgs = diag_frac

         if (sqrtqt>0.0 .AND. sqrtw2>0.0) then
            rwqt = (1.-0.5)*wqwsec/(sqrtqt*sqrtw2)
!            rwqt = (wqwsec)/(sqrtqt*sqrtw2)
!            rwqt = max(-1.,min(1.,pdf_rwqt))
         else
            rwqt = 0.0
         end if
         if (sqrtthl>0.0 .AND. sqrtw2>0.0) then
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

!                          + ((lstarn1/cp)-thv(i,j,k))*0.5*(wqp_sec(i,j,kd)+wqp_sec(i,j,ku))

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
         QCl         , &
         QAl         , &
         QCi         , &
         QAi         , &
         TE          , &
         CF          , &
         AF          , &
         NL          , &
         NI          , &
         WHL         , &
         WQT         , &
!         wqtfac      , &
!         whlfac      , &
         HL2         , &
         QT2         , &
         HLQT        , & 
         W3          , &
         W2          , &
         MFQT3       , &
         MFHL3       , &
         MF_FRC      , &
         PDF_A,      &  ! can remove these after development
         PDFITERS,   &
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
         USE_BERGERON )

      real, intent(in)    :: DT,ALPHA,PL,ZL
      integer, intent(in) :: PDFSHAPE
      real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF,PDF_A
      real, intent(in)    :: NL,NI,CNVFRC,SRF_TYPE
      real, intent(in)    :: WHL,WQT,HL2,QT2,HLQT,W3,W2,MF_FRC,MFQT3,MFHL3
#ifdef PDFDIAG
      real, intent(out)   :: PDF_SIGW1, PDF_SIGW2, PDF_W1, PDF_W2, &
                             PDF_SIGHL1, PDF_SIGHL2, PDF_HL1, PDF_HL2, &
                             PDF_SIGQT1, PDF_SIGQT2, PDF_QT1, PDF_QT2, &
                             PDF_RHLQT,  PDF_RWHL, PDF_RWQT
#endif
      real, intent(out)   :: WTHV2, WQL
      real, intent(out)   :: PDFITERS
      logical, intent(in) :: needs_preexisting, USE_BERGERON

      ! internal arrays
      real :: QCO, QVO, CFO, QAO, TAU,HL
      real :: QT, QMX, QMN, DQ, sigmaqt1, sigmaqt2

      real :: TEO,QSx,DQsx,QS,DQs

      real :: TEp, QSp, CFp, QVp, QCp
      real :: TEn, QSn, CFn, QVn, QCn

      real :: QCx, QVx, CFx, QAx, QC, QA, fQi
      real :: dQAi, dQAl, dQCi, dQCl, Nfac, NLv, NIv 

!      real :: fQip

      real :: tmpARR
      real :: ALHX, DQCALL
      ! internal scalars
      integer :: N, nmax

      QC = QCl + QCi
      QA = QAl + QAi
      QT  =  QC  + QA + QV  !Total water after microphysics
      tmpARR = 0.0
      nmax =  20
      QAx = 0.0

      if ( AF < 1.0 )  tmpARR = 1./(1.-AF)

      TEo = TE

      DQSx  = GEOS_DQSAT( TE, PL, QSAT=QSx )
      CFx = CF*tmpARR
      QCx = QC*tmpARR
      QVx = ( QV - QSx*AF )*tmpARR

      if ( AF >= 1.0 )  QVx = QSx*1.e-4 
      if ( AF >  0.0 )  QAx = QA/AF

      QT  = QCx + QVx

      TEp = TEo
      QSn = QSx
      TEn = TEo
      CFn = CFx
      QVn = QVx
      QCn = QCx
      DQS = DQSx

      fQi = ice_fraction( TE, CNVFRC,SRF_TYPE )
      ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

      HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn

      do n=1,nmax

         QVp = QVn
         QCp = QCn
         CFp = CFn
         TEp = TEn

         if(PDFSHAPE.lt.2) then

            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn

         elseif(PDFSHAPE.eq.2) then  ! triangular
            ! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
            ! for triangular, skewed r : sigmaqt1 < sigmaqt2
            ! try: skewed right below 500 mb
            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn
!         elseif(pdffrac .eq. 3) then ! single gaussian

         elseif(PDFSHAPE .eq. 4) then !lognormal (sigma is dimmensionless)
            sigmaqt1 =  max(ALPHA/sqrt(3.0), 0.001)
         endif

         if (PDFSHAPE.lt.5) then
           call pdffrac(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,CFn)
           call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)
         elseif (PDFSHAPE.eq.5) then

            ! Update the liquid water static energy
            ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS
            HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn

           call partition_dblgss(DT/nmax,           &
                                 TEn,          &
                                 QVn,          &
                                 QCn,          &
                                 0.0,          & ! assume OMEGA=0
                                 ZL,           &
                                 PL*100.,      &
                                 QT,           &
                                 HL,          &
                                 WHL,         &
                                 WQT,         &
                                 HL2,         &
                                 QT2,         &
                                 HLQT,        & 
                                 W3,           &
                                 W2,           &
                                 MFQT3,        &
                                 MFHL3,        &
                                 MF_FRC,       &
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
           CF  = CFn * ( 1.-AF)
           Nfac = 100.*PL*R_AIR/TEp !density times conversion factor
           NLv = NL/Nfac
           NIv = NI/Nfac
           call Bergeron_Partition( &         !Microphysically-based partitions the new condensate
                 DT               , &
                 PL               , &
                 TEp              , &
                 QT               , &
                 QCi              , &
                 QAi              , &
                 QCl              , &
                 QAl              , &
                 CF               , &
                 AF               , &
                 NLv              , &
                 NIv              , &
                 DQCALL           , &
                 fQi              , & 
                 CNVFRC,SRF_TYPE  , &
                 needs_preexisting)
         ELSE
           fQi = ice_fraction( TEp, CNVFRC,SRF_TYPE )
         ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! These lines represent adjustments
         ! to anvil condensate due to the 
         ! assumption of a stationary TOTAL 
         ! water PDF subject to a varying 
         ! QSAT value during the iteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         if ( AF > 0. ) then
            QAo = QAx  ! + QSx - QS 
         else
            QAo = 0.
         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

         if(PDFSHAPE.eq.1) then 
            QCn = QCp + ( QCn - QCp ) / ( 1. - (CFn * (ALPHA-1.) - (QCn/QSn))*DQS*ALHX/MAPL_CP)             
         elseif(PDFSHAPE.eq.2 .or. PDFSHAPE.eq.5) then
            ! This next line needs correcting - need proper d(del qc)/dT derivative for triangular
            ! for now, just use relaxation of 1/2.
            if (n.ne.nmax) QCn = QCp + ( QCn - QCp ) *0.5
         endif

         QVn = QVp - (QCn - QCp)
         TEn = TEp + (1.0-fQi)*(MAPL_ALHL/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
               +      fQi* (MAPL_ALHS/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

         PDFITERS = n
         if (abs(Ten - Tep) .lt. 0.00001) exit 

         DQS  = GEOS_DQSAT( TEn, PL, QSAT=QSn )

      enddo ! qsat iteration

      CFo = CFn
      CF = CFn
      QCo = QCn
      QVo = QVn
      TEo = TEn

      ! Update prognostic variables.  Deal with special case of AF=1
      ! Temporary variables QCo, QAo become updated grid means.
      if ( AF < 1.0 ) then
         CF  = CFo * ( 1.-AF)
         QCo = QCo * ( 1.-AF)
         QAo = QAo *   AF  
      else

         ! Special case AF=1, i.e., box filled with anvil. 
         !   - Note: no guarantee QV_box > QS_box
         CF  = 0.          ! Remove any other cloud
         QAo = QA  + QC    ! Add any LS condensate to anvil type
         QCo = 0.          ! Remove same from LS   
         QT  = QAo + QV    ! Total water
         ! Now set anvil condensate to any excess of total water 
         ! over QSx (saturation value at top)
         QAo = MAX( QT - QSx, 0. )
      end if

      ! Now take {\em New} condensate and partition into ice and liquid
      ! taking care to keep both >=0 separately. New condensate can be
      ! less than old, so $\Delta$ can be < 0.

      dQCl = 0.0
      dQCi = 0.0
      dQAl = 0.0
      dQAi = 0.0

      !large scale   

      QCx   = QCo - QC
      if  (QCx .lt. 0.0) then  !net evaporation. Water evaporates first
         dQCl = max(QCx, -QCl)   
         dQCi = max(QCx - dQCl, -QCi)
      else
         dQCl  = (1.0-fQi)*QCx
         dQCi  =    fQi  * QCx
      end if

      !Anvil   
      QAx   = QAo - QA

      if  (QAx .lt. 0.0) then  !net evaporation. Water evaporates first
         dQAl = max(QAx, -QAl)   
         dQAi = max(QAx - dQAl, -QAi)
      else            
         dQAl  =  (1.0-fQi)*QAx
         dQAi  = QAx*fQi
      end if

      ! Clean-up cloud if fractions are too small
      if ( AF < 1.e-5 ) then
         dQAi = -QAi
         dQAl = -QAl
      end if
      if ( CF < 1.e-5 ) then
         dQCi = -QCi
         dQCl = -QCl
      end if

      QAi    = QAi + dQAi
      QAl    = QAl + dQAl
      QCi    = QCi + dQCi
      QCl    = QCl + dQCl
      QV     = QV  - ( dQAi+dQCi+dQAl+dQCl) 

      TE  = TE + (MAPL_ALHL*( dQAi+dQCi+dQAl+dQCl)+MAPL_ALHF*(dQAi+dQCi))/ MAPL_CP

      ! We need to take care of situations where QS moves past QA
      ! during QSAT iteration. This should be only when QA/AF is small
      ! to begin with. Effect is to make QAo negative. So, we 
      ! "evaporate" offending QA's
      !
      ! We get rid of anvil fraction also, although strictly
      ! speaking, PDF-wise, we should not do this.
      if ( QAo <= 0. ) then
         QV  = QV + QAi + QAl
         TE  = TE - (MAPL_ALHS/MAPL_CP)*QAi - (MAPL_ALHL/MAPL_CP)*QAl
         QAi = 0.
         QAl = 0.
         AF  = 0.  
      end if

   end subroutine hystpdf

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
      DEP=0.0 
      QI = QILS + QICN !neccesary because NI is for convective and large scale 
      QL = QLLS +QLCN
      QTOT=QI+QL
      FQA = 0.0
      if (QTOT .gt. 0.0) FQA = (QICN+QILS)/QTOT
      NIX= (1.0-FQA)*NI

      DQALL=DQALL/DTIME                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      CFALL= min(CF+AF, 1.0)
      TC=TE-MAPL_TICE
      fQI_0 = fQI

      !Completelely glaciated cloud:
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
                  fQi  =   ice_fraction( TE, CNVFRC, SRF_TYPE )
              end if                      
              return 
         end if 
         
         QVINC=  QV 
         QSLIQ  = GEOS_QsatLQU( TE, PL*100.0 , DQ=DQSL )
         QSICE  = GEOS_QsatICE( TE, PL*100.0 , DQ=DQSI )
         QVINC =MIN(QVINC, QSLIQ) !limit to below water saturation 

         ! Calculate deposition onto preexisting ice 

         DIFF=(0.211*1013.25/(PL+0.1))*(((TE+0.1)/MAPL_TICE)**1.94)*1e-4  !From Seinfeld and Pandis 2006
         DENAIR=PL*100.0/MAPL_RGAS/TE
         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
         LHcorr = ( 1.0 + DQSI*MAPL_ALHS/MAPL_CP) !must be ice deposition

         if  ((NIX .gt. 1.0) .and. (QILS .gt. 1.0e-10)) then 
            DC=max((QILS/(NIX*DENICE*MAPL_PI))**(0.333), 20.0e-6) !Assumme monodisperse size dsitribution 
         else
            DC = 20.0e-6
         end if

         TEFF= NIX*DENAIR*2.0*MAPL_PI*DIFF*DC/LHcorr ! 1/Dep time scale 

         DEP=0.0
         if ((TEFF .gt. 0.0) .and. (QILS .gt. 1.0e-14)) then 
            AUX =max(min(DTIME*TEFF, 20.0), 0.0)
            DEP=(QVINC-QSICE)*(1.0-EXP(-AUX))/DTIME
         end if
         DEP=MAX(DEP, -QILS/DTIME) !only existing ice can be sublimated

         DQI = 0.0
         DQL = 0.0
         FQI=0.0
         !Partition DQALL accounting for Bergeron-Findensen process
         if  (DQALL .ge. 0.0) then !net condensation. Note: do not allow bergeron with QLCN
            if (DEP .gt. 0.0) then 
               DQI = min(DEP, DQALL + QLLS/DTIME)
               DQL = DQALL - DQI
            else
               DQL=DQALL ! could happen because the PDF allows condensation in subsaturated conditions
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
      real  ::  taufrz
      integer :: K
      ! freeze liquid first
      if ( TE <= MAPL_TICE ) then
         fQi  = ice_fraction( TE, CNVFRC, SRFTYPE )
         taufrz = 450.
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
          GAMMA850 =  (MAPL_GRAV/MAPL_CP)*(1.0-GAMMA850)

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
      real, parameter:: pi = 3.1415926536
      real, intent(in):: q_ice, temp
      integer idx_rei
      real corr, reice, deice
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
      corr = temp - int(temp)
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
      make_IceNumber = Q_ice * lambda*lambda*lambda / (PI*Ice_density)

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
      real, parameter:: am_r = 3.1415927*1000./6.
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

end module GEOSmoist_Process_Library
