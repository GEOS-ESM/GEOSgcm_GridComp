! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

module GEOSmoist_Process_Library

  use MAPL

  implicit none
  private

  character(len=ESMF_MAXSTR)              :: IAm="GEOSmoist_Process_Library"
  integer                                 :: STATUS

  interface ICE_FRACTION
    module procedure ICE_FRACTION_3D
    module procedure ICE_FRACTION_2D
    module procedure ICE_FRACTION_1D
    module procedure ICE_FRACTION_SC
  end interface ICE_FRACTION

 ! parameters
  real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
  real, parameter :: K_COND  =  2.4e-2    ! J m**-1 s**-1 K**-1
  real, parameter :: DIFFU   =  2.2e-5    ! m**2 s**-1
  real, parameter :: RHO_I   =  916.8     ! Density of ice crystal in kg/m^3
  real, parameter :: RHO_W   = 1000.0     ! Density of liquid water in kg/m^3
  real, parameter :: r13  = 1./3.
  real, parameter :: be   = r13 - 0.11
  real, parameter :: aewc = 0.13*(3./(4.*MAPL_PI*RHO_W*1.e3))**r13
  real, parameter :: aeic = 0.13*(3./(4.*MAPL_PI*RHO_I*1.e3))**r13

  public :: ICE_FRACTION, EVAP3, SUBL3, LDRADIUS4, BUOYANCY, RADCOUPLE, FIX_UP_CLOUDS

  contains

  function ICE_FRACTION_3D (TEMP) RESULT(ICEFRCT)
      real, intent(in) :: TEMP(:,:,:)
      real :: ICEFRCT(size(TEMP,1),size(TEMP,2),size(TEMP,3))
      integer :: i,j,l
      do l=1,size(TEMP,3)
      do j=1,size(TEMP,2)
      do i=1,size(TEMP,1)
        ICEFRCT(i,j,l) = ICE_FRACTION_SC(TEMP(i,j,l))
      enddo
      enddo
      enddo
  end function ICE_FRACTION_3D

  function ICE_FRACTION_2D (TEMP) RESULT(ICEFRCT)
      real, intent(in) :: TEMP(:,:)
      real :: ICEFRCT(size(TEMP,1),size(TEMP,2))
      integer :: i,j
      do j=1,size(TEMP,2)
      do i=1,size(TEMP,1)
        ICEFRCT(i,j) = ICE_FRACTION_SC(TEMP(i,j))
      enddo
      enddo
  end function ICE_FRACTION_2D

  function ICE_FRACTION_1D (TEMP) RESULT(ICEFRCT)
      real, intent(in) :: TEMP(:)
      real :: ICEFRCT(size(TEMP))
      integer :: i
      do i=1,size(TEMP)
        ICEFRCT(i) = ICE_FRACTION_SC(TEMP(i))
      enddo
  end function ICE_FRACTION_1D

  function ICE_FRACTION_SC (TEMP) RESULT(ICEFRCT)
      real, intent(in) :: TEMP
      real             :: ICEFRCT
      real             :: tc, ptc

      ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384) 
      tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
      ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
      ICEFRCT = 1.0 - (1.0/(1.0 + exp(-1*ptc)))

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

      real :: ES,RADIUS,K1,K2,TEFF,QCm,EVAP,RHx,QC  !,QS

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
      
      if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
         TEFF   =   ( RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
      else
         TEFF   = 0.0 ! -999.
      end if

      EVAP = a_eff*QL*DT*TEFF
      EVAP = MIN( EVAP , QL  )

      QC=QL+QI
      if (QC > 0.) then
         F    = F * ( QC - EVAP ) / QC
      end if

      QV   = QV   + EVAP
      QL   = QL   - EVAP
      TE   = TE   - (MAPL_ALHL/MAPL_CP)*EVAP

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

      real :: ES,RADIUS,K1,K2,TEFF,QCm,SUBL,RHx,QC !, QS

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

      K1 = (MAPL_ALHL**2) * RHO_W / ( K_COND*MAPL_RVAP*(TE**2))

      K2 = MAPL_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

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
      
      if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
         TEFF   =   ( RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
      else
         TEFF   = 0.0 ! -999.
      end if

      SUBL = a_eff*QI*DT*TEFF
      SUBL = MIN( SUBL , QI  )

      QC=QL+QI
      if (QC > 0.) then
         F    = F * ( QC - SUBL ) / QC
      end if

      QV   = QV   + SUBL
      QI   = QI   - SUBL
      TE   = TE   - (MAPL_ALHS/MAPL_CP)*SUBL

   end subroutine SUBL3

   function LDRADIUS4(PL,TE,QC,NNL,NNI,ITYPE) RESULT(RADIUS)

       REAL   , INTENT(IN) :: TE,PL,QC,NNL,NNI
       INTEGER, INTENT(IN) :: ITYPE
       REAL  :: RADIUS
       INTEGER, PARAMETER  :: CLOUD=1, ICE=2
       REAL :: NNX,RHO,BB,WC

       !- air density (kg/m^3)
       RHO = 100.*PL / (MAPL_RGAS*TE )
       IF(ITYPE == CLOUD) THEN

       !- liquid cloud effective radius ----- 
          !- [liu&daum, 2000 and 2005. liu et al 2008]
          !- liquid water content
          WC = 1.e3*RHO*QC  !g/m3
          !- cloud drop number concentration #/m3
          !- from the aerosol model + ....
          NNX = MAX(NNL,1.e3)
          !- radius in meters
          RADIUS = MIN(60.e-6,MAX(2.5e-6,aewc * (WC/NNX)**be))

       ELSEIF(ITYPE == ICE) THEN

       !- ice cloud effective radius ----- 
        !- ice water content
         WC = 1.e3*RHO*QC  !g/m3
        !------ice cloud effective radius ----- [klaus wyser, 1998]
         if(TE>MAPL_TICE .or. QC <=0.) then
            BB = -2.
         else
            BB = -2. + log10(WC/50.)*(1.e-3*(MAPL_TICE-TE)**1.5)
         endif
         BB     = MIN((MAX(BB,-6.)),-2.)
         RADIUS = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
         RADIUS = RADIUS * 1.e-6 !- convert to meter

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

    BUOY(:,:,LM) =  T(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,LM) + (MAPL_ALHL/MAPL_CP)*Q(:,:,LM)

    do L=LM-1,1,-1
       BUOY(:,:,L) = BUOY(:,:,LM) - (T(:,:,L) + (MAPL_GRAV/MAPL_CP)*ZLO(:,:,L) + (MAPL_ALHL/MAPL_CP)*QS(:,:,L))
       BUOY(:,:,L) = MAPL_GRAV*BUOY(:,:,L) / ( (1.+ (MAPL_ALHL/MAPL_CP)*DQS(:,:,L))*T(:,:,L) )
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
         TE  = TE - (MAPL_ALHL/MAPL_CP)*QLA - (MAPL_ALHS/MAPL_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if

      ! Fix if LS cloud fraction too small
      if ( CF < 1.E-5 ) then
         QV = QV + QLC + QIC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC - (MAPL_ALHS/MAPL_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if
      
      ! LS LIQUID too small
      if ( QLC  < 1.E-8 ) then
         QV = QV + QLC 
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC
         QLC = 0.
      end if
      ! LS ICE too small
      if ( QIC  < 1.E-8 ) then
         QV = QV + QIC 
         TE = TE - (MAPL_ALHS/MAPL_CP)*QIC
         QIC = 0.
      end if

      ! Anvil LIQUID too small
      if ( QLA  < 1.E-8 ) then
         QV = QV + QLA 
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLA
         QLA = 0.
      end if
      ! Anvil ICE too small
      if ( QIA  < 1.E-8 ) then
         QV = QV + QIA 
         TE = TE - (MAPL_ALHS/MAPL_CP)*QIA
         QIA = 0.
      end if

      ! Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
      if ( ( QLA + QIA ) < 1.E-8 ) then
         QV = QV + QLA + QIA
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLA - (MAPL_ALHS/MAPL_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if
      ! Ditto if LS cloud LIQUID+ICE too small
      if ( ( QLC + QIC ) < 1.E-8 ) then
         QV = QV + QLC + QIC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC - (MAPL_ALHS/MAPL_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if

   end subroutine FIX_UP_CLOUDS

end module GEOSmoist_Process_Library
