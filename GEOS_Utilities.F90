
!  $Id$

!BOP

! !MODULE: GEOS_Utils -- A Module to containing computational utilities

  module GEOS_UtilsMod

! !USES:

  use MAPL_ConstantsMod
!   use pFlogger
!   use ieee_arithmetic

! #include "qsatlqu.code"
! #include "qsatice.code"
! #include "esatlqu.code"
! #include "esatice.code"
! #include "trisolve.code"
! #include "trilu.code"

  implicit none
  private


! !PUBLIC MEMBER FUNCTIONS:

  public GEOS_QsatSet

  public GEOS_QsatLQU
  public GEOS_QsatICE
  public GEOS_Qsat
  public GEOS_DQsat

  public ESINIT_v2

!   public GEOS_TRILU
!   public GEOS_TRISOLVE

!EOP

  interface GEOS_QsatICE
     module procedure QSATICE0
     module procedure QSATICE1
     module procedure QSATICE2
     module procedure QSATICE3
  end interface

  interface GEOS_QsatLQU
     module procedure QSATLQU0
     module procedure QSATLQU1
     module procedure QSATLQU2
     module procedure QSATLQU3
  end interface

  interface GEOS_DQsat
     module procedure DQSAT0
     module procedure DQSAT1
     module procedure DQSAT2
     module procedure DQSAT3
  end interface

  interface GEOS_Qsat
     module procedure QSAT0
     module procedure QSAT1
     module procedure QSAT2
     module procedure QSAT3
  end interface

!   interface GEOS_TRILU
!      module procedure GEOS_TRILU1
!      module procedure GEOS_TRILU2
!      module procedure GEOS_TRILU3
!   end interface

!   interface GEOS_TRISOLVE
!      module procedure GEOS_TRISOLVE1
!      module procedure GEOS_TRISOLVE2
!      module procedure GEOS_TRISOLVE3
!   end interface


      real,    parameter :: ESFAC = MAPL_H2OMW/MAPL_AIRMW
      real,    parameter :: MAX_MIXING_RATIO = 1.  
      real,    parameter :: ZEROC   = MAPL_TICE

      real,    parameter :: TMINTBL    =  150.0
      real,    parameter :: TMAXTBL    =  333.0
      integer, parameter :: DEGSUBS    =  100
      real,    parameter :: ERFAC      = (DEGSUBS/ESFAC)
      real,    parameter :: DELTA_T    =  1.0 / DEGSUBS
      integer, parameter :: TABLESIZE  =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
      real,    parameter :: TMIX       = -20.

      logical, save      :: UTBL       = .true.
      integer, save      :: TYPE       =  1

      logical, save      :: FIRST      = .true.

      real,    save      :: ESTFRZ
      real,    save      :: ESTLQU

      real,    save      :: ESTBLE(TABLESIZE)
      real,    save      :: ESTBLW(TABLESIZE)
      real,    save      :: ESTBLX(TABLESIZE)

      real,    parameter :: TMINSTR = -95.
      real,    parameter :: TSTARR1 = -75.
      real,    parameter :: TSTARR2 = -65.
      real,    parameter :: TSTARR3 = -50.
      real,    parameter :: TSTARR4 = -40.
      real,    parameter :: TMAXSTR = +60.

      real*8,  parameter :: B6 = 6.136820929E-11*100.0
      real*8,  parameter :: B5 = 2.034080948E-8 *100.0
      real*8,  parameter :: B4 = 3.031240396E-6 *100.0
      real*8,  parameter :: B3 = 2.650648471E-4 *100.0
      real*8,  parameter :: B2 = 1.428945805E-2 *100.0
      real*8,  parameter :: B1 = 4.436518521E-1 *100.0
      real*8,  parameter :: B0 = 6.107799961E+0 *100.0
      real*8,  parameter :: BI6= 1.838826904E-10*100.0
      real*8,  parameter :: BI5= 4.838803174E-8 *100.0
      real*8,  parameter :: BI4= 5.824720280E-6 *100.0
      real*8,  parameter :: BI3= 4.176223716E-4 *100.0
      real*8,  parameter :: BI2= 1.886013408E-2 *100.0
      real*8,  parameter :: BI1= 5.034698970E-1 *100.0
      real*8,  parameter :: BI0= 6.109177956E+0 *100.0
      real*8,  parameter :: S16= 0.516000335E-11*100.0
      real*8,  parameter :: S15= 0.276961083E-8 *100.0
      real*8,  parameter :: S14= 0.623439266E-6 *100.0
      real*8,  parameter :: S13= 0.754129933E-4 *100.0
      real*8,  parameter :: S12= 0.517609116E-2 *100.0
      real*8,  parameter :: S11= 0.191372282E+0 *100.0
      real*8,  parameter :: S10= 0.298152339E+1 *100.0
      real*8,  parameter :: S26= 0.314296723E-10*100.0
      real*8,  parameter :: S25= 0.132243858E-7 *100.0
      real*8,  parameter :: S24= 0.236279781E-5 *100.0
      real*8,  parameter :: S23= 0.230325039E-3 *100.0
      real*8,  parameter :: S22= 0.129690326E-1 *100.0
      real*8,  parameter :: S21= 0.401390832E+0 *100.0
      real*8,  parameter :: S20= 0.535098336E+1 *100.0


      real*8, parameter  :: DI(0:3)=(/ 57518.5606E08, 2.01889049, 3.56654, 20.947031 /)
      real*8, parameter  :: CI(0:3)=(/ 9.550426, -5723.265, 3.53068, -.00728332 /)
      real*8, parameter  :: DL(6)=(/  -7.902980, 5.02808, -1.3816, 11.344, 8.1328, -3.49149 /)
      real*8, parameter  :: LOGPS = 3.005714898  ! log10(1013.246)
      real*8, parameter  :: TS = 373.16
      real*8, parameter  :: CL(0:9)=(/54.842763, -6763.22, -4.21000, .000367, &
                       .0415, 218.8,  53.878000, -1331.22, -9.44523, .014025  /)


      real, save :: TMINLQU    =  ZEROC - 40.0
      real, save :: TMINICE    =  ZEROC + TMINSTR

      ! class(Logger), pointer :: lgr
      logical :: debugIsEnabled
!$acc declare create(UTBL, ESTBLW, FIRST, ESTBLE, ESTFRZ, TMINLQU, TMINICE, ESTBLX, ESTLQU, TYPE)
  contains

!BOPI

! !IROUTINE: GEOS_QsatLqu Computes saturation specific humidity over
!            liquid water.
! !IROUTINE: GEOS_QsatIce Computes saturation specific humidity over
!            frozen water.

! !INTERFACE:

!    function GEOS_QsatLqu(TL,PL,DQ) result(QS)
!    function GEOS_QsatIce(TL,PL,DQ) result(QS)
!
! Overloads:
!
!      real,                               intent(IN)               :: TL
!      logical,                  optional, intent(IN)               :: PL
!      real,                     optional, intent(OUT)              :: DQ
!      real                                                         :: QS
!
!      real,                               intent(IN)               :: TL(:)
!      logical,                  optional, intent(IN)               :: PL(:)
!      real,                     optional, intent(OUT)              :: DQ(:)
!      real, dimension(size(TL,1))                                  :: QS
!
!      real,                               intent(IN)               :: TL(:,:)
!      logical,                  optional, intent(IN)               :: PL(:,:)
!      real,                     optional, intent(OUT)              :: DQ(:,:)
!      real, dimension(size(TL,1),size(TL,2))                       :: QS
!
!      real,                               intent(IN)               :: TL(:,:,:)
!      logical,                  optional, intent(IN)               :: PL(:,:,:)
!      real,                     optional, intent(OUT)              :: DQ(:,:,:)
!      real, dimension(size(TL,1),size(TL,2),size(TL,3))            :: QS
!
!

! !DESCRIPTION:  Uses various formulations of the saturation
!                vapor pressure to compute the saturation specific 
!    humidity and, optionally, its derivative with respect to temperature
!    for temperature TL and pressure PL. If PL is not present
!    it returns the saturation vapor pressure and, optionally, its derivative. 
!
!    All pressures are in Pascals and all temperatures in Kelvins.
!
!    The choice of saturation vapor pressure formulation is controlled by  GEOS_QsatSet.
!    Three choices are currently supported: The CAM formulation,
!    Murphy and Koop (2005, QJRMS), and the Staar formulation from NSIPP-1.
!    The default is Starr. All three are valid up to 333K. Above the 
!    freezing point, GEOS_QsatIce returns values at the freezing point.
!    Murphy and Koop is valid down to 150K, for both liquid and ice.
!    The other two are valid down to 178K for ice and 233K for super-cooled liquid. 
!
!    Another choice is whether to use the exact formulation
!    or a table look-up. This can also be controlled with GEOS_QsatSet.
!    The default is to do a table look-up. The tables are generated
!    at 0.1C intervals, controlled by parameter DEGSUBS=10.
! 
!    
!EOPI


       function QSATLQU0(TL,PL,DQ) result(QS)
!$acc routine seq
         real,              intent(IN) :: TL
         real, optional,    intent(IN) :: PL
         real, optional,    intent(OUT):: DQ
         real    :: QS

         real    :: TI
         real    :: DD
         real    :: TT
         real    :: DDQ
         integer :: IT
#define TX TL
#define PX PL
#define EX QS
#define DX DQ
#include "qsatlqu.code"
#undef  DX
#undef  TX
#undef  EX
#undef  PX
         return
       end function QSATLQU0

       function QSATLQU1(TL,PL,DQ) result(QS)
         real,              intent(IN) :: TL(:)
         real, optional,    intent(IN) :: PL(:)
         real, optional,    intent(OUT):: DQ(:)
         real    :: QS(SIZE(TL,1))
         integer :: I
         real    :: TI
         real    :: TT
         real    :: DDQ
         real    :: DD
         integer :: IT
         do I=1,size(TL,1)
#define TX TL(I)
#define PX PL(I)
#define EX QS(I)
#define DX DQ(I)
#include "qsatlqu.code"
#undef  DX
#undef  TX
#undef  PX
#undef  EX
         end do
       end function QSATLQU1

       function QSATLQU2(TL,PL,DQ) result(QS)
         real,              intent(IN) :: TL(:,:)
         real, optional,    intent(IN) :: PL(:,:)
         real, optional,    intent(OUT):: DQ(:,:)
         real    :: QS(SIZE(TL,1),SIZE(TL,2))
         integer :: I, J
         real    :: TI
         real    :: TT
         real    :: DDQ
         real    :: DD
         integer :: IT
         do J=1,size(TL,2)
            do I=1,size(TL,1)
#define TX TL(I,J)
#define PX PL(I,J)
#define EX QS(I,J)
#define DX DQ(I,J)
#include "qsatlqu.code"
#undef  DX
#undef  TX
#undef  PX
#undef  EX
            end do
         end do
       end function QSATLQU2

       function QSATLQU3(TL,PL,DQ) result(QS)
         real,              intent(IN) :: TL(:,:,:)
         real, optional,    intent(IN) :: PL(:,:,:)
         real, optional,    intent(OUT):: DQ(:,:,:)
         real    :: QS(SIZE(TL,1),SIZE(TL,2),SIZE(TL,3))
         integer :: I, J, K
         real    :: TI
         real    :: TT
         real    :: DDQ
         real    :: DD
         integer :: IT
         do K=1,size(TL,3)
            do J=1,size(TL,2)
               do I=1,size(TL,1)
#define TX TL(I,J,K)
#define PX PL(I,J,K)
#define EX QS(I,J,K)
#define DX DQ(I,J,K)
#include "qsatlqu.code"
#undef  DX
#undef  TX
#undef  PX
#undef  EX

               end do
            end do
         end do
       end function QSATLQU3



       function QSATICE0(TL,PL,DQ) result(QS)
!$acc routine seq
         real,              intent(IN) :: TL
         real, optional,    intent(IN) :: PL
         real, optional,    intent(OUT):: DQ
         real    :: QS

         real    :: TI,W
         real    :: DD
         real    :: TT
         real    :: DDQ
         integer :: IT
#define TX TL
#define PX PL
#define EX QS
#define DX DQ
#include "qsatice.code"
#undef  DX
#undef  TX
#undef  EX
#undef  PX
         return
       end function QSATICE0

       function QSATICE1(TL,PL,DQ) result(QS)
         real,              intent(IN) :: TL(:)
         real, optional,    intent(IN) :: PL(:)
         real, optional,    intent(OUT):: DQ(:)
         real    :: QS(SIZE(TL,1))
         integer :: I
         real    :: TI,W  
         real    :: TT
         real    :: DDQ
         real    :: DD
         integer :: IT
         do I=1,size(TL,1)
#define TX TL(I)
#define PX PL(I)
#define EX QS(I)
#define DX DQ(I)
#include "qsatice.code"
#undef  DX
#undef  TX
#undef  PX
#undef  EX
         end do
       end function QSATICE1

       function QSATICE2(TL,PL,DQ) result(QS)
         real,              intent(IN) :: TL(:,:)
         real, optional,    intent(IN) :: PL(:,:)
         real, optional,    intent(OUT):: DQ(:,:)
         real    :: QS(SIZE(TL,1),SIZE(TL,2))
         integer :: I, J
         real    :: TI,W  
         real    :: TT
         real    :: DDQ
         real    :: DD
         integer :: IT
         do J=1,size(TL,2)
            do I=1,size(TL,1)
#define TX TL(I,J)
#define PX PL(I,J)
#define EX QS(I,J)
#define DX DQ(I,J)
#include "qsatice.code"
#undef  DX
#undef  TX
#undef  PX
#undef  EX
            end do
         end do
       end function QSATICE2

       function QSATICE3(TL,PL,DQ) result(QS)
         real,              intent(IN) :: TL(:,:,:)
         real, optional,    intent(IN) :: PL(:,:,:)
         real, optional,    intent(OUT):: DQ(:,:,:)
         real    :: QS(SIZE(TL,1),SIZE(TL,2),SIZE(TL,3))
         integer :: I, J, K
         real    :: TI,W  
         real    :: TT
         real    :: DDQ
         real    :: DD
         integer :: IT
         do K=1,size(TL,3)
            do J=1,size(TL,2)
               do I=1,size(TL,1)
#define TX TL(I,J,K)
#define PX PL(I,J,K)
#define EX QS(I,J,K)
#define DX DQ(I,J,K)
#include "qsatice.code"
#undef  DX
#undef  TX
#undef  PX
#undef  EX

               end do
            end do
         end do
       end function QSATICE3



!==============================================
!==============================================

!  Traditional Qsat and Dqsat (these are deprecated)

!==============================================
!==============================================

!BOPI

! !IROUTINE: GEOS_Qsat -- Computes satuation specific humidity.

! !INTERFACE:

!    function GEOS_Qsat(TL,PL,RAMP,PASCALS,DQSAT) result(QSAT)
!
! Overloads:
!
!      real,                      intent(IN)                        :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real,                     optional, intent(OUT)              :: DQSAT
!      real                                                         :: QSAT
!
!      real, dimension(:),        intent(IN)                        :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real,                     optional, intent(OUT)              :: DQSAT(:)
!      real, dimension(size(PL,1))                                  :: QSAT
!
!      real, dimension(:,:),      intent(IN)                        :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real,                     optional, intent(OUT)              :: DQSAT(:,:)
!      real, dimension(size(PL,1),size(PL,2))                       :: QSAT
!
!      real, dimension(:,:,:),    intent(IN)                        :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real,                     optional, intent(OUT)              :: DQSAT(:,:,:)
!      real, dimension(size(PL,1),size(PL,2),size(PL,3))            :: QSAT
!

! !DESCRIPTION:  Uses various formulations of the saturation
!                vapor pressure to compute the saturation specific 
!    humidity for temperature TL and pressure PL.
!
!    For temperatures <= TMIX (-20C)
!    the calculation is done over ice; for temperatures >= ZEROC (0C) the calculation
!    is done over liquid water; and in between these values,
!    it interpolates linearly between the two.
!
!    The optional RAMP is the width of this
!    ice/water ramp (i.e., TMIX = ZEROC-RAMP); its default is 20.
!
!    If PASCALS is true, PL is
!    assumed to be in Pa; if false or not present, it is assumed to be in mb.
!
!    The choice of saturation vapor pressure formulation is a compile-time
!    option. Three choices are currently supported: The CAM formulation,
!    Murphy and Koop (2005, QJRMS), and Staars formulation from NSIPP-1.
!
!    Another compile time choice is whether to use the exact formulation
!    or a table look-up.
!    If UTBL is true, tabled values of the saturation vapor pressures
!    are used. These tables are automatically generated at a 0.1K resolution
!    for whatever vapor pressure formulation is being used.
! 
!    
!EOPI
    
       
  function QSAT0(TL,PL,RAMP,PASCALS,DQSAT) result(QSAT)
!$acc routine seq
    real,   intent(IN) :: TL, PL
    logical, optional, intent(IN) :: PASCALS
    real,    optional, intent(IN) :: RAMP
    real,    optional, intent(OUT):: DQSAT
    real    :: QSAT

    real    :: URAMP, DD, QQ, TI, DQ, PP
    integer :: IT

   !  if (debugIsEnabled) then
   !    if (ieee_is_nan(TL)) call lgr%warning(' QSAT0: TL contains NaN')
   !    if (ieee_is_nan(PL)) call lgr%warning(' QSAT0: PL contains NaN')
   !  end if

    if(present(RAMP)) then
       URAMP = -abs(RAMP)
    else
       URAMP = TMIX
    end if

    if(present(PASCALS)) then
       if(PASCALS) then
          PP = PL
       else
          PP = PL*100.
       end if
    else
       PP = PL*100.
    end if

    if((URAMP==TMIX .OR. URAMP==0.) .and. UTBL) then
       if(FIRST) then
          FIRST = .false.
         !  call ESINIT
         !  call LOGGER_INIT
       end if

       if    (TL<=TMINTBL) then
          TI = TMINTBL
       elseif(TL>=TMAXTBL-.001) then
          TI = TMAXTBL-.001
       else
          TI = TL
       end if

       TI = (TI - TMINTBL)*DEGSUBS+1
       IT = int(TI)

       if(URAMP==TMIX) then
          DQ    = ESTBLX(IT+1) - ESTBLX(IT)
          QSAT  = (TI-IT)*DQ + ESTBLX(IT)
       else
          DQ    = ESTBLE(IT+1) - ESTBLE(IT)
          QSAT  = (TI-IT)*DQ + ESTBLE(IT)
       endif

       if(present(DQSAT)) DQSAT = DQ*DEGSUBS

       if(PP <= QSAT) then
          QSAT = MAX_MIXING_RATIO
          if(present(DQSAT)) DQSAT = 0.0
       else
          DD = 1.0/(PP - (1.0-ESFAC)*QSAT)
          QSAT = ESFAC*QSAT*DD
          if(present(DQSAT)) DQSAT = ESFAC*DQSAT*PP*(DD*DD)
       end if

    else

       if(FIRST) then
          FIRST = .false.
         !  call LOGGER_INIT
       end if

       TI = TL - ZEROC

       if    (TI <= URAMP) then
          QSAT  =  QSATICE0(TL,PP,DQ=DQSAT)
       elseif(TI >= 0.0  ) then
          QSAT  =  QSATLQU0(TL,PP,DQ=DQSAT)
       else
          QSAT  =  QSATICE0(TL,PP,DQ=DQSAT)
          QQ    =  QSATLQU0(TL,PP,DQ=DQ   )
          TI    =  TI/URAMP
          QSAT  =  TI*(QSAT - QQ) +  QQ
          if(PRESENT(DQSAT)) DQSAT = TI*(DQSAT-DQ) + DQ
       end if

    end if

  end function QSAT0

    function QSAT1(TL,PL,RAMP,PASCALS,DQSAT) result(QSAT)
      real,              intent(IN) :: TL(:), PL(:)
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: DQSAT(:)
      real :: QSAT(size(TL,1))
      integer :: I

      ! if (debugIsEnabled) then
      !    if (any(ieee_is_nan(TL))) call lgr%warning(' QSAT1: TL contains NaN')
      !    if (any(ieee_is_nan(PL))) call lgr%warning(' QSAT1: PL contains NaN')
      ! end if

      do I=1,SIZE(TL,1)
         if (present(DQSAT)) then
            QSAT(I) = QSAT0(TL(I),PL(I),RAMP,PASCALS,DQSAT(I))
         else
            QSAT(I) = QSAT0(TL(I),PL(I),RAMP,PASCALS)
         end if
      end do
    end function QSAT1

    function QSAT2(TL,PL,RAMP,PASCALS,DQSAT) result(QSAT)
      real,              intent(IN) :: TL(:,:), PL(:,:)
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: DQSAT(:,:)
      real :: QSAT(size(TL,1),size(TL,2))
      integer :: I, J

      ! if (debugIsEnabled) then
      !    if (any(ieee_is_nan(TL))) call lgr%warning(' QSAT2: TL contains NaN')
      !    if (any(ieee_is_nan(PL))) call lgr%warning(' QSAT2: PL contains NaN')
      ! end if

      do J=1,SIZE(TL,2)
         do I=1,SIZE(TL,1)
            if (present(DQSAT)) then
               QSAT(I,J) = QSAT0(TL(I,J),PL(I,J),RAMP,PASCALS,DQSAT(I,J))
            else
               QSAT(I,J) = QSAT0(TL(I,J),PL(I,J),RAMP,PASCALS)
            end if
         end do
      end do
    end function QSAT2

    function QSAT3(TL,PL,RAMP,PASCALS,DQSAT) result(QSAT)
      real,              intent(IN) :: TL(:,:,:), PL(:,:,:)
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: DQSAT(:,:,:)
      real :: QSAT(size(TL,1),size(TL,2),size(TL,3))
      integer :: I, J, K

      ! if (debugIsEnabled) then
      !    if (any(ieee_is_nan(TL))) call lgr%warning(' QSAT3: TL contains NaN')
      !    if (any(ieee_is_nan(PL))) call lgr%warning(' QSAT3: PL contains NaN')
      ! end if

      do K=1,SIZE(TL,3)
         do J=1,SIZE(TL,2)
            do I=1,SIZE(TL,1)
               if (present(DQSAT)) then
                  QSAT(I,J,K) = QSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS,DQSAT(I,J,K))
               else
                  QSAT(I,J,K) = QSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS)
               end if
            end do
         end do
      end do
    end function QSAT3

!=======================================================================================

!BOPI

! !IROUTINE: GEOS_DQsat -- Computes derivative satuation specific humidity wrt temperature.

! !INTERFACE:

!    function GEOS_DQsat(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
!
! Overloads:
!
!      real,                               intent(IN)               :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real,                     optional, intent(OUT)              :: QSAT
!      real                                                         :: DQSAT
!
!      real, dimension(:),                 intent(IN)               :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real, dimension(:),       optional, intent(OUT)              :: QSAT
!      real, dimension(size(PL,1))                                  :: DQSAT
!
!      real, dimension(:,:),               intent(IN)               :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real, dimension(:,:),     optional, intent(OUT)              :: QSAT
!      real, dimension(size(PL,1),size(PL,2))                       :: DQSAT
!
!      real, dimension(:,:,:),             intent(IN)               :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real, dimension(:,:,:),   optional, intent(OUT)              :: QSAT
!      real, dimension(size(PL,1),size(PL,2),size(PL,3))            :: DQSAT
!
!      real, dimension(:,:,:,:),           intent(IN)               :: TL, PL
!      logical,                  optional, intent(IN)               :: PASCALS
!      real,                     optional, intent(IN)               :: RAMP
!      real, dimension(:,:,:,:), optional, intent(OUT)              :: QSAT
!      real, dimension(size(PL,1),size(PL,2),size(PL,3),size(PL,4)) :: DQSAT

! !DESCRIPTION:  Differentiates the approximations used
!                by GEOS_Qsat with respect to temperature,
!    using the same scheme to handle ice. Arguments are as in 
!    GEOS_Qsat, with the addition of QSAT, which is the saturation specific
!    humidity. This is for economy, in case both qsat and dqsat are 
!    required.
!                

!EOPI

    
    function DQSAT0(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
      real,   intent(IN) :: TL, PL
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: QSAT
      real    :: DQSAT
      real    :: URAMP, TT, DD, DQQ, QQ, TI, DQI, QI, PP
      integer :: IT

      ! if (debugIsEnabled) then
      !    if (ieee_is_nan(TL)) call lgr%warning('DQSAT0: TL contains NaN')
      !    if (ieee_is_nan(PL)) call lgr%warning('DQSAT0: PL contains NaN')
      ! end if

      if(present(RAMP)) then
         URAMP = -abs(RAMP)
      else
         URAMP = TMIX
      end if

      if(present(PASCALS)) then
         if(PASCALS) then
            PP = PL
         else
            PP = PL*100.
         end if
      else
         PP = PL*100.
      end if

      if((URAMP==TMIX .OR. URAMP==0.) .and. UTBL) then

         if(FIRST) then
            FIRST = .false.
            call ESINIT
            ! call LOGGER_INIT
         end if

         if    (TL<=TMINTBL) then
            TI = TMINTBL
         elseif(TL>=TMAXTBL-.001) then
            TI = TMAXTBL-.001
         else
            TI = TL
         end if

         TT = (TI - TMINTBL)*DEGSUBS+1
         IT = int(TT)

         if(URAMP==TMIX) then
            DQQ =  ESTBLX(IT+1) - ESTBLX(IT)
            QQ  =  (TT-IT)*DQQ + ESTBLX(IT)
         else
            DQQ =  ESTBLE(IT+1) - ESTBLE(IT)
            QQ  =  (TT-IT)*DQQ + ESTBLE(IT)
         endif

         if(PP <= QQ) then
            if(present(QSAT)) QSAT = MAX_MIXING_RATIO
            DQSAT = 0.0
         else
            DD = 1.0/(PP - (1.0-ESFAC)*QQ)
            if(present(QSAT)) QSAT = ESFAC*QQ*DD
            DQSAT = (ESFAC*DEGSUBS)*DQQ*PP*(DD*DD)
         end if

      else

         if(FIRST) then
            FIRST = .false.
            ! call LOGGER_INIT
         end if

         TI = TL - ZEROC

         if    (TI <= URAMP) then
            QQ  = QSATICE0(TL,PP,DQ=DQSAT)
            if(present(QSAT)) QSAT  = QQ
         elseif(TI >= 0.0  ) then
            QQ  = QSATLQU0(TL,PP,DQ=DQSAT)
            if(present(QSAT)) QSAT  = QQ
         else
            QQ  = QSATLQU0(TL,PP,DQ=DQQ)
            QI  = QSATICE0(TL,PP,DQ=DQI)
            TI  = TI/URAMP
            DQSAT = TI*(DQI - DQQ) + DQQ
            if(present(QSAT)) QSAT  = TI*(QI - QQ) +  QQ
         end if

      end if

    end function DQSAT0
        
    function DQSAT1(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
      real,              intent(IN) :: TL(:), PL(:)
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: QSAT(:)
      real :: DQSAT(size(TL,1))
      integer :: I

      ! if (debugIsEnabled) then
      !    if (any(ieee_is_nan(TL))) call lgr%warning('DQSAT1: TL contains NaN')
      !    if (any(ieee_is_nan(PL))) call lgr%warning('DQSAT1: PL contains NaN')
      ! end if

      do I=1,SIZE(TL,1)
         if (present(QSAT)) then
            DQSAT(I) = DQSAT0(TL(I),PL(I),RAMP,PASCALS,QSAT(I))
         else
            DQSAT(I) = DQSAT0(TL(I),PL(I),RAMP,PASCALS)
         endif
      end do
    end function DQSAT1

    function DQSAT2(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
      real,              intent(IN) :: TL(:,:), PL(:,:)
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: QSAT(:,:)
      real :: DQSAT(size(TL,1),size(TL,2))
      integer :: I, J

      ! if (debugIsEnabled) then
      !    if (any(ieee_is_nan(TL))) call lgr%warning('DQSAT2: TL contains NaN')
      !    if (any(ieee_is_nan(PL))) call lgr%warning('DQSAT2: PL contains NaN')
      ! end if

      do J=1,SIZE(TL,2)
         do I=1,SIZE(TL,1)
            if (present(QSAT)) then
               DQSAT(I,J) = DQSAT0(TL(I,J),PL(I,J),RAMP,PASCALS,QSAT(I,J))
            else
               DQSAT(I,J) = DQSAT0(TL(I,J),PL(I,J),RAMP,PASCALS)
            end if
         end do
      end do
    end function DQSAT2

    function DQSAT3(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
      real,              intent(IN) :: TL(:,:,:), PL(:,:,:)
      logical, optional, intent(IN) :: PASCALS
      real,    optional, intent(IN) :: RAMP
      real,    optional, intent(OUT):: QSAT(:,:,:)
      real :: DQSAT(size(TL,1),size(TL,2),size(TL,3))
      integer :: I, J, K

      ! if (debugIsEnabled) then
      !    if (any(ieee_is_nan(TL))) call lgr%warning('DQSAT3: TL contains NaN')
      !    if (any(ieee_is_nan(PL))) call lgr%warning('DQSAT3: PL contains NaN')
      ! end if

      do K=1,SIZE(TL,3)
         do J=1,SIZE(TL,2)
            do I=1,SIZE(TL,1)
               if (present(QSAT)) then
                  DQSAT(I,J,K) = DQSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS,QSAT(I,J,K))
               else
                  DQSAT(I,J,K) = DQSAT0(TL(I,J,K),PL(I,J,K),RAMP,PASCALS)
               end if
            end do
         end do
      end do
    end function DQSAT3

!==============================================

!BOPI

! !IROUTINE: GEOS_QsatSet -- Sets behavior of GEOS_QsatLqu an GEOS_QsatIce

! !INTERFACE:

       subroutine GEOS_QsatSet(USETABLE,FORMULATION)
         logical, optional, intent(IN) :: USETABLE
         integer, optional, intent(IN) :: FORMULATION

! !DESCRIPTION: GEOS_QsatSet can be used to modify 
!  the behavior of GEOS_QsatLqu an GEOS_QsatIce 
!  from its default setting.

!  If {\tt \bf USETABLE} is true, tabled values of the saturation vapor pressures are used.
!  These tables are automatically generated at a 0.1K resolution for whatever
!  vapor pressure formulation is being used. The default is to use the table.

!  {\tt \bf FORMULATION} sets the saturation vapor pressure function.
!  Three formulations of saturation vapor pressure are supported: 
!  the Starr code that was in NSIPP-1 (FORMULATION==1), the formulation in  CAM 
!  (FORMULATION==2), and Murphy and Koop (2005, QJRMS) (FORMULATION==3).
!  The default is FORMULATION=1.

!  If appropriate, GEOS_QsatSet also initializes the tables. If GEOS_QsatSet is
!  not called and tables are required, they will be initialized the first time
!  a Qsat function is called.

!EOPI
       
         if(present(UseTable   )) UTBL = UseTable
         if(present(Formulation)) TYPE = max(min(Formulation,3),1)
         
         if(TYPE==3)  then ! Murphy and Koop (2005, QJRMS)
            TMINICE    =  max(TMINTBL,110.)
            TMINLQU    =  max(TMINTBL,123.)
         else
            TMINLQU    =  ZEROC - 40.0
            TMINICE    =  ZEROC + TMINSTR
         endif

         if(UTBL) then
            call ESINIT
            ! call LOGGER_INIT
         end if

         return
       end subroutine GEOS_QsatSet

! !=======================================================================================

        subroutine ESINIT
!$acc routine seq
! Saturation vapor pressure table initialization. This is invoked if UTBL is true 
! on the first call to any qsat routine or whenever GEOS_QsatSet is called 
! N.B.--Tables are in Pa
 
          integer :: I
          real    :: T
          logical :: UT

          UT = UTBL
          UTBL=.false.

          do I=1,TABLESIZE

            T = (I-1)*DELTA_T + TMINTBL

            ESTBLW(I) = QSATLQU0(T)

            if(T>ZEROC) then
               ESTBLE(I) = ESTBLW(I)
            else
               ESTBLE(I) = QSATICE0(T)
            end if

            T = T-ZEROC

            if(T>=TMIX .and. T<0.0) then
               ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
            else
               ESTBLX(I) = ESTBLE(I)
            end if

         end do

         ESTFRZ = QSATLQU0(ZEROC  )
         ESTLQU = QSATLQU0(TMINLQU)

         UTBL = UT
!!$acc update device(ESTFRZ, ESTLQU, UTBL, ESTBLX, ESTBLE, ESTBLW)
       end subroutine ESINIT

      subroutine ESINIT_
         ! Saturation vapor pressure table initialization. This is invoked if UTBL is true 
         ! on the first call to any qsat routine or whenever GEOS_QsatSet is called 
         ! N.B.--Tables are in Pa
          
         integer :: I
         real    :: T
         logical :: UT

         UT = UTBL
         UTBL=.false.

         do I=1,TABLESIZE

            T = (I-1)*DELTA_T + TMINTBL

            ESTBLW(I) = QSATLQU0(T)

            if(T>ZEROC) then
               ESTBLE(I) = ESTBLW(I)
            else
               ESTBLE(I) = QSATICE0(T)
            end if

            T = T-ZEROC

            if(T>=TMIX .and. T<0.0) then
               ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
            else
               ESTBLX(I) = ESTBLE(I)
            end if

         end do

         ESTFRZ = QSATLQU0(ZEROC  )
         ESTLQU = QSATLQU0(TMINLQU)

         UTBL = UT

      end subroutine ESINIT_

      subroutine ESINIT_v2
         if(FIRST) then
            FIRST = .false.
            call ESINIT_
         endif
      !$acc update device(FIRST)
      ! print*, 'ESTFRZ = ', ESTFRZ
      ! print*, 'ESTLQU = ', ESTLQU
      ! print*, "UBTL = ", UTBL
      ! print*, "sum(ESTBLX) = ", sum(ESTBLX)
      ! print*, "sum(ESTBLE) = ", sum(ESTBLE)
      ! print*, "sum(ESTBLW) = ", sum(ESTBLW)
      !$acc update device(FIRST, ESTFRZ, ESTLQU, UTBL, ESTBLX, ESTBLE, ESTBLW)
      ! !$acc update host(ESTFRZ, ESTLQU, UTBL, ESTBLX, ESTBLE, ESTBLW)
      ! print*, 'ESTFRZ = ', ESTFRZ
      ! print*, 'ESTLQU = ', ESTLQU
      ! print*, "UBTL = ", UTBL
      ! print*, "sum(ESTBLX) = ", sum(ESTBLX)
      ! print*, "sum(ESTBLE) = ", sum(ESTBLE)
      ! print*, "sum(ESTBLW) = ", sum(ESTBLW)
      !call exit(1)
      end subroutine
!        subroutine LOGGER_INIT

!           implicit none

!           lgr => logging%get_logger('SHARED.GMAOSHARED.GEOSSHARED.QSAT')
!           debugIsEnabled = lgr%isEnabledFor(DEBUG)

!        end subroutine LOGGER_INIT











!*************************************************************************
!*************************************************************************

!  Tridiagonal solvers

!*************************************************************************
!*************************************************************************

!BOP

! !IROUTINE:  VTRISOLVE -- Solves for tridiagonal system that has been decomposed by VTRILU


! !INTERFACE:

!  subroutine GEOS_TRISOLVE ( A,B,C,Y,YG )

! !ARGUMENTS:

!    real, dimension([:,[:,]] :),  intent(IN   ) ::  A, B, C
!    real, dimension([:,[:,]] :),  intent(INOUT) ::  Y
!    real, dimension([:,[:,]] :),  intent(IN   ) ::  YG

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



!BOP

! !IROUTINE:  VTRILU --  Does LU decomposition of tridiagonal matrix.

! !INTERFACE:

!  subroutine GEOS_TRILU  ( A,B,C )

! !ARGUMENTS:

!    real, dimension ([:,[:,]] :), intent(IN   ) ::  C
!    real, dimension ([:,[:,]] :), intent(INOUT) ::  A, B

! !DESCRIPTION: {\tt VTRILU} performs an $LU$ decomposition on
! a tridiagonal matrix $M=LU$.

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

! On input, A, B, and C contain, $a_k$, $b_k$, and $c_k$
! the lower, main, and upper diagonals of the matrix, respectively.
! On output, B contains $1/\hat{b}_k$, the inverse of the main diagonal of $U$,
! and A contains $\hat{a}_k$,
! the lower diagonal of $L$. C contains the upper diagonal of the original matrix and of $U$.

! The new diagonals $\hat{a}_k$ and $\hat{b}_k$ are:
! $$
! \begin{array}{rcl}
! \hat{b}_1 & = & b_1, \\
! \hat{a}_k & = & \makebox[2 in][l]{$a_k / \hat{b}_{k-1}$,}  k=2, K, \\
! \hat{b}_k & = & \makebox[2 in][l]{$b_k - c_{k-1} \hat{a}_k$,} k=2, K. 
! \end{array}
! $$
! EOP



! #define DIMS
!   subroutine GEOS_TRILU1 ( A,B,C )
! #include "trilu.code"
!   end subroutine GEOS_TRILU1

!   subroutine GEOS_TRISOLVE1 ( A,B,C,Y,YG )
! #include "trisolve.code"
!   end subroutine GEOS_TRISOLVE1
! #undef DIMS

! #define DIMS :,
!   subroutine GEOS_TRILU2 ( A,B,C )
! #include "trilu.code"
!   end subroutine GEOS_TRILU2

!   subroutine GEOS_TRISOLVE2 ( A,B,C,Y,YG )
! #include "trisolve.code"
!   end subroutine GEOS_TRISOLVE2
! #undef DIMS

! #define DIMS :,:,
!   subroutine GEOS_TRILU3 ( A,B,C )
! #include "trilu.code"
!   end subroutine GEOS_TRILU3

!   subroutine GEOS_TRISOLVE3 ( A,B,C,Y,YG )
! #include "trisolve.code"
!   end subroutine GEOS_TRISOLVE3
! #undef DIMS



  end module GEOS_UtilsMod

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”

! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at

! http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.