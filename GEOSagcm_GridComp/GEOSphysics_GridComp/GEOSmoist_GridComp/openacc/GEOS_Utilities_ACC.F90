module GEOS_Utils_ACC

  use MAPL_ConstantsMod
  use ieee_arithmetic

  implicit none
  private

  public GEOS_DQSAT_NOACC

  interface GEOS_DQSAT_NOACC
     module procedure DQSAT0_NOACC
     module procedure DQSAT3_NOACC
  end interface GEOS_DQSAT_NOACC

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

  logical            :: UTBL       = .true.
  integer            :: TYPE       =  1

  logical            :: FIRST      = .true.

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
  real*8, parameter  :: &
       CL(0:9)=(/54.842763, -6763.22, -4.21000, .000367, &
       .0415, 218.8,  53.878000, -1331.22, -9.44523, .014025  /)

  real       :: TMINLQU    =  ZEROC - 40.0
  real       :: TMINICE    =  ZEROC + TMINSTR

contains

  function QSATLQU0_NOACC(TL,PL,DQ) result(QS)

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
#include "qsatlqu_noacc.code"
#undef  DX
#undef  TX
#undef  EX
#undef  PX
    return

  end function QSATLQU0_NOACC

  function QSATICE0_NOACC(TL,PL,DQ) result(QS)

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
#include "qsatice_noacc.code"
#undef  DX
#undef  TX
#undef  EX
#undef  PX
    return

  end function QSATICE0_NOACC

  function DQSAT0_NOACC(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)

    real,   intent(IN) :: TL, PL
    logical, optional, intent(IN) :: PASCALS
    real,    optional, intent(IN) :: RAMP
    real,    optional, intent(OUT):: QSAT
    real    :: DQSAT
    real    :: URAMP, TT, DD, DQQ, QQ, TI, DQI, QI, PP
    integer :: IT

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
          call LOGGER_INIT
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
          call LOGGER_INIT
       end if

       TI = TL - ZEROC

       if    (TI <= URAMP) then
          QQ  = QSATICE0_NOACC(TL,PP,DQ=DQSAT)
          if(present(QSAT)) QSAT  = QQ
       elseif(TI >= 0.0  ) then
          QQ  = QSATLQU0_NOACC(TL,PP,DQ=DQSAT)
          if(present(QSAT)) QSAT  = QQ
       else
          QQ  = QSATLQU0_NOACC(TL,PP,DQ=DQQ)
          QI  = QSATICE0_NOACC(TL,PP,DQ=DQI)
          TI  = TI/URAMP
          DQSAT = TI*(DQI - DQQ) + DQQ
          if(present(QSAT)) QSAT  = TI*(QI - QQ) +  QQ
       end if

    end if

  end function DQSAT0_NOACC

  function DQSAT3_NOACC(TL,PL,RAMP,PASCALS,QSAT) result(DQSAT)
    real,              intent(IN) :: TL(:,:,:), PL(:,:,:)
    logical, optional, intent(IN) :: PASCALS
    real,    optional, intent(IN) :: RAMP
    real,    optional, intent(OUT):: QSAT(:,:,:)
    real :: DQSAT(size(TL,1),size(TL,2),size(TL,3))
    integer :: I, J, K

    do K=1,SIZE(TL,3)
       do J=1,SIZE(TL,2)
          do I=1,SIZE(TL,1)
             if (present(QSAT)) then
                DQSAT(I,J,K) = DQSAT0_NOACC(TL(I,J,K),PL(I,J,K),RAMP,PASCALS,QSAT(I,J,K))
             else
                DQSAT(I,J,K) = DQSAT0_NOACC(TL(I,J,K),PL(I,J,K),RAMP,PASCALS)
             end if
          end do
       end do
    end do
  end function DQSAT3_NOACC

  subroutine ESINIT

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

       ESTBLW(I) = QSATLQU0_NOACC(T)

       if(T>ZEROC) then
          ESTBLE(I) = ESTBLW(I)
       else
          ESTBLE(I) = QSATICE0_NOACC(T)
       end if

       T = T-ZEROC

       if(T>=TMIX .and. T<0.0) then
          ESTBLX(I) = ( T/TMIX )*( ESTBLE(I) - ESTBLW(I) ) + ESTBLW(I)
       else
          ESTBLX(I) = ESTBLE(I)
       end if

    end do

    ESTFRZ = QSATLQU0_NOACC(ZEROC  )
    ESTLQU = QSATLQU0_NOACC(TMINLQU)

    UTBL = UT

  end subroutine ESINIT

end module GEOS_Utils_ACC

