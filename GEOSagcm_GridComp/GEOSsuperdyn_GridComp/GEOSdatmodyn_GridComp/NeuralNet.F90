! $Id$

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_ErrLog.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Singcol -- A Module to drive single column model with  
!  profile data.

! !INTERFACE:

module NeuralNetMod

! !USES:

  use ESMF
  use GEOS_Mod

  implicit none
  private


  PUBLIC NN_SCMI
  PUBLIC NN_INTRFC
  PUBLIC NN_INIT
  PUBLIC T_NN_XRNG
  PUBLIC T_NN_WGT

  PUBLIC QTESTDATA
  PUBLIC TTESTDATA

  integer, parameter :: NPREDS = 5

  TYPE T_NN_XRNG
     REAL*8, pointer,   DIMENSION(:,:)  :: UDQDX
     REAL*8, pointer,   DIMENSION(:,:)  :: UDTDX
     REAL*8, pointer,   DIMENSION(:,:)  :: VDQDY
     REAL*8, pointer,   DIMENSION(:,:)  :: VDTDY 
     REAL*8, pointer,   DIMENSION(:,:)  :: OMEGA
  end TYPE T_NN_XRNG

  TYPE T_NN_WGT
     REAL*8, pointer,   DIMENSION(:)    :: UDQDX
     REAL*8, pointer,   DIMENSION(:)    :: UDTDX
     REAL*8, pointer,   DIMENSION(:)    :: VDQDY
     REAL*8, pointer,   DIMENSION(:)    :: VDTDY 
     REAL*8, pointer,   DIMENSION(:)    :: OMEGA
  end TYPE T_NN_WGT

  TYPE T_NN_RW
     CHARACTER(LEN=64)              :: NAME
     REAL*8, pointer, DIMENSION(:,:)  :: XRNG
     REAL*8, pointer, DIMENSION(:,:)  :: WGT
  end TYPE T_NN_RW


  real, dimension(72), parameter :: TTESTDATA = (/     &
 294.87,  293.44,  292.05,  290.77,  289.61,  288.55,  &
 288.66,  288.46,  287.67,  286.66,  285.76,  284.59,  &
 283.64,  283.18,  282.28,  280.77,  279.33,  277.50,  &
 276.06,  274.21,  272.62,  270.47,  266.86,  262.80,  &
 258.72,  254.30,  249.73,  242.97,  235.04,  226.41,  &
 217.72,  208.57,  200.88,  195.36,  192.56,  193.07,  &
 196.76,  205.30,  202.73,  211.92,  210.59,  217.66,  &
 215.98,  220.69,  220.97,  224.92,  224.84,  229.17,  &
 230.10,  236.04,  237.89,  242.14,  243.45,  248.42,  &
 252.54,  258.70,  262.55,  265.93,  265.38,  263.04,  &
 260.77,  256.62,  253.78,  248.20,  243.48,  237.78,  &
 231.85,  226.67,  220.86,  217.02,  209.41,  210.31   /)

  real, dimension(72), parameter :: QTESTDATA = (/                                                           &
     0.012145169,      0.011775635,      0.011550725,      0.011353433,      0.011130421,      0.010783325,  &
     0.010073253,     0.0097698579,     0.0095104435,     0.0093158633,     0.0089346385,     0.0084399348,  &
    0.0073577804,     0.0049268133,     0.0033021194,     0.0023144188,     0.0017462018,     0.0013301875,  &
   0.00091290235,    0.00049789826,    0.00041268408,    0.00035149863,    0.00033418552,    0.00023334385,  &
   0.00015293497,    0.00016023224,    0.00016299008,    9.0956491e-05,    3.0488196e-05,    1.4280661e-05,  &
   1.1828073e-05,    9.7573311e-06,    2.7460148e-06,    1.5195064e-06,    1.2196456e-06,    1.5783821e-06,  &
   1.6860496e-06,    1.9997979e-06,    2.0628661e-06,    2.0594791e-06,    1.8036558e-06,    1.6057670e-06,  &
   1.5282549e-06,    1.5742063e-06,    1.6001330e-06,    1.6366392e-06,    1.6382760e-06,    1.6478992e-06,  &
   1.6628369e-06,    1.6739592e-06,    1.6752618e-06,    1.6873224e-06,    1.6996736e-06,    1.7280622e-06,  &
   1.7446847e-06,    1.7510930e-06,    1.7531524e-06,    1.7431888e-06,    1.7264066e-06,    1.7308356e-06,  &
   1.7606135e-06,    1.7770291e-06,    1.7914211e-06,    1.7947577e-06,    1.8051473e-06,    1.8164109e-06,  &
   1.8268433e-06,    1.8353438e-06,    1.8427180e-06,    1.8470054e-06,    1.8472193e-06,    1.8468419e-06   /)





! integer, parameter :: LM = 72      ! model level Bottom -> Top
! integer, parameter :: nx = 40      ! layer number to calculate
! integer, parameter :: NHID = 12    ! hiden nodes
! integer, parameter :: INPUT= 80    ! input number
! integer, parameter :: NOUT = 40    ! output number
! integer, parameter :: NUNIT= 20    ! weight data unit to link

contains

subroutine NN_INIT( NHID, INPUT, NOUT, wgt, xrng )

 integer, intent(in)              :: NHID, NOUT, INPUT
 type( T_NN_WGT )        :: wgt
 type( T_NN_XRNG )       :: xrng


 character(len=120), dimension(NPREDS)   :: nfile

 integer  ::   I_UDTDX, I_UDQDX, I_VDQDY, I_VDTDY, I_OMEGA
 integer  ::   NUNIT, NDU,NWT,INPUTom

 I_UDTDX = 1
 I_UDQDX = 2
 I_VDTDY = 3
 I_VDQDY = 4
 I_OMEGA = 5


 nfile( I_UDTDX ) = 'NN_UDTDX.dat'
 nfile( I_UDQDX ) = 'NN_UDQDX.dat'
 nfile( I_VDTDY ) = 'NN_VDTDY.dat'
 nfile( I_VDQDY ) = 'NN_VDQDY.dat'
 nfile( I_OMEGA ) = 'NN_OMEGA.dat'


 NDU = INPUT + NOUT
 NWT = NHID*(INPUT+1) + (NHID + 1)*NOUT

  allocate(  WGT%UDQDX( NWT )   )
  allocate(  WGT%UDTDX( NWT )   )
  allocate(  WGT%VDQDY( NWT )   )
  allocate(  WGT%VDTDY( NWT )   )

  allocate(  XRNG%UDQDX( NDU, 3 )   )
  allocate(  XRNG%UDTDX( NDU, 3 )   )
  allocate(  XRNG%VDQDY( NDU, 3 )   )
  allocate(  XRNG%VDTDY( NDU, 3 )   )

!! BDC added PS and TS to inputs for OMEGA
 INPUTom = INPUT+2
 NDU     = INPUTom + NOUT
 NWT     = NHID*(INPUTom+1) + (NHID + 1)*NOUT
  allocate(  WGT%OMEGA( NWT )   )
  allocate(  XRNG%OMEGA( NDU, 3 )   )

 nunit=20


 open(nunit,file=trim(nfile(I_UDQDX) ),status = 'old')
 call read_weight(  WGT%UDQDX ,  XRNG%UDQDX , nunit)
 close(nunit)

 open(nunit,file=trim(nfile(I_UDTDX) ),status = 'old')
 call read_weight(  WGT%UDTDX ,  XRNG%UDTDX , nunit)
 close(nunit)

 open(nunit,file=trim(nfile(I_VDQDY) ),status = 'old')
 call read_weight(  WGT%VDQDY ,  XRNG%VDQDY , nunit)
 close(nunit)

 open(nunit,file=trim(nfile(I_VDTDY) ),status = 'old')
 call read_weight(  WGT%VDTDY ,  XRNG%VDTDY , nunit)
 close(nunit)

 

end subroutine NN_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine NN_INTRFC(  U, V, T, Q, P, TS, PS, &
                       XRNG,  WGT,    &
                       NHID, NOUT, INPUT, &
                       UdQdx, VdQdy,  &
                       UdTdx, VdTdy,  &
                       omega        )  

integer, intent(in)              :: NHID, NOUT, INPUT

type( T_NN_WGT ),  intent(in)   :: wgt
type( T_NN_XRNG ), intent(in)   :: xrng
real, dimension(:), intent(in ) :: U, V, T, Q
real, dimension(:), intent(in)  :: P
real,               intent(in)  :: PS, TS
real, dimension(:), intent(out) :: UdQdx, VdQdy, UdTdx, VdTdy, omega 


   integer,parameter:: R=6.370E+6
   real :: PI,Xlat,DX,DY

  PI = 4.0*ATAN(1.0)/180.
  xLAT = 20.0
  DX = 2.*1.25 * cos(xlat*PI)*PI*R
  DY = 2.*1.00 * PI*R




CALL CAL_NN(t, q, udtdx, xrng%udtdx, wgt%udtdx, NHID, NOUT, INPUT )
CALL CAL_NN(t, q, udqdx, xrng%udqdx, wgt%udqdx, NHID, NOUT, INPUT )
CALL CAL_NN(t, q, vdtdy, xrng%vdtdy, wgt%vdtdy, NHID, NOUT, INPUT )
CALL CAL_NN(t, q, vdqdy, xrng%vdqdy, wgt%vdqdy, NHID, NOUT, INPUT )



udtdx=udtdx/DX
udqdx=udqdx/DX
vdtdy=vdtdy/DY
vdqdy=vdqdy/DY


end subroutine NN_INTRFC

subroutine CAL_NN(t,q,out,xrange,w, NHID, NOUT, INPUT )
 
 integer, intent(in)              :: NHID, NOUT, INPUT

 real,   intent(in)  , dimension(:)    :: t, q    ! profile inputs
 real*8, intent(in)  , dimension(:)    :: w       ! Neural net weights
 real*8, intent(in)  , dimension(:,:)  :: xrange  ! Neural Net ranges
 real,   intent(out) , dimension(:)    :: out     ! outputs
  
 integer :: NDU,I,J,nwt,LM


 real*8, allocatable, dimension(:)     :: X_in, X_OUT
 real*8, allocatable, dimension(:)     :: X

 real, allocatable, dimension(:) ::twk,qwk,x1,x2,hrt ! working arrays
 character(len=120):: nfile


 LM  = SIZE(T)
 OUT = 0.0*T

 NDU = INPUT + NOUT
 NWT = NHID*(INPUT+1) + (NHID + 1)*NOUT

 allocate(twk(NOUT))
 allocate(qwk(NOUT))
 allocate( x1(NOUT))
 allocate( x2(NOUT))
 allocate(hrt(NOUT))

 allocate(x_in (0:INPUT))
 allocate(x_out(nout))
 allocate(X(NDU))



 x (1:40) = dble(log10(t(1:40)))
 x(41:80) = dble(q(1:40))

 call scaling(x,xrange,input,NDU)

 x_in(1:INPUT) = X(1:INPUT)

 call net(x_in,x_out,NHID,INPUT,nout,NWT,W)
 call unscaling(x_out,xrange,ndu,nout)

 out(1:40) = x_out(1:40)

 deallocate(x )
 deallocate(x_in)
 deallocate(x_out)
 deallocate(twk )
 deallocate(qwk )
 deallocate( x1 )
 deallocate( x2 )
 deallocate(hrt )

 end subroutine CAL_NN
  

subroutine net(x_in,x_out,NHID,INPUT,nout,nwt,W)

integer :: nhid, input, nout,nwt

real*8::x_in(0:input),x_out(nout)
real*8::W(nwt),B, SUM
real*8::G(0:NHID),GPRM(NHID)

integer :: is, im,k,i

B = 1.0D0         !SLOP
GPRM = 0.0D0
G(1:NHID) = 0.0D0
  X_in(0) = 1.0D0
  G(0)    = 1.0D0
IS=0
do k = 1, NHID
sum = 0.0D0
do I = 0, INPUT
IS = IS + 1
SUM = SUM + W(IS)*X_IN(I)
END DO
CALL FNET(SUM,G(K),B)
END DO
DO IM = 1, NOUT
SUM = 0.0D0
DO I = 0, NHID
IS = IS + 1
SUM = SUM + W(IS)*G(I)
END DO
X_OUT(im) = SUM             ! Ilinear = 1
!OUTPRM(im) = SUM
!call FNET(sum,X_OUT(im),B)
END DO

END subroutine net

SUBROUTINE FNET(SUM, Z, B)
real*8 :: SUM,Z,B
real*8 :: E, A

E = 2.718281828D0
A = B
IF(SUM.GE.0.0D0) Z = (1.0D0 - E**(-A*SUM))/(1.0D0 + E**(-A*SUM))
IF(SUM.LT.0.0D0) Z =-(1.0D0 - E**( A*SUM))/(1.0D0 + E**( A*SUM))

END SUBROUTINE FNET

subroutine scaling(x,TRANGE,input,NDU)
! ndu=input + output
integer :: input, ndu
real*8::x(input),TRANGE(NDU,3),srange
real*8,parameter::scal_min=-0.9D0,scal_max=0.9D0

integer :: io

srange = scal_max - scal_min
do io = 1, input  
!TRANGE(IO,1) = min(TRANGE(IO,1),X(IO))
!TRANGE(IO,2) = max(TRANGE(IO,2),X(IO))
!TRANGE(IO,3) = TRANGE(IO,2) - TRANGE(IO,1)
x(IO) = scal_min + ((x(IO)-TRANGE(IO,1))/TRANGE(IO,3))*srange
end do


end subroutine scaling

subroutine unscaling(x,TRANGE,ndu,nout)
real*8::x(nout),TRANGE(NDU,3),srange,T(3)
real*8,parameter::scal_min=-0.9D0,scal_max=0.9D0
integer :: ndu, nout

integer :: i,m

srange = scal_max - scal_min
do i = 1,nout
m = ndu-nout + i
!T(1) = min(TRANGE(m,1),X(I))
!T(2) = max(TRANGE(m,2),X(I))
!T(3) = T(2) - T(1)
T(1) = TRANGE(m,1)
T(2) = TRANGE(m,2)
T(3) = T(2) - T(1)
x(i) = ((x(i)-scal_min)*T(3))/srange + T(1)
end do

end subroutine unscaling

subroutine read_WEIGHT(w,TRANGE,iunit)
integer, intent(in)   :: iunit

!!!   real*8, dimension(:) ::W(nwt),TRANGE(NDU,3)
real*8, dimension(:)   ::W
real*8, dimension(:,:) ::TRANGE

integer :: nwt, ndu
integer :: i,j


nwt=size(W,1)

ndu=size(TRANGE,1)

do i = 1, nwt 
  read(iunit,*) w(i)
end do

!w=34.

do i = 1, NDU
  read(iunit,*) (TRANGE(i,j),j=1,3)
end do
close(iunit)

!trange=43.

end subroutine read_WEIGHT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ugly but necessary routine to initial myriad things needed for SCM 
! run.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine NN_SCMI(                                      &
                     nlayr,                              &
                     PREF_MODEL_E,                       &
                     time,                               &
                     yy,                                 &
                     mo,                                 &
                     dd,                                 &
                     hh,                                 &
                     mm,                                 & 
                     PCP_OBS,                            &
                     TS_AIR,                             &
                     TG_SOIL,                            &
                     TSKIN,                              &
                     QSFCAIR,                            &
                     QSKIN,                              &
                     PSFC,                               &
                     tt,                                 &
                     qq,                                 &
                     uu,                                 &
                     vv,                                 &
                     T_H_adv,                            &
                     T_V_adv,                            &
                     Q_H_adv,                            &
                     Q_V_adv,                            &
                     Q1,                                 &
                     Q2,                                 & 
                     OMG,                                & 
                     P_MODEL_E                           )

      INTEGER, parameter                           :: NT=1

      INTEGER,                       INTENT(IN   ) :: nlayr
      REAL,  DIMENSION(0:NLAYR),     INTENT(IN   ) :: PREF_MODEL_E

! OUTPUT
! Singli-layer:
      real, DIMENSION(nt       ), INTENT(  OUT) :: time       !Calenday day
      real, DIMENSION(nt       ), INTENT(  OUT) ::  yy        !Year
      real, DIMENSION(nt       ), INTENT(  OUT) ::  mo        !Month
      real, DIMENSION(nt       ), INTENT(  OUT) ::  dd        !Day
      real, DIMENSION(nt       ), INTENT(  OUT) ::  hh        !Hour
      real, DIMENSION(nt       ), INTENT(  OUT) ::  mm        !Minutes
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: PCP_OBS
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TS_AIR 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TG_SOIL 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TSKIN 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: QSKIN 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: QSFCAIR 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: PSFC

! Multiple-layer

      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: tt
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: qq
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: uu
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: vv
!     REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: ww  
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: T_H_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: T_V_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q_H_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q_V_adv
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q1
      REAL, DIMENSION(NT, nlayr), INTENT(  OUT) :: Q2
      REAL, DIMENSION(NT, 0:nlayr), INTENT(INOUT) :: OMG
      REAL, DIMENSION(NT, 0:nlayr), INTENT(  OUT) :: P_MODEL_E

      REAL, DIMENSION(NLAYR)       :: THETA00, TEMP, Q, PLO, PKO

      real :: r1,r2,r3,r4,r5,r6 ,lat0,lon0,hour,omeg0
      integer :: i0,i1,L,LM,iyyyy,imo,idd,ihh,imm,iss,I,UNIT,n
      REAL    :: SFCAIR_TEMP_0, SST_0 , OMG_0 


       P_MODEL_E(1,:) = PREF_MODEL_E(:)
       PLO   = ( P_MODEL_E(1,1:NLAYR) +  P_MODEL_E(1,0:NLAYR-1) )/2.0
       PKO   = ( PLO / MAPL_P00 )**MAPL_KAPPA
 
       THETA00= SFCAIR_TEMP_0

       THETA00 = THETA00 + (100000.- PLO)*30./90000.

       TEMP   = THETA00*PKO

       where( TEMP<190.)
         TEMP=190.
       endwhere

       Q      = 0.8*GEOS_QSAT (TEMP, PLO/100.)

       TT(1,:)= TEMP
       QQ(1,:)= Q
       

       iyyyy =  1900
       imo   =  7
       idd   =  1



       yy    =  1.0*iyyyy
       mo    =  1.0*imo   
       dd    =  1.0*idd
       hh    =  0.00
       mm    =  0.00
      

       time       =  -999.9
       PCP_OBS    = -999.9
       TS_AIR     = THETA00*( ( P_MODEL_E(1,NLAYR) / MAPL_P00 )**MAPL_KAPPA )
       TG_SOIL    = -999.9 
       PSFC       = P_MODEL_E(1,NLAYR)
       TSKIN      = SST_0
       QSFCAIR    = QQ(1,NLAYR)
       QSKIN      = GEOS_QSAT (TSKIN , PSFC/100.)

      uu = 0.0 
      vv = 0.0
      T_H_adv =0.
      T_V_adv =0.
      Q_H_adv =0.
      Q_V_adv =0.
      Q1      =0.
      Q2      = 0.
      OMEG0   = OMG_0 * 100./86400.

      do n=1,nt
           omg(n,:) = ( p_model_e(n,:)-psfc(n))*( 20000.-p_model_e(n,:) ) &
                       /( ( (psfc(n)-20000.)/2)**2 )
      end do

      where(omg < 0. )
        omg=0.
      endwhere

      omg = OMEG0 *  omg 

      do n=1,nlayr
         vv(:,n) = ((plo(n) - psfc(:) )/10000.)**2
      enddo

      vv=min(vv,5.0) +3.0
      uu=5.0


	write(*,*) 'Done!'

end subroutine NN_SCMI




end module NeuralNetMod
