MODULE math_routines

implicit none
private
public :: gammln, betai, betacf, rsq, prob_den_func, nbins,    &
     N_RANDOM_YEARS, shuffle, OPTIMIZ,  N_PARAMS !,init_MPI
real, save :: U(97), CS, CD, CM
integer, save :: I97, J97
integer, parameter :: N_PARAMS=3
integer, parameter :: nbins = 40, N_RANDOM_YEARS = 41

! initialize to non-MPI values
  include 'mpif.h'	
!integer,public   :: myid=0, numprocs=1, mpierr, mpistatus(MPI_STATUS_SIZE)  
!logical, public  :: root_proc=.true.

contains

!
! ----------------------------------------------------------------
!

subroutine Shuffle(a)  
  integer, intent(inout) :: a(:)  
  integer :: i, randpos, temp  
  real :: r   

  call init_random_seed ()

  do i = size(a), 2, -1    
     call random_number(r)    
     randpos = int(r * i) + 1    
     temp = a(randpos)    
     a(randpos) = a(i)    
     a(i) = temp 
   end do
 end subroutine Shuffle

!
! -------------------------------------------
!

 SUBROUTINE init_random_seed()
   INTEGER :: i, n, clock
   INTEGER, DIMENSION(:), ALLOCATABLE :: seed
   
   CALL RANDOM_SEED(size = n)
   ALLOCATE(seed(n))
   
   CALL SYSTEM_CLOCK(COUNT=clock)
   
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   CALL RANDOM_SEED(PUT = seed)
   
   DEALLOCATE(seed)
 
END SUBROUTINE init_random_seed

!
! ---------------------------------------------------------------------
!

  SUBROUTINE random_rain(eta)

!======================================================================
! NPAC $Id: math_routines.F90,v 1.1 2018/05/10 15:05:29 smahanam Exp $
!======================================================================
! _____________________________________________________________________
!                                                                     |
! Test program to generate a vector of Gaussian deviates using        |
! the Box-Muller method and the system-supplied uniform RNG.          |
! K.A.Hawick, 15 July 1994.                                           |
!                                                                     |
! H W Yau.  23rd of August, 1996.                                     |
! Northeast Parallel Architectures Center.                            |
! Syracuse University.                                                |
!_____________________________________________________________________|
!                                                                     |
! Edit record:                                                        |
! --/Jul/1996: M.McMahon -- ALIGN statement and EXTRINSIC statements  |
!              PROCESSORS statement as well.                          |
! 23/Aug/1996: HWY.  Cleaning up.                                     |
!              Removed superfluous definition of module.              |
! 26/Aug/1996: Added `IMPLICIT NONE'.                                 |
!              Fixed HPF directives.  Cleaned up interface block to   |
!              gaussian_vector().                                     |
!              F90 version.                                           |
!_____________________________________________________________________|
!

      IMPLICIT NONE

      INTEGER, PARAMETER :: n = N_RANDOM_YEARS
      REAL ti, tf,tcalc,tcom,tend,start
      REAL, DIMENSION(1:n) :: a
      REAL, DIMENSION(1:N)::ETA
!
!  The vector of Guassian deviates
      INTEGER, DIMENSION(1:n) :: j
!
! Histograms to check the distribution
      INTEGER, DIMENSION(1:nbins) :: bins
      INTEGER i,less, more,nb
      INTEGER :: nover = 0
!
!      INTERFACE
!         SUBROUTINE timer(return_time, initial_time)
!         REAL, INTENT(IN)  :: initial_time
!         REAL, INTENT(OUT) :: return_time
!         END SUBROUTINE
!      END INTERFACE
!
!      INTERFACE
!         SUBROUTINE gaussian_vector(x,n,tcalc2,tcomm2)
!         REAL, DIMENSION(:) :: x
!         REAL,INTENT(INOUT) :: tcalc2,tcomm2
!         INTEGER,INTENT(IN) :: n
!         REAL ETA
!         END SUBROUTINE gaussian_vector
!      END INTERFACE
!______________________________________________________________________
!
! Executable Code.
!______________________________________________________________________
!	
! start_timer
      CALL timer(start,0.0)
      tcom  = 0.0
      tcalc = 0.0
      CALL timer(ti,0.0)
      bins(1:nbins) = 0
      CALL timer(tend,ti)
      tcalc = tcalc + tend
!
! Ask for n deviate
      CALL gaussian_vector( a, n,tend,tcom )
      tcalc = tcalc + tend
      CALL timer(ti,0.0)
!F95  FORALL(i=1:n)j(i) = a(i) * real(nbins/7) + nbins/2 + 1
      DO i=1,n
         j(i) = a(i) * real(nbins/7) + nbins/2 + 1
      END DO
      CALL timer(tend,ti)
      tcalc = tcalc + tend
!	
      CALL timer(ti,0.0)
      DO i=1,n
         ETA(I)=(REAL(J(I))-10.-1.)/5.
      ENDDO
!      write(*,*)eta
!      write(*,*)sum(j)/300.
      DO nb = 1,nbins
         bins(nb) = count(j.EQ.nb,1)
!CCCCCCCC write(6,*) nb, bins(nb)
      ENDDO
      CALL timer(tend,ti)
      tcom = tcom + tend
!
! Stop timer.
      CALL timer(tf,start)
!
! Correction for serial execution.
      tcalc = tcalc + tcom
      tcom = 0.0
!      WRITE(6,6001) 0,n,
!     1              tcom,tcalc,(tf-tcom-tcalc),tf
!
! Write histogram
      DO nb=1,nbins
!        WRITE(6,*) nb,(REAL(NB)-10.-1.)/5., bins(nb)
      ENDDO
!
!      STOP
      RETURN
!
 6001 FORMAT('Number of Processors = ',I4/  &
            'Problem size = ',I6/           &
            'Communications = ',F9.3/       &
            'Compute        = ',F9.3/       &
            'Others         = ',F9.3/       &
            'Total time     = ',F9.3)
!
    END SUBROUTINE random_rain
!
! ------------------------------------------------------
!

      SUBROUTINE gaussian_vector( x, n, tcalc2,tcomm2 )
        
        IMPLICIT NONE
        !
        ! Box-Muller Method
        ! See Knuth, Vol 2, 2nd Edn, PP 117
        INTEGER,INTENT(IN) :: n
        REAL, DIMENSION(:) :: x
        REAL,INTENT(INOUT) :: tcalc2,tcomm2
        REAL, DIMENSION(:), ALLOCATABLE :: v1, v2, r, f
        REAL ti,cal1,cal2,com1,com2 
        LOGICAL, DIMENSION(:), ALLOCATABLE :: mask
        !
        ! Accept/reject efficiency is about 1.27
        integer m, i, np
        !______________________________________________________________________
        !
        ! Executable code.
        !______________________________________________________________________
        !
        CALL timer(ti,0.0)
        np = n
        !
        ! avoid fluctuation problems
        IF( np .lt. 10 ) np = 10
        ALLOCATE( v1(np), v2(np), r(np), f(np), mask(np) )
        !
        ! Generate two deviates:
        CALL random_number( v1 )
        CALL random_number( v2 )
        v1 = 2.0 * v1 -1.0
        v2 = 2.0 * v2 -1.0
        r = v1**2 + v2**2
        !
        ! are they in the unit circle?
        mask = r < 1.0
        r = merge( r, 0.5, mask )
        f = sqrt( -2.0 * log(r) / r)
        v1 = v1 * f
        v2 = v2 * f
        !
        ! since pack is a new intrinsic, here is serial code to show
        ! what is being done.
        !      m = 0
        !      do i=1,np
        !        if( mask(i) )then
        !          m = m + 1
        !          r(m) = v1(i)
        !          f(m) = v2(i)
        !        endif
        !      enddo 
        CALL timer(cal1,ti)
        
        CALL timer(ti,0.0)
        !
        ! we now have 2 * m deviates
        m = count( mask )
        CALL timer(com1,ti)
        
        CALL timer(ti,0.0)
        !
        ! and to save space, we'll store them in r and f
        !      r = pack( v1, mask )
        !      f = pack( v2, mask )
        !
        ! We now get performance at the expense of memory.
        r = v1
        f = v2
        !
        ! Statistically, this should not happen for large n
        IF( 2*m .lt. n )THEN
           WRITE(6,*) 'Not enough deviates! Got: ', 2*m, ', Needed: ', n
           !        WRITE(6,*) 'Increase accept reject allowance in',
           !                   ' xgaussian_vector'
           STOP
        ENDIF
        !
        ! use the two result vectors to patch up enough as	
        IF( m .LT. n )THEN
           x(1:m)=r(1:m)
           CALL timer(cal2,ti)
           CALL timer(ti,0.0)
           x(m+1:n) = f(1:n-m)
           CALL timer(com2,ti)
        ELSE
           x(1:n) = r(1:n)
        ENDIF
        
        tcalc2 = cal1 + cal2
        tcomm2 = com1 + com2
        DEALLOCATE( v1, v2, r, f, mask )
        
        RETURN
      END SUBROUTINE gaussian_vector
!
! ------------------------------------------------------------------- 
!
      SUBROUTINE timer(return_time, initial_time)
        implicit none
        REAL, INTENT(IN)  :: initial_time
        REAL, INTENT(OUT) :: return_time
        INTEGER finish,rate
        CALL system_clock( COUNT=finish,COUNT_RATE=rate)
        return_time = FLOAT(finish) / FLOAT(rate) - initial_time
        RETURN
      END SUBROUTINE timer
!
! ------------------------------------------------------------------- 
!

subroutine prob_den_func (ndata, datain, cdf, bins, &
     mean, std, skew, pdf, lwval, upval)

implicit none
integer, intent (in) :: ndata
real, dimension(ndata), intent (in) :: datain
real, intent(out), dimension(nbins) :: cdf, bins
real, optional, intent(out), dimension(nbins) :: pdf
real, optional, intent (out) :: mean, std, skew  
real, optional               :: lwval, upval
real :: lw,up,db,var1,var2, variance
integer :: i,n

lw = minval (datain)
up = maxval (datain)

if(present(upval)) then
   up = upval
   lw = lwval
endif

db = (up - lw)/real(nbins)
cdf = 0.
if(present(pdf)) pdf =0.

do n = 1, nbins
   bins (n) = lw + real(n)*db - db/2.
   do i = 1,ndata
      if(datain(i) <= bins (n) + db/2.) cdf(n) = cdf(n) + 1.
      if(present(pdf)) then
         if(n==1) then
            if((datain(i) >=  bins (n) - db/2.).and.   &
                 (datain(i) <= bins (n) + db/2.))      &
                 pdf(n) = pdf(n) + 1.
         else
            if((datain(i) >   bins (n) - db/2.).and.   &
                 (datain(i) <= bins (n) + db/2.))      &
                 pdf(n) = pdf(n) + 1.
         endif
      endif
   end do
end do

cdf = cdf/real(ndata)

if(present(pdf)) pdf = pdf/real(ndata)

if(present (mean)) mean = sum(datain)/real(ndata)

if(present (std)) then

   var1 = 0.
   mean = sum(datain)/real(ndata)

   do i = 1,ndata
      var1 = var1 + (datain(i) - mean)*(datain(i) - mean)
   end do

   std = sqrt (var1/real(ndata - 1))

   if(present (skew)) then
      var2 = 0.
      do i = 1,ndata
         var2 = var2 + ((datain(i) - mean)/std)* &
             ((datain(i) - mean)/std)*           &
             ((datain(i) - mean)/std) 
      end do

      skew = var2/real(ndata - 1)

   endif

endif

END subroutine prob_den_func

!
! -------------------------------------------------------
!

FUNCTION betai(a,b,x)
REAL betai,a,b,x
REAL bt
!external gammln

if(x.lt.0..or.x.gt.1.)print *, 'bad argument x in betai',x
if(x.lt.0..or.x.gt.1.)stop
if(x.eq.0..or.x.eq.1.)then
   bt=0.
else 
   bt=exp(gammln(a+b)-gammln(a)-gammln(b) &
    +a*log(x)+b*log(1.-x))
endif

if(x.lt.(a+1.)/(a+b+2.))then 
   betai=bt*betacf(a,b,x)/a
   return
else
   betai=1.-bt*betacf(b,a,1.-x)/b 
   return
endif
END FUNCTION betai
!
! -------------------------------------------------------
!
FUNCTION betacf(a,b,x)
INTEGER MAXIT
REAL betacf,a,b,x,EPS,FPMIN
PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
INTEGER m,m2
REAL aa,c,d,del,h,qab,qam,qap

qab=a+b 
qap=a+1. 
qam=a-1.
c=1. 
d=1.-qab*x/qap

if(abs(d).lt.FPMIN)d=FPMIN
d=1./d
h=d
do m=1,MAXIT
   m2=2*m
   aa=m*(b-m)*x/((qam+m2)*(a+m2))
   d=1.+aa*d 
   if(abs(d).lt.FPMIN)d=FPMIN
   c=1.+aa/c
   if(abs(c).lt.FPMIN)c=FPMIN
   d=1./d
   h=h*d*c
   aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
   d=1.+aa*d 
   if(abs(d).lt.FPMIN)d=FPMIN
   c=1.+aa/c
   if(abs(c).lt.FPMIN)c=FPMIN
   d=1./d
   del=d*c
   h=h*del
   if(abs(del-1.).lt.EPS)exit 
enddo
 betacf=h
return
END FUNCTION betacf
!
! --------------------------------------------------------------
!
FUNCTION gammln(xx)
REAL gammln,xx
INTEGER j
DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)

SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,          &
 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
 -.5395239384953d-5,2.5066282746310005d0/
x=xx
y=x
tmp=x+5.5d0
tmp=(x+0.5d0)*log(tmp)-tmp
ser=1.000000000190015d0
do  j=1,6
   y=y+1.d0
   ser=ser+cof(j)/y
enddo
gammln=tmp+log(stp*ser/x)
return
END FUNCTION gammln

!
! ------------------------------------------------
!

SUBROUTINE RSQ (NDATA, X, Y, R2, RMSE, limits, slope, intercept)

implicit none
integer, intent (in) :: NDATA
real, dimension (ndata), intent(in) :: x,y
real, optional, dimension(4), intent (in) :: limits
real, optional, intent(out) :: slope, intercept, RMSE
real, intent(out) :: r2

integer ::n
real ::  ic
real :: sumx, sumy,sumxy, sumx2,sumy2,sig2x,sig2y,sig2xy,error

ic  =0. 
sumx=0.
sumy=0.
sumxy=0.
sumy2=0.
sumx2=0.
error=0.
do n = 1,ndata

   if(present (limits)) then
     if((x(n) > limits(1)).and.(x(n) < limits(2)).and. &
        (y(n) > limits(3)).and.(y(n) < limits(4))) then
        ic = ic + 1.
        sumx = sumx + x(n)
        sumy = sumy + y(n)
        sumxy= sumxy + x(n)*y(n)
        sumy2= sumy2 + y(n)*y(n)
        sumx2= sumx2 + x(n)*x(n)
        error= error + (y(n) - x(n))*(y(n) - x(n))
     endif
   else 
      ic = ic + 1.
      sumx = sumx + x(n)
      sumy = sumy + y(n)
      sumxy= sumxy + x(n)*y(n)
      sumy2= sumy2 + y(n)*y(n)
      sumx2= sumx2 + x(n)*x(n)
      error= error + (y(n) - x(n))*(y(n) - x(n))
   endif

end do

if (present(intercept)) intercept = -9999.
if (present(slope))     slope     = -9999.
if (present(rmse))      rmse      = -9999.
r2 = -9999.

if(ic /= 0) then
   if(ic*sumx2 /= sumx*sumx) then
      
      if (present(intercept)) intercept = &
           (sumy*sumx2 - sumx*sumxy) / (ic*sumx2 - sumx*sumx)
      if (present(slope))     slope     = &
           (ic*sumxy - sumx*sumy) / (ic*sumx2 - sumx*sumx)
      if (present(rmse)) rmse = sqrt(error/real(ic))
   endif
   
   sumx =sumx/ic
   sumy =sumy/ic
   sumxy=sumxy/ic
   sumy2=sumy2/ic
   sumx2=sumx2/ic
   sig2x=sumx2-sumx*sumx
   sig2y=sumy2-sumy*sumy
   sig2xy=(sumxy-sumx*sumy)*(sumxy-sumx*sumy)
   
   r2 = sig2xy/(sig2x*sig2y + 1.e-20)
endif

END SUBROUTINE RSQ

!
! ------------------------------------------
!
     SUBROUTINE OPTIMIZ(NDATA,wet,eff,X)

     implicit none

     integer, PARAMETER :: N = N_PARAMS, NEPS = 4
     integer, intent (in) :: ndata
     REAL (kind = 8) ::  LB(N), UB(N), X(N), XOPT(N), CON(N), VM(N), &
                       FSTAR(NEPS), XP(N), T, EPS, RT, FOPT,        &
                       EFF(NDATA),WET(NDATA)

      INTEGER  NACP(N), NS, NT, NFCNEV, IER, ISEED1, ISEED2, &
               MAXEVL, IPRINT, NACC, NOBDS, I                   

      LOGICAL  MAX

!  Set underflows to zero on IBM mainframes.
!      CALL XUFLOW(0)

      MAX = .false.
      EPS = 1.0D-6
      RT =  .5
      ISEED1 = 1
      ISEED2 = 2
      NS = 20
      NT = 5
      MAXEVL = 100000
      IPRINT = 0
      
      DO I=1,N
         LB(I)=0.0001
         UB(I)=100.
         CON(I) = 2.0
      END DO
      UB(3)= 1.
      LB(3)= 0.1
!
!  Set input values of the input/output parameters.
!
      T = 5.0
      DO I = 1, N
         VM(I) = 1.0
      END DO

!      WRITE(*,1000) N, MAX, T, RT, EPS, NS, NT, NEPS, MAXEVL, IPRINT, &
!                   ISEED1, ISEED2
!
!      CALL PRTVEC(X,N,'STARTING VALUES')
!      CALL PRTVEC(VM,N,'INITIAL STEP LENGTH')
!      CALL PRTVEC(LB,N,'LOWER BOUND')
!      CALL PRTVEC(UB,N,'UPPER BOUND')
!      CALL PRTVEC(C,N,'C VECTOR')
!      WRITE(*,'(/,''  ****   END OF DRIVER ROUTINE OUTPUT   ****''     &
!               /,''  ****   BEFORE CALL TO SA.             ****'')') 

      CALL SA(X,MAX,RT,EPS,NS,NT,NEPS,MAXEVL,LB,UB,CON,IPRINT,ISEED1, &
             ISEED2,T,VM,XOPT,FOPT,NACC,NFCNEV,NOBDS,IER,             &
             FSTAR,XP,NACP,NDATA,EFF,WET)

!      WRITE(*,'(/,''  ****   RESULTS AFTER SA   ****   '')') 
!      CALL PRTVEC(XOPT,N,'SOLUTION')
!      CALL PRTVEC(VM,N,'FINAL STEP LENGTH')
!      WRITE(*,1001) FOPT, NFCNEV, NACC, NOBDS, T, IER

1000  FORMAT(/,' SIMULATED ANNEALING EXAMPLE',/,                        &
            /,' NUMBER OF PARAMETERS: ',I3,'   MAXIMAZATION: ',L5,      &
            /,' INITIAL TEMP: ', G8.2, '   RT: ',G8.2, '   EPS: ',G8.2, &
            /,' NS: ',I3, '   NT: ',I2, '   NEPS: ',I2,                 &
            /,' MAXEVL: ',I10, '   IPRINT: ',I1, '   ISEED1: ',I4,      &
            '   ISEED2: ',I4)                                           
1001  FORMAT(/,' OPTIMAL FUNCTION VALUE: ',G20.13              &
            /,' NUMBER OF FUNCTION EVALUATIONS:     ',I10,     &
            /,' NUMBER OF ACCEPTED EVALUATIONS:     ',I10,     &
            /,' NUMBER OF OUT OF BOUND EVALUATIONS: ',I10,     &
            /,' FINAL TEMP: ', G20.13,'  IER: ', I3)           
                                                               
      RETURN
 END SUBROUTINE OPTIMIZ

!
! ---------------------------------------------------------------
!

      SUBROUTINE SA(X,MAX,RT,EPS,NS,NT,NEPS,MAXEVL,LB,UB,CON,IPRINT,    &
                   ISEED1,ISEED2,T,VM,XOPT,FOPT,NACC,NFCNEV,NOBDS,IER, &
                   FSTAR,XP,NACP,NDATA,EFF,WET)

!  Type all external variables.

      INTEGER,PARAMETER  :: N=N_PARAMS
      integer :: ndata
      REAL (KIND = 8)  X(N), LB(N), UB(N), CON(N), VM(N), FSTAR(N),  &
                       XOPT(N), XP(N), T, EPS, RT, FOPT
      REAL (KIND = 8) EFF(NDATA),WET(NDATA) 
      INTEGER  NACP(N), NS, NT, NEPS, NACC, MAXEVL, IPRINT,       &
              NOBDS, IER, NFCNEV, ISEED1, ISEED2
      LOGICAL  MAX

!  Type all internal variables.
      REAL (KIND = 8)  F, FP, P, PP, RATIO
      INTEGER  NUP, NDOWN, NREJ, NNEW, LNOBDS, H, I, J, M
      LOGICAL  QUIT

!  Type all functions.
!      REAL (KIND = 8)  EXPREP

!  Initialize the random number generator RANMAR.
      CALL RMARIN(ISEED1,ISEED2)

!  Set initial values.
      NACC = 0
      NOBDS = 0
      NFCNEV = 0
      IER = 99

      DO I = 1, N
         XOPT(I) = X(I)
         NACP(I) = 0
      END DO

      DO I = 1, NEPS
         FSTAR(I) = 1.0D+20
      END DO

!  If the initial temperature is not positive, notify the user and 
!  return to the calling routine. 
      IF (T .LE. 0.0) THEN
         WRITE(*,'(/,''  THE INITIAL TEMPERATURE IS NOT POSITIVE. '' &
                  /,''  RESET THE VARIABLE T. ''/)')
         IER = 3
         RETURN
      END IF

!  If the initial value is out of bounds, notify the user and return
!  to the calling routine.
      DO I = 1, N
         IF ((X(I) .GT. UB(I)) .OR. (X(I) .LT. LB(I))) THEN
            CALL PRT1
            IER = 2
            RETURN
         END IF
      END DO

!  Evaluate the function with input X and return value as F.
      CALL FCN(X,F,NDATA,EFF,WET)

!  If the function is to be minimized, switch the sign of the function.
!  Note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.
      IF(.NOT. MAX) F = -F
      NFCNEV = NFCNEV + 1
      FOPT = F
      FSTAR(1) = F
      IF(IPRINT .GE. 1) CALL PRT2(MAX,N,X,F)
 
!  Start the main loop. Note that it terminates if (i) the algorithm
!  succesfully optimizes the function or (ii) there are too many
!  function evaluations (more than MAXEVL).
100   NUP = 0
      NREJ = 0
      NNEW = 0
      NDOWN = 0
      LNOBDS = 0

      DO M = 1, NT
         DO J = 1, NS
            DO H = 1, N

!  Generate XP, the trial value of X. Note use of VM to choose XP.
               DO I = 1, N
                  IF (I .EQ. H) THEN
                     XP(I) = X(I) + (RANMAR()*2.- 1.) * VM(I)
                  ELSE
                     XP(I) = X(I)
                  END IF

!  If XP is out of bounds, select a point in bounds for the trial.
                  IF((XP(I) .LT. LB(I)) .OR. (XP(I) .GT. UB(I))) THEN
                    XP(I) = LB(I) + (UB(I) - LB(I))*RANMAR()
                    LNOBDS = LNOBDS + 1
                    NOBDS = NOBDS + 1
                    IF(IPRINT .GE. 3) CALL PRT3(MAX,N,XP,X,FP,F)
                  END IF
               END DO
!  Evaluate the function with the trial point XP and return as FP.
      CALL FCN(XP,FP,NDATA,EFF,WET)
               IF(.NOT. MAX) FP = -FP
               NFCNEV = NFCNEV + 1
               IF(IPRINT .GE. 3) CALL PRT4(MAX,N,XP,X,FP,F)

!  If too many function evaluations occur, terminate the algorithm.
               IF(NFCNEV .GE. MAXEVL) THEN
                  CALL PRT5
                  IF (.NOT. MAX) FOPT = -FOPT
                  IER = 1
                  RETURN
               END IF

!  Accept the new point if the function value increases.
               IF(FP .GE. F) THEN
                  IF(IPRINT .GE. 3) THEN
                     WRITE(*,'(''  POINT ACCEPTED'')')
                  END IF
                  DO I = 1, N
                     X(I) = XP(I)
                  END DO
                  F = FP
                  NACC = NACC + 1
                  NACP(H) = NACP(H) + 1
                  NUP = NUP + 1
                  
!  If greater than any other point, record as new optimum.
                  IF (FP .GT. FOPT) THEN
                     IF(IPRINT .GE. 3) THEN
                        WRITE(*,'(''  NEW OPTIMUM'')')
                     END IF
                     DO I = 1, N
                        XOPT(I) = XP(I)
                     END DO
                     FOPT = FP
                     NNEW = NNEW + 1
                  END IF

                  !  If the point is lower, use the Metropolis criteria to decide on
                  !  acceptance or rejection.
               ELSE
                  P = EXPREP((FP - F)/T)
                  PP = RANMAR()
                  IF (PP .LT. P) THEN
                     IF(IPRINT .GE. 3) CALL PRT6(MAX)
                     DO I = 1, N
                        X(I) = XP(I)
                     END DO
                     F = FP
                     NACC = NACC + 1
                     NACP(H) = NACP(H) + 1
                     NDOWN = NDOWN + 1
                  ELSE
                     NREJ = NREJ + 1
                     IF(IPRINT .GE. 3) CALL PRT7(MAX)
                  END IF
               END IF

            END DO
         END DO

!  Adjust VM so that approximately half of all evaluations are accepted.
         DO I = 1, N
            RATIO = DFLOAT(NACP(I)) /DFLOAT(NS)
            IF (RATIO .GT. .6) THEN
               VM(I) = VM(I)*(1. + CON(I)*(RATIO - .6)/.4)
            ELSE IF (RATIO .LT. .4) THEN
               VM(I) = VM(I)/(1. + CON(I)*((.4 - RATIO)/.4))
            END IF
            IF (VM(I) .GT. (UB(I)-LB(I))) THEN
               VM(I) = UB(I) - LB(I)
            END IF
         END DO

         IF(IPRINT .GE. 2) THEN
            CALL PRT8(N,VM,XOPT,X)
         END IF

         DO I = 1, N
            NACP(I) = 0
         END DO

      END DO

      IF(IPRINT .GE. 1) THEN
         CALL PRT9(MAX,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
      END IF

!  Check termination criteria.
      QUIT = .FALSE.
      FSTAR(1) = F
      IF ((FOPT - FSTAR(1)) .LE. EPS) QUIT = .TRUE.
      DO I = 1, NEPS
         IF (ABS(F - FSTAR(I)) .GT. EPS) QUIT = .FALSE.
      END DO

!  Terminate SA if appropriate.
      IF (QUIT) THEN
         DO I = 1, N
            X(I) = XOPT(I)
         END DO
         IER = 0
         IF (.NOT. MAX) FOPT = -FOPT
         IF(IPRINT .GE. 1) CALL PRT10
         RETURN
      END IF

!  If termination criteria is not met, prepare for another loop.
      T = RT*T
      DO I = NEPS, 2, -1
         FSTAR(I) = FSTAR(I-1)
      END DO
      F = FOPT
      DO I = 1, N
         X(I) = XOPT(I)
      END DO

!  Loop again.
      GO TO 100

    END SUBROUTINE SA
!
! ---------------------------------------------------------
!
    
  SUBROUTINE FCN(X,F,NDATA,EFF,WET)
    implicit none
    integer, parameter :: N = N_PARAMS 
    REAL (KIND = 8) X(n),F
    INTEGER, intent (in) :: ndata
    INTEGER Nd
    REAL (KIND = 8) EFF(NDATA),WET(NDATA) 
    real :: a,b,xv, r2, rmse
    real , dimension(ndata) :: yval, xval
    
    a = X(1)
    b = X(2)
    
    do nd = 1, ndata
       xv = wet(nd)
       xval(nd) = eff(nd)
       yval(nd) = X(3)*betai(a,b,xv)
    end do

!   if(maxval(yval) > 1.) print *,maxval(yval), minval(yval)

    call rsq (ndata, xval, yval, r2, rmse)
    F = rmse
!    if(maxval(yval) > 1.) F = 1.d10
!    F = r2
!    if(F > 1.) F = -1.d-10
    
    RETURN
  END SUBROUTINE FCN
  
!
! ---------------------------------------------------------------
!
      FUNCTION  EXPREP(RDUM)
        implicit none
        REAL (KIND = 8)  RDUM, EXPREP
        
        IF (RDUM .GT. 174.) THEN
           EXPREP = 3.69D+75
        ELSE IF (RDUM .LT. -180.) THEN
           EXPREP = 0.0
        ELSE
           EXPREP = EXP(RDUM)
        END IF
        
        RETURN
      END FUNCTION EXPREP
!
! ---------------------------------------------
!
      subroutine RMARIN(IJ,KL)

        implicit none
        integer, intent (in) :: ij,kl
        integer :: i,j,k,l,m, ii,jj
        real :: t,s

!        real U(97), C, CD, CM
!        integer I97, J97
!        common /raset1/ U, C, CD, CM, I97, J97

        if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.  &
             KL .lt. 0  .or.  KL .gt. 30081 ) then
           print '(A)', ' The first random number seed must have a value between 0 and 31328'
           print '(A)',' The second seed must have a value between 0 and 30081'
           stop
        endif

        i = mod(IJ/177, 177) + 2
        j = mod(IJ    , 177) + 2
        k = mod(KL/169, 178) + 1
        l = mod(KL,     169)

        do ii = 1, 97
           s = 0.0
           t = 0.5
           do jj = 1, 24
              m = mod(mod(i*j, 179)*k, 179)
              i = j
              j = k
              k = m
              l = mod(53*l+1, 169)
              if (mod(l*m, 64) .ge. 32) then
                 s = s + t
              endif
              t = 0.5 * t
           end do
           U(ii) = s
        end do
        
        
        CS = 362436.0 / 16777216.0
        CD = 7654321.0 / 16777216.0
        CM = 16777213.0 /16777216.0
        I97 = 97
        J97 = 33
        return
      end subroutine RMARIN
     
!
! ---------------------------------------
!         
  
      real function ranmar()
        real :: uni
!      real U(97), C, CD, CM
!      integer I97, J97
!      common /raset1/ U, C, CD, CM, I97, J97
         uni = U(I97) - U(J97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         U(I97) = uni
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         CS = CS - CD
         if( CS .lt. 0.0 ) CS = CS + CM
         uni = uni - CS
         if( uni .lt. 0.0 ) uni = uni + 1.0
         RANMAR = uni
      return
    END function ranmar

    SUBROUTINE PRT1
      implicit none
      WRITE(*,'(/,''  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS ''  &
                /,''  (LB AND UB). EXECUTION TERMINATED WITHOUT ANY''  &
                /,''  OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  ''  &
                /,''  LB(I) .LT. X(I) .LT. UB(I), I = 1, N. ''/)')

      RETURN
    END SUBROUTINE PRT1

    SUBROUTINE PRT2(MAX,N,X,F)
      implicit none
      REAL (KIND = 8)  X(N), F
      INTEGER  N
      LOGICAL  MAX

      WRITE(*,'(''  '')')
      CALL PRTVEC(X,N,'INITIAL X')
      IF (MAX) THEN
         WRITE(*,'(''  INITIAL F: '',/, G25.18)') F
      ELSE
         WRITE(*,'(''  INITIAL F: '',/, G25.18)') -F
      END IF
      
      RETURN
    END SUBROUTINE PRT2

    SUBROUTINE PRT3(MAX,N,XP,X,FP,F)
      implicit none
      REAL (KIND = 8)  XP(N), X(N), FP, F
      INTEGER  N
      LOGICAL  MAX
      
      WRITE(*,'(''  '')')
      CALL PRTVEC(X,N,'CURRENT X')
      IF (MAX) THEN
         WRITE(*,'(''  CURRENT F: '',G25.18)') F
      ELSE
         WRITE(*,'(''  CURRENT F: '',G25.18)') -F
      END IF
      CALL PRTVEC(XP,N,'TRIAL X')
      WRITE(*,'(''  POINT REJECTED SINCE OUT OF BOUNDS'')')

      RETURN
    END SUBROUTINE PRT3

    SUBROUTINE PRT4(MAX,N,XP,X,FP,F)
      implicit none
      REAL (KIND = 8)  XP(N), X(N), FP, F
      INTEGER  N
      LOGICAL  MAX
      
      WRITE(*,'(''  '')')
      CALL PRTVEC(X,N,'CURRENT X')
      IF (MAX) THEN
         WRITE(*,'(''  CURRENT F: '',G25.18)') F
         CALL PRTVEC(XP,N,'TRIAL X')
         WRITE(*,'(''  RESULTING F: '',G25.18)') FP
      ELSE
         WRITE(*,'(''  CURRENT F: '',G25.18)') -F
         CALL PRTVEC(XP,N,'TRIAL X')
         WRITE(*,'(''  RESULTING F: '',G25.18)') -FP
      END IF
      
      RETURN
    END SUBROUTINE PRT4
    
    SUBROUTINE PRT5
      implicit none
      WRITE(*,'(/,''  TOO MANY FUNCTION EVALUATIONS; CONSIDER '' &
          /,''  INCREASING MAXEVL OR EPS, OR DECREASING ''           &
          /,''  NT OR RT. THESE RESULTS ARE LIKELY TO BE ''          &
          /,''  POOR.'',/)')

      RETURN
    END SUBROUTINE PRT5

    SUBROUTINE PRT6(MAX)
      implicit none
      LOGICAL  MAX
      
      IF (MAX) THEN
         WRITE(*,'(''  THOUGH LOWER, POINT ACCEPTED'')')
      ELSE
         WRITE(*,'(''  THOUGH HIGHER, POINT ACCEPTED'')')
      END IF
      
      RETURN
    END SUBROUTINE PRT6
    
    SUBROUTINE PRT7(MAX)
      implicit none
      LOGICAL  MAX
      
      IF (MAX) THEN
         WRITE(*,'(''  LOWER POINT REJECTED'')')
      ELSE
         WRITE(*,'(''  HIGHER POINT REJECTED'')')
      END IF
      
      RETURN
    END SUBROUTINE PRT7
    
    SUBROUTINE PRT8(N,VM,XOPT,X)
      implicit none
      REAL (KIND = 8)  VM(N), XOPT(N), X(N)
      INTEGER  N
      
      WRITE(*,'(/,'' INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT'',/)')
      CALL PRTVEC(VM,N,'NEW STEP LENGTH (VM)')
      CALL PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
      CALL PRTVEC(X,N,'CURRENT X')
      WRITE(*,'('' '')')
      
      RETURN
    END SUBROUTINE PRT8
    
    SUBROUTINE PRT9(MAX,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
      implicit none
      REAL (KIND = 8)  XOPT(N), VM(N), T, FOPT
      INTEGER  N, NUP, NDOWN, NREJ, LNOBDS, NNEW, TOTMOV
      LOGICAL  MAX
      
      TOTMOV = NUP + NDOWN + NREJ
      
      WRITE(*,'(/, '' INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION'',/)')
      WRITE(*,'(''  CURRENT TEMPERATURE:            '',G12.5)') T
      IF (MAX) THEN
         WRITE(*,'(''  MAX FUNCTION VALUE SO FAR:  '',G25.18)') FOPT
         WRITE(*,'(''  TOTAL MOVES:                '',I8)') TOTMOV
         WRITE(*,'(''     UPHILL:                  '',I8)') NUP
         WRITE(*,'(''     ACCEPTED DOWNHILL:       '',I8)') NDOWN
         WRITE(*,'(''     REJECTED DOWNHILL:       '',I8)') NREJ
         WRITE(*,'(''  OUT OF BOUNDS TRIALS:       '',I8)') LNOBDS
         WRITE(*,'(''  NEW MAXIMA THIS TEMPERATURE:'',I8)') NNEW
      ELSE
         WRITE(*,'(''  MIN FUNCTION VALUE SO FAR:  '',G25.18)') -FOPT
         WRITE(*,'(''  TOTAL MOVES:                '',I8)') TOTMOV
         WRITE(*,'(''     DOWNHILL:                '',I8)')  NUP
         WRITE(*,'(''     ACCEPTED UPHILL:         '',I8)')  NDOWN
         WRITE(*,'(''     REJECTED UPHILL:         '',I8)')  NREJ
         WRITE(*,'(''  TRIALS OUT OF BOUNDS:       '',I8)')  LNOBDS
         WRITE(*,'(''  NEW MINIMA THIS TEMPERATURE:'',I8)')  NNEW
      END IF
      CALL PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
      CALL PRTVEC(VM,N,'STEP LENGTH (VM)')
      WRITE(*,'('' '')')
      
      RETURN
    END SUBROUTINE PRT9
    
    SUBROUTINE PRT10
      implicit none
      WRITE(*,'(/,''  SA ACHIEVED TERMINATION CRITERIA. IER = 0. '',/)')
      
      RETURN
    END SUBROUTINE PRT10
    
    SUBROUTINE PRTVEC(VECTOR,NCOLS,NAME)
      implicit none
      INTEGER NCOLS, LL,I,J, LINES
      REAL (KIND = 8) VECTOR(NCOLS)
      CHARACTER *(*) NAME
      
      WRITE(*,1001) NAME
      
      IF (NCOLS .GT. 10) THEN
         LINES = INT(NCOLS/10.)
         
         DO I = 1, LINES
            LL = 10*(I - 1)
            WRITE(*,1000) (VECTOR(J),J = 1+LL, 10+LL)
         END DO
            
         WRITE(*,1000) (VECTOR(J),J = 11+LL, NCOLS)
      ELSE
         WRITE(*,1000) (VECTOR(J),J = 1, NCOLS)
      END IF
      
1000  FORMAT( 10(G12.5,1X))
1001  FORMAT(/,25X,A)
      
      RETURN

 END SUBROUTINE PRTVEC

 ! *****************************************************************************
 ! 
 ! subroutine init_MPI()
 !   
 !   ! initialize MPI
 !   
 !   call MPI_INIT(mpierr)
 !   
 !   call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
 !   call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpierr )
 !
 !   if (myid .ne. 0)  root_proc = .false.
 !   
!!    call init_MPI_types()
 !   
 !   write (*,*) "MPI process ", myid, " of ", numprocs, " is alive"    
 !   write (*,*) "MPI process ", myid, ": root_proc=", root_proc
!
!  end subroutine init_MPI


END MODULE math_routines
