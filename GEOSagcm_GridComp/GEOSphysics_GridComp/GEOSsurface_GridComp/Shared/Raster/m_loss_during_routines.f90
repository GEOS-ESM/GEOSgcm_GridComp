MODULE loss_during_routines

  implicit none

  private
  
  public tridag, GMRES, DDOTL, DAXPYL, ope
  
  contains 
    
  SUBROUTINE tridag(a,b,c,r,u,n)
    IMPLICIT NONE
    integer, intent(in)  :: n
    real a(n),b(n),c(n),r(n),u(n)
    real :: gam(N)
    integer :: j
    real :: bet
!     if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
         if(bet.eq.0.)pause 'tridag failed'
         u(j)=(r(j)-a(j)*u(j-1))/bet
      end do
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do
      return
    END SUBROUTINE tridag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE  GMRES (N, IM, RHS, SOL, SS, EPS, MAXITS, IOUT, &
           c1,c2,c3)

!C*************************************************************
!c
!c     This is gmrd.f (double precision), the original Saad version
!c     of GMRES, with trivial changes in output format.
!c
!C GMRES ALGORITHM . SIMPLE VERSION .  (MAY 23, 1985)
!C PARAMETER LIST:
!C N     == SIZE OF PROBLEM
!C IM    == SIZE OF KRYLOV SUBSPACE:  SHOULD NOT EXCEED 50 IN THIS
!C          VERSION (CAN BE RESET IN CODE. LOOKING AT COMMENT BELOW)
!C RHS   == RIGHT HAND SIDE
!C SOL   == INITIAL GUESS ON INPUT, APPROXIMATE SOLUTION ON OUTPUT
!C SS    == WORK SPACE OF SIZE N X (IM+1)
!C EPS   == TOLERANCE FOR STOPPING CRITERION. PROCESS IS STOPPED
!C          AS SOON AS ( ||.|| IS THE EUCLIDEAN NORM):
!C          || CURRENT RESIDUAL||/||INITIAL RESIDUAL|| <= EPS
!C          ON OUTPUT, EPS = ::FINAL RESIDUAL::/::INITIAL RESIDUAL::
!C MAXITS== MAXIMUM NUMBER OF ITERATIONS ALLOWED ON INPUT
!C          ON OUTPUT, MAXITS = TOTAL NUMBER OF ITERATIONS
!C IOUT  == OUTPUT UNIT NUMBER NUMBER FOR PRINTING INTERMEDIATE RESULTS
!C          IF (IOUT .LE. 0) NO STATISTICS ARE PRINTED.
!C ----------------------------------------------------------------
!C SUBROUTINES USED =
!C OPE(N,X,Y)  ==  MATRIX BY VECTOR MULTIPLICATION DELIVERS Y=AX, GIVEN X.
!C DDOTL       == DOT PRODUCT FUNCTION.  REPLACE BY BLAS ROUTINE DDOT
!C DAXPYL      == Y <-- Y+AX  ROUTINE. REPLACE BY BLAS ROUTINE DAXPY
!C*************************************************************
!c      IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
      integer, intent (in) ::  n, im, maxits, iout
      real, intent(in) :: c1(n),c2(n),c3(n)
      integer i, j, k, k1, ii, i1, n1, its
      REAL SS(N,IM+1), RHS(N), SOL(N)
      real HH(51,50), C(50), S(50), RS(51)
      real eps, eps1, gam, mcheps, ro, t

!C-------------------------------------------------------------
!C ARNOLDI SIZE SHOULD NOT EXCEED 50 IN THIS VERSION..
!C TO RESET MODIFY SIZES OF HH, C, S, RS       ----------------
!C-------------------------------------------------------------
      mcheps = 1.d-30
      N1 = N + 1
      ITS = 0
!C-------------------------------------------------------------
!C **  OUTER LOOP STARTS HERE..
!C-------------- COMPUTE INITIAL RESIDUAL VECTOR --------------
 10     CONTINUE
       CALL ope (N, SOL, SS,c1,c2,c3)

        DO J=1,N
           SS(J,1) = RHS(J) - SS(J,1)
        ENDDO

!C-------------------------------------------------------------

         RO = SQRT( DDOTL(N, SS,SS) )
         IF (RO .EQ. 0.0D0) RETURN
         DO J=1, N
            SS(J,1) = SS(J,1)/RO
         ENDDO
         IF (ITS .EQ. 0) EPS1=EPS*RO

!C ** INITIALIZE 1-ST TERM  OF RHS OF HESSENBERG SYSTEM..
         RS(1) = RO
         I = 0
 4       I=I+1
         ITS = ITS + 1
         I1 = I + 1
         CALL ope (N, SS(1,I), SS(1,I1),c1,c2,c3)

!C-----------------------------------------
!C  MODIFIED GRAM - SCHMIDT...
!C-----------------------------------------
        DO J=1, I
           T = DDOTL(N, SS(1,J),SS(1,I1))
           HH(J,I) = T
           CALL DAXPYL(N, -T, SS(1,J), SS(1,I1))
        ENDDO
        T = SQRT(DDOTL(N, SS(1,I1), SS(1,I1)))
        HH(I1,I) = T
        DO K=1,N
           SS(K,I1) = SS(K,I1) / T
        ENDDO

!C--------DONE WITH MODIFIED GRAM SCHIMD AND ARNOLDI STEP..
!C NOW  UPDATE FACTORIZATION OF HH
!C---------------------------------------------------------
        IF (I .EQ. 1) GOTO 121
!C-------- PERFROM PREVIOUS TRANSFORMATIONS  ON I-TH COLUMN OF H
        DO K=2,I
           K1 = K-1
           T = HH(K1,I)
           HH(K1,I) = C(K1)*T + S(K1)*HH(K,I)
           HH(K,I) = -S(K1)*T + C(K1)*HH(K,I)
        ENDDO
 121    GAM = SQRT(HH(I,I)**2 + HH(I1,I)**2)
        IF (GAM .EQ. 0.0D0) GAM = MCHEPS
!C-----------#  DETERMINE NEXT PLANE ROTATION  #-------------------
        C(I) = HH(I,I)/GAM
        S(I) = HH(I1,I)/GAM
        RS(I1) = -S(I)*RS(I)
        RS(I) =  C(I)*RS(I)
!C---DETERMINE RESIDUAL NORM AND TEST FOR CONVERGENCE-
        HH(I,I) = C(I)*HH(I,I) + S(I)*HH(I1,I)
        RO = ABS(RS(I1))
!c        IF (IOUT .GT. 0)
!c     *          WRITE(IOUT, 199) ITS, RO
        IF ( (I .LT. IM) .AND. (RO .GT. EPS1) .and. &
            (its .lt. maxits) )  GOTO 4
!C
!C NOW COMPUTE SOLUTION. FIRST SOLVE UPPER TRIANGULAR SYSTEM.
!C
        RS(I) = RS(I)/HH(I,I)
        DO II=2,I
           K=I-II+1
           K1 = K+1
           T=RS(K)
           DO J=K1,I
              T = T-HH(K,J)*RS(J)
           ENDDO
           RS(K) = T/HH(K,K)
        ENDDO
!C DONE WITH BACK SUBSTITUTION..
!C NOW FORM LINEAR COMBINATION TO GET SOLUTION
        DO J=1, I
           T = RS(J)
           CALL DAXPYL(N, T, SS(1,J), SOL)
        ENDDO
!C RESTART OUTER LOOP  WHEN NECESSARY
        IF (RO .GT. EPS1 .AND. ITS .LT. MAXITS) GOTO 10
 199    format(1x, i5, 3x, 1pe12.4 )
      EPS = (RO / EPS1) * EPS
!crdk      MAXITS = ITS
!cjcm
!crdk      WRITE(IOUT, 199) ITS, RO
!cjcm
      RETURN
!C------------------------------- END OF GMRES ----------------------
    END SUBROUTINE GMRES

!C------------------------------- BEG OF DDOTL ----------------------
    FUNCTION DDOTL(N,RVA,RVB)
      IMPLICIT NONE
      INTEGER I,N
      REAL RVA(N),RVB(N),DDOTL
      DDOTL = 0.0
      DO I=1,N
         DDOTL=DDOTL+RVA(I)*RVB(I)
      ENDDO
      RETURN
    END FUNCTION DDOTL
!C------------------------------- END OF DDOTL ----------------------
!C------------------------------- BEG OF DAXPYL ----------------------
      SUBROUTINE DAXPYL(N,ALPHA,X,Y)
      IMPLICIT NONE
      INTEGER I,N
      REAL X(N),Y(N),ALPHA
      DO I=1,N
         Y(I)=Y(I)+ALPHA*X(I)
      ENDDO     
      RETURN
    END SUBROUTINE DAXPYL

!C------------------------------- END OF DAXPYL ----------------------

    subroutine ope(n,sol,ss,c1,c2,c3)
      
      IMPLICIT NONE
      integer, intent(in) :: n
      real, intent (out) :: ss(n)
      real,intent(in) :: c1(n),c2(n),c3(n),sol(n)
      integer :: nx
      
      do nx=2,n-1
         ss(nx)=c1(nx)*sol(nx+1)+c2(nx)*sol(nx)+c3(nx)*sol(nx-1)
      enddo
      
      ss(1)=c1(1)*sol(2)+c2(1)*sol(1)
      ss(n)=c2(n)*sol(n)+c3(n)*sol(n-1)
      
      return
    end subroutine ope

  END MODULE loss_during_routines
