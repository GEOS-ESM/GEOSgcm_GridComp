!  $Id$
!
!-----------------------------
module ppm

IMPLICIT NONE
PRIVATE

PUBLIC UPSTAD3
PUBLIC UPSTAD1
PUBLIC UPVERT1

CONTAINS

SUBROUTINE upvert1( U , Q , X, DT, I1 , IN  )

IMPLICIT NONE

INTEGER, INTENT(IN) :: I1, IN
REAL, INTENT(IN)    :: U(I1:IN+1)
REAL, INTENT(IN)    :: X(I1:IN)
REAL, INTENT(IN)    :: DT
REAL, INTENT(INOUT) :: Q(I1:IN)

REAL :: qnew(I1:IN)
integer :: I
REAL :: QQ(I1-1:IN+1), XX(I1-1:IN+1)
real :: damp

qq(I1:IN) = q(I1:IN)
xx(I1:IN) = x(I1:IN)

qq(i1-1) = q(i1) + (qq(i1+1)-qq(i1))
xx(i1-1) = xx(i1) - (xx(i1+1)-xx(i1))

qq(iN+1) = q(iN)
xx(iN+1) = xx(iN) +   (xx(iN)-xx(iN-1))

do i=i1,iN
   damp = min ( 1., max ( 1.e-3, 0.0111 * xx(i) - 0.11 ) )
   if ( u(i) <= 0. ) then
      qnew(i) = qq(i) - ( u(i+1) * (qq(i+1)-qq(i)) ) * dt*damp/(xx(i+1)-xx(i))
   else
      qnew(i) = qq(i) - ( u(i)   * (qq(i-1)-qq(i)) ) * dt*damp/(xx(i-1)-xx(i))
   endif
end do 

do i=i1,iN
   q(i) = qnew(i)
end do

END SUBROUTINE upvert1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE upstad3( U , Q , X, DT, I1 , IN  )

IMPLICIT NONE

INTEGER, INTENT(IN)         :: I1, IN

REAL, INTENT(IN)            :: U(I1:IN+1)
REAL, INTENT(IN)            :: X(I1:IN)
REAL, INTENT(IN)            :: DT
REAL, INTENT(INOUT)         :: Q(I1:IN)

REAL :: QQ(I1-2:IN+2), XX(I1-2:IN+2)

REAL :: XDP, QDP  ,QDP1  ,QDP2  ,QDP3  ,QDP4, sumq 

integer :: I

qq(I1:IN) = q(I1:IN)
xx(I1:IN) = x(I1:IN)

qq(i1-2:i1-1) = q(i1)  ! q(iN-1:iN)
xx(i1-1) = xx(i1) -   (xx(i1+1)-xx(i1))
xx(i1-2) = xx(i1) - 2*(xx(i1+1)-xx(i1))

qq(iN+1:iN+2) = q(iN)  ! q(i1:i1+1)
xx(iN+1) = xx(iN) +   (xx(iN)-xx(iN-1))
xx(iN+2) = xx(iN) + 2*(xx(iN)-xx(iN-1))

do i=i1,iN

   xdp = xx(i) - u(i)*dt

   if ( u(i) >= 0. )  then 

      qdp  =                                                                              &
               qq(i) * ( xdp - xx(i-2) )*( xdp - xx(i-1) )*( xdp - xx(i+1) )              &
                  /  ( ( xx(i) - xx(i-2) )*( xx(i) - xx(i-1) )*( xx(i) - xx(i+1) ) )      &

            +  qq(i-1) * ( xdp - xx(i-2) )*( xdp - xx(i) )*( xdp - xx(i+1) )              &
                  /  ( ( xx(i-1) - xx(i-2) )*( xx(i-1) - xx(i) )*( xx(i-1) - xx(i+1) ) )  &

            +  qq(i+1) * ( xdp - xx(i-2) )*( xdp - xx(i) )*( xdp - xx(i-1) )              &
                  /  ( ( xx(i+1) - xx(i-2) )*( xx(i+1) - xx(i) )*( xx(i+1) - xx(i-1) ) )  &

            +  qq(i-2) * ( xdp - xx(i-1) )*( xdp - xx(i) )*( xdp - xx(i+1) )              &
                  /  ( ( xx(i-2) - xx(i-1) )*( xx(i-2) - xx(i) )*( xx(i-2) - xx(i+1) ) )     

      qdp1  = qq(i) * ( xdp - xx(i-2) )*( xdp - xx(i-1) )*( xdp - xx(i+1) )              &
                  /  ( ( xx(i) - xx(i-2) )*( xx(i) - xx(i-1) )*( xx(i) - xx(i+1) ) ) 
      qdp2  = qq(i-1) * ( xdp - xx(i-2) )*( xdp - xx(i) )*( xdp - xx(i+1) )              &
                  /  ( ( xx(i-1) - xx(i-2) )*( xx(i-1) - xx(i) )*( xx(i-1) - xx(i+1) ) ) 
      qdp3  = qq(i+1) * ( xdp - xx(i-2) )*( xdp - xx(i) )*( xdp - xx(i-1) )              &
                  /  ( ( xx(i+1) - xx(i-2) )*( xx(i+1) - xx(i) )*( xx(i+1) - xx(i-1) ) )
      qdp4  = qq(i-2) * ( xdp - xx(i-1) )*( xdp - xx(i) )*( xdp - xx(i+1) )              &
                  /  ( ( xx(i-2) - xx(i-1) )*( xx(i-2) - xx(i) )*( xx(i-2) - xx(i+1) ) )     
      sumq=qdp1+qdp2+qdp3+qdp4

   else 

      qdp1  =                                                                             &
               qq(i) * ( xdp - xx(i+2) )*( xdp - xx(i-1) )*( xdp - xx(i+1) )              &
                  / ( ( xx(i) - xx(i+2) )*( xx(i) - xx(i-1) )*( xx(i) - xx(i+1) ) )         ! &

      qdp2 =0.       +  qq(i-1) * ( xdp - xx(i+2) )*( xdp - xx(i) )*( xdp - xx(i+1) )     &
                  / ( ( xx(i-1) - xx(i+2) )*( xx(i-1) - xx(i) )*( xx(i-1) - xx(i+1) ) )     ! &

      qdp3 =0.          +  qq(i+1) * ( xdp - xx(i+2) )*( xdp - xx(i) )*( xdp - xx(i-1) )  &
                  / ( ( xx(i+1) - xx(i+2) )*( xx(i+1) - xx(i) )*( xx(i+1) - xx(i-1) ) )     ! &

      qdp4 =0.          +  qq(i+2) * ( xdp - xx(i-1) )*( xdp - xx(i) )*( xdp - xx(i+1) )  &
                  / ( ( xx(i+2) - xx(i-1) )*( xx(i+2) - xx(i) )*( xx(i+2) - xx(i+1) ) )       


      qdp=qdp1+qdp2+qdp3+qdp4


   endif

   q(i) = qdp     

end do 

END SUBROUTINE upstad3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE upstad1( U , Q , X, DT, I1 , IN  )

IMPLICIT NONE

INTEGER, INTENT(IN)         :: I1, IN

REAL, INTENT(IN)            :: U(I1:IN+1)
REAL, INTENT(IN)            :: X(I1:IN)
REAL, INTENT(IN)            :: DT
REAL, INTENT(INOUT)         :: Q(I1:IN)

REAL :: QQ(I1-1:IN+1), XX(I1-1:IN+1)

REAL :: XDP, QDP  ,QDP1  ,QDP2  ,QDP3  ,QDP4 

integer :: I

qq(I1:IN) = q(I1:IN)
xx(I1:IN) = x(I1:IN)

qq(i1-1) = q(i1)  ! q(iN-1:iN)
xx(i1-1) = xx(i1) - (xx(i1+1)-xx(i1))

qq(iN+1) = q(iN)  ! q(i1:i1+1)
xx(iN+1) = xx(iN) + (xx(iN)-xx(iN-1))

do i=i1,iN

   xdp = xx(i) - u(i)*dt

   if ( u(i) >= 0. )  then 

       qdp  =  qq(i) * ( xdp - xx(i-1) ) / ( xx(i) - xx(i-1) )   &
          +  qq(i-1) * ( xdp - xx(i)   ) / ( xx(i-1) - xx(i) )      

     else 

       qdp  =  qq(i) * ( xdp - xx(i+1) ) / ( xx(i) - xx(i+1) )         &
          +  qq(i+1) * ( xdp - xx(i)   ) / ( xx(i+1) - xx(i) )    

      endif

      q(i) = qdp     

end do 

END SUBROUTINE upstad1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ppm
