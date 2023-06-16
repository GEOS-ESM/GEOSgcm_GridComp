module Fillq2zeroMod

  implicit none

  private
  public :: fillq2zero

contains

  subroutine FILLQ2ZERO(Q, MASS, FILLQ)

    ! New algorithm to fill the negative q values in a mass conserving way.
    ! Conservation of TPW was checked. Donifan Barahona
    ! Updated from FILLQ2ZERO, avoids the usage of scalars

    real, dimension(:,:,:),   intent(inout)  :: Q
    real, dimension(:,:,:),   intent(in)     :: MASS
    real, dimension(:,:),     intent(  out)  :: FILLQ
    real                                     :: TPW1, TPW2, TPWC
    integer                                  :: IM,JM,LM, l, i, j

    IM = SIZE( Q, 1 )
    JM = SIZE( Q, 2 )
    LM = SIZE( Q, 3 )

!$acc parallel loop gang vector collapse(2) private(TPW2, TPWC,TPW1)
    do j = 1,JM
       do i = 1, IM
          TPW2 = 0.0
          TPWC = 0.0
          TPW1 = 0.0
          do l = 1,LM
             TPW1 = TPW1 + Q(I,J,L)*MASS(I,J,L)
             IF(Q(I,J,L) < 0.0) Q(I,J,L) = 0.0
          enddo

          do l = 1,LM
             TPW2 = TPW2 + Q(I,J,L)*MASS(I,J,L)
          enddo
          if (TPW2 > 0.0) TPWC = (TPW2-TPW1)/TPW2

          do l = 1,LM
             Q(I,J,L) = Q(I,J,L) * (1.0 - TPWC)
          enddo
          FILLQ(I,J) = TPW2-TPW1
       enddo
    enddo
!$acc end parallel loop

  end subroutine FILLQ2ZERO

end module Fillq2zeroMod
