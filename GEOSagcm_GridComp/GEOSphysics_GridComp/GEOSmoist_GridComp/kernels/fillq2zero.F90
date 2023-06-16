module Fillq2zeroMod

  implicit none

  private
  public :: fillq2zero

contains

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

end module Fillq2zeroMod
