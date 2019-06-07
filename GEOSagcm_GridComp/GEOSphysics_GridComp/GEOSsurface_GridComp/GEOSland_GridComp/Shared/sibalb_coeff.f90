      module sibalb_coeff
      implicit none

      contains

      FUNCTION COEFFSIB(TABLE, NTABL, LAI ,DX, DY)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTABL, LAI

      REAL, INTENT(IN) :: DX, DY
      REAL, INTENT(IN), DIMENSION(NTABL,2) :: TABLE
      REAL COEFFSIB

      COEFFSIB = (TABLE(LAI,  1)                                                  &
             + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX)             &
             + (TABLE(LAI+1,1)                                                 &
             + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX

      RETURN
      END FUNCTION COEFFSIB     



    end module sibalb_coeff
