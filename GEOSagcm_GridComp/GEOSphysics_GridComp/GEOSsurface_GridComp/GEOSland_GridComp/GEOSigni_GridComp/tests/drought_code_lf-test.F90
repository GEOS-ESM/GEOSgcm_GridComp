!
! Test the implementation of the monthly day length adjustment 
! factors for Drought Code (Lf) in the CFFWI module.
!
!
! Anton Darmenov, NASA, 2023
!


!
! Buld executable: 
!    $ gfortran -g -Wall -o cffwi_dc_lf-test cffwi.F90 drought_code_lf-test.F90
!

program cffwi_dc_le_test

    use cffwi, only: drought_code_lf, &
                     CFFWI_REFERENCE_LATITUDE


    implicit none

    ! local
    integer :: month
    

    real, dimension(12), parameter :: &
        DAY_LENGTH_FACTOR_NH = (/ -1.6, -1.6, -1.6, 0.9,  3.8,  5.8, &
                                   6.4,  5.0,  2.4, 0.4, -1.6, -1.6 /)

    real, dimension(12), parameter :: &
        DAY_LENGTH_FACTOR_SH = (/  6.4,  5.0,  2.4, 0.4, -1.6, -1.6, &
                                  -1.6, -1.6, -1.6, 0.9,  3.8,  5.8 /)


    write (*, '(A)') "Testing the monthly day length adjustment factors for Drought Code (Lf):"
    write (*, '(A)') ""

    write (*, '(A)') "Northern hemisphere:"
    do month = 1, 12
        write (*, '(A, I2, F6.2, F6.2)') "  month = ", month, &
            DAY_LENGTH_FACTOR_NH(month), &
            drought_code_lf(CFFWI_REFERENCE_LATITUDE, month)
    end do

    write (*, '(A)') ""

    write (*, '(A)') "Southern hemisphere:"
    do month = 1, 12
        write (*, '(A, I2, F6.2, F6.2)') "  month = ", month, &
            DAY_LENGTH_FACTOR_SH(month), &
            drought_code_lf(-CFFWI_REFERENCE_LATITUDE, month)
    end do

end program
