module Wnet

use omp_lib 

  implicit none
  private
 
  integer, parameter :: NF = 14  ! Number of features
  integer, parameter :: NL = 5   ! Number of hidden layers
  integer, parameter :: NODES = 128

  real, allocatable :: W_hidden(:, :, :) ! Shape: [NODES, NODES, NL]
  real, allocatable :: b_hidden(:, :)    ! Shape: [NODES, NL]
  real, allocatable :: W_input(:, :)     ! Shape: [NF, NODES]
  real, allocatable :: b_input(:)        ! Shape: [NODES]
  real, allocatable :: W_output(:)       ! Shape: [NODES]
  real :: b_output

  ! Module-level mean and standard deviation
  real :: mean(NF) = (/ &
      243.9, 0.6, 6.3, 0.013, 0.0002, 5.04, 21.8, 0.002, 9.75e-7, 7.87e-6, 0.6, 5.04, 21.8, 0.002/)
  real :: stddev(NF) = (/30.3, 0.42, 16.1, 7.9, 0.05, 20.6, 20.8, 0.0036, 7.09e-6, 2.7e-5, 0.42, 20.6, 20.8, 0.0036/)
  
  public :: load_wnet_weights, Wnet_forward_pass    

!========== Parameterization of SigmaW ====================
!This module calculates the subgrid scale standard deviation in vertical velocity using the Wnet model according to:
! Barahona, D., Breen, K. H., Kalesse-Los, H., & Röttenbacher, J. (2024). 
!Deep Learning Parameterization of Vertical Wind Velocity Variability via Constrained Adversarial Training. 
!Artificial Intelligence for the Earth Systems, 3(1), e230025. doi: 10.1175/AIES-D-23-0025.1
!
!Wnet is called using run_Wnet(input, output, Ns) where Ns is the number of samples.
!The input to Wnet consists of a 14-dimensional vector [Ns, 14] including the Richardson number (Ri, dimensionless), total scalar
!diffusivity for momentum (KM, in m2 s-1), the 3-dimensional wind velocity (U, V and
!W, in m s-1), the water vapor, liquid and ice mass mixing ratios (QV, QL and QI
!in Kg Kg-1), air density (ÏAIRD, in Kg m-3), and temperature (T in K). For each sample they must be in the order:
!['T', 'AIRD', 'U', 'V', 'W', 'KM', 'RI', 'QV', 'QI', 'QL', AIRD_sfc, 'KM_sfc', 'RI_sfc', 'QV_sfc']  
!The subscript _sfc indicates surface values.
!
!Output (shape [Ns]) is the standard deviation in W (m s-1)
! 
!This version has been optimized for parallel performance (use the -fopenmp -O3 flags).
!
!Developed by Donifan Barahona donifan.o.barahona@nasa.gov

contains

!========================================
!weight loading
!========================================
subroutine load_wnet_weights()
    implicit none
    integer :: i, j, k
    character(len=10) :: header
    logical :: file_exists

    inquire(file='Wnet_weights.txt', exist=file_exists)
    if (.not. file_exists) then
        print *, "Error: File 'Wnet_weights.txt' not found!"
        stop
    end if

    open(10, file='Wnet_weights.txt', status='old')

    if (.not. allocated(W_input)) allocate(W_input(NF, NODES))
    if (.not. allocated(b_input)) allocate(b_input(NODES))
    if (.not. allocated(W_hidden)) allocate(W_hidden(NODES, NODES, NL))
    if (.not. allocated(b_hidden)) allocate(b_hidden(NODES, NL))
    if (.not. allocated(W_output)) allocate(W_output(NODES))

    read(10, *) header, i, j
    read(10, *) ((W_input(i, j), j = 1, NODES), i = 1, NF)
    read(10, *) header, i
    read(10, *) (b_input(i), i = 1, NODES)

    do k = 1, NL
        read(10, *) header, i, j
        read(10, *) ((W_hidden(i, j, k), j = 1, NODES), i = 1, NODES)
        read(10, *) header, i
        read(10, *) (b_hidden(i, k), i = 1, NODES)
    end do

    read(10, *) header, i
    read(10, *) (W_output(i), i = 1, NODES)
    read(10, *) header, i
    read(10, *) b_output

    close(10)
end subroutine load_wnet_weights

!========================================
! Leaky ReLU 
!========================================
function leaky_relu(x, alpha) result(y)
    real, intent(in) :: x, alpha
    real :: y
    y = max(x, alpha * x) !more efficient than if
end function leaky_relu

!========================================
! Standardize Input Data
!========================================
subroutine standardize_input(input, Ns)
    integer, intent(in) :: Ns
    real, intent(inout) :: input(Ns, NF)

    integer :: i, j
    !$OMP PARALLEL DO PRIVATE(i, j) SCHEDULE(static)
    do i = 1, Ns
        do j = 1, NF
            if (stddev(j) /= 0.0) then
                input(i, j) = (input(i, j) - mean(j)) / stddev(j)
            else
                input(i, j) = 0.0
            end if
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine standardize_input

!========================================
! Forward Pass
!========================================
subroutine Wnet_forward_pass(input, output, Ns)
    implicit none

    integer, intent(in) :: Ns
    real, intent(in) :: input(Ns, NF)      ! Uses default REAL (typically single precision)
    real, intent(out) :: output(Ns)

    integer :: i, j, k
    real, allocatable, save :: layer_output(:,:), temp_output(:,:)
    real :: alpha

    ! Use BLAS SGEMM (Single Precision)
    external SGEMM
    real :: one, zero
    parameter (one = 1.0, zero = 0.0)  ! No need for `_sp`, follows default precision

    ! Allocate once (performance boost)
    if (.not. allocated(layer_output)) allocate(layer_output(Ns, NODES))
    if (.not. allocated(temp_output)) allocate(temp_output(Ns, NODES))

    ! First layer computation using SGEMM
    call SGEMM('N', 'N', Ns, NODES, NF, one, input, Ns, W_input, NF, zero, temp_output, Ns)

    alpha = 0.2
    !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(static)
    do i = 1, Ns
        do j = 1, NODES
            layer_output(i, j) = leaky_relu(temp_output(i, j) + b_input(j), alpha)
        end do
    end do
    !$OMP END PARALLEL DO

    ! Hidden layers (looping over NL)
    do k = 1, NL
        ! Perform matrix multiplication using SGEMM
        call SGEMM('N', 'N', Ns, NODES, NODES, one, layer_output, Ns, W_hidden(:, :, k), NODES, zero, temp_output, Ns)

        if (k == NL) alpha = 0.1

        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(static)
        do i = 1, Ns
            do j = 1, NODES
                layer_output(i, j) = leaky_relu(temp_output(i, j) + b_hidden(j, k), alpha)
            end do
        end do
        !$OMP END PARALLEL DO
    end do

    ! Final layer (dot product, parallelized)
    !$OMP PARALLEL DO SCHEDULE(static)
    do i = 1, Ns
        output(i) = dot_product(layer_output(i, :), W_output) + b_output
    end do
    !$OMP END PARALLEL DO

end subroutine Wnet_forward_pass


!!!!!!!!!!!!!!!!!
subroutine run_Wnet(input, output, Ns)
    integer, intent(in) :: Ns
    real, intent(inout) :: input(Ns, NF)
    real, intent(out) :: output(Ns)

    call load_wnet_weights()
    call standardize_input(input, Ns)
    call Wnet_forward_pass(input, output, Ns)
end subroutine run_Wnet

end module Wnet
