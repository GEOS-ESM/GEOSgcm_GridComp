program vtrilu_standalone

    use MAPL_ConstantsMod

    implicit none

    real, dimension(:,:,:), allocatable :: C
    real(kind=MAPL_R8), dimension(:,:,:), allocatable :: A, B
    
    real, dimension(:,:,:), allocatable :: A_inout, B_inout
    real, dimension(:,:,:), allocatable :: A_ref, B_ref

    character*100 :: dirName, rank_str
    integer :: IM, JM, LM, fileID
    real :: t1, t2

    dirName = './c90-90-90-72'
    rank_str = '0'

    ! C12 Data Set
    if (dirName == './c12-6-6-72') then
        IM = 6
        JM = 6
        LM = 72
    else if (dirName == './c90-90-90-72') then
        IM = 90
        JM = 90
        LM = 72
    else if(dirName == './c180-180-180-72') then
        IM = 180
        JM = 180
        LM = 72
    endif

    allocate(A(IM, JM, LM))
    allocate(B(IM, JM, LM))
    allocate(A_inout(IM, JM, LM))
    allocate(B_inout(IM, JM, LM))
    allocate(A_ref(IM, JM, LM))
    allocate(B_ref(IM, JM, LM))
    allocate(C(IM, JM, LM))

    ! *** Test AKS and BKS ***

    open(newunit=fileID, file=trim(dirName) // "/AKS_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) A_inout
    close(fileID)
    write(*,*) 'sum(abs(AKS)) input = ', sum(abs(A_inout))

    open(newunit=fileID, file=trim(dirName) // "/BKS_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) B_inout
    close(fileID)
    write(*,*) 'sum(abs(BKS)) input = ', sum(abs(B_inout))

    open(newunit=fileID, file=trim(dirName) // "/CKS_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) C
    close(fileID)
    write(*,*) 'sum(abs(CKS)) = ', sum(abs(C))



    A = A_inout
    B = B_inout
!$acc data copy(A, B, C)
    call cpu_time(t1)
    call VTRILU(A, B, C)
    call cpu_time(t2)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/AKS_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) A_ref
    close(fileID)
    write(*,*) 'sum(abs(AKS_ref)) input = ', sum(abs(A_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKS_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) B_ref
    close(fileID)
    write(*,*) 'sum(abs(BKS_ref)) input = ', sum(abs(B_ref))

    A_inout = A
    B_inout = B

    write(*,*) 'SUM DIFF AKS = ', sum(A_inout - A_REF)
    write(*,*) 'SUM DIFF BKS = ', sum(B_inout - B_REF)
    write(*,*) 'Execution Time = ', t2 - t1
    write(*,*) '**************************************'

    ! *** Test AKQ and BKQ ***

    open(newunit=fileID, file=trim(dirName) // "/AKQ_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) A_inout
    close(fileID)
    write(*,*) 'sum(abs(AKQ)) input = ', sum(abs(A_inout))

    open(newunit=fileID, file=trim(dirName) // "/BKQ_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) B_inout
    close(fileID)
    write(*,*) 'sum(abs(BKQ)) input = ', sum(abs(B_inout))

    open(newunit=fileID, file=trim(dirName) // "/CKQ_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) C
    close(fileID)
    write(*,*) 'sum(abs(CKQ)) = ', sum(abs(C))

    A = A_inout
    B = B_inout
!$acc data copy(A, B, C)
    call cpu_time(t1)
    call VTRILU(A, B, C)
    call cpu_time(t2)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/AKQ_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) A_ref
    close(fileID)
    write(*,*) 'sum(abs(AKQ_ref)) input = ', sum(abs(A_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKQ_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) B_ref
    close(fileID)
    write(*,*) 'sum(abs(BKQ_ref)) input = ', sum(abs(B_ref))

    A_inout = A
    B_inout = B

    write(*,*) 'SUM DIFF AKQ = ', sum(A_inout - A_REF)
    write(*,*) 'SUM DIFF BKQ = ', sum(B_inout - B_REF)
    write(*,*) 'Execution Time = ', t2 - t1
    write(*,*) '**************************************'

    ! *** Test AKV and BKV ***

    open(newunit=fileID, file=trim(dirName) // "/AKV_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) A_inout
    close(fileID)
    write(*,*) 'sum(abs(AKV)) input = ', sum(abs(A_inout))

    open(newunit=fileID, file=trim(dirName) // "/BKV_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) B_inout
    close(fileID)
    write(*,*) 'sum(abs(BKV)) input = ', sum(abs(B_inout))

    open(newunit=fileID, file=trim(dirName) // "/CKV_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) C
    close(fileID)
    write(*,*) 'sum(abs(CKV)) = ', sum(abs(C))

    A = A_inout
    B = B_inout
!$acc data copy(A, B, C)
    call cpu_time(t1)
    call VTRILU(A, B, C)
    call cpu_time(t2)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/AKV_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) A_ref
    close(fileID)
    write(*,*) 'sum(abs(AKV_ref)) input = ', sum(abs(A_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKV_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) B_ref
    close(fileID)
    write(*,*) 'sum(abs(BKV_ref)) input = ', sum(abs(B_ref))

    A_inout = A
    B_inout = B

    write(*,*) 'SUM DIFF AKV = ', sum(A_inout - A_REF)
    write(*,*) 'SUM DIFF BKV = ', sum(B_inout - B_REF)
    write(*,*) 'Execution Time = ', t2 - t1
    write(*,*) '**************************************'

    ! *** Test AKSS and BKSS ***

    open(newunit=fileID, file=trim(dirName) // "/AKSS_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) A_inout
    close(fileID)
    write(*,*) 'sum(abs(AKSS)) input = ', sum(abs(A_inout))

    open(newunit=fileID, file=trim(dirName) // "/BKSS_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) B_inout
    close(fileID)
    write(*,*) 'sum(abs(BKSS)) input = ', sum(abs(B_inout))

    open(newunit=fileID, file=trim(dirName) // "/CKSS_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) C
    close(fileID)
    write(*,*) 'sum(abs(CKSS)) = ', sum(abs(C))

    A = A_inout
    B = B_inout
!$acc data copy(A, B, C)
    call cpu_time(t1)
    call VTRILU(A, B, C)
    call cpu_time(t2)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/AKSS_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) A_ref
    close(fileID)
    write(*,*) 'sum(abs(AKSS_ref)) input = ', sum(abs(A_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKSS_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) B_ref
    close(fileID)
    write(*,*) 'sum(abs(BKSS_ref)) input = ', sum(abs(B_ref))

    A_inout = A
    B_inout = B

    write(*,*) 'SUM DIFF AKSS = ', sum(A_inout - A_REF)
    write(*,*) 'SUM DIFF BKSS = ', sum(B_inout - B_REF)
    write(*,*) 'Execution Time = ', t2 - t1
    write(*,*) '**************************************'

    ! *** Test AKQQ and BKQQ ***

    open(newunit=fileID, file=trim(dirName) // "/AKQQ_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) A_inout
    close(fileID)
    write(*,*) 'sum(abs(AKQQ)) input = ', sum(abs(A_inout))

    open(newunit=fileID, file=trim(dirName) // "/BKQQ_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) B_inout
    close(fileID)
    write(*,*) 'sum(abs(BKQQ)) input = ', sum(abs(B_inout))

    open(newunit=fileID, file=trim(dirName) // "/CKQQ_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) C
    close(fileID)
    write(*,*) 'sum(abs(CKQQ)) = ', sum(abs(C))

    A = A_inout
    B = B_inout
!$acc data copy(A, B, C)
    call cpu_time(t1)
    call VTRILU(A, B, C)
    call cpu_time(t2)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/AKQQ_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) A_ref
    close(fileID)
    write(*,*) 'sum(abs(AKQQ_ref)) input = ', sum(abs(A_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKQQ_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) B_ref
    close(fileID)
    write(*,*) 'sum(abs(BKQQ_ref)) input = ', sum(abs(B_ref))

    A_inout = A
    B_inout = B

    write(*,*) 'SUM DIFF AKQQ = ', sum(A_inout - A_REF)
    write(*,*) 'SUM DIFF BKQQ = ', sum(B_inout - B_REF)
    write(*,*) 'Execution Time = ', t2 - t1
    write(*,*) '**************************************'

    ! *** Test AKUU and BKUU ***

    open(newunit=fileID, file=trim(dirName) // "/AKUU_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) A_inout
    close(fileID)
    write(*,*) 'sum(abs(AKUU)) input = ', sum(abs(A_inout))

    open(newunit=fileID, file=trim(dirName) // "/BKUU_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) B_inout
    close(fileID)
    write(*,*) 'sum(abs(BKUU)) input = ', sum(abs(B_inout))

    open(newunit=fileID, file=trim(dirName) // "/CKUU_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) C
    close(fileID)
    write(*,*) 'sum(abs(CKUU)) = ', sum(abs(C))

    A = A_inout
    B = B_inout
!$acc data copy(A, B, C)
    call cpu_time(t1)
    call VTRILU(A, B, C)
    call cpu_time(t2)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/AKUU_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) A_ref
    close(fileID)
    write(*,*) 'sum(abs(AKUU_ref)) input = ', sum(abs(A_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKUU_" // trim(rank_str) // ".out", &
    status='old', form="unformatted", action="read")
    read(fileID) B_ref
    close(fileID)
    write(*,*) 'sum(abs(BKUU_ref)) input = ', sum(abs(B_ref))

    A_inout = A
    B_inout = B

    write(*,*) 'SUM DIFF AKUU = ', sum(A_inout - A_REF)
    write(*,*) 'SUM DIFF BKUU = ', sum(B_inout - B_REF)
    write(*,*) 'Execution Time = ', t2 - t1
    write(*,*) '**************************************'

    contains
    ! !IROUTINE:  VTRILU --  Does LU decomposition of tridiagonal matrix.

    ! !INTERFACE:

    subroutine VTRILU(A,B,C)

        ! !ARGUMENTS:
        
        real,               dimension(:,:,:), intent(IN   ) :: C
        real(kind=MAPL_R8), dimension(:,:,:), intent(INOUT) :: A, B
        
        ! !DESCRIPTION: {\tt VTRILU} performs an $LU$ decomposition on
        ! a tridiagonal matrix $M=LU$.
        !
        ! $$
        ! M = \left( \begin{array}{ccccccc}
        !      b_1 & c_1 & & & & & \\
        !      a_2 & b_2 & c_2 & & & &  \\
        !      &  \cdot& \cdot & \cdot & & &  \\
        !      & & \cdot& \cdot & \cdot & &  \\
        !      &&  & \cdot& \cdot & \cdot &  \\
        !      &&&& a_{K-1} & b_{K-1} & c_{K-1}   \\
        !      &&&&& a_{K} & b_{K}
        !    \end{array} \right)
        ! $$
        !
        !
        ! $$
        ! \begin{array}{lr}
        ! L = \left( \begin{array}{ccccccc}
        !      1 &&&&&& \\
        !      \hat{a}_2 & 1 & &&&&  \\
        !      &  \cdot& \cdot &  & & &  \\
        !      & & \cdot& \cdot &  &&  \\
        !      &&  & \cdot& \cdot &  &  \\
        !      &&&& \hat{a}_{K-1} & 1 &   \\
        !      &&&&& \hat{a}_{K} & 1
        !    \end{array} \right)
        ! &
        ! U = \left( \begin{array}{ccccccc}
        !      \hat{b}_1 & c_1 &&&&& \\
        !       & \hat{b}_2 & c_2 &&&&  \\
        !      &  & \cdot & \cdot & & &  \\
        !      & & & \cdot & \cdot &&  \\
        !      &&  & & \cdot & \cdot &  \\
        !      &&&&  & \hat{b}_{K-1} & c_{K-1}   \\
        !      &&&&&  & \hat{b}_{K}
        !    \end{array} \right)
        ! \end{array}
        ! $$
        !
        !
        ! On input, A, B, and C contain, $a_k$, $b_k$, and $c_k$
        ! the lower, main, and upper diagonals of the matrix, respectively.
        ! On output, B contains $1/\hat{b}_k$, the inverse of the main diagonal of $U$,
        ! and A contains $\hat{a}_k$,
        ! the lower diagonal of $L$. C contains the upper diagonal of the original matrix and of $U$.
        !
        ! The new diagonals $\hat{a}_k$ and $\hat{b}_k$ are:
        ! $$
        ! \begin{array}{rcl}
        ! \hat{b}_1 & = & b_1, \\
        ! \hat{a}_k & = & \makebox[2 in][l]{$a_k / \hat{b}_{k-1}$,}  k=2, K, \\
        ! \hat{b}_k & = & \makebox[2 in][l]{$b_k - c_{k-1} \hat{a}_k$,} k=2, K. 
        ! \end{array}
        ! $$
        !EOP
        
        integer :: IM, JM, LM, I, J, L

        LM = size(C,3)
        JM = size(C,2)
        IM = size(C,1)
!$acc data present(A, B, C)
        ! B(:,:,1) = 1. / B(:,:,1)
!$acc parallel loop gang vector collapse(2)
        do J = 1, JM
            do I = 1, IM
                B(I,J,1) = 1. / B(I, J, 1)
            enddo
        enddo
!$acc end parallel loop

        ! do L = 2,LM
        !     A(:,:,L) = A(:,:,L) * B(:,:,L-1)
        !     B(:,:,L) = 1. / ( B(:,:,L) - C(:,:,L-1) * A(:,:,L) )
        ! end do
!$acc parallel
!$acc loop seq
        do L = 2, LM
!$acc loop gang vector collapse(2)
            do J = 1, JM
                do I = 1, IM
                    A(I,J,L) = A(I,J,L) * B(I,J,L-1)
                    B(I,J,L) = 1. / ( B(I,J,L) - C(I,J,L-1) * A(I,J,L) )
                enddo
            enddo
        enddo
!$acc end parallel
!$acc end data
    end subroutine VTRILU

end program

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”

! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at

! http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.