program beljaars_standalone

    implicit none

    character*100 :: dirName, rank_str
    integer :: IM, JM, LM, fileID
    real    :: DT, LAMBDA_B, C_B
    real, dimension(:,:,:), allocatable :: U, V, Z, PLE, BKV, BKUU, FKV
    real, dimension(:,:,:), allocatable :: BKV_ref, BKUU_ref, FKV_ref
    real, dimension(:,:),   allocatable :: VARFLT

    real :: ts, te 

    dirName = './c90-90-90-72'
    rank_str = '0'

    ! C12 Data Set
    if (dirName == './c12-6-6-72') then
        IM = 6
        JM = 6
        LM = 72
        DT = 900.0
    else if (dirName == './c90-90-90-72') then
        IM = 90
        JM = 90
        LM = 72
        DT = 450.0
    else if(dirName == './c180-180-180-72') then
        IM = 180
        JM = 180
        LM = 72
        DT = 450.0
    endif

    LAMBDA_B = 1500.0
    C_B      = 6.000e-7

    allocate(U       (IM, JM, LM))
    allocate(V       (IM, JM, LM))
    allocate(Z       (IM, JM, LM))
    allocate(VARFLT  (IM, JM))
    allocate(PLE     (IM, JM, 0:LM))
    allocate(BKV     (IM, JM, LM))
    allocate(BKUU    (IM, JM, LM))
    allocate(FKV     (IM, JM, LM))
    allocate(BKV_ref (IM, JM, LM))
    allocate(BKUU_ref(IM, JM, LM))
    allocate(FKV_ref (IM, JM, LM))

    open(newunit=fileID, file=trim(dirName) // "/U_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) U
    close(fileID)
    write(*,*) 'sum(abs(U)) = ', sum(abs(U))

    open(newunit=fileID, file=trim(dirName) // "/V_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) V
    close(fileID)
    write(*,*) 'sum(abs(V)) = ', sum(abs(V))

    open(newunit=fileID, file=trim(dirName) // "/Z_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) Z
    close(fileID)
    write(*,*) 'sum(abs(Z)) = ', sum(abs(Z))

    open(newunit=fileID, file=trim(dirName) // "/VARFLT_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) VARFLT
    close(fileID)
    write(*,*) 'sum(abs(VARFLT)) = ', sum(abs(VARFLT))

    open(newunit=fileID, file=trim(dirName) // "/PLE_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) PLE
    close(fileID)
    write(*,*) 'sum(abs(PLE)) = ', sum(abs(PLE))

    open(newunit=fileID, file=trim(dirName) // "/BKV_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) BKV
    close(fileID)
    write(*,*) 'sum(abs(BKV)) = ', sum(abs(BKV))

    open(newunit=fileID, file=trim(dirName) // "/BKUU_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) BKUU
    close(fileID)
    write(*,*) 'sum(abs(BKUU)) = ', sum(abs(BKUU))

    open(newunit=fileID, file=trim(dirName) // "/FKV_" // trim(rank_str) // ".in", &
         status='old', form="unformatted", action="read")
    read(fileID) FKV
    close(fileID)
    write(*,*) 'sum(abs(FKV)) = ', sum(abs(FKV))

    open(newunit=fileID, file=trim(dirName) // "/BKV_" // trim(rank_str) // ".out", &
         status='old', form="unformatted", action="read")
    read(fileID) BKV_ref
    close(fileID)
    write(*,*) 'sum(abs(BKV_ref)) = ', sum(abs(BKV_ref))

    open(newunit=fileID, file=trim(dirName) // "/BKUU_" // trim(rank_str) // ".out", &
         status='old', form="unformatted", action="read")
    read(fileID) BKUU_ref
    close(fileID)
    write(*,*) 'sum(abs(BKUU_ref)) = ', sum(abs(BKUU_ref))

    open(newunit=fileID, file=trim(dirName) // "/FKV_" // trim(rank_str) // ".out", &
         status='old', form="unformatted", action="read")
    read(fileID) FKV_ref
    close(fileID)
    write(*,*) 'sum(abs(FKV_ref)) = ', sum(abs(FKV_ref))

!$acc enter data copyin(U, V, Z, VARFLT, PLE) &
!$acc            copyin(BKV, BKUU) &
!$acc            create(FKV)

    call cpu_time(ts)

    call BELJAARS(IM, JM, LM, DT, &
                    LAMBDA_B, C_B,  &
                    U, V, Z,        &
                    VARFLT, PLE,    &
                    BKV, BKUU, FKV  )

    call cpu_time(te)

!$acc exit data copyout(BKV, BKUU, FKV)

    write(*,*) 'Sum diff abs BKV = ',  sum(BKV  - BKV_ref),   sum(abs(BKV)),  sum(abs(BKV_ref))
    write(*,*) 'Sum diff abs BKUU = ', sum(BKUU  - BKUU_ref), sum(abs(BKUU)), sum(abs(BKUU_ref))
    write(*,*) 'Sum diff abs FKV = ',  sum(FKV  - FKV_ref),   sum(abs(FKV)),  sum(abs(FKV_ref))
    print*,'****'
    print*,'BELJAARS Run Time = ', te - ts
    print*,'****'
    
    contains

    subroutine BELJAARS(IM, JM, LM, DT, &
        LAMBDA_B, C_B,  &
        U, V, Z,        &
        VARFLT, PLE,    &
        BKV, BKVV, FKV  )

    !BOP
    !
    !   Orographic drag follows  Beljaars (2003):
    !   $$
    !   \frac{\partial}{\partial z}\frac{\tau}{\rho} = \frac{C_B}{\lambda_B} |U(z)| U(z) 
    !          e^{-\tilde{z}^\frac{3}{2}}\tilde{z}^{-1.2},
    !   $$
    !   where $z$ is the height above the surface in meters, 
    !   $\tilde{z}=\frac{z}{\lambda_B}$, $\tau$ is the orographic stress at $z$,
    !   $\rho$ is the air density, $U(z)$ is the wind velocity, and $\lambda_B$ is a vertical length scale.
    !   Beljaars uses $\lambda_B = 1500$m, for which the non-dimensional parameter $C_B = 2.5101471 \times 10^{-8}$.
    !   These are the default values, but both can be modified from the configuration. To avoid underflow.
    !   the tendency is set to zero once $\tilde{z}$ exceeds 4 (i.e., 6 km from the surface for default values). 
    !
    !EOP

        integer, intent(IN   )                    :: IM,JM,LM
        real,    intent(IN   )                    :: DT
        real,    intent(IN   )                    :: LAMBDA_B
        real,    intent(IN   )                    :: C_B

        real,    intent(IN   ), dimension(:,:,: ) :: U
        real,    intent(IN   ), dimension(:,:,: ) :: V
        real,    intent(IN   ), dimension(:,:,: ) :: Z
        real,    intent(IN   ), dimension(:,:   ) :: VARFLT
        real,    intent(IN   ), dimension(:,:,0:) :: PLE

        real,    intent(INOUT), dimension(:,:,: ) :: BKV,BKVV

        real,    intent(  OUT), dimension(:,:,: ) :: FKV

        integer :: I,J,L
        real    :: FKV_temp

!$acc data present(U, V, Z, VARFLT, PLE, BKV, BKVV, FKV)

!$acc parallel
!$acc loop gang vector collapse(3) private(FKV_temp)
        do L = 1,LM
            do J = 1,JM
                do I = 1,IM
                    FKV(I,J,L) = 0.0

                    if (Z(I,J,L) < 4.0*LAMBDA_B) then
                        FKV_temp = Z(I,J,L)*(1.0/LAMBDA_B)
                        FKV_temp = VARFLT(I,J) * exp(-FKV_temp*sqrt(FKV_temp))*(FKV_temp**(-1.2))
                        FKV_temp = (C_B/LAMBDA_B)*min( sqrt(U(I,J,L)**2+V(I,J,L)**2),5.0 )*FKV_temp

                        BKV(I,J,L) = BKV(I,J,L) + DT*FKV_temp
                        BKVV(I,J,L) = BKVV(I,J,L) + DT*FKV_temp
                        FKV(I,J,L) = FKV_temp * (PLE(I,J,L)-PLE(I,J,L-1))
                    end if
                end do
            end do 
        end do 

!$acc end parallel

!$acc end data

    end subroutine BELJAARS
end program beljaars_standalone

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