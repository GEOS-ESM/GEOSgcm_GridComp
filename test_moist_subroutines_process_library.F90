module test_process_library_subroutines

    use Process_Library_standalone
    use timing_module
    use GEOS_UtilsMod

    implicit none

    public test_BUOYANCY2

    private

    real, dimension(:,:,:), allocatable :: MASS, Q, Q_comp
    real, dimension(:,:,:), allocatable :: T, QST3, DQST3, DZET, ZL0, BYNCY, BYNCY_comp, PLmb, PLEmb
    real, dimension(:,:),   allocatable :: SBCAPE, SBCIN, SBCAPE_comp, SBCIN_comp
    real, dimension(:,:),   allocatable :: LFC, LNB, LFC_comp, LNB_comp
    real, dimension(:,:),   pointer     :: MLCAPE, MUCAPE, MLCIN, MUCIN
    real, dimension(:,:),   allocatable :: MLCAPE_comp, MUCAPE_comp, MLCIN_comp, MUCIN_comp

    contains

    subroutine test_BUOYANCY2(IM, JM, LM, dirName, rank_str)

        integer :: IM, JM, LM, fileID, iter
        character*100 :: dirName, rank_str

        allocate(T(IM, JM, LM))
        allocate(Q(IM, JM, LM))
        allocate(QST3(IM, JM, LM))
        allocate(DQST3(IM, JM, LM))
        allocate(DZET(IM, JM, LM))
        allocate(ZL0(IM, JM, LM))
        allocate(PLMB(IM, JM, LM))
        allocate(PLEmb(IM, JM, 0:LM))
        allocate(SBCAPE(IM, JM))
        allocate(SBCAPE_comp(IM, JM))
        allocate(MLCAPE(IM, JM))
        allocate(MLCAPE_comp(IM, JM))
        allocate(MUCAPE(IM, JM))
        allocate(MUCAPE_comp(IM, JM))
        allocate(SBCIN(IM, JM))
        allocate(SBCIN_comp(IM, JM))
        allocate(MLCIN(IM, JM))
        allocate(MLCIN_comp(IM, JM))
        allocate(MUCIN(IM, JM))
        allocate(MUCIN_comp(IM, JM))
        allocate(BYNCY(IM, JM, LM))
        allocate(BYNCY_comp(IM, JM, LM))
        allocate(LFC(IM, JM))
        allocate(LNB(IM, JM))
        allocate(LFC_comp(IM, JM))
        allocate(LNB_comp(IM, JM))

        open(newunit=fileID, file=trim(dirName) // "/T_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) T
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) Q
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/QST3_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) QST3
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/DQST3_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) DQST3
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/DZET_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) DZET
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ZL0_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ZL0
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/PLmb_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) PLmb
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/PLEmb_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) PLEmb
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/SBCAPE_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) SBCAPE
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MLCAPE_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) MLCAPE
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MUCAPE_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) MUCAPE
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/SBCIN_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) SBCIN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MLCIN_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) MLCIN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MUCIN_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) MUCIN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/LFC_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) LFC
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/LNB_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) LNB
        close(fileID)

!$acc data copyin(T, Q, QST3, DQST3, DZET, ZL0, PLmb, PLEmb), &
!$acc      copyout(BYNCY) &
!$acc      copy(SBCAPE, MLCAPE, MUCAPE, SBCIN, MLCIN, MUCIN, LFC, LNB)

        ! This is to initialize arrays in GEOS_UtilsMod ahead of time
        call ESINIT_v2

        print*,'Testing Buoyancy2 Subroutine'

        do iter = 1,1
            call start_timing()
            call BUOYANCY2( IM, JM, LM, T, Q, QST3, DQST3, DZET, ZL0, PLmb, PLEmb(:,:,LM), SBCAPE, MLCAPE, MUCAPE, SBCIN, MLCIN, MUCIN, BYNCY, LFC, LNB )
            call end_timing()
            call print_timing()
        enddo
!$acc end data
        open(newunit=fileID, file=trim(dirName) // "/BYNCY_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) BYNCY_comp
        close(fileID)

        print*,'sum(BYNCY - BYNCY_comp) = ', sum(BYNCY - BYNCY_comp)
        print*,'sum(BYNCY) = ', sum(BYNCY)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/SBCAPE_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) SBCAPE_comp
        close(fileID)

        print*,'sum(SBCAPE - SBCAPE_comp) = ', sum(SBCAPE - SBCAPE_comp)
        print*,'sum(SBCAPE) = ', sum(SBCAPE)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/SBCIN_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) SBCIN_comp
        close(fileID)

        print*,'sum(SBCIN - SBCIN_comp) = ', sum(SBCIN - SBCIN_comp)
        print*,'sum(SBCIN) = ', sum(SBCIN)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/LFC_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) LFC_comp
        close(fileID)

        print*,'sum(LFC - LFC_comp) = ', sum(LFC - LFC_comp)
        print*,'sum(LFC) = ', sum(LFC)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/LNB_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) LNB_comp
        close(fileID)

        print*,'sum(LNB - LNB_comp) = ', sum(LNB - LNB_comp)
        print*,'sum(LNB) = ', sum(LNB)

        print*,'*****'
        
        open(newunit=fileID, file=trim(dirName) // "/MLCAPE_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) MLCAPE_comp
        close(fileID)

        print*,'sum(MLCAPE - MLCAPE_comp) = ', sum(MLCAPE - MLCAPE_comp)
        print*,'sum(MLCAPE) = ', sum(MLCAPE)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/MUCAPE_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) MUCAPE_comp
        close(fileID)

        print*,'sum(MUCAPE - MUCAPE_comp) = ', sum(MUCAPE - MUCAPE_comp)
        print*,'sum(MUCAPE) = ', sum(MUCAPE)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/MLCIN_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) MLCIN_comp
        close(fileID)

        print*,'sum(MLCIN - MLCIN_comp) = ', sum(MLCIN - MLCIN_comp)
        print*,'sum(MLCIN) = ', sum(MLCIN)

        print*,'*****'

        open(newunit=fileID, file=trim(dirName) // "/MUCIN_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) MUCIN_comp
        close(fileID)

        print*,'sum(MUCIN - MUCIN_comp) = ', sum(MUCIN - MUCIN_comp)
        print*,'sum(MUCIN) = ', sum(MUCIN)
        
    end subroutine
end module

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