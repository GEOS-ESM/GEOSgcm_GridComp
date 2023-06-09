module test_aer_activation_subroutine

    use moist_subroutines_aer_activation

    implicit none

    public test_aer_activation

    private

    real, dimension(:,:,:), &
          allocatable :: Q, T, PLmb, PLO, PLE, ZL0, ZLE0, QLCN, &
                         QICN, QLLS, QILS, TKE, TMP3D, NACTL, NACTI, &
                         NWFA, NACTL_comp, NACTI_comp, NWFA_comp
                    
    real, dimension(:,:), allocatable :: SH, EVAP, KPBL, FRLAND

    type(AerProps), dimension(:,:,:), allocatable :: AeroProps

    integer :: AERO

    logical :: USE_AERO_BUFFER
    real :: tStart, tEnd, time
!$acc declare create(AeroProps,AeroProps%num,AeroProps%dpg,AeroProps%sig,AeroProps%den, &
!$acc                AeroProps%kap,AeroProps%fdust,AeroProps%fsoot,AeroProps%forg)
    contains

    subroutine test_aer_activation(IM, JM, LM, dirName, rank_str)

        integer :: IM, JM, LM, fileID
        character*100 :: dirName, rank_str
        real :: CCN_OCN, CCN_LND

        allocate(Q   (IM, JM, LM))
        allocate(T   (IM, JM, LM))
        allocate(PLmb(IM, JM, LM))
        allocate(PLO (IM, JM, LM))
        allocate(PLE (IM, JM, 0:LM))
        allocate(ZL0 (IM, JM, LM))
        allocate(ZLE0(IM, JM, 0:LM))
        allocate(QLCN(IM, JM, LM))
        allocate(QICN(IM, JM, LM))
        allocate(QLLS(IM, JM, LM))
        allocate(QILS(IM, JM, LM))
        allocate(TKE(IM, JM, LM))
        allocate(TMP3D(IM, JM, LM))
        allocate(NACTL(IM, JM, LM))
        allocate(NACTI(IM, JM, LM))
        allocate(NWFA (IM, JM, LM))

        allocate(NACTL_comp(IM, JM, LM))
        allocate(NACTI_comp(IM, JM, LM))
        allocate(NWFA_comp (IM, JM, LM))

        allocate(SH    (IM, JM))
        allocate(EVAP  (IM, JM))
        allocate(KPBL  (IM, JM))
        allocate(FRLAND(IM, JM))
        
        allocate(AeroProps(IM, JM, LM))

        open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) Q
        close(fileID)
        ! write(*,*) 'sum(abs(Q)) = ', sum(abs(Q))

        open(newunit=fileID, file=trim(dirName) // "/T_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) T
        close(fileID)
        ! write(*,*) 'sum(abs(T)) = ', sum(abs(T))

        open(newunit=fileID, file=trim(dirName) // "/PLmb_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) PLmb
        close(fileID)
        ! write(*,*) 'sum(abs(PLmb)) = ', sum(abs(PLmb))

        open(newunit=fileID, file=trim(dirName) // "/PLE_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) PLE
        close(fileID)
        ! write(*,*) 'sum(abs(PLE)) = ', sum(abs(PLE))

        open(newunit=fileID, file=trim(dirName) // "/ZL0_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ZL0
        close(fileID)
        ! write(*,*) 'sum(abs(ZL0)) = ', sum(abs(ZL0))

        open(newunit=fileID, file=trim(dirName) // "/ZLE0_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ZLE0
        close(fileID)
        ! write(*,*) 'sum(abs(ZLE0)) = ', sum(abs(ZLE0))

        open(newunit=fileID, file=trim(dirName) // "/QLCN_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) QLCN
        close(fileID)
        ! write(*,*) 'sum(abs(QLCN)) = ', sum(abs(QLCN))

        open(newunit=fileID, file=trim(dirName) // "/QICN_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) QICN
        close(fileID)
        ! write(*,*) 'sum(abs(QICN)) = ', sum(abs(QICN))

        open(newunit=fileID, file=trim(dirName) // "/QLLS_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) QLLS
        close(fileID)
        ! write(*,*) 'sum(abs(QLLS)) = ', sum(abs(QLLS))

        open(newunit=fileID, file=trim(dirName) // "/QILS_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) QILS
        close(fileID)
        ! write(*,*) 'sum(abs(QILS)) = ', sum(abs(QILS))

        open(newunit=fileID, file=trim(dirName) // "/SH_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) SH
        close(fileID)
        ! write(*,*) 'sum(abs(SH)) = ', sum(abs(SH))

        open(newunit=fileID, file=trim(dirName) // "/EVAP_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) EVAP
        close(fileID)
        ! write(*,*) 'sum(abs(EVAP)) = ', sum(abs(EVAP))

        open(newunit=fileID, file=trim(dirName) // "/KPBL_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) KPBL
        close(fileID)
        ! write(*,*) 'sum(abs(KPBL)) = ', sum(abs(KPBL))

        open(newunit=fileID, file=trim(dirName) // "/TKE_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) TKE
        close(fileID)
        
        open(newunit=fileID, file=trim(dirName) // "/TMP3D_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) TMP3D
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/FRLAND_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) FRLAND
        close(fileID)
        ! write(*,*) 'sum(abs(FRLAND)) = ', sum(abs(FRLAND))

        open(newunit=fileID, file=trim(dirName) // "/USE_AERO_BUFFER_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) USE_AERO_BUFFER
        close(fileID)
        ! write(*,*) 'USE_AERO_BUFFER = ', USE_AERO_BUFFER

        open(newunit=fileID, file=trim(dirName) // "/CCN_LND_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) CCN_LND
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/CCN_OCN_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) CCN_OCN
        close(fileID)

        AERO = 7
        PLO = PLmb*100.0

        print*,'Testing Aer_Activation Subroutine'

!$acc data copy(AeroProps,AeroProps%num,AeroProps%dpg,AeroProps%sig,AeroProps%kap,AeroProps%den,AeroProps%fdust,&
!$acc           AeroProps%fsoot,AeroProps%forg) &
!$acc      copyin(Q,T,PLO,PLE,ZL0,ZLE0,QLCN,QICN,QLLS,QILS,SH,EVAP,KPBL,TKE,TMP3D,FRLAND) &
!$acc      copyout(NACTL, NACTI, NWFA)

        call cpu_time(tStart)
        call Aer_Activation(IM, JM, LM, Q, T, PLO, PLE, ZL0, ZLE0, QLCN, QICN, QLLS, QILS, &
                            SH, EVAP, KPBL, TKE, TMP3D, FRLAND, USE_AERO_BUFFER, &
                            AeroProps, AERO, NACTL, NACTI, NWFA, CCN_LND*1.e6, CCN_OCN*1.e6,dirName, rank_str)
        call cpu_time(tEnd)

!$acc end data

        time = tEnd - tStart
        print*, 'Total Aer_Activation Subroutine Time = ', time

        open(newunit=fileID, file=trim(dirName) // "/NACTL_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) NACTL_comp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/NACTI_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) NACTI_comp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/NWFA_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) NWFA_comp
        close(fileID)
        
        print*,'NACTL Sum Diff = ', sum(NACTL-NACTL_comp)
        print*,'sum(NACTL) = ', sum(NACTL)
        print*,'sum(NACTL_comp) = ', sum(NACTL_comp)
        print*,'***'
        print*,'NACTI Sum Diff = ', sum(NACTI-NACTI_comp)
        print*,'sum(NACTI) = ', sum(NACTI)
        print*,'sum(NACTI_comp) = ', sum(NACTI_comp)
        print*,'***'
        print*,'NWFA Sum Diff = ', sum(NWFA - NWFA_comp)
        print*,'sum(NWFA) = ', sum(NWFA)
        print*,'sum(NWFA_comp) = ', sum(NWFA_comp)
        write(*,*) '*******************************************************************'

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
