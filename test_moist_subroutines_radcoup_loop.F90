module test_radcoup_subroutines

    use radcoup_loop
    use timing_module

    implicit none

    public test_radcoup_loop

    private

    contains

    subroutine test_radcoup_loop(IM, JM, LM, dirName, rank_str)

        implicit none

        integer, intent(IN) :: IM, JM, LM

        character*100, intent(IN) :: dirName, rank_str

        integer :: fileID

        logical :: do_qa

        real :: FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI

        real, dimension(:,:,:), allocatable :: PLmb, QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, &
            Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, &
            RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, &
            RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, RHX

        real, dimension(:,:,:), allocatable :: Q_ref, T_ref, QLLS_ref, QILS_ref, CLLS_ref, &
            QLCN_ref, QICN_ref, CLCN_ref, RAD_QV_ref, RAD_QL_ref, RAD_QI_ref, RAD_QR_ref, &
            RAD_QS_ref, RAD_QG_ref, RAD_CF_ref, CLDREFFL_ref, CLDREFFI_ref, RHX_ref

        allocate(PLmb(IM, JM, LM))
        allocate(QRAIN(IM, JM, LM))
        allocate(QSNOW(IM, JM, LM))
        allocate(QGRAUPEL(IM, JM, LM))
        allocate(NACTL(IM, JM, LM))
        allocate(NACTI(IM, JM, LM))
        allocate(Q(IM, JM, LM))
        allocate(T(IM, JM, LM))
        allocate(QLLS(IM, JM, LM))
        allocate(QILS(IM, JM, LM))
        allocate(CLLS(IM, JM, LM))
        allocate(QLCN(IM, JM, LM))
        allocate(QICN(IM, JM, LM))
        allocate(CLCN(IM, JM, LM))
        allocate(RAD_QV(IM, JM, LM))
        allocate(RAD_QL(IM, JM, LM))
        allocate(RAD_QI(IM, JM, LM))
        allocate(RAD_QR(IM, JM, LM))
        allocate(RAD_QS(IM, JM, LM))
        allocate(RAD_QG(IM, JM, LM))
        allocate(RAD_CF(IM, JM, LM))
        allocate(CLDREFFL(IM, JM, LM))
        allocate(CLDREFFI(IM, JM, LM))
        allocate(RHX(IM, JM, LM))
        allocate(Q_ref(IM, JM, LM))
        allocate(T_ref(IM, JM, LM))
        allocate(QLLS_ref(IM, JM, LM))
        allocate(QILS_ref(IM, JM, LM))
        allocate(CLLS_ref(IM, JM, LM))
        allocate(QLCN_ref(IM, JM, LM))
        allocate(QICN_ref(IM, JM, LM))
        allocate(CLCN_ref(IM, JM, LM))
        allocate(RAD_QV_ref(IM, JM, LM))
        allocate(RAD_QL_ref(IM, JM, LM))
        allocate(RAD_QI_ref(IM, JM, LM))
        allocate(RAD_QR_ref(IM, JM, LM))
        allocate(RAD_QS_ref(IM, JM, LM))
        allocate(RAD_QG_ref(IM, JM, LM))
        allocate(RAD_CF_ref(IM, JM, LM))
        allocate(CLDREFFL_ref(IM, JM, LM))
        allocate(CLDREFFI_ref(IM, JM, LM))
        allocate(RHX_ref(IM, JM, LM))

        open(newunit=fileID, file=trim(dirName) // '/PLmb_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLmb
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(PLmb) = ', sum(PLmb)

        open(newunit=fileID, file=trim(dirName) // '/QRAIN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QRAIN
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QRAIN) = ', sum(QRAIN)

        open(newunit=fileID, file=trim(dirName) // '/QSNOW_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QSNOW
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QSNOW) = ', sum(QSNOW)

        open(newunit=fileID, file=trim(dirName) // '/QGRAUPEL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QGRAUPEL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QGRAUPEL) = ', sum(QGRAUPEL)

        open(newunit=fileID, file=trim(dirName) // '/NACTL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(NACTL) = ', sum(NACTL)

        open(newunit=fileID, file=trim(dirName) // '/NACTI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTI
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(NACTI) = ', sum(NACTI)

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) Q
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(Q) = ', sum(Q)

        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) T
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(T) = ', sum(T)

        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLLS
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QLLS) = ', sum(QLLS)

        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QILS
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QILS) = ', sum(QILS)

        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLLS
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(CLLS) = ', sum(CLLS)

        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLCN
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QLCN) = ', sum(QLCN)

        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QICN
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(QICN) = ', sum(QICN)

        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLCN
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(CLCN) = ', sum(CLCN)

        open(newunit=fileID, file=trim(dirName) // '/FAC_RL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) FAC_RL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(FAC_RL) = ', FAC_RL

        open(newunit=fileID, file=trim(dirName) // '/MIN_RL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MIN_RL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(MIN_RL) = ', MIN_RL

        open(newunit=fileID, file=trim(dirName) // '/MAX_RL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MAX_RL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(MAX_RL) = ', MAX_RL

        open(newunit=fileID, file=trim(dirName) // '/FAC_RI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) FAC_RI
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(FAC_RI) = ', FAC_RI

        open(newunit=fileID, file=trim(dirName) // '/MIN_RI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MIN_RI
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(MIN_RI) = ', MIN_RI

        open(newunit=fileID, file=trim(dirName) // '/MAX_RI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MAX_RI
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(MAX_RI) = ', MAX_RI

        open(newunit=fileID, file=trim(dirName) // '/do_qa_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_qa
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_qa) = ', do_qa

        open(newunit=fileID, file=trim(dirName) // '/RHX_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RHX
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(RHX) = ', sum(RHX)

        call start_timing()

        call radcoup_loop_standalone(IM, JM, LM, FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI, do_qa, &
            PLmb, QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, &
            Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, &
            RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, &
            RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, RHX)

        call end_timing()

        call print_timing()

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) Q_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) T_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QLLS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QILS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLLS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QLCN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QICN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLCN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QV_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QV_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QI_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QI_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QR_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QR_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QG_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QG_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_CF_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_CF_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLDREFFL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLDREFFL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLDREFFI_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLDREFFI_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RHX_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RHX_ref
        close(fileID)

        print*,'Compare sum(diff(Q)) = ',sum(Q_ref - Q)
        print*,'Compare sum(Q) = ',sum(Q)
        print*,'Compare sum(Q_ref) = ',sum(Q_ref)
        print*,'***'
        print*,'Compare sum(diff(T)) = ',sum(T_ref - T)
        print*,'Compare sum(T) = ',sum(T)
        print*,'Compare sum(T_ref) = ',sum(T_ref)
        print*,'***'
        print*,'Compare sum(diff(QLLS)) = ',sum(QLLS_ref - QLLS)
        print*,'Compare sum(QLLS) = ',sum(QLLS)
        print*,'Compare sum(QLLS_ref) = ',sum(QLLS_ref)
        print*,'***'
        print*,'Compare sum(diff(QILS)) = ',sum(QILS_ref - QILS)
        print*,'Compare sum(QILS) = ',sum(QILS)
        print*,'Compare sum(QILS_ref) = ',sum(QILS_ref)
        print*,'***'
        print*,'Compare sum(diff(CLLS)) = ',sum(CLLS_ref - CLLS)
        print*,'Compare sum(CLLS) = ',sum(CLLS)
        print*,'Compare sum(CLLS_ref) = ',sum(CLLS_ref)
        print*,'***'
        print*,'Compare sum(diff(QLCN)) = ',sum(QLCN_ref - QLCN)
        print*,'Compare sum(QLCN) = ',sum(QLCN)
        print*,'Compare sum(QLCN_ref) = ',sum(QLCN_ref)
        print*,'***'
        print*,'Compare sum(diff(QICN)) = ',sum(QICN_ref - QICN)
        print*,'Compare sum(QICN) = ',sum(QICN)
        print*,'Compare sum(QICN_ref) = ',sum(QICN_ref)
        print*,'***'
        print*,'Compare sum(diff(CLCN)) = ',sum(CLCN_ref - CLCN)
        print*,'Compare sum(CLCN) = ',sum(CLCN)
        print*,'Compare sum(CLCN_ref) = ',sum(CLCN_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QV)) = ',sum(RAD_QV_ref - RAD_QV)
        print*,'Compare sum(RAD_QV) = ',sum(RAD_QV)
        print*,'Compare sum(RAD_QV_ref) = ',sum(RAD_QV_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QL)) = ',sum(RAD_QL_ref - RAD_QL)
        print*,'Compare sum(RAD_QL) = ',sum(RAD_QL)
        print*,'Compare sum(RAD_QL_ref) = ',sum(RAD_QL_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QI)) = ',sum(RAD_QI_ref - RAD_QI)
        print*,'Compare sum(RAD_QI) = ',sum(RAD_QI)
        print*,'Compare sum(RAD_QI_ref) = ',sum(RAD_QI_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QR)) = ',sum(RAD_QR_ref - RAD_QR)
        print*,'Compare sum(RAD_QR) = ',sum(RAD_QR)
        print*,'Compare sum(RAD_QR_ref) = ',sum(RAD_QR_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QS)) = ',sum(RAD_QS_ref - RAD_QS)
        print*,'Compare sum(RAD_QS) = ',sum(RAD_QS)
        print*,'Compare sum(RAD_QS_ref) = ',sum(RAD_QS_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QG)) = ',sum(RAD_QG_ref - RAD_QG)
        print*,'Compare sum(RAD_QG) = ',sum(RAD_QG)
        print*,'Compare sum(RAD_QG_ref) = ',sum(RAD_QG_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_CF)) = ',sum(RAD_CF_ref - RAD_CF)
        print*,'Compare sum(RAD_CF) = ',sum(RAD_CF)
        print*,'Compare sum(RAD_CF_ref) = ',sum(RAD_CF_ref)
        print*,'***'
        print*,'Compare sum(diff(CLDREFFL)) = ',sum(CLDREFFL_ref - CLDREFFL)
        print*,'Compare sum(CLDREFFL) = ',sum(CLDREFFL)
        print*,'Compare sum(CLDREFFL_ref) = ',sum(CLDREFFL_ref)
        print*,'***'
        print*,'Compare sum(diff(CLDREFFI)) = ',sum(CLDREFFI_ref - CLDREFFI)
        print*,'Compare sum(CLDREFFI) = ',sum(CLDREFFI)
        print*,'Compare sum(CLDREFFI_ref) = ',sum(CLDREFFI_ref)
        print*,'***'
        print*,'Compare sum(diff(RHX)) = ',sum(RHX_ref - RHX)
        print*,'Compare sum(RHX) = ',sum(RHX)
        print*,'Compare sum(RHX_ref) = ',sum(RHX_ref)
        print*,'***'

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