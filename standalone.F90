program update_moments_standalone

    use shoc

    implicit none

    character*100 :: dirName, rank_str
    integer :: IM, JM, LM, fileID
    real :: DT, HL2tune, qt2tune, Hlqt2tune, qt2scale, qt3_tscale

    real, dimension(:,:),   allocatable :: SH, EVAP
    real, dimension(:,:,:), allocatable :: Z, KH, TKESHOC, ISOTROPY, QT, HL, EDMF_FRC
    real, dimension(:,:,:), allocatable :: MFQT2, MFQT3, MFHL2, MFHL3, MFW2, MFW3
    real, dimension(:,:,:), allocatable :: MFWQT, MFWHL, MFHLQT, QT2, QT3, QT2_ref, QT3_ref
    real, dimension(:,:,:), allocatable :: HL2, HL2_ref, HL3, HL3_ref, W2, W2_ref, W3, W3_ref
    real, dimension(:,:,:), allocatable :: Wqt, Wqt_ref, Whl, Whl_ref, Hlqt, Hlqt_ref

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
    endif

    HL2tune = 0.30
    qt2tune = 2.00
    Hlqt2tune = 1.0
    qt2scale = 2000.0
    qt3_tscale = 3600.0

    allocate(SH(IM,JM))
    allocate(EVAP(IM,JM))
    allocate(Z(IM,JM,LM))
    allocate(KH(IM,JM,0:LM))
    allocate(TKESHOC(IM,JM,LM))
    allocate(ISOTROPY(IM,JM,LM))
    allocate(QT(IM,JM,LM))
    allocate(HL(IM,JM,LM))
    allocate(EDMF_FRC(IM,JM,LM))
    allocate(MFQT2(IM,JM,LM))
    allocate(MFQT3(IM,JM,LM))
    allocate(MFHL2(IM,JM,LM))
    allocate(MFHL3(IM,JM,LM))
    allocate(MFW2(IM,JM,LM))
    allocate(MFW3(IM,JM,LM))
    allocate(MFWQT(IM,JM,LM))
    allocate(MFWHL(IM,JM,LM))
    allocate(MFHLQT(IM,JM,LM))
    allocate(QT2(IM,JM,LM))
    allocate(QT2_ref(IM,JM,LM))
    allocate(QT3(IM,JM,LM))
    allocate(QT3_ref(IM,JM,LM))
    allocate(HL2(IM,JM,LM))
    allocate(HL2_ref(IM,JM,LM))
    allocate(HL3(IM,JM,LM))
    allocate(HL3_ref(IM,JM,LM))
    allocate(W2(IM,JM,LM))
    allocate(W2_ref(IM,JM,LM))
    allocate(W3(IM,JM,LM))
    allocate(W3_ref(IM,JM,LM))
    allocate(Wqt(IM,JM,LM))
    allocate(Wqt_ref(IM,JM,LM))
    allocate(Whl(IM,JM,LM))
    allocate(Whl_ref(IM,JM,LM))
    allocate(Hlqt(IM,JM,LM))
    allocate(Hlqt_ref(IM,JM,LM))

!   Reading Input Arrays

    open(newunit=fileID, file=trim(dirName) // "/SH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) SH
    close(fileID)
    !write(*,*) 'sum(SH) = ', sum(SH)

    open(newunit=fileID, file=trim(dirName) // "/EVAP_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) EVAP
    close(fileID)
    !write(*,*) 'sum(EVAP) = ', sum(EVAP)

    open(newunit=fileID, file=trim(dirName) // "/Z_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Z
    close(fileID)
    !write(*,*) 'sum(Z) = ', sum(Z)

    open(newunit=fileID, file=trim(dirName) // "/KH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) KH
    close(fileID)
    !write(*,*) 'sum(KH) = ', sum(KH)

    open(newunit=fileID, file=trim(dirName) // "/TKESHOC_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TKESHOC
    close(fileID)
    !write(*,*) 'sum(TKESHOC) = ', sum(TKESHOC)

    open(newunit=fileID, file=trim(dirName) // "/ISOTROPY_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) ISOTROPY
    close(fileID)
    !write(*,*) 'sum(ISOTROPY) = ', sum(ISOTROPY)

    open(newunit=fileID, file=trim(dirName) // "/QT_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QT
    close(fileID)
    !write(*,*) 'sum(QT) = ', sum(QT)

    open(newunit=fileID, file=trim(dirName) // "/HL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) HL
    close(fileID)
    !write(*,*) 'sum(HL) = ', sum(HL)

    open(newunit=fileID, file=trim(dirName) // "/EDMF_FRC_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) EDMF_FRC
    close(fileID)
    !write(*,*) 'sum(EDMF_FRC) = ', sum(EDMF_FRC)

    open(newunit=fileID, file=trim(dirName) // "/MFQT2_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFQT2
    close(fileID)
    !write(*,*) 'sum(MFQT2) = ', sum(MFQT2)

    open(newunit=fileID, file=trim(dirName) // "/MFQT3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFQT3
    close(fileID)
    !write(*,*) 'sum(MFQT3) = ', sum(MFQT3)

    open(newunit=fileID, file=trim(dirName) // "/MFHL2_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFHL2
    close(fileID)
    !write(*,*) 'sum(MFHL2) = ', sum(MFHL2)

    open(newunit=fileID, file=trim(dirName) // "/MFHL3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFHL3
    close(fileID)
    !write(*,*) 'sum(MFHL3) = ', sum(MFHL3)

    open(newunit=fileID, file=trim(dirName) // "/MFW2_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFW2
    close(fileID)
    !write(*,*) 'sum(MFW2) = ', sum(MFW2)

    open(newunit=fileID, file=trim(dirName) // "/MFW3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFW3
    close(fileID)
    !write(*,*) 'sum(MFW3) = ', sum(MFW3)

    open(newunit=fileID, file=trim(dirName) // "/MFWQT_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFWQT
    close(fileID)
    !write(*,*) 'sum(MFWQT) = ', sum(MFWQT)

    open(newunit=fileID, file=trim(dirName) // "/MFWHL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFWHL
    close(fileID)
    !write(*,*) 'sum(MFWHL) = ', sum(MFWHL)

    open(newunit=fileID, file=trim(dirName) // "/MFHLQT_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MFHLQT
    close(fileID)
    !write(*,*) 'sum(MFHLQT) = ', sum(MFHLQT)

    open(newunit=fileID, file=trim(dirName) // "/QT2_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QT2
    close(fileID)
    !write(*,*) 'sum(QT2) = ', sum(QT2)

    open(newunit=fileID, file=trim(dirName) // "/QT3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QT3
    close(fileID)
    !write(*,*) 'sum(QT3) = ', sum(QT3)

! Reading output reference solution arrays

    open(newunit=fileID, file=trim(dirName) // "/QT2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) QT2_ref
    close(fileID)
    !write(*,*) 'sum(QT2_ref) = ', sum(QT2_ref)

    open(newunit=fileID, file=trim(dirName) // "/QT3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) QT3_ref
    close(fileID)
    !write(*,*) 'sum(QT3_ref) = ', sum(QT3_ref)

    open(newunit=fileID, file=trim(dirName) // "/HL2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) HL2_ref
    close(fileID)
    !write(*,*) 'sum(HL2_ref) = ', sum(HL2_ref)

    open(newunit=fileID, file=trim(dirName) // "/HL3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) HL3_ref
    close(fileID)
    !write(*,*) 'sum(HL3_ref) = ', sum(HL3_ref)

    open(newunit=fileID, file=trim(dirName) // "/W2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) W2_ref
    close(fileID)
    !write(*,*) 'sum(W2_ref) = ', sum(W2_ref)

    open(newunit=fileID, file=trim(dirName) // "/W3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) W3_ref
    close(fileID)
    !write(*,*) 'sum(W3_ref) = ', sum(W3_ref)

    open(newunit=fileID, file=trim(dirName) // "/Wqt_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) Wqt_ref
    close(fileID)
    !write(*,*) 'sum(Wqt_ref) = ', sum(Wqt_ref)

    open(newunit=fileID, file=trim(dirName) // "/Whl_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) Whl_ref
    close(fileID)
    !write(*,*) 'sum(Whl_ref) = ', sum(Whl_ref)

    open(newunit=fileID, file=trim(dirName) // "/Hlqt_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) Hlqt_ref
    close(fileID)
    !write(*,*) 'sum(Hlqt_ref) = ', sum(Hlqt_ref)

!$acc enter data copyin(SH, EVAP, Z, KH, TKESHOC, ISOTROPY, QT, HL, EDMF_FRC, MFQT2, MFQT3) &
!$acc            copyin(MFHL2, MFHL3, MFW2, MFW3, MFWQT, MFWHL, MFHLQT, qt2, qt3) &
!$acc            create(hl2, hl3, w2, w3, wqt, whl, hlqt)

    call update_moments(IM, JM, LM, DT, &
                          SH,             &  ! in
                          EVAP,           &
                          Z,              &
                          KH,             &
                          TKESHOC,        &
                          ISOTROPY,       &
                          QT,             &
                          HL,             &
                          EDMF_FRC,       &
                          MFQT2,          &
                          MFQT3,          &
                          MFHL2,          &
                          MFHL3,          &
                          MFW2,           &
                          MFW3,           &
                          MFWQT,          &
                          MFWHL,          &
                          MFHLQT,         &
                          qt2,            &  ! inout
                          qt3,            &
                          hl2,            &  ! out
                          hl3,            &
                          w2,             &
                          w3,             &
                          wqt,            &
                          whl,            &
                          hlqt,           &
                          hl2tune,        &  ! tuning parameters
                          qt2tune,        &
                          hlqt2tune,      &
                          qt2scale,       &
                          qt3_tscale )
!$acc exit data copyout(qt2, qt3, hl2, hl3, w2, w3, wqt, whl, hlqt)

    write(*,*) 'Sum diff abs QT2 = ', sum(QT2  - QT2_ref),  sum(abs(QT2)),  sum(abs(QT2_ref))
    write(*,*) 'Sum diff abs QT3 = ', sum(QT3  - QT3_ref),  sum(abs(QT3)),  sum(abs(QT3_ref))
    write(*,*) 'Sum diff abs HL2 = ', sum(HL2  - HL2_ref),  sum(abs(HL2)),  sum(abs(HL2_ref))
    write(*,*) 'Sum diff abs HL3 = ', sum(HL3  - HL3_ref),  sum(abs(HL3)),  sum(abs(HL3_ref))
    write(*,*) 'Sum diff abs W2  = ', sum(W2   - W2_ref),   sum(abs(W2)),   sum(abs(W2_ref))
    write(*,*) 'Sum diff abs W3  = ', sum(W3   - W3_ref),   sum(abs(W3)),   sum(abs(W3_ref))
    write(*,*) 'Sum diff abs WQT = ', sum(WQT  - Wqt_ref),  sum(abs(WQT)),  sum(abs(Wqt_ref))
    write(*,*) 'Sum diff abs WHL = ', sum(WHL  - Whl_ref),  sum(abs(WHL)),  sum(abs(Whl_ref))
    write(*,*) 'Sum diff abs HLQT =', sum(HLQT - Hlqt_ref), sum(abs(HLQT)), sum(abs(Hlqt_ref))

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