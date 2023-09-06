program test_UWSC

    use uwshcu
    use tracer_module
    implicit none

    integer   :: IM, JM, LM, num_args, USE_TRACER_TRANSP_UW, fileID
    real      :: UW_DT
    integer(kind=8) :: t_start, t_end
    real(kind=8) :: rate

    real, dimension(:,:,:), allocatable :: PL, ZL0, PK, PLE, ZLE0, PKE, DP
    real, dimension(:,:,:), allocatable :: U, V, Q, QLTOT, QITOT, T, TKE
    real, dimension(:,:),   allocatable :: KPBL_SC, SH, EVAP, CNPCPRATE, FRLAND, CUSH, CUSH_ref

    real, dimension(:,:,:), allocatable :: UMF_SC, DCM_SC, DQVDT_SC, DQLDT_SC, DQIDT_SC
    real, dimension(:,:,:), allocatable :: DTDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC
    real, dimension(:,:,:), allocatable :: DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC
    real, dimension(:,:,:), allocatable :: QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC
    real, dimension(:,:,:), allocatable :: SC_NDROP, SC_NICE, QTFLX_SC, SLFLX_SC
    real, dimension(:,:,:), allocatable :: UFLX_SC, VFLX_SC
    real, dimension(:,:),   allocatable :: TPERT_SC, QPERT_SC

    real, dimension(:,:,:), allocatable :: UMF_SC_ref, DCM_SC_ref, DQVDT_SC_ref, DQLDT_SC_ref, DQIDT_SC_ref
    real, dimension(:,:,:), allocatable :: DTDT_SC_ref, DUDT_SC_ref, DVDT_SC_ref, DQRDT_SC_ref
    real, dimension(:,:,:), allocatable :: DQSDT_SC_ref, CUFRC_SC_ref, ENTR_SC_ref, DETR_SC_ref
    real, dimension(:,:,:), allocatable :: QLDET_SC_ref, QIDET_SC_ref, QLSUB_SC_ref, QISUB_SC_ref
    real, dimension(:,:,:), allocatable :: SC_NDROP_ref, SC_NICE_ref, QTFLX_SC_ref, SLFLX_SC_ref
    real, dimension(:,:,:), allocatable :: UFLX_SC_ref, VFLX_SC_ref
    real, dimension(:,:),   allocatable :: TPERT_SC_ref, QPERT_SC_ref

    character *100 :: BUFFER, dirName, rank_str

    logical :: print_compare = .true.

    num_args = command_argument_count()

    if(num_args.ne.2) then
        print*, 'Missing arguments : <executable> <data directory> <trim(rank_str)>'
        call exit(1)
    else
        call get_command_argument(1, dirName)
        call get_command_argument(2, rank_str)
    endif

    ! *** These are the assumed dimensions of the arrays ***
    IM = 180
    JM = 180
    LM = 72

    write(*,*) 'IM = ', IM
    write(*,*) 'JM = ', JM
    write(*,*) 'LM = ', LM
    write(*,*) 'dirName = ', trim(dirName)

    allocate(PL(IM, JM, LM))
    allocate(ZL0(IM, JM, LM))
    allocate(PK(IM, JM, LM))
    allocate(PLE(IM, JM, 0:LM))
    allocate(ZLE0(IM, JM, 0:LM))
    allocate(PKE(IM, JM, 0:LM))
    allocate(DP(IM, JM, LM))
    allocate(U(IM, JM, LM))
    allocate(V(IM, JM, LM))
    allocate(Q(IM, JM, LM))
    allocate(QLTOT(IM, JM, LM))
    allocate(QITOT(IM, JM, LM))
    allocate(T(IM, JM, LM))
    allocate(TKE(IM, JM, 0:LM))
    allocate(KPBL_SC(IM, JM))
    allocate(SH(IM, JM))
    allocate(EVAP(IM, JM))
    allocate(CNPCPRATE(IM, JM))
    allocate(FRLAND(IM, JM))
    allocate(CUSH(IM, JM))

    allocate(UMF_SC(IM, JM, 0:LM))
    allocate(DCM_SC(IM, JM, LM))
    allocate(DQVDT_SC(IM, JM, LM))
    allocate(DQLDT_SC(IM, JM, LM))
    allocate(DQIDT_SC(IM, JM, LM))
    allocate(DTDT_SC(IM, JM, LM))
    allocate(DUDT_SC(IM, JM, LM))
    allocate(DVDT_SC(IM, JM, LM))
    allocate(DQRDT_SC(IM, JM, LM))
    allocate(DQSDT_SC(IM, JM, LM))
    allocate(CUFRC_SC(IM, JM, LM))
    allocate(ENTR_SC(IM, JM, LM))
    allocate(DETR_SC(IM, JM, LM))
    allocate(QLDET_SC(IM, JM, LM))
    allocate(QIDET_SC(IM, JM, LM))
    allocate(QLSUB_SC(IM, JM, LM))
    allocate(QISUB_SC(IM, JM, LM))
    allocate(SC_NDROP(IM, JM, LM))
    allocate(SC_NICE(IM, JM, LM))
    allocate(TPERT_SC(IM, JM))
    allocate(QPERT_SC(IM, JM))
    allocate(QTFLX_SC(IM, JM, 0:LM))
    allocate(SLFLX_SC(IM, JM, 0:LM))
    allocate(UFLX_SC(IM, JM, 0:LM))
    allocate(VFLX_SC(IM, JM, 0:LM))

    allocate(CUSH_ref(IM, JM))
    allocate(UMF_SC_ref(IM, JM, 0:LM))
    allocate(DCM_SC_ref(IM, JM, LM))
    allocate(DQVDT_SC_ref(IM, JM, LM))
    allocate(DQLDT_SC_ref(IM, JM, LM))
    allocate(DQIDT_SC_ref(IM, JM, LM))
    allocate(DTDT_SC_ref(IM, JM, LM))
    allocate(DUDT_SC_ref(IM, JM, LM))
    allocate(DVDT_SC_ref(IM, JM, LM))
    allocate(DQRDT_SC_ref(IM, JM, LM))
    allocate(DQSDT_SC_ref(IM, JM, LM))
    allocate(CUFRC_SC_ref(IM, JM, LM))
    allocate(ENTR_SC_ref(IM, JM, LM))
    allocate(DETR_SC_ref(IM, JM, LM))
    allocate(QLDET_SC_ref(IM, JM, LM))
    allocate(QIDET_SC_ref(IM, JM, LM))
    allocate(QLSUB_SC_ref(IM, JM, LM))
    allocate(QISUB_SC_ref(IM, JM, LM))
    allocate(SC_NDROP_ref(IM, JM, LM))
    allocate(SC_NICE_ref(IM, JM, LM))
    allocate(TPERT_SC_ref(IM, JM))
    allocate(QPERT_SC_ref(IM, JM))
    allocate(QTFLX_SC_ref(IM, JM, 0:LM))
    allocate(SLFLX_SC_ref(IM, JM, 0:LM))
    allocate(UFLX_SC_ref(IM, JM, 0:LM))
    allocate(VFLX_SC_ref(IM, JM, 0:LM))

    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_niter_xc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%niter_xc
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_iter_cin_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%iter_cin
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_use_CINcin_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%use_CINcin
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_use_self_detrain_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%use_self_detrain
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_use_momenflx_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%use_momenflx
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_use_cumpenent_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%use_cumpenent
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_scverbose_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%scverbose
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_windsrcavg_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%windsrcavg
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rpen_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rpen
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rle_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rle
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rkm_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rkm
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_mixscale_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%mixscale
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_detrhgt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%detrhgt
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rkfre_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rkfre
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rmaxfrac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rmaxfrac
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_mumin1_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%mumin1
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rbuoy_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rbuoy
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rdrag_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rdrag
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_epsvarw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%epsvarw
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_PGFc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%PGFc
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_criqc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%criqc
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_frc_rasn_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%frc_rasn
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_kevp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%kevp
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_rdrop_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%rdrop
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_thlsrc_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%thlsrc_fac
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_qtsrc_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%qtsrc_fac
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_qtsrchgt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%qtsrchgt
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // '/SHLWPARAMS_cridist_opt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SHLWPARAMS%cridist_opt
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/UW_DT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) UW_DT
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(UW_DT) = ', size(UW_DT)
    
    open(newunit=fileID, file=trim(dirName) // '/PL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) PL
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(PL) = ', size(PL)
    
    open(newunit=fileID, file=trim(dirName) // '/ZL0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) ZL0
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(ZL0) = ', size(ZL0)
    
    open(newunit=fileID, file=trim(dirName) // '/PK_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) PK
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(PK) = ', size(PK)
    
    open(newunit=fileID, file=trim(dirName) // '/PLE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) PLE
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(PLE) = ', size(PLE)
    
    open(newunit=fileID, file=trim(dirName) // '/ZLE0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) ZLE0
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(ZLE0) = ', size(ZLE0)
    
    open(newunit=fileID, file=trim(dirName) // '/PKE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) PKE
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(PKE) = ', size(PKE)
    
    open(newunit=fileID, file=trim(dirName) // '/DP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) DP
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(DP) = ', size(DP)
    
    open(newunit=fileID, file=trim(dirName) // '/U_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) U
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(U) = ', size(U)
    
    open(newunit=fileID, file=trim(dirName) // '/V_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) V
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(V) = ', size(V)
    
    open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) Q
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(Q) = ', size(Q)
    
    open(newunit=fileID, file=trim(dirName) // '/QLTOT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) QLTOT
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(QLTOT) = ', size(QLTOT)
    
    open(newunit=fileID, file=trim(dirName) // '/QITOT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) QITOT
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(QITOT) = ', size(QITOT)
    
    open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) T
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(T) = ', size(T)
    
    open(newunit=fileID, file=trim(dirName) // '/TKE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) TKE
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(TKE) = ', size(TKE)
    
    open(newunit=fileID, file=trim(dirName) // '/KPBL_SC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) KPBL_SC
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(KPBL_SC) = ', size(KPBL_SC)
    
    open(newunit=fileID, file=trim(dirName) // '/SH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) SH
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(SH) = ', size(SH)
    
    open(newunit=fileID, file=trim(dirName) // '/EVAP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) EVAP
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(EVAP) = ', size(EVAP)
    
    open(newunit=fileID, file=trim(dirName) // '/CNPCPRATE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) CNPCPRATE
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(CNPCPRATE) = ', size(CNPCPRATE)
    
    open(newunit=fileID, file=trim(dirName) // '/FRLAND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) FRLAND
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(FRLAND) = ', size(FRLAND)
    
    open(newunit=fileID, file=trim(dirName) // '/USE_TRACER_TRANSP_UW_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) USE_TRACER_TRANSP_UW
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(USE_TRACER_TRANSP_UW) = ', size(USE_TRACER_TRANSP_UW)
    
    open(newunit=fileID, file=trim(dirName) // '/CUSH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
    read(fileID) CUSH
    close(fileID)
    ! print*, 'Rank ', trim(rank_str),': In size(CUSH) = ', size(CUSH)

    call read_tracers(IM, JM, LM, dirName, rank_str)

    write(*,*) 'Calling compute_uwshcu_inv...'
!$acc data copyin(pmid0_inv,zmid0_inv,exnmid0_inv,pifc0_inv,zifc0_inv,&
!$acc      exnifc0_inv,dp0_inv,u0_inv,v0_inv,qv0_inv,ql0_inv,qi0_inv,&
!$acc      th0_inv,tke_inv,kpbl_inv)&
!$acc      copy(cush,tr0_inv) &
!$acc      copyout(umf_inv,qvten_inv,qlten_inv,qiten_inv,thten_inv,uten_inv,&
!$acc      vten_inv,qrten_inv,qsten_inv,cufrc_inv,fer_inv,fdr_inv, &
!$acc      qldet_inv,qidet_inv,qlsub_inv,qisub_inv,ndrop_inv,nice_inv)
    call system_clock(t_start,rate)

    call compute_uwshcu_inv(IM*JM, LM, UW_DT,           & ! IN
            PL, ZL0, PK, PLE, ZLE0, PKE, DP,              &
            U, V, Q, QLTOT, QITOT, T, TKE, KPBL_SC,       &
            SH, EVAP, CNPCPRATE, FRLAND,                  &
            CUSH,                                         & ! INOUT
            UMF_SC, DCM_SC, DQVDT_SC, DQLDT_SC, DQIDT_SC, & ! OUT
            DTDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC,          &
            DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC,         &
            QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC,       &
            SC_NDROP, SC_NICE, TPERT_SC, QPERT_SC,        &
            QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC,         &
#ifdef UWDIAG 
            QCU_SC, QLU_SC,                               & ! DIAG ONLY 
            QIU_SC, CBMF_SC, SHL_DQCDT, CNT_SC, CNB_SC,   &
            CIN_SC, PLCL_SC, PLFC_SC, PINV_SC, PREL_SC,   &
            PBUP_SC, WLCL_SC, QTSRC_SC, THLSRC_SC,        &
            THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, WUP_SC,     &
            QTUP_SC, THLUP_SC, THVUP_SC, UUP_SC, VUP_SC,  &
            XC_SC,                                        &
#endif 
            USE_TRACER_TRANSP_UW)
    call system_clock(t_end)
!$acc end data

    open(newunit=fileID, file=trim(dirName) // '/CUSH_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) CUSH_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/UMF_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) UMF_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DCM_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DCM_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DQVDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DQVDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DQLDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DQLDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DQIDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DQIDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DTDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DTDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DUDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DUDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DVDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DVDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DQRDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DQRDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DQSDT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DQSDT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/CUFRC_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) CUFRC_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/ENTR_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) ENTR_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/DETR_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) DETR_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/QLDET_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) QLDET_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/QIDET_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) QIDET_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/QLSUB_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) QLSUB_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/QISUB_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) QISUB_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/SC_NDROP_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) SC_NDROP_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/SC_NICE_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) SC_NICE_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/TPERT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) TPERT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/QPERT_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) QPERT_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/QTFLX_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) QTFLX_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/SLFLX_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) SLFLX_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/UFLX_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) UFLX_SC_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // '/VFLX_SC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
    read(fileID) VFLX_SC_ref
    close(fileID)

    if(print_compare.eq..true.) then
    
        print*, 'DIFF : sum(CUSH - CUSH_ref) = ', sum(CUSH - CUSH_ref)
        print*, 'sum(CUSH) = ', sum(CUSH)
        print*, 'sum(CUSH_ref) = ', sum(CUSH_ref)
        
        print*, 'DIFF : sum(UMF_SC - UMF_SC_ref) = ', sum(UMF_SC - UMF_SC_ref)
        print*, 'sum(UMF_SC) = ', sum(UMF_SC)
        print*, 'sum(UMF_SC_ref) = ', sum(UMF_SC_ref)
        
        print*, 'DIFF : sum(DCM_SC - DCM_SC_ref) = ', sum(DCM_SC - DCM_SC_ref)
        print*, 'sum(DCM_SC) = ', sum(DCM_SC)
        print*, 'sum(DCM_SC_ref) = ', sum(DCM_SC_ref)
        
        print*, 'DIFF : sum(DQVDT_SC - DQVDT_SC_ref) = ', sum(DQVDT_SC - DQVDT_SC_ref)
        print*, 'sum(DQVDT_SC) = ', sum(DQVDT_SC)
        print*, 'sum(DQVDT_SC_ref) = ', sum(DQVDT_SC_ref)
        
        print*, 'DIFF : sum(DQLDT_SC - DQLDT_SC_ref) = ', sum(DQLDT_SC - DQLDT_SC_ref)
        print*, 'sum(DQLDT_SC) = ', sum(DQLDT_SC)
        print*, 'sum(DQLDT_SC_ref) = ', sum(DQLDT_SC_ref)
        
        print*, 'DIFF : sum(DQIDT_SC - DQIDT_SC_ref) = ', sum(DQIDT_SC - DQIDT_SC_ref)
        print*, 'sum(DQIDT_SC) = ', sum(DQIDT_SC)
        print*, 'sum(DQIDT_SC_ref) = ', sum(DQIDT_SC_ref)
        
        print*, 'DIFF : sum(DTDT_SC - DTDT_SC_ref) = ', sum(DTDT_SC - DTDT_SC_ref)
        print*, 'sum(DTDT_SC) = ', sum(DTDT_SC)
        print*, 'sum(DTDT_SC_ref) = ', sum(DTDT_SC_ref)
        
        print*, 'DIFF : sum(DUDT_SC - DUDT_SC_ref) = ', sum(DUDT_SC - DUDT_SC_ref)
        print*, 'sum(DUDT_SC) = ', sum(DUDT_SC)
        print*, 'sum(DUDT_SC_ref) = ', sum(DUDT_SC_ref)
        
        print*, 'DIFF : sum(DVDT_SC - DVDT_SC_ref) = ', sum(DVDT_SC - DVDT_SC_ref)
        print*, 'sum(DVDT_SC) = ', sum(DVDT_SC)
        print*, 'sum(DVDT_SC_ref) = ', sum(DVDT_SC_ref)
        
        print*, 'DIFF : sum(DQRDT_SC - DQRDT_SC_ref) = ', sum(DQRDT_SC - DQRDT_SC_ref)
        print*, 'sum(DQRDT_SC) = ', sum(DQRDT_SC)
        print*, 'sum(DQRDT_SC_ref) = ', sum(DQRDT_SC_ref)
        
        print*, 'DIFF : sum(DQSDT_SC - DQSDT_SC_ref) = ', sum(DQSDT_SC - DQSDT_SC_ref)
        print*, 'sum(DQSDT_SC) = ', sum(DQSDT_SC)
        print*, 'sum(DQSDT_SC_ref) = ', sum(DQSDT_SC_ref)
        
        print*, 'DIFF : sum(CUFRC_SC - CUFRC_SC_ref) = ', sum(CUFRC_SC - CUFRC_SC_ref)
        print*, 'sum(CUFRC_SC) = ', sum(CUFRC_SC)
        print*, 'sum(CUFRC_SC_ref) = ', sum(CUFRC_SC_ref)
        
        print*, 'DIFF : sum(ENTR_SC - ENTR_SC_ref) = ', sum(ENTR_SC - ENTR_SC_ref)
        print*, 'sum(ENTR_SC) = ', sum(ENTR_SC)
        print*, 'sum(ENTR_SC_ref) = ', sum(ENTR_SC_ref)
        
        print*, 'DIFF : sum(DETR_SC - DETR_SC_ref) = ', sum(DETR_SC - DETR_SC_ref)
        print*, 'sum(DETR_SC) = ', sum(DETR_SC)
        print*, 'sum(DETR_SC_ref) = ', sum(DETR_SC_ref)
        
        print*, 'DIFF : sum(QLDET_SC - QLDET_SC_ref) = ', sum(QLDET_SC - QLDET_SC_ref)
        print*, 'sum(QLDET_SC) = ', sum(QLDET_SC)
        print*, 'sum(QLDET_SC_ref) = ', sum(QLDET_SC_ref)
        
        print*, 'DIFF : sum(QIDET_SC - QIDET_SC_ref) = ', sum(QIDET_SC - QIDET_SC_ref)
        print*, 'sum(QIDET_SC) = ', sum(QIDET_SC)
        print*, 'sum(QIDET_SC_ref) = ', sum(QIDET_SC_ref)
        
        print*, 'DIFF : sum(QLSUB_SC - QLSUB_SC_ref) = ', sum(QLSUB_SC - QLSUB_SC_ref)
        print*, 'sum(QLSUB_SC) = ', sum(QLSUB_SC)
        print*, 'sum(QLSUB_SC_ref) = ', sum(QLSUB_SC_ref)
        
        print*, 'DIFF : sum(QISUB_SC - QISUB_SC_ref) = ', sum(QISUB_SC - QISUB_SC_ref)
        print*, 'sum(QISUB_SC) = ', sum(QISUB_SC)
        print*, 'sum(QISUB_SC_ref) = ', sum(QISUB_SC_ref)
        
        print*, 'DIFF : sum(SC_NDROP - SC_NDROP_ref) = ', sum(SC_NDROP - SC_NDROP_ref)
        print*, 'sum(SC_NDROP) = ', sum(SC_NDROP)
        print*, 'sum(SC_NDROP_ref) = ', sum(SC_NDROP_ref)
        
        print*, 'DIFF : sum(SC_NICE - SC_NICE_ref) = ', sum(SC_NICE - SC_NICE_ref)
        print*, 'sum(SC_NICE) = ', sum(SC_NICE)
        print*, 'sum(SC_NICE_ref) = ', sum(SC_NICE_ref)
        
        print*, 'DIFF : sum(TPERT_SC - TPERT_SC_ref) = ', sum(TPERT_SC - TPERT_SC_ref)
        print*, 'sum(TPERT_SC) = ', sum(TPERT_SC)
        print*, 'sum(TPERT_SC_ref) = ', sum(TPERT_SC_ref)
        
        print*, 'DIFF : sum(QPERT_SC - QPERT_SC_ref) = ', sum(QPERT_SC - QPERT_SC_ref)
        print*, 'sum(QPERT_SC) = ', sum(QPERT_SC)
        print*, 'sum(QPERT_SC_ref) = ', sum(QPERT_SC_ref)
        
        print*, 'DIFF : sum(QTFLX_SC - QTFLX_SC_ref) = ', sum(QTFLX_SC - QTFLX_SC_ref)
        print*, 'sum(QTFLX_SC) = ', sum(QTFLX_SC)
        print*, 'sum(QTFLX_SC_ref) = ', sum(QTFLX_SC_ref)
        
        print*, 'DIFF : sum(SLFLX_SC - SLFLX_SC_ref) = ', sum(SLFLX_SC - SLFLX_SC_ref)
        print*, 'sum(SLFLX_SC) = ', sum(SLFLX_SC)
        print*, 'sum(SLFLX_SC_ref) = ', sum(SLFLX_SC_ref)
        
        print*, 'DIFF : sum(UFLX_SC - UFLX_SC_ref) = ', sum(UFLX_SC - UFLX_SC_ref)
        print*, 'sum(UFLX_SC) = ', sum(UFLX_SC)
        print*, 'sum(UFLX_SC_ref) = ', sum(UFLX_SC_ref)
        
        print*, 'DIFF : sum(VFLX_SC - VFLX_SC_ref) = ', sum(VFLX_SC - VFLX_SC_ref)
        print*, 'sum(VFLX_SC) = ', sum(VFLX_SC)
        print*, 'sum(VFLX_SC_ref) = ', sum(VFLX_SC_ref)

    endif
    print*,'Elapsed Time = ', (t_end - t_start)/rate

end program
