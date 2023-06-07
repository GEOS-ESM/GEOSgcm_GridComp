module test_GEOS5_subroutines

    use moist_subroutines_GEOS5
    use ConvPar_GF_SharedParams
    use timing_module

    implicit none

    public test_cup_gf_GEOS5, test_cup_gf_sh

    private
    
    character*2   :: i_string

    integer :: itf, ktf, its, ite, kts, kte, mtp, mxp, mzp, use_excess, ens4

    real :: dt, cum_entr_rate_shal, cum_entr_rate_deep

    integer, dimension(:), allocatable :: kpbli, ierr4d, ierr4d_ref, jmin4d, jmin4d_ref, klcl4d, klcl4d_ref, &
        k224d, k224d_ref, kbcon4d, kbcon4d_ref, ktop4d, ktop4d_ref, kstabi4d, kstabi4d_ref, &
        kstabm4d, kstabm4d_ref, ierr4d_mid, ierr4d_deep

    real, dimension(:), allocatable :: dx2d, h_sfc_flux, h_sfc_flux_ref, le_sfc_flux, le_sfc_flux_ref, & 
        tsur, tsur_ref, psur, psur_ref, ter11, ter11_ref, xlandi, xlandi_ref, ztexec, zqexec, ccn, ccn_ref, &
        zws, zws_ref, xmb4d, xmb4d_ref, edt4d, edt4d_ref, pwav4d, pwav4d_ref, stochastic_sig, xlons, xlons_ref, &
        xlats, xlats_ref, cprr4d, cprr4d_ref, &
        sigma4d, sigma4d_ref, col_sat, cum_ztexec, cum_ztexec_ref, cum_zqexec, cum_zqexec_ref, &
        AA0, AA0_ref, AA1, AA1_ref, AA2, AA2_ref, AA3, AA3_ref, AA1_BL, AA1_BL_ref, AA1_CIN, AA1_CIN_ref, TAU_BL, TAU_BL_ref, &
        TAU_EC, TAU_EC_ref

    real, dimension(:,:), allocatable :: rhoi, dhdt, zo, dm2d, temp_old, qv_old, temp_new_sh, qv_new_sh, &
        po, us, vs, outt, outt_ref, outq, outq_ref, outqc, outqc_ref, outu, outu_ref, outv, outv_ref, &
        pcup5d, pcup5d_ref, up_massentr5d, up_massentr5d_ref, up_massdetr5d, up_massdetr5d_ref, &
        dd_massentr5d, dd_massentr5d_ref, dd_massdetr5d, dd_massdetr5d_ref, &
        zup5d, zup5d_ref, zdn5d, zdn5d_ref, prup5d, prup5d_ref, &
        prdn5d, prdn5d_ref, clwup5d, clwup5d_ref, tup5d, tup5d_ref, &
        conv_cld_fr5d, conv_cld_fr5d_ref, temp_new_dp, qv_new_dp, &
        temp_new_BL, temp_new_BL_ref, qv_new_BL, qv_new_BL_ref, revsu_gf_2d, revsu_gf_2d_ref, &
        prfil_gf_2d, prfil_gf_2d_ref

    real, dimension(:,:,:), allocatable :: se_chem, out_chem, out_chem_ref, omeg, omeg_ref

    contains 

    subroutine test_cup_gf_sh(IM, JM, LM, dirName, rank_str)
        integer :: IM, JM, LM, fileID
        character*100 :: dirName, rank_str

        print*, 'Testing cup_gf_sh'

        itf = IM
        ktf = LM - 1
        its = 1
        ite = IM
        kts = 1
        kte = LM
        if (dirName(1:10) == './c24_data') then
            mtp = 39
            mxp = IM
            mzp = LM
        else if (dirName(1:10) == './c90_data') then
            mtp = 39
            mxp = IM
            mzp = LM
        else if (dirName(1:11) == './c180_data') then
            mtp = 39
            mxp = IM
            mzp = LM
        endif

        allocate(dx2d       (its:ite))
        allocate(kpbli      (its:ite))
        allocate(h_sfc_flux (its:ite))
        allocate(le_sfc_flux(its:ite))
        allocate(tsur       (its:ite))
        allocate(psur       (its:ite))
        allocate(ter11      (its:ite))
        allocate(xlandi     (its:ite))
        allocate(ztexec     (its:ite))
        allocate(zqexec     (its:ite))
        allocate(ccn        (its:ite))
        allocate(zws        (its:ite))

        allocate(rhoi       (its:ite, kts:kte))
        allocate(dhdt       (its:ite, kts:kte))
        allocate(zo         (its:ite, kts:kte))
        allocate(dm2d       (its:ite, kts:kte))
        allocate(temp_old   (its:ite, kts:kte))
        allocate(qv_old     (its:ite, kts:kte))
        allocate(temp_new_sh(its:ite, kts:kte))
        allocate(qv_new_sh  (its:ite, kts:kte))
        allocate(po         (its:ite, kts:kte))
        allocate(us         (its:ite, kts:kte))
        allocate(vs         (its:ite, kts:kte))
        allocate(outt       (its:ite, kts:kte))
        allocate(outt_ref   (its:ite, kts:kte))
        allocate(outq       (its:ite, kts:kte))
        allocate(outq_ref   (its:ite, kts:kte))
        allocate(outqc      (its:ite, kts:kte))
        allocate(outqc_ref  (its:ite, kts:kte))
        allocate(outu       (its:ite, kts:kte))
        allocate(outu_ref   (its:ite, kts:kte))
        allocate(outv       (its:ite, kts:kte))
        allocate(outv_ref   (its:ite, kts:kte))

        allocate(se_chem     (mtp, its:ite, kts:kte))
        allocate(out_chem    (mtp, its:ite, kts:kte))
        allocate(out_chem_ref(mtp, its:ite, kts:kte))

        allocate(ierr4d      (mxp))
        allocate(ierr4d_ref  (mxp))
        allocate(jmin4d      (mxp))
        allocate(jmin4d_ref  (mxp))
        allocate(klcl4d      (mxp))
        allocate(klcl4d_ref  (mxp))
        allocate(k224d       (mxp))
        allocate(k224d_ref   (mxp))
        allocate(kbcon4d     (mxp))
        allocate(kbcon4d_ref (mxp))
        allocate(ktop4d      (mxp))
        allocate(ktop4d_ref  (mxp))
        allocate(kstabi4d    (mxp))
        allocate(kstabi4d_ref(mxp))
        allocate(kstabm4d    (mxp))
        allocate(kstabm4d_ref(mxp))
        allocate(xmb4d       (mxp))
        allocate(xmb4d_ref   (mxp))
        allocate(edt4d       (mxp))
        allocate(edt4d_ref   (mxp))
        allocate(pwav4d      (mxp))
        allocate(pwav4d_ref  (mxp))

        allocate(pcup5d           (mxp, mzp))
        allocate(pcup5d_ref       (mxp, mzp))
        allocate(up_massentr5d    (mxp, mzp))
        allocate(up_massentr5d_ref(mxp, mzp))
        allocate(up_massdetr5d    (mxp, mzp))
        allocate(up_massdetr5d_ref(mxp, mzp))
        allocate(dd_massentr5d    (mxp, mzp))
        allocate(dd_massentr5d_ref(mxp, mzp))
        allocate(dd_massdetr5d    (mxp, mzp))
        allocate(dd_massdetr5d_ref(mxp, mzp))
        allocate(zup5d            (mxp, mzp))
        allocate(zup5d_ref        (mxp, mzp))
        allocate(zdn5d            (mxp, mzp))
        allocate(zdn5d_ref        (mxp, mzp))
        allocate(prup5d           (mxp, mzp))
        allocate(prup5d_ref       (mxp, mzp))
        allocate(prdn5d           (mxp, mzp))
        allocate(prdn5d_ref       (mxp, mzp))
        allocate(clwup5d          (mxp, mzp))
        allocate(clwup5d_ref      (mxp, mzp))
        allocate(tup5d            (mxp, mzp))
        allocate(tup5d_ref        (mxp, mzp))
        allocate(conv_cld_fr5d    (mxp, mzp))
        allocate(conv_cld_fr5d_ref(mxp, mzp))

        open(newunit=fileID, file=trim(dirName) // "/use_excess_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) use_excess
        ! print*,'use_excess = ', use_excess
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dt_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dt
        ! print*,'dt = ', dt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/cum_entr_rate_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) cum_entr_rate_shal
        ! print*,'cum_entr_rate_shal = ', cum_entr_rate_shal
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dx2d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dx2d
        ! print*,'sum(dx2d) = ', sum(dx2d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/kpbli_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) kpbli
        ! print*,'sum(kpbli) ', sum(kpbli)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/h_sfc_flux_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) h_sfc_flux
        ! print*,'sum(h_sfc_flux) = ', sum(h_sfc_flux)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/le_sfc_flux_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) le_sfc_flux
        ! print*,'sum(le_sfc_flux) = ', sum(le_sfc_flux)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/tsur_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) tsur
        ! print*,'sum(tsur) = ', sum(tsur)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/psur_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) psur
        ! print*,'sum(psur) = ', sum(psur)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ter11_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ter11
        ! print*,'sum(ter11) = ', sum(ter11)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/xlandi_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) xlandi
        ! print*,'sum(xlandi) = ', sum(xlandi)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ztexec_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ztexec
        ! print*,'sum(ztexec) = ', sum(ztexec)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zqexec_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zqexec
        ! print*,'sum(zqexec) = ', sum(zqexec)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ccn_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ccn
        ! print*,'sum(ccn) = ', sum(ccn)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/rhoi_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) rhoi
        ! print*,'sum(rhoi) = ', sum(rhoi)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dhdt_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dhdt
        ! print*,'sum(dhdt) = ', sum(dhdt)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zws_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zws
        ! print*,'sum(zws) = ', sum(zws)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zo_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zo
        ! print*,'sum(zo) = ', sum(zo)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dm2d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dm2d
        ! print*,'sum(dm2d) = ', sum(dm2d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/temp_old_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) temp_old
        ! print*,'sum(temp_old) = ', sum(temp_old)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/qv_old_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) qv_old
        ! print*,'sum(qv_old) = ', sum(qv_old)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/temp_new_sh_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) temp_new_sh
        ! print*,'sum(temp_new_sh) = ', sum(temp_new_sh)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/qv_new_sh_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) qv_new_sh
        ! print*,'sum(qv_new_sh) = ', sum(qv_new_sh)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/po_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) po
        ! print*,'sum(po) = ', sum(po)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/us_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) us
        ! print*,'sum(us) = ', sum(us)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/vs_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) vs
        ! print*,'sum(vs) = ', sum(vs)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/se_chem_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) se_chem
        ! print*,'sum(se_chem) = ', sum(se_chem)
        close(fileID)

        call start_timing()

        CALL CUP_gf_sh(itf,ktf,its,ite, kts,kte, mtp        &
                      ,cumulus_type  (shal) &
                      ,closure_choice(shal) &
                    !   , cum_entr_rate(shal) &
                      ,cum_entr_rate_shal &
                    !   ,use_excess    (shal) &
                      , use_excess &
                      !input data
                      ,dt                                   &
                    !   ,dx2d(:,j)                            &
                      ,dx2d(:)                              &
                      ,kpbli                                &
                      , h_sfc_flux                          &
                      ,le_sfc_flux                          &
                      ,tsur                                 &
                      ,psur                                 &
                      ,ter11                                &
                      ,xlandi                               &
                      ,ztexec                               &
                      ,zqexec                               &
                      ,ccn                                  &
                      ,rhoi                                 &
                      ,dhdt                                 &
                      ,zws                                  &
                      ,zo                                   &
                      ,dm2d                                 &
                      ,temp_old                             &
                      ,qv_old                               &
                      ,temp_new_sh                          &
                      ,qv_new_sh                            &
                      ,po                                   &
                      ,us                                   &
                      ,vs                                   &
                      ,se_chem                              &
                      !- output data
                    !   ,outt                   (:,:,shal)    &
                    !   ,outq                   (:,:,shal)    &
                    !   ,outqc                  (:,:,shal)    &
                    !   ,outu                   (:,:,shal)    &
                    !   ,outv                   (:,:,shal)    &
                    !   ,out_chem             (:,:,:,shal)    &
                      ,outt                   (:,:)    &
                      ,outq                   (:,:)    &
                      ,outqc                  (:,:)    &
                      ,outu                   (:,:)    &
                      ,outv                   (:,:)    &
                      ,out_chem             (:,:,:)    &
                      !- for shallow convective transport
                    !   ,ierr4d               (:,j,shal)      &
                    !   ,jmin4d               (:,j,shal)      &
                    !   ,klcl4d               (:,j,shal)      &
                    !   ,k224d                (:,j,shal)      &
                    !   ,kbcon4d              (:,j,shal)      &
                    !   ,ktop4d               (:,j,shal)      &
                    !   ,kstabi4d             (:,j,shal)      &
                    !   ,kstabm4d             (:,j,shal)      &
                    !   ,xmb4d                (:,j,shal)      &
                    !   ,edt4d                (:,j,shal)      &
                    !   ,pwav4d               (:,j,shal)      &
                    !   ,pcup5d             (:,j,:,shal)      &
                    !   ,up_massentr5d      (:,j,:,shal)      &
                    !   ,up_massdetr5d      (:,j,:,shal)      &
                    !   ,dd_massentr5d      (:,j,:,shal)      &
                    !   ,dd_massdetr5d      (:,j,:,shal)      &
                    !   ,zup5d              (:,j,:,shal)      &
                    !   ,zdn5d              (:,j,:,shal)      &
                    !   ,prup5d             (:,j,:,shal)      &
                    !   ,prdn5d             (:,j,:,shal)      &
                    !   ,clwup5d            (:,j,:,shal)      &
                    !   ,tup5d              (:,j,:,shal)      &
                    !   ,conv_cld_fr5d      (:,j,:,shal)      &
                      ,ierr4d               (:)      &
                      ,jmin4d               (:)      &
                      ,klcl4d               (:)      &
                      ,k224d                (:)      &
                      ,kbcon4d              (:)      &
                      ,ktop4d               (:)      &
                      ,kstabi4d             (:)      &
                      ,kstabm4d             (:)      &
                      ,xmb4d                (:)      &
                      ,edt4d                (:)      &
                      ,pwav4d               (:)      &
                    !   ,pcup5d             (:,j,:,shal)      &
                    !   ,up_massentr5d      (:,j,:,shal)      &
                    !   ,up_massdetr5d      (:,j,:,shal)      &
                    !   ,dd_massentr5d      (:,j,:,shal)      &
                    !   ,dd_massdetr5d      (:,j,:,shal)      &
                    !   ,zup5d              (:,j,:,shal)      &
                    !   ,zdn5d              (:,j,:,shal)      &
                    !   ,prup5d             (:,j,:,shal)      &
                    !   ,prdn5d             (:,j,:,shal)      &
                    !   ,clwup5d            (:,j,:,shal)      &
                    !   ,tup5d              (:,j,:,shal)      &
                    !   ,conv_cld_fr5d      (:,j,:,shal)      &
                      ,pcup5d             (:,:)      &
                      ,up_massentr5d      (:,:)      &
                      ,up_massdetr5d      (:,:)      &
                      ,dd_massentr5d      (:,:)      &
                      ,dd_massdetr5d      (:,:)      &
                      ,zup5d              (:,:)      &
                      ,zdn5d              (:,:)      &
                      ,prup5d             (:,:)      &
                      ,prdn5d             (:,:)      &
                      ,clwup5d            (:,:)      &
                      ,tup5d              (:,:)      &
                      ,conv_cld_fr5d      (:,:)      &
                      )

        call end_timing()

        call print_timing()

        open(newunit=fileID, file=trim(dirName) // "/outt_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outt_ref
        ! print*,'sum(outt_ref) = ', sum(outt_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outt)) = ', sum(outt_ref-outt)
        ! print*,'sum(outt_ref) = ', sum(outt_ref)
        ! print*,'sum(outt) = ', sum(outt)

        open(newunit=fileID, file=trim(dirName) // "/outq_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outq_ref
        ! print*,'sum(outq_ref) = ', sum(outq_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outq)) = ', sum(outq_ref-outq)
        ! print*,'sum(outq_ref) = ', sum(outq_ref)
        ! print*,'sum(outq) = ', sum(outq)

        open(newunit=fileID, file=trim(dirName) // "/outqc_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outqc_ref
        ! print*,'sum(outqc_ref) = ', sum(outqc_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outqc)) = ', sum(outqc_ref-outqc)
        ! print*,'sum(outqc_ref) = ', sum(outqc_ref)
        ! print*,'sum(outqc) = ', sum(outqc)

        open(newunit=fileID, file=trim(dirName) // "/outu_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outu_ref
        ! print*,'sum(outu_ref) = ', sum(outu_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outu)) = ', sum(outu_ref-outu)
        ! print*,'sum(outu_ref) = ', sum(outu_ref)
        ! print*,'sum(outu)) = ', sum(outu)


        open(newunit=fileID, file=trim(dirName) // "/outv_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outv_ref
        ! print*,'sum(outv_ref) = ', sum(outv_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outv)) = ', sum(outv_ref-outv)
        ! print*,'sum(outv_ref) = ', sum(outv_ref)
        ! print*,'sum(outv) = ', sum(outv)

        open(newunit=fileID, file=trim(dirName) // "/out_chem_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) out_chem_ref
        ! print*,'sum(out_chem_ref) = ', sum(out_chem_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(out_chem)) = ', sum(out_chem_ref-out_chem)
        ! print*,'sum(out_chem_ref) = ', sum(out_chem_ref)
        ! print*,'sum(out_chem) = ', sum(out_chem)

        open(newunit=fileID, file=trim(dirName) // "/ierr4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) ierr4d_ref
        ! print*,'sum(ierr4d_ref) = ', sum(ierr4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(ierr4d)) = ', sum(ierr4d_ref-ierr4d)
        ! print*,'sum(ierr4d_ref) = ', sum(ierr4d_ref)
        ! print*,'sum(ierr4d) = ', sum(ierr4d)

        open(newunit=fileID, file=trim(dirName) // "/jmin4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) jmin4d_ref
        ! print*,'sum(jmin4d_ref) = ', sum(jmin4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(jmin4d)) = ', sum(jmin4d_ref-jmin4d)
        ! print*,'sum(jmin4d_ref) = ', sum(jmin4d_ref)
        ! print*,'sum(jmin4d) = ', sum(jmin4d)

        open(newunit=fileID, file=trim(dirName) // "/klcl4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) klcl4d_ref
        ! print*,'sum(klcl4d_ref) = ', sum(klcl4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(klcl4d)) = ', sum(klcl4d_ref-klcl4d)
        ! print*,'sum(klcl4d_ref) = ', sum(klcl4d_ref)
        ! print*,'sum(klcl4d) = ', sum(klcl4d)

        open(newunit=fileID, file=trim(dirName) // "/k224d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) k224d_ref
        ! print*,'sum(k224d_ref) = ', sum(k224d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(k224d)) = ', sum(k224d_ref-k224d)
        ! print*,'sum(k224d_ref) = ', sum(k224d_ref)
        ! print*,'sum(k224d) = ', sum(k224d)

        open(newunit=fileID, file=trim(dirName) // "/kbcon4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) kbcon4d_ref
        ! print*,'sum(kbcon4d_ref) = ', sum(kbcon4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(kbcon4d)) = ', sum(kbcon4d_ref-kbcon4d)
        ! print*,'sum(kbcon4d_ref) = ', sum(kbcon4d_ref)
        ! print*,'sum(kbcon4d) = ', sum(kbcon4d)

        open(newunit=fileID, file=trim(dirName) // "/ktop4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) ktop4d_ref
        ! print*,'sum(ktop4d_ref) = ', sum(ktop4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(ktop4d)) = ', sum(ktop4d_ref-ktop4d)
        ! print*,'sum(ktop4d_ref) = ', sum(ktop4d_ref)
        ! print*,'sum(ktop4d) = ', sum(ktop4d)

        open(newunit=fileID, file=trim(dirName) // "/kstabi4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) kstabi4d_ref
        ! print*,'sum(kstabi4d_ref) = ', sum(kstabi4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(kstabi4d)) = ', sum(kstabi4d_ref-kstabi4d)
        ! print*,'sum(kstabi4d_ref) = ', sum(kstabi4d_ref)
        ! print*,'sum(kstabi4d) = ', sum(kstabi4d)

        open(newunit=fileID, file=trim(dirName) // "/kstabm4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) kstabm4d_ref
        ! print*,'sum(kstabm4d_ref) = ', sum(kstabm4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(kstabm4d)) = ', sum(kstabm4d_ref-kstabm4d)
        ! print*,'sum(kstabm4d_ref) = ', sum(kstabm4d_ref)
        ! print*,'sum(kstabm4d) = ', sum(kstabm4d)

        open(newunit=fileID, file=trim(dirName) // "/xmb4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) xmb4d_ref
        ! print*,'sum(xmb4d_ref) = ', sum(xmb4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(xmb4d)) = ', sum(xmb4d_ref-xmb4d)
        ! print*,'sum(xmb4d_ref) = ', sum(xmb4d_ref)
        ! print*,'sum(xmb4d) = ', sum(xmb4d)

        open(newunit=fileID, file=trim(dirName) // "/edt4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) edt4d_ref
        ! print*,'sum(edt4d_ref) = ', sum(edt4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(edt4d)) = ', sum(edt4d_ref-edt4d)
        ! print*,'sum(edt4d_ref) = ', sum(edt4d_ref)
        ! print*,'sum(edt4d) = ', sum(edt4d)


        open(newunit=fileID, file=trim(dirName) // "/pwav4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) pwav4d_ref
        ! print*,'sum(pwav4d_ref) = ', sum(pwav4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(pwav4d)) = ', sum(pwav4d_ref-pwav4d)
        ! print*,'sum(pwav4d_ref) = ', sum(pwav4d_ref)
        ! print*,'sum(pwav4d) = ', sum(pwav4d)

        open(newunit=fileID, file=trim(dirName) // "/pcup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) pcup5d_ref
        ! print*,'sum(pcup5d_ref) = ', sum(pcup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(pcup5d)) = ', sum(pcup5d_ref-pcup5d)
        ! print*,'sum(pcup5d_ref) = ', sum(pcup5d_ref)
        ! print*,'sum(pcup5d) = ', sum(pcup5d)

        open(newunit=fileID, file=trim(dirName) // "/up_massentr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) up_massentr5d_ref
        ! print*,'sum(up_massentr5d_ref) = ', sum(up_massentr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(up_massentr5d)) = ', sum(up_massentr5d_ref-up_massentr5d)
        ! print*,'sum(up_massentr5d_ref) = ', sum(up_massentr5d_ref)
        ! print*,'sum(up_massentr5d) = ', sum(up_massentr5d)

        open(newunit=fileID, file=trim(dirName) // "/up_massdetr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) up_massdetr5d_ref
        ! print*,'sum(up_massdetr5d_ref) = ', sum(up_massdetr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(up_massdetr5d)) = ', sum(up_massdetr5d_ref-up_massdetr5d)
        ! print*,'sum(up_massdetr5d_ref) = ', sum(up_massdetr5d_ref)
        ! print*,'sum(up_massdetr5d) = ', sum(up_massdetr5d)

        open(newunit=fileID, file=trim(dirName) // "/dd_massentr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) dd_massentr5d_ref
        ! print*,'sum(dd_massentr5d_ref) = ', sum(dd_massentr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(dd_massentr5d)) = ', sum(dd_massentr5d_ref-dd_massentr5d)
        ! print*,'sum(dd_massentr5d_ref) = ', sum(dd_massentr5d_ref)
        ! print*,'sum(dd_massentr5d) = ', sum(dd_massentr5d)

        open(newunit=fileID, file=trim(dirName) // "/dd_massdetr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) dd_massdetr5d_ref
        ! print*,'sum(dd_massdetr5d_ref) = ', sum(dd_massdetr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(dd_massdetr5d)) = ', sum(dd_massdetr5d_ref-dd_massdetr5d)
        ! print*,'sum(dd_massdetr5d_ref) = ', sum(dd_massdetr5d_ref)
        ! print*,'sum(dd_massdetr5d) = ', sum(dd_massdetr5d)

        open(newunit=fileID, file=trim(dirName) // "/zup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) zup5d_ref
        ! print*,'sum(zup5d_ref) = ', sum(zup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(zup5d)) = ', sum(zup5d_ref-zup5d)
        ! print*,'sum(zup5d_ref) = ', sum(zup5d_ref)
        ! print*,'sum(zup5d) = ', sum(zup5d)

        open(newunit=fileID, file=trim(dirName) // "/zdn5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) zdn5d_ref
        ! print*,'sum(zdn5d_ref) = ', sum(zdn5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(zdn5d)) = ', sum(zdn5d_ref-zdn5d)
        ! print*,'sum(zdn5d_ref) = ', sum(zdn5d_ref)
        ! print*,'sum(zdn5d) = ', sum(zdn5d)

        open(newunit=fileID, file=trim(dirName) // "/prup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) prup5d_ref
        ! print*,'sum(prup5d_ref) = ', sum(prup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(prup5d)) = ', sum(prup5d_ref-prup5d)
        ! print*,'sum(prup5d_ref) = ', sum(prup5d_ref)
        ! print*,'sum(prup5d) = ', sum(prup5d)

        open(newunit=fileID, file=trim(dirName) // "/prdn5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) prdn5d_ref
        ! print*,'sum(prdn5d_ref) = ', sum(prdn5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(prdn5d)) = ', sum(prdn5d_ref-prdn5d)
        ! print*,'sum(prdn5d_ref) = ', sum(prdn5d_ref)
        ! print*,'sum(prdn5d) = ', sum(prdn5d)

        open(newunit=fileID, file=trim(dirName) // "/clwup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) clwup5d_ref
        ! print*,'sum(clwup5d_ref) = ', sum(clwup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(clwup5d)) = ', sum(clwup5d_ref-clwup5d)
        ! print*,'sum(clwup5d_ref) = ', sum(clwup5d_ref)
        ! print*,'sum(clwup5d) = ', sum(clwup5d)

        open(newunit=fileID, file=trim(dirName) // "/tup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) tup5d_ref
        ! print*,'sum(tup5d_ref) = ', sum(tup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(tup5d)) = ', sum(tup5d_ref-tup5d)
        ! print*,'sum(tup5d_ref) = ', sum(tup5d_ref)
        ! print*,'sum(tup5d) = ', sum(tup5d)

        open(newunit=fileID, file=trim(dirName) // "/conv_cld_fr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) conv_cld_fr5d_ref
        ! print*,'sum(conv_cld_fr5d_ref) = ', sum(conv_cld_fr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(conv_cld_fr5d)) = ', sum(conv_cld_fr5d_ref-conv_cld_fr5d)
        ! print*,'sum(conv_cld_fr5d_ref) = ', sum(conv_cld_fr5d_ref)
        ! print*,'sum(conv_cld_fr5d) = ', sum(conv_cld_fr5d)

    end subroutine

    subroutine test_cup_gf_GEOS5(IM, JM, LM, dirName, rank_str)
        integer :: IM, JM, LM, fileID, ii
        character*100 :: dirName, rank_str

        print*, 'Testing cup_gf (ver GEOS5)'

        ens4 = 1

        itf = IM
        ktf = LM - 1
        its = 1
        ite = IM
        kts = 1
        kte = LM
        if (dirName(1:10) == './c24_data') then
            mtp = 39
            mxp = IM
            mzp = LM
        else if (dirName(1:10) == './c90_data') then
            mtp = 39
            mxp = IM
            mzp = LM
        else if (dirName(1:11) == './c180_data') then
            mtp = 39
            mxp = IM
            mzp = LM
        endif

        ! Size of CNV_Tracers seem to match to mpt.  Note that this may change later
        allocate(CNV_Tracers(mtp))

        do ii = 1,mtp
            if(ii.ge.10) then
                write(i_string,'(i2)') ii
            else
                write(i_string,'(i1)') ii
            endif

            allocate(CNV_Tracers(ii)%Q(IM, JM, LM))
            open(newunit=fileID, file=trim(dirName) // "/CNV_Tracers_Q_" // trim(i_string) // "_" // trim(rank_str) // ".in", &
                 status='old', form="unformatted", action="read")
            read(fileID) CNV_Tracers(ii)%Q
            ! print*,'Rank = ', trim(rank_str), ': Index = ', ii, 'sum(CNV_Tracers%Q) = ', sum(CNV_Tracers(ii)%Q)
            close(fileID)

            open(newunit=fileID, file=trim(dirName) // "/CNV_Tracers_fscav_" // trim(i_string) // "_" // trim(rank_str) // ".in", &
                 status='old', form="unformatted", action="read")
            read(fileID) CNV_Tracers(ii)%fscav
            ! print*,'Rank = ', trim(rank_str), ': Index = ', ii, 'sum(CNV_Tracers%fscav) = ', CNV_Tracers(ii)%fscav
            close(fileID)

            open(newunit=fileID, file=trim(dirName) // "/CNV_Tracers_Vect_Hcts_" // trim(i_string) // "_" // trim(rank_str) // ".in", &
                 status='old', form="unformatted", action="read")
            read(fileID) CNV_Tracers(ii)%Vect_Hcts
            ! print*,'Rank = ', trim(rank_str), ': Index = ', ii, 'sum(CNV_Tracers%Vect_Hcts) = ', sum(CNV_Tracers(ii)%Vect_Hcts)
            close(fileID)

            open(newunit=fileID, file=trim(dirName) // "/CNV_Tracers_QNAME_" // trim(i_string) // "_" // trim(rank_str) // ".in", &
                 status='old', form="unformatted", action="read")
            read(fileID) CNV_Tracers(ii)%QNAME
            ! print*,'Rank = ', trim(rank_str), ': Index = ', ii, 'CNV_Tracers%QNAME = ', CNV_Tracers(ii)%QNAME
            close(fileID)

            open(newunit=fileID, file=trim(dirName) // "/CNV_Tracers_CNAME_" // trim(i_string) // "_" // trim(rank_str) // ".in", &
                 status='old', form="unformatted", action="read")
            read(fileID) CNV_Tracers(ii)%CNAME
            ! print*,'Rank = ', trim(rank_str), ': Index = ', ii, 'CNV_Tracers%CNAME = ', CNV_Tracers(ii)%CNAME
            close(fileID)
        enddo

        allocate(dx2d           (its:ite))
        allocate(stochastic_sig (its:ite))
        allocate(kpbli      (its:ite))
        allocate(ztexec     (its:ite))
        allocate(zqexec     (its:ite))
        allocate(ccn        (its:ite))
        allocate(ter11      (its:ite))
        allocate(h_sfc_flux (its:ite))
        allocate(le_sfc_flux(its:ite))
        allocate(xlons      (its:ite))
        allocate(xlats      (its:ite))
        allocate(xlandi     (its:ite))
        allocate(tsur       (its:ite))
        allocate(psur       (its:ite))
        allocate(zws        (its:ite))

        allocate(rhoi       (its:ite, kts:kte))
        allocate(temp_old   (its:ite, kts:kte))
        allocate(qv_old     (its:ite, kts:kte))
        allocate(temp_new_dp (its:ite, kts:kte))
        allocate(qv_new_dp (its:ite, kts:kte))
        allocate(temp_new_BL (its:ite, kts:kte))
        allocate(qv_new_BL (its:ite, kts:kte))
        allocate(zo         (its:ite, kts:kte))
        allocate(po         (its:ite, kts:kte))
        allocate(us         (its:ite, kts:kte))
        allocate(vs         (its:ite, kts:kte))
        allocate(dm2d       (its:ite, kts:kte))
        allocate(dhdt       (its:ite, kts:kte))
        allocate(revsu_gf_2d(its:ite, kts:kte))
        allocate(revsu_gf_2d_ref(its:ite, kts:kte))
        allocate(prfil_gf_2d(its:ite, kts:kte))
        allocate(prfil_gf_2d_ref(its:ite, kts:kte))
        allocate(outt       (its:ite, kts:kte))
        allocate(outt_ref   (its:ite, kts:kte))
        allocate(outq       (its:ite, kts:kte))
        allocate(outq_ref   (its:ite, kts:kte))
        allocate(outqc      (its:ite, kts:kte))
        allocate(outqc_ref  (its:ite, kts:kte))
        allocate(outu       (its:ite, kts:kte))
        allocate(outu_ref   (its:ite, kts:kte))
        allocate(outv       (its:ite, kts:kte))
        allocate(outv_ref   (its:ite, kts:kte))
        allocate(out_chem    (mtp, its:ite, kts:kte))
        allocate(out_chem_ref(mtp, its:ite, kts:kte))
        
        allocate(omeg(its:ite, kts:kte, 1:ens4))
        allocate(se_chem     (mtp, its:ite, kts:kte))

        allocate(ierr4d_mid      (mxp))
        allocate(ierr4d_deep      (mxp))
        allocate(jmin4d      (mxp))
        allocate(jmin4d_ref  (mxp))
        allocate(klcl4d      (mxp))
        allocate(klcl4d_ref  (mxp))
        allocate(k224d       (mxp))
        allocate(k224d_ref   (mxp))
        allocate(kbcon4d     (mxp))
        allocate(kbcon4d_ref (mxp))
        allocate(ktop4d      (mxp))
        allocate(ktop4d_ref  (mxp))
        allocate(kstabi4d    (mxp))
        allocate(kstabi4d_ref(mxp))
        allocate(kstabm4d    (mxp))
        allocate(kstabm4d_ref(mxp))
        allocate(cprr4d    (mxp))
        allocate(cprr4d_ref    (mxp))
        allocate(xmb4d       (mxp))
        allocate(xmb4d_ref   (mxp))
        allocate(edt4d       (mxp))
        allocate(edt4d_ref   (mxp))
        allocate(pwav4d      (mxp))
        allocate(pwav4d_ref  (mxp))
        allocate(sigma4d     (mxp))
        allocate(sigma4d_ref     (mxp))
        allocate(ierr4d_ref  (mxp))
        allocate(AA0(mxp))
        allocate(AA1(mxp))
        allocate(AA2(mxp))
        allocate(AA3(mxp))
        allocate(AA1_BL(mxp))
        allocate(AA1_CIN(mxp))
        allocate(TAU_BL(mxp))
        allocate(TAU_EC(mxp))


        allocate(pcup5d           (mxp, mzp))
        allocate(pcup5d_ref       (mxp, mzp))
        allocate(up_massentr5d    (mxp, mzp))
        allocate(up_massentr5d_ref(mxp, mzp))
        allocate(up_massdetr5d    (mxp, mzp))
        allocate(up_massdetr5d_ref(mxp, mzp))
        allocate(dd_massentr5d    (mxp, mzp))
        allocate(dd_massentr5d_ref(mxp, mzp))
        allocate(dd_massdetr5d    (mxp, mzp))
        allocate(dd_massdetr5d_ref(mxp, mzp))
        allocate(zup5d            (mxp, mzp))
        allocate(zup5d_ref        (mxp, mzp))
        allocate(zdn5d            (mxp, mzp))
        allocate(zdn5d_ref        (mxp, mzp))
        allocate(prup5d           (mxp, mzp))
        allocate(prup5d_ref       (mxp, mzp))
        allocate(prdn5d           (mxp, mzp))
        allocate(prdn5d_ref       (mxp, mzp))
        allocate(clwup5d          (mxp, mzp))
        allocate(clwup5d_ref      (mxp, mzp))
        allocate(tup5d            (mxp, mzp))
        allocate(tup5d_ref        (mxp, mzp))
        allocate(conv_cld_fr5d    (mxp, mzp))
        allocate(conv_cld_fr5d_ref(mxp, mzp))

        open(newunit=fileID, file=trim(dirName) // "/use_excess_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) use_excess
        ! print*,'use_excess = ', use_excess
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/cum_entr_rate_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) cum_entr_rate_deep
        ! print*,'cum_entr_rate_shal = ', cum_entr_rate_shal
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dx2d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dx2d
        ! print*,'sum(dx2d) = ', sum(dx2d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/stochastic_sig_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) stochastic_sig
        ! print*,'sum(stochastic_sig) = ', sum(stochastic_sig)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dt_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dt
        ! print*,'dt = ', dt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/kpbli_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) kpbli
        ! print*,'sum(kpbli) ', sum(kpbli)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ztexec_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ztexec
        ! print*,'sum(ztexec) = ', sum(ztexec)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zqexec_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zqexec
        ! print*,'sum(zqexec) = ', sum(zqexec)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ccn_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ccn
        ! print*,'sum(ccn) = ', sum(ccn)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/rhoi_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) rhoi
        ! print*,'sum(rhoi) = ', sum(rhoi)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/omeg_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) omeg
        ! print*,'sum(omeg) = ', sum(omeg)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/temp_old_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) temp_old
        ! print*,'sum(temp_old) = ', sum(temp_old)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/qv_old_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) qv_old
        ! print*,'sum(qv_old) = ', sum(qv_old)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ter11_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ter11
        ! print*,'sum(ter11) = ', sum(ter11)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/h_sfc_flux_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) h_sfc_flux
        ! print*,'sum(h_sfc_flux) = ', sum(h_sfc_flux)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/le_sfc_flux_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) le_sfc_flux
        ! print*,'sum(le_sfc_flux) = ', sum(le_sfc_flux)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/xlons_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) xlons
        ! print*,'sum(xlons) = ', sum(xlons)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/xlats_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) xlats
        ! print*,'sum(xlats) = ', sum(xlats)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/xlandi_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) xlandi
        ! print*,'sum(xlandi) = ', sum(xlandi)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/temp_new_dp_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) temp_new_dp
        ! print*,'sum(temp_new_dp) = ', sum(temp_new_dp)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/qv_new_dp_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) qv_new_dp
        ! print*,'sum(qv_new_dp) = ', sum(qv_new_dp)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/temp_new_BL_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) temp_new_BL
        ! print*,'sum(temp_new_BL) = ', sum(temp_new_BL)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/qv_new_BL_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) qv_new_BL
        ! print*,'sum(qv_new_BL) = ', sum(qv_new_BL)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zo_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zo
        ! print*,'sum(zo) = ', sum(zo)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/po_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) po
        ! print*,'sum(po) = ', sum(po)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/tsur_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) tsur
        ! print*,'sum(tsur) = ', sum(tsur)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/psur_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) psur
        ! print*,'sum(psur) = ', sum(psur)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/us_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) us
        ! print*,'sum(us) = ', sum(us)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/vs_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) vs
        ! print*,'sum(vs) = ', sum(vs)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dm2d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dm2d
        ! print*,'sum(dm2d) = ', sum(dm2d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/se_chem_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) se_chem
        ! print*,'sum(se_chem) = ', sum(se_chem)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zws_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zws
        ! print*,'sum(zws) = ', sum(zws)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dhdt_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dhdt
        ! print*,'sum(dhdt) = ', sum(dhdt)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ierr4d_mid_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ierr4d_mid
        ! print*,'sum(ierr4d_mid) = ', sum(ierr4d_mid)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ierr4d_deep_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ierr4d_deep
        ! print*,'sum(ierr4d_deep) = ', sum(ierr4d_deep)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/jmin4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) jmin4d
        ! print*,'sum(jmin4d) = ', sum(jmin4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/klcl4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) klcl4d
        ! print*,'sum(klcl4d) = ', sum(klcl4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/k224d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) k224d
        ! print*,'sum(k224d) = ', sum(k224d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/kbcon4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) kbcon4d
        ! print*,'sum(kbcon4d) = ', sum(kbcon4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/ktop4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ktop4d
        ! print*,'sum(ktop4d) = ', sum(ktop4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/kstabi4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) kstabi4d
        ! print*,'sum(kstabi4d) = ', sum(kstabi4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/kstabm4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) kstabm4d
        ! print*,'sum(kstabm4d) = ', sum(kstabm4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/cprr4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) cprr4d
        ! print*,'sum(cprr4d) = ', sum(cprr4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/xmb4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) xmb4d
        ! print*,'sum(xmb4d) = ', sum(xmb4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/edt4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) edt4d
        ! print*,'sum(edt4d) = ', sum(edt4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/pwav4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) pwav4d
        ! print*,'sum(pwav4d) = ', sum(pwav4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/sigma4d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) sigma4d
        ! print*,'sum(sigma4d) = ', sum(sigma4d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/pcup5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) pcup5d
        ! print*,'sum(pcup5d) = ', sum(pcup5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/up_massentr5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) up_massentr5d
        ! print*,'sum(up_massentr5d) = ', sum(up_massentr5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/up_massdetr5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) up_massdetr5d
        ! print*,'sum(up_massdetr5d) = ', sum(up_massdetr5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dd_massentr5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dd_massentr5d
        ! print*,'sum(dd_massentr5d) = ', sum(dd_massentr5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/dd_massdetr5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) dd_massdetr5d
        ! print*,'sum(dd_massdetr5d) = ', sum(dd_massdetr5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zup5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zup5d
        ! print*,'sum(zup5d) = ', sum(zup5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/zdn5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) zdn5d
        ! print*,'sum(zdn5d) = ', sum(zdn5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/prup5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) prup5d
        ! print*,'sum(prup5d) = ', sum(prup5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/clwup5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) clwup5d
        ! print*,'sum(clwup5d) = ', sum(clwup5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/tup5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) tup5d
        ! print*,'sum(tup5d) = ', sum(tup5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/conv_cld_fr5d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) conv_cld_fr5d
        ! print*,'sum(conv_cld_fr5d) = ', sum(conv_cld_fr5d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/revsu_gf_2d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) revsu_gf_2d
        ! print*,'sum(revsu_gf_2d) = ', sum(revsu_gf_2d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/prfil_gf_2d_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) prfil_gf_2d
        ! print*,'sum(prfil_gf_2d) = ', sum(prfil_gf_2d)
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/c1_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) c1
        ! print*,'c1 = ', c1
        close(fileID)

        call start_timing()

        CALL CUP_gf_GEOS5(its,ite,kts,kte, itf,ktf , mtp    &
                  ,cumulus_type   (deep)             &
                  ,closure_choice (deep)             &
                !   ,cum_entr_rate  (deep)             &
                  ,cum_entr_rate_deep                &
                !   ,use_excess     (deep)             &
                  ,use_excess                        &
                  ! input data
                !   ,dx2d(:,j)                         &
                  ,dx2d(:)                            &
                !   ,stochastic_sig(:,j)               &
                  ,stochastic_sig(:)                 &
                  ,dt                                &
                  ,kpbli                             &
                  ,ztexec                            &
                  ,zqexec                            &
                  ,ccn                               &
                  ,rhoi                              &
                  ,omeg                              &
                  ,temp_old                          &
                  ,qv_old                            &
                  ,ter11                             &
                  , h_sfc_flux                       &
                  ,le_sfc_flux                       &
                  ,xlons                             &
                  ,xlats                             &
                  ,xlandi                            &
                  ,temp_new_dp                       &
                  ,qv_new_dp                         &
                  ,temp_new_BL                       &
                  ,qv_new_BL                         &
                  ,zo                                &
                  ,po                                &
                  ,tsur                              &
                  ,psur                              &
                  ,us                                &
                  ,vs                                &
		          ,dm2d                              &
                  ,se_chem                           &
                  ,zws                               &
                  ,dhdt                              &
                !   ,ierr4d               (:,j, mid)   &
                  ,ierr4d_mid               (:)      &
                  !output data
                !   ,outt                 (:,:,deep)   &
                !   ,outq                 (:,:,deep)   &
                !   ,outqc                (:,:,deep)   &
                !   ,outu                 (:,:,deep)   &,j
                !   ,outv                 (:,:,deep)   &
                !   ,out_chem           (:,:,:,deep)   &
                  ,outt                 (:,:)        &
                  ,outq                 (:,:)        &
                  ,outqc                (:,:)        &
                  ,outu                 (:,:)        &
                  ,outv                 (:,:)        &
                  ,out_chem           (:,:,:)        &
                  !
                !   ,ierr4d               (:,j,deep)   &
                !   ,jmin4d               (:,j,deep)   &
                !   ,klcl4d               (:,j,deep)   &
                !   ,k224d                (:,j,deep)   &
                !   ,kbcon4d              (:,j,deep)   &
                !   ,ktop4d               (:,j,deep)   &
                !   ,kstabi4d             (:,j,deep)   &
                !   ,kstabm4d             (:,j,deep)   &
                !   ,cprr4d               (:,j,deep)   &
                !   ,xmb4d                (:,j,deep)   &
                !   ,edt4d                (:,j,deep)   &
                !   ,pwav4d               (:,j,deep)   &
                !   ,sigma4d              (:,j,deep)   &
                !   ,pcup5d             (:,j,:,deep)   &
                !   ,up_massentr5d      (:,j,:,deep)   &
                !   ,up_massdetr5d      (:,j,:,deep)   &
                !   ,dd_massentr5d      (:,j,:,deep)   &
                !   ,dd_massdetr5d      (:,j,:,deep)   &
                !   ,zup5d              (:,j,:,deep)   &
                !   ,zdn5d              (:,j,:,deep)   &
                !   ,prup5d             (:,j,:,deep)   &
                !   ,prdn5d             (:,j,:,deep)   &
                !   ,clwup5d            (:,j,:,deep)   &
                !   ,tup5d              (:,j,:,deep)   &
                !   ,conv_cld_fr5d      (:,j,:,deep)   &
                  ,ierr4d_deep          (:)   &
                  ,jmin4d               (:)   &
                  ,klcl4d               (:)   &
                  ,k224d                (:)   &
                  ,kbcon4d              (:)   &
                  ,ktop4d               (:)   &
                  ,kstabi4d             (:)   &
                  ,kstabm4d             (:)   &
                  ,cprr4d               (:)   &
                  ,xmb4d                (:)   &
                  ,edt4d                (:)   &
                  ,pwav4d               (:)   &
                  ,sigma4d              (:)   &
                  ,pcup5d             (:,:)   &
                  ,up_massentr5d      (:,:)   &
                  ,up_massdetr5d      (:,:)   &
                  ,dd_massentr5d      (:,:)   &
                  ,dd_massdetr5d      (:,:)   &
                  ,zup5d              (:,:)   &
                  ,zdn5d              (:,:)   &
                  ,prup5d             (:,:)   &
                  ,prdn5d             (:,:)   &
                  ,clwup5d            (:,:)   &
                  ,tup5d              (:,:)   &
                  ,conv_cld_fr5d      (:,:)   &
                  !-- for debug/diag
                  ,AA0(:),AA1(:),AA2(:),AA3(:),AA1_BL(:),AA1_CIN(:),TAU_BL(:),TAU_EC(:) &
                  ,revsu_gf_2d                      &
                  ,prfil_gf_2d                      &
                  )

        call end_timing()

        call print_timing()

        open(newunit=fileID, file=trim(dirName) // "/outt_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outt_ref
        ! print*,'sum(outt_ref) = ', sum(outt_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outt)) = ', sum(outt_ref-outt)
        ! print*,'sum(outt_ref) = ', sum(outt_ref)
        ! print*,'sum(outt) = ', sum(outt)

        open(newunit=fileID, file=trim(dirName) // "/outq_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outq_ref
        ! print*,'sum(outq_ref) = ', sum(outq_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outq)) = ', sum(outq_ref-outq)
        ! print*,'sum(outq_ref) = ', sum(outq_ref)
        ! print*,'sum(outq) = ', sum(outq)

        open(newunit=fileID, file=trim(dirName) // "/outqc_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outqc_ref
        ! print*,'sum(outqc_ref) = ', sum(outqc_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outqc)) = ', sum(outqc_ref-outqc)
        ! print*,'sum(outqc_ref) = ', sum(outqc_ref)
        ! print*,'sum(outqc) = ', sum(outqc)

        open(newunit=fileID, file=trim(dirName) // "/outu_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outu_ref
        ! print*,'sum(outu_ref) = ', sum(outu_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outu)) = ', sum(outu_ref-outu)
        ! print*,'sum(outu_ref) = ', sum(outu_ref)
        ! print*,'sum(outu) = ', sum(outu)

        open(newunit=fileID, file=trim(dirName) // "/outv_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) outv_ref
        ! print*,'sum(outv_ref) = ', sum(outv_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(outv)) = ', sum(outv_ref-outv)
        ! print*,'sum(outv_ref) = ', sum(outv_ref)
        ! print*,'sum(outv) = ', sum(outv)

        open(newunit=fileID, file=trim(dirName) // "/out_chem_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) out_chem_ref
        ! print*,'sum(out_chem_ref) = ', sum(out_chem_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(out_chem_ref)) = ', sum(out_chem_ref-out_chem)
        ! print*,'sum(out_chem_ref) = ', sum(out_chem_ref)
        ! print*,'sum(out_chem) = ', sum(out_chem)
        print*,'Note : out_chem_ref reference will not match up to computed'

        open(newunit=fileID, file=trim(dirName) // "/ierr4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) ierr4d_ref
        ! print*,'sum(ierr4d_ref) = ', sum(ierr4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(ierr4d)) = ', sum(ierr4d_ref-ierr4d_deep)
        ! print*,'sum(ierr4d_ref) = ', sum(ierr4d_ref)
        ! print*,'sum(ierr4d) = ', sum(ierr4d_deep)

        open(newunit=fileID, file=trim(dirName) // "/jmin4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) jmin4d_ref
        ! print*,'sum(jmin4d_ref) = ', sum(jmin4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(jmin4d)) = ', sum(jmin4d_ref-jmin4d)
        ! print*,'sum(jmin4d_ref) = ', sum(jmin4d_ref)
        ! print*,'sum(jmin4d) = ', sum(jmin4d)

        open(newunit=fileID, file=trim(dirName) // "/klcl4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) klcl4d_ref
        ! print*,'sum(klcl4d_ref) = ', sum(klcl4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(klcl4d)) = ', sum(klcl4d_ref-klcl4d)
        ! print*,'sum(klcl4d_ref) = ', sum(klcl4d_ref)
        ! print*,'sum(klcl4d) = ', sum(klcl4d)

        open(newunit=fileID, file=trim(dirName) // "/k224d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) k224d_ref
        ! print*,'sum(k224d_ref) = ', sum(k224d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(k224d)) = ', sum(k224d_ref-k224d)
        ! print*,'sum(k224d_ref) = ', sum(k224d_ref)
        ! print*,'sum(k224d) = ', sum(k224d)

        open(newunit=fileID, file=trim(dirName) // "/kbcon4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) kbcon4d_ref
        ! print*,'sum(kbcon4d_ref) = ', sum(kbcon4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(kbcon4d)) = ', sum(kbcon4d_ref-kbcon4d)
        ! print*,'sum(kbcon4d_ref) = ', sum(kbcon4d_ref)
        ! print*,'sum(kbcon4d) = ', sum(kbcon4d)

        open(newunit=fileID, file=trim(dirName) // "/ktop4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) ktop4d_ref
        ! print*,'sum(ktop4d_ref) = ', sum(ktop4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(ktop4d)) = ', sum(ktop4d_ref-ktop4d)
        ! print*,'sum(ktop4d_ref) = ', sum(ktop4d_ref)
        ! print*,'sum(ktop4d) = ', sum(ktop4d)

        open(newunit=fileID, file=trim(dirName) // "/kstabi4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) kstabi4d_ref
        ! print*,'sum(kstabi4d_ref) = ', sum(kstabi4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(kstabi4d)) = ', sum(kstabi4d_ref-kstabi4d)
        ! print*,'sum(kstabi4d_ref) = ', sum(kstabi4d_ref)
        ! print*,'sum(kstabi4d) = ', sum(kstabi4d)

        open(newunit=fileID, file=trim(dirName) // "/kstabm4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) kstabm4d_ref
        ! print*,'sum(kstabm4d_ref) = ', sum(kstabm4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(kstabm4d)) = ', sum(kstabm4d_ref-kstabm4d)
        ! print*,'sum(kstabm4d_ref) = ', sum(kstabm4d_ref)
        ! print*,'sum(kstabm4d) = ', sum(kstabm4d)

        open(newunit=fileID, file=trim(dirName) // "/cprr4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) cprr4d_ref
        ! print*,'sum(cprr4d_ref) = ', sum(cprr4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(cprr4d)) = ', sum(cprr4d_ref-cprr4d)
        ! print*,'sum(cprr4d_ref) = ', sum(cprr4d_ref)
        ! print*,'sum(cprr4d) = ', sum(cprr4d)

        open(newunit=fileID, file=trim(dirName) // "/xmb4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) xmb4d_ref
        ! print*,'sum(xmb4d_ref) = ', sum(xmb4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(xmb4d)) = ', sum(xmb4d_ref-xmb4d)
        ! print*,'sum(xmb4d_ref) = ', sum(xmb4d_ref)
        ! print*,'sum(xmb4d) = ', sum(xmb4d)

        open(newunit=fileID, file=trim(dirName) // "/edt4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) edt4d_ref
        ! print*,'sum(edt4d_ref) = ', sum(edt4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(edt4d)) = ', sum(edt4d_ref-edt4d)
        ! print*,'sum(edt4d_ref) = ', sum(edt4d_ref)
        ! print*,'sum(edt4d) = ', sum(edt4d)

        open(newunit=fileID, file=trim(dirName) // "/pwav4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) pwav4d_ref
        ! print*,'sum(pwav4d_ref) = ', sum(pwav4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(pwav4d)) = ', sum(pwav4d_ref-pwav4d)
        ! print*,'sum(pwav4d_ref) = ', sum(pwav4d_ref)
        ! print*,'sum(pwav4d) = ', sum(pwav4d)

        open(newunit=fileID, file=trim(dirName) // "/sigma4d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) sigma4d_ref
        ! print*,'sum(sigma4d_ref) = ', sum(sigma4d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(sigma4d)) = ', sum(sigma4d_ref-sigma4d)
        ! print*,'sum(sigma4d_ref) = ', sum(sigma4d_ref)
        ! print*,'sum(sigma4d) = ', sum(sigma4d)

        open(newunit=fileID, file=trim(dirName) // "/pcup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) pcup5d_ref
        ! print*,'sum(pcup5d_ref) = ', sum(pcup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(pcup5d)) = ', sum(pcup5d_ref-pcup5d)
        ! print*,'sum(pcup5d_ref) = ', sum(pcup5d_ref)
        ! print*,'sum(pcup5d) = ', sum(pcup5d)

        open(newunit=fileID, file=trim(dirName) // "/up_massentr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) up_massentr5d_ref
        ! print*,'sum(up_massentr5d_ref) = ', sum(up_massentr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(up_massentr5d)) = ', sum(up_massentr5d_ref-up_massentr5d)
        ! print*,'sum(up_massentr5d_ref) = ', sum(up_massentr5d_ref)
        ! print*,'sum(up_massentr5d) = ', sum(up_massentr5d)

        open(newunit=fileID, file=trim(dirName) // "/up_massdetr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) up_massdetr5d_ref
        ! print*,'sum(up_massdetr5d_ref) = ', sum(up_massdetr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(up_massdetr5d)) = ', sum(up_massdetr5d_ref-up_massdetr5d)
        ! print*,'sum(up_massdetr5d_ref) = ', sum(up_massdetr5d_ref)
        ! print*,'sum(up_massdetr5d) = ', sum(up_massdetr5d)

        open(newunit=fileID, file=trim(dirName) // "/dd_massentr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) dd_massentr5d_ref
        ! print*,'sum(dd_massentr5d_ref) = ', sum(dd_massentr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(dd_massentr5d)) = ', sum(dd_massentr5d_ref-dd_massentr5d)
        ! print*,'sum(dd_massentr5d_ref) = ', sum(dd_massentr5d_ref)
        ! print*,'sum(dd_massentr5d) = ', sum(dd_massentr5d)

        open(newunit=fileID, file=trim(dirName) // "/dd_massdetr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) dd_massdetr5d_ref
        ! print*,'sum(dd_massdetr5d_ref) = ', sum(dd_massdetr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(dd_massdetr5d)) = ', sum(dd_massdetr5d_ref-dd_massdetr5d)
        ! print*,'sum(dd_massdetr5d_ref) = ', sum(dd_massdetr5d_ref)
        ! print*,'sum(dd_massdetr5d) = ', sum(dd_massdetr5d)

        open(newunit=fileID, file=trim(dirName) // "/zup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) zup5d_ref
        ! print*,'sum(zup5d_ref) = ', sum(zup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(zup5d)) = ', sum(zup5d_ref-zup5d)
        ! print*,'sum(zup5d_ref) = ', sum(zup5d_ref)
        ! print*,'sum(zup5d) = ', sum(zup5d)

        open(newunit=fileID, file=trim(dirName) // "/zdn5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) zdn5d_ref
        ! print*,'sum(zdn5d_ref) = ', sum(zdn5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(zdn5d)) = ', sum(zdn5d_ref-zdn5d)
        ! print*,'sum(zdn5d_ref) = ', sum(zdn5d_ref)
        ! print*,'sum(zdn5d) = ', sum(zdn5d)

        open(newunit=fileID, file=trim(dirName) // "/prup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) prup5d_ref
        ! print*,'sum(prup5d_ref) = ', sum(prup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(prup5d)) = ', sum(prup5d_ref-prup5d)
        ! print*,'sum(prup5d_ref) = ', sum(prup5d_ref)
        ! print*,'sum(prup5d) = ', sum(prup5d)

        open(newunit=fileID, file=trim(dirName) // "/prdn5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) prdn5d_ref
        ! print*,'sum(prdn5d_ref) = ', sum(prdn5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(prdn5d)) = ', sum(prdn5d_ref-prdn5d)
        ! print*,'sum(prdn5d_ref) = ', sum(prdn5d_ref)
        ! print*,'sum(prdn5d) = ', sum(prdn5d)

        open(newunit=fileID, file=trim(dirName) // "/clwup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) clwup5d_ref
        ! print*,'sum(clwup5d_ref) = ', sum(clwup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(clwup5d)) = ', sum(clwup5d_ref-clwup5d)
        ! print*,'sum(clwup5d_ref) = ', sum(clwup5d_ref)
        ! print*,'sum(clwup5d) = ', sum(clwup5d)

        open(newunit=fileID, file=trim(dirName) // "/tup5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) tup5d_ref
        ! print*,'sum(tup5d_ref) = ', sum(tup5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(tup5d)) = ', sum(tup5d_ref-tup5d)
        ! print*,'sum(tup5d_ref) = ', sum(tup5d_ref)
        ! print*,'sum(tup5d) = ', sum(tup5d)

        open(newunit=fileID, file=trim(dirName) // "/conv_cld_fr5d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) conv_cld_fr5d_ref
        ! print*,'sum(conv_cld_fr5d_ref) = ', sum(conv_cld_fr5d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(conv_cld_fr5d)) = ', sum(conv_cld_fr5d_ref-conv_cld_fr5d)
        ! print*,'sum(conv_cld_fr5d_ref) = ', sum(conv_cld_fr5d_ref)
        ! print*,'sum(conv_cld_fr5d) = ', sum(conv_cld_fr5d)

        open(newunit=fileID, file=trim(dirName) // "/revsu_gf_2d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) revsu_gf_2d_ref
        ! print*,'sum(revsu_gf_2d_ref) = ', sum(revsu_gf_2d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(revsu_gf_2d)) = ', sum(revsu_gf_2d_ref-revsu_gf_2d)
        ! print*,'sum(revsu_gf_2d_ref) = ', sum(revsu_gf_2d_ref)
        ! print*,'sum(revsu_gf_2d) = ', sum(revsu_gf_2d)

        open(newunit=fileID, file=trim(dirName) // "/prfil_gf_2d_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) prfil_gf_2d_ref
        ! print*,'sum(prfil_gf_2d_ref) = ', sum(prfil_gf_2d_ref)
        close(fileID)

        print*,'***'
        print*,'sum(diff(prfil_gf_2d)) = ', sum(prfil_gf_2d_ref-prfil_gf_2d)
        ! print*,'sum(prfil_gf_2d_ref) = ', sum(prfil_gf_2d_ref)
        ! print*,'sum(prfil_gf_2d) = ', sum(prfil_gf_2d)


        ! Encoding results for CI
        print*, '#CI#VAR|outt#DIFF|',sum(outt_ref - outt)
        print*, '#CI#VAR|outt#NEW|',sum(outt)
        print*, '#CI#VAR|outt#REF|',sum(outt_ref)
        print*, '#CI#VAR|outt#THRSH|',sum(outt_ref)*1.0e-9

        print*, '#CI#VAR|outq#DIFF|',sum(outq_ref - outq)
        print*, '#CI#VAR|outq#NEW|',sum(outq)
        print*, '#CI#VAR|outq#REF|',sum(outq_ref)
        print*, '#CI#VAR|outq#THRSH|',sum(outq_ref)*1.0e-9

        print*, '#CI#VAR|outqc#DIFF|',sum(outqc_ref - outqc)
        print*, '#CI#VAR|outqc#NEW|',sum(outqc)
        print*, '#CI#VAR|outqc#REF|',sum(outqc_ref)
        print*, '#CI#VAR|outqc#THRSH|',sum(outqc_ref)*1.0e-9

        print*, '#CI#VAR|outu#DIFF|',sum(outu_ref - outu)
        print*, '#CI#VAR|outu#NEW|',sum(outu)
        print*, '#CI#VAR|outu#REF|',sum(outu_ref)
        print*, '#CI#VAR|outu#THRSH|',sum(outu_ref)*1.0e-9

        print*, '#CI#VAR|outv#DIFF|',sum(outv_ref - outv)
        print*, '#CI#VAR|outv#NEW|',sum(outv)
        print*, '#CI#VAR|outv#REF|',sum(outv_ref)
        print*, '#CI#VAR|outv#THRSH|',sum(outv_ref)*1.0e-9

        print*, '#CI#VAR|out_chem#DIFF|',sum(out_chem_ref - out_chem)
        print*, '#CI#VAR|out_chem#NEW|',sum(out_chem)
        print*, '#CI#VAR|out_chem#REF|',sum(out_chem_ref)
        print*, '#CI#VAR|out_chem#THRSH|',sum(out_chem_ref)*1.0e-9

        print*, '#CI#VAR|ierr4d_deep#DIFF|',sum(ierr4d_ref - ierr4d_deep)
        print*, '#CI#VAR|ierr4d_deep#NEW|',sum(ierr4d_deep)
        print*, '#CI#VAR|ierr4d_deep#REF|',sum(ierr4d_ref)
        print*, '#CI#VAR|ierr4d_deep#THRSH|',sum(ierr4d_ref)*1.0e-9

        print*, '#CI#VAR|jmin4d#DIFF|',sum(jmin4d_ref - jmin4d)
        print*, '#CI#VAR|jmin4d#NEW|',sum(jmin4d)
        print*, '#CI#VAR|jmin4d#REF|',sum(jmin4d_ref)
        print*, '#CI#VAR|jmin4d#THRSH|',sum(jmin4d_ref)*1.0e-9

        print*, '#CI#VAR|klcl4d#DIFF|',sum(klcl4d_ref - klcl4d)
        print*, '#CI#VAR|klcl4d#NEW|',sum(klcl4d)
        print*, '#CI#VAR|klcl4d#REF|',sum(klcl4d_ref)
        print*, '#CI#VAR|klcl4d#THRSH|',sum(klcl4d_ref)*1.0e-9

        print*, '#CI#VAR|k224d#DIFF|',sum(k224d_ref - k224d)
        print*, '#CI#VAR|k224d#NEW|',sum(k224d)
        print*, '#CI#VAR|k224d#REF|',sum(k224d_ref)
        print*, '#CI#VAR|k224d#THRSH|',sum(k224d_ref)*1.0e-9

        print*, '#CI#VAR|kbcon4d#DIFF|',sum(kbcon4d_ref - kbcon4d)
        print*, '#CI#VAR|kbcon4d#NEW|',sum(kbcon4d)
        print*, '#CI#VAR|kbcon4d#REF|',sum(kbcon4d_ref)
        print*, '#CI#VAR|kbcon4d#THRSH|',sum(kbcon4d_ref)*1.0e-9

        print*, '#CI#VAR|ktop4d#DIFF|',sum(ktop4d_ref - ktop4d)
        print*, '#CI#VAR|ktop4d#NEW|',sum(ktop4d)
        print*, '#CI#VAR|ktop4d#REF|',sum(ktop4d_ref)
        print*, '#CI#VAR|ktop4d#THRSH|',sum(ktop4d_ref)*1.0e-9

        print*, '#CI#VAR|kstabi4d#DIFF|',sum(kstabi4d_ref - kstabi4d)
        print*, '#CI#VAR|kstabi4d#NEW|',sum(kstabi4d)
        print*, '#CI#VAR|kstabi4d#REF|',sum(kstabi4d_ref)
        print*, '#CI#VAR|kstabi4d#THRSH|',sum(kstabi4d_ref)*1.0e-9

        print*, '#CI#VAR|kstabm4d#DIFF|',sum(kstabm4d_ref - kstabm4d)
        print*, '#CI#VAR|kstabm4d#NEW|',sum(kstabm4d)
        print*, '#CI#VAR|kstabm4d#REF|',sum(kstabm4d_ref)
        print*, '#CI#VAR|kstabm4d#THRSH|',sum(kstabm4d_ref)*1.0e-9

        print*, '#CI#VAR|cprr4d#DIFF|',sum(cprr4d_ref - cprr4d)
        print*, '#CI#VAR|cprr4d#NEW|',sum(cprr4d)
        print*, '#CI#VAR|cprr4d#REF|',sum(cprr4d_ref)
        print*, '#CI#VAR|cprr4d#THRSH|',sum(cprr4d_ref)*1.0e-9

        print*, '#CI#VAR|xmb4d#DIFF|',sum(xmb4d_ref - xmb4d)
        print*, '#CI#VAR|xmb4d#NEW|',sum(xmb4d)
        print*, '#CI#VAR|xmb4d#REF|',sum(xmb4d_ref)
        print*, '#CI#VAR|xmb4d#THRSH|',sum(xmb4d_ref)*1.0e-9

        print*, '#CI#VAR|edt4d#DIFF|',sum(edt4d_ref - edt4d)
        print*, '#CI#VAR|edt4d#NEW|',sum(edt4d)
        print*, '#CI#VAR|edt4d#REF|',sum(edt4d_ref)
        print*, '#CI#VAR|edt4d#THRSH|',sum(edt4d_ref)*1.0e-9

        print*, '#CI#VAR|pwav4d#DIFF|',sum(pwav4d_ref - pwav4d)
        print*, '#CI#VAR|pwav4d#NEW|',sum(pwav4d)
        print*, '#CI#VAR|pwav4d#REF|',sum(pwav4d_ref)
        print*, '#CI#VAR|pwav4d#THRSH|',sum(pwav4d_ref)*1.0e-9

        print*, '#CI#VAR|sigma4d#DIFF|',sum(sigma4d_ref - sigma4d)
        print*, '#CI#VAR|sigma4d#NEW|',sum(sigma4d)
        print*, '#CI#VAR|sigma4d#REF|',sum(sigma4d_ref)
        print*, '#CI#VAR|sigma4d#THRSH|',sum(sigma4d_ref)*1.0e-9

        print*, '#CI#VAR|pcup5d#DIFF|',sum(pcup5d_ref - pcup5d)
        print*, '#CI#VAR|pcup5d#NEW|',sum(pcup5d)
        print*, '#CI#VAR|pcup5d#REF|',sum(pcup5d_ref)
        print*, '#CI#VAR|pcup5d#THRSH|',sum(pcup5d_ref)*1.0e-9

        print*, '#CI#VAR|up_massentr5d#DIFF|',sum(up_massentr5d_ref - up_massentr5d)
        print*, '#CI#VAR|up_massentr5d#NEW|',sum(up_massentr5d)
        print*, '#CI#VAR|up_massentr5d#REF|',sum(up_massentr5d_ref)
        print*, '#CI#VAR|up_massentr5d#THRSH|',sum(up_massentr5d_ref)*1.0e-9

        print*, '#CI#VAR|dd_massentr5d#DIFF|',sum(dd_massentr5d_ref - dd_massentr5d)
        print*, '#CI#VAR|dd_massentr5d#NEW|',sum(dd_massentr5d)
        print*, '#CI#VAR|dd_massentr5d#REF|',sum(dd_massentr5d_ref)
        print*, '#CI#VAR|dd_massentr5d#THRSH|',sum(dd_massentr5d_ref)*1.0e-9

        print*, '#CI#VAR|dd_massdetr5d#DIFF|',sum(dd_massdetr5d_ref - dd_massdetr5d)
        print*, '#CI#VAR|dd_massdetr5d#NEW|',sum(dd_massdetr5d)
        print*, '#CI#VAR|dd_massdetr5d#REF|',sum(dd_massdetr5d_ref)
        print*, '#CI#VAR|dd_massdetr5d#THRSH|',sum(dd_massdetr5d_ref)*1.0e-9

        print*, '#CI#VAR|zup5d#DIFF|',sum(zup5d_ref - zup5d)
        print*, '#CI#VAR|zup5d#NEW|',sum(zup5d)
        print*, '#CI#VAR|zup5d#REF|',sum(zup5d_ref)
        print*, '#CI#VAR|zup5d#THRSH|',sum(zup5d_ref)*1.0e-9

        print*, '#CI#VAR|zdn5d#DIFF|',sum(zdn5d_ref - zdn5d)
        print*, '#CI#VAR|zdn5d#NEW|',sum(zdn5d)
        print*, '#CI#VAR|zdn5d#REF|',sum(zdn5d_ref)
        print*, '#CI#VAR|zdn5d#THRSH|',sum(zdn5d_ref)*1.0e-9

        print*, '#CI#VAR|prup5d#DIFF|',sum(prup5d_ref - prup5d)
        print*, '#CI#VAR|prup5d#NEW|',sum(prup5d)
        print*, '#CI#VAR|prup5d#REF|',sum(prup5d_ref)
        print*, '#CI#VAR|prup5d#THRSH|',sum(prup5d_ref)*1.0e-9

        print*, '#CI#VAR|prdn5d#DIFF|',sum(prdn5d_ref - prdn5d)
        print*, '#CI#VAR|prdn5d#NEW|',sum(prdn5d)
        print*, '#CI#VAR|prdn5d#REF|',sum(prdn5d_ref)
        print*, '#CI#VAR|prdn5d#THRSH|',sum(prdn5d_ref)*1.0e-9

        print*, '#CI#VAR|clwup5d#DIFF|',sum(clwup5d_ref - clwup5d)
        print*, '#CI#VAR|clwup5d#NEW|',sum(clwup5d)
        print*, '#CI#VAR|clwup5d#REF|',sum(clwup5d_ref)
        print*, '#CI#VAR|clwup5d#THRSH|',sum(clwup5d_ref)*1.0e-9

        print*, '#CI#VAR|tup5d#DIFF|',sum(tup5d_ref - tup5d)
        print*, '#CI#VAR|tup5d#NEW|',sum(tup5d)
        print*, '#CI#VAR|tup5d#REF|',sum(tup5d_ref)
        print*, '#CI#VAR|tup5d#THRSH|',sum(tup5d_ref)*1.0e-9

        print*, '#CI#VAR|conv_cld_fr5d#DIFF|',sum(conv_cld_fr5d_ref - conv_cld_fr5d)
        print*, '#CI#VAR|conv_cld_fr5d#NEW|',sum(conv_cld_fr5d)
        print*, '#CI#VAR|conv_cld_fr5d#REF|',sum(conv_cld_fr5d_ref)
        print*, '#CI#VAR|conv_cld_fr5d#THRSH|',sum(conv_cld_fr5d_ref)*1.0e-9

        print*, '#CI#VAR|revsu_gf_2d#DIFF|',sum(revsu_gf_2d_ref - revsu_gf_2d)
        print*, '#CI#VAR|revsu_gf_2d#NEW|',sum(revsu_gf_2d)
        print*, '#CI#VAR|revsu_gf_2d#REF|',sum(revsu_gf_2d_ref)
        print*, '#CI#VAR|revsu_gf_2d#THRSH|',sum(revsu_gf_2d_ref)*1.0e-9

        print*, '#CI#VAR|prfil_gf_2d#DIFF|',sum(prfil_gf_2d_ref - prfil_gf_2d)
        print*, '#CI#VAR|prfil_gf_2d#NEW|',sum(prfil_gf_2d)
        print*, '#CI#VAR|prfil_gf_2d#REF|',sum(prfil_gf_2d_ref)
        print*, '#CI#VAR|prfil_gf_2d#THRSH|',sum(prfil_gf_2d_ref)*1.0e-9

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