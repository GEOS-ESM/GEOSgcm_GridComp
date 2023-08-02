module test_GF2020_subroutines
    use moist_subroutines_GF2020
    use ConvPar_GF_SharedParams
    use timing_module

    implicit none

    public test_cup_gf_GF2020

    private

    integer :: ens4 = 1

    integer :: cum_use_excess

    real :: dt, cum_entr_rate_shal, cum_entr_rate_deep

    integer, dimension(:), allocatable :: kpbli, kpbli_ref, ierr4d, ierr4d_ref, jmin4d, jmin4d_ref, klcl4d, klcl4d_ref, &
        k224d, k224d_ref, kbcon4d, kbcon4d_ref, ktop4d, ktop4d_ref, kstabi4d, kstabi4d_ref, &
        kstabm4d, kstabm4d_ref, ierr4d_mid, ierr4d_deep, last_ierr, last_ierr_ref

    real, dimension(:), allocatable :: dx2d, h_sfc_flux, h_sfc_flux_ref, le_sfc_flux, le_sfc_flux_ref, & 
        tsur, tsur_ref, psur, psur_ref, ter11, ter11_ref, xlandi, xlandi_ref, ztexec, zqexec, ccn, ccn_ref, &
        zws, zws_ref, xmb4d, xmb4d_ref, edt4d, edt4d_ref, pwav4d, pwav4d_ref, stochastic_sig, xlons, xlons_ref, &
        xlats, xlats_ref, cprr4d, cprr4d_ref, &
        sigma4d, sigma4d_ref, col_sat, cum_ztexec, cum_ztexec_ref, cum_zqexec, cum_zqexec_ref, &
        AA0, AA0_ref, AA1, AA1_ref, AA2, AA2_ref, AA3, AA3_ref, AA1_BL, AA1_BL_ref, AA1_CIN, AA1_CIN_ref, TAU_BL, TAU_BL_ref, &
        TAU_EC, TAU_EC_ref, lightn_dens, lightn_dens_ref, var2d, var2d_ref, &
        cnvfrc, srftype

    real, dimension(:,:), allocatable :: rhoi, rhoi_ref, dhdt, dhdt_ref, zo, zo_ref, dm2d, dm2d_ref, temp_old, temp_old_ref, qv_old, qv_old_ref, &
        temp_new_sh, qv_new_sh, &
        po, po_ref, us, us_ref, vs, vs_ref, outt, outt_ref, outq, outq_ref, outqc, outqc_ref, outu, outu_ref, outv, outv_ref, &
        pcup5d, pcup5d_ref, up_massentr5d, up_massentr5d_ref, up_massdetr5d, up_massdetr5d_ref, &
        dd_massentr5d, dd_massentr5d_ref, dd_massdetr5d, dd_massdetr5d_ref, &
        zup5d, zup5d_ref, zdn5d, zdn5d_ref, prup5d, prup5d_ref, &
        prdn5d, prdn5d_ref, clwup5d, clwup5d_ref, tup5d, tup5d_ref, &
        conv_cld_fr5d, conv_cld_fr5d_ref, temp_new_dp, qv_new_dp, &
        temp_new_BL, temp_new_BL_ref, qv_new_BL, qv_new_BL_ref, revsu_gf_2d, revsu_gf_2d_ref, &
        prfil_gf_2d, prfil_gf_2d_ref, temp_new, temp_new_ref, qv_new, qv_new_ref, temp_new_ADV, temp_new_ADV_ref, &
        qv_new_ADV, qv_new_ADV_ref, buoy_exc2d, buoy_exc2d_ref, outnliq, outnliq_ref, outnice, outnice_ref, outbuoy, outbuoy_ref, &
        var3d_agf_2d, var3d_agf_2d_ref, var3d_bgf_2d, var3d_bgf_2d_ref, Tpert_2d, Tpert_2d_ref

    real, dimension(:,:,:), allocatable :: se_chem, se_chem_ref, out_chem, out_chem_ref, omeg, omeg_ref, &
        mpqi, mpqi_ref, mpql, mpql_ref, mpcf, mpcf_ref, &
        outmpqi, outmpqi_ref, outmpql, outmpql_ref, outmpcf, outmpcf_ref

    contains

    subroutine test_cup_gf_GF2020(IM, JM, LM, dirName, rank_str)

        integer :: IM, JM, LM, fileID, itf, ktf, its, ite, kts, kte, mtp, mxp, mzp, ii, nmp, plume
        character*100 :: dirName, rank_str
        character*2   :: i_string

        itf = IM
        ktf = LM - 1
        its = 1
        ite = IM
        kts = 1
        kte = LM
        ! if (dirName(1:10) == './c24_data') then
        !     mtp = 40
        !     mxp = IM
        !     mzp = LM
        ! else if (dirName(1:10) == './c90_data') then
        !     mtp = 40
        !     mxp = IM
        !     mzp = LM
        ! else if (dirName(1:11) == './c180_data') then
            mtp = 40
            mxp = IM
            mzp = LM
        ! endif

        print*, 'Testing CUP_gf (ver GF2020)'

        allocate(dx2d       (its:ite))
        allocate(stochastic_sig (its:ite))
        allocate(col_sat (its:ite))
        allocate(kpbli      (its:ite))
        allocate(kpbli_ref      (its:ite))
        allocate(cum_ztexec      (its:ite))
        allocate(cum_ztexec_ref      (its:ite))
        allocate(cum_zqexec      (its:ite))
        allocate(cum_zqexec_ref      (its:ite))
        allocate(ccn        (its:ite))
        allocate(ccn_ref        (its:ite))
        allocate(ter11      (its:ite))
        allocate(ter11_ref      (its:ite))
        allocate(h_sfc_flux (its:ite))
        allocate(h_sfc_flux_ref (its:ite))
        allocate(le_sfc_flux(its:ite))
        allocate(le_sfc_flux_ref(its:ite))
        allocate(xlons      (its:ite))
        allocate(xlons_ref      (its:ite))
        allocate(xlats      (its:ite))
        allocate(xlats_ref      (its:ite))
        allocate(xlandi     (its:ite))
        allocate(xlandi_ref     (its:ite))
        allocate(cnvfrc     (its:ite))
        allocate(srftype     (its:ite))
        allocate(tsur       (its:ite))
        allocate(tsur_ref       (its:ite))
        allocate(psur       (its:ite))
        allocate(psur_ref       (its:ite))
        allocate(zws        (its:ite))
        allocate(zws_ref        (its:ite))
        allocate(last_ierr  (its:ite))
        allocate(last_ierr_ref  (its:ite))
        allocate(lightn_dens (its:ite))
        allocate(lightn_dens_ref (its:ite))

        allocate(ierr4d      (mxp))
        allocate(ierr4d_ref      (mxp))
        allocate(jmin4d      (mxp))
        allocate(jmin4d_ref      (mxp))
        allocate(klcl4d      (mxp))
        allocate(klcl4d_ref      (mxp))
        allocate(k224d       (mxp))
        allocate(k224d_ref       (mxp))
        allocate(kbcon4d     (mxp))
        allocate(kbcon4d_ref     (mxp))
        allocate(ktop4d      (mxp))
        allocate(ktop4d_ref      (mxp))
        allocate(kstabi4d    (mxp))
        allocate(kstabi4d_ref    (mxp))
        allocate(kstabm4d    (mxp))
        allocate(kstabm4d_ref    (mxp))
        allocate(cprr4d    (mxp))
        allocate(cprr4d_ref    (mxp))
        allocate(xmb4d       (mxp))
        allocate(xmb4d_ref       (mxp))
        allocate(edt4d       (mxp))
        allocate(edt4d_ref       (mxp))
        allocate(pwav4d      (mxp))
        allocate(pwav4d_ref      (mxp))
        allocate(AA0(mxp))
        allocate(AA1(mxp))
        allocate(AA2(mxp))
        allocate(AA3(mxp))
        allocate(AA1_BL(mxp))
        allocate(AA1_CIN(mxp))
        allocate(TAU_BL(mxp))
        allocate(TAU_EC(mxp))
        allocate(AA0_ref(mxp))
        allocate(AA1_ref(mxp))
        allocate(AA2_ref(mxp))
        allocate(AA3_ref(mxp))
        allocate(AA1_BL_ref(mxp))
        allocate(AA1_CIN_ref(mxp))
        allocate(TAU_BL_ref(mxp))
        allocate(TAU_EC_ref(mxp))
        allocate(sigma4d     (mxp))
        allocate(sigma4d_ref     (mxp))
        allocate(var2d(mxp))
        allocate(var2d_ref(mxp))

        allocate(pcup5d           (mxp, mzp))
        allocate(up_massentr5d    (mxp, mzp))
        allocate(up_massdetr5d    (mxp, mzp))
        allocate(dd_massentr5d    (mxp, mzp))
        allocate(dd_massdetr5d    (mxp, mzp))
        allocate(zup5d            (mxp, mzp))
        allocate(zdn5d            (mxp, mzp))
        allocate(prup5d           (mxp, mzp))
        allocate(prdn5d           (mxp, mzp))
        allocate(clwup5d          (mxp, mzp))
        allocate(tup5d            (mxp, mzp))
        allocate(conv_cld_fr5d    (mxp, mzp))
        allocate(pcup5d_ref           (mxp, mzp))
        allocate(up_massentr5d_ref    (mxp, mzp))
        allocate(up_massdetr5d_ref    (mxp, mzp))
        allocate(dd_massentr5d_ref    (mxp, mzp))
        allocate(dd_massdetr5d_ref    (mxp, mzp))
        allocate(zup5d_ref            (mxp, mzp))
        allocate(zdn5d_ref            (mxp, mzp))
        allocate(prup5d_ref           (mxp, mzp))
        allocate(prdn5d_ref           (mxp, mzp))
        allocate(clwup5d_ref          (mxp, mzp))
        allocate(tup5d_ref            (mxp, mzp))
        allocate(conv_cld_fr5d_ref    (mxp, mzp))

        allocate(rhoi       (its:ite, kts:kte))
        allocate(rhoi_ref       (its:ite, kts:kte))
        allocate(temp_old   (its:ite, kts:kte))
        allocate(temp_old_ref   (its:ite, kts:kte))
        allocate(qv_old     (its:ite, kts:kte))
        allocate(qv_old_ref     (its:ite, kts:kte))
        allocate(temp_new   (its:ite, kts:kte))
        allocate(temp_new_ref   (its:ite, kts:kte))
        allocate(qv_new     (its:ite, kts:kte))
        allocate(qv_new_ref     (its:ite, kts:kte))
        allocate(temp_new_BL (its:ite, kts:kte))
        allocate(temp_new_BL_ref (its:ite, kts:kte))
        allocate(qv_new_BL (its:ite, kts:kte))
        allocate(qv_new_BL_ref (its:ite, kts:kte))
        allocate(temp_new_ADV (its:ite, kts:kte))
        allocate(temp_new_ADV_ref (its:ite, kts:kte))
        allocate(qv_new_ADV (its:ite, kts:kte))
        allocate(qv_new_ADV_ref (its:ite, kts:kte))
        allocate(zo         (its:ite, kts:kte))
        allocate(zo_ref         (its:ite, kts:kte))
        allocate(po         (its:ite, kts:kte))
        allocate(po_ref         (its:ite, kts:kte))
        allocate(us         (its:ite, kts:kte))
        allocate(us_ref         (its:ite, kts:kte))
        allocate(vs         (its:ite, kts:kte))
        allocate(vs_ref         (its:ite, kts:kte))
        allocate(dm2d       (its:ite, kts:kte))
        allocate(dm2d_ref       (its:ite, kts:kte))
        allocate(dhdt       (its:ite, kts:kte))
        allocate(dhdt_ref       (its:ite, kts:kte))
        allocate(buoy_exc2d       (its:ite, kts:kte))
        allocate(buoy_exc2d_ref       (its:ite, kts:kte))
        allocate(outt       (its:ite, kts:kte))
        allocate(outt_ref       (its:ite, kts:kte))
        allocate(outq       (its:ite, kts:kte))
        allocate(outq_ref       (its:ite, kts:kte))
        allocate(outqc      (its:ite, kts:kte))
        allocate(outqc_ref      (its:ite, kts:kte))
        allocate(outu       (its:ite, kts:kte))
        allocate(outu_ref       (its:ite, kts:kte))
        allocate(outv       (its:ite, kts:kte))
        allocate(outv_ref       (its:ite, kts:kte))
        allocate(outnliq    (its:ite, kts:kte))
        allocate(outnliq_ref    (its:ite, kts:kte))
        allocate(outnice    (its:ite, kts:kte))
        allocate(outnice_ref    (its:ite, kts:kte))
        allocate(outbuoy    (its:ite, kts:kte))
        allocate(outbuoy_ref    (its:ite, kts:kte))
        allocate(revsu_gf_2d(its:ite, kts:kte))
        allocate(prfil_gf_2d(its:ite, kts:kte))
        allocate(var3d_agf_2d(its:ite, kts:kte))
        allocate(var3d_bgf_2d(its:ite, kts:kte))
        allocate(Tpert_2d(its:ite, kts:kte))
        allocate(revsu_gf_2d_ref(its:ite, kts:kte))
        allocate(prfil_gf_2d_ref(its:ite, kts:kte))
        allocate(var3d_agf_2d_ref(its:ite, kts:kte))
        allocate(var3d_bgf_2d_ref(its:ite, kts:kte))
        allocate(Tpert_2d_ref(its:ite, kts:kte))

        allocate(omeg(its:ite, kts:kte, 1:ens4))
        allocate(omeg_ref(its:ite, kts:kte, 1:ens4))
        allocate(se_chem     (mtp, its:ite, kts:kte))
        allocate(se_chem_ref     (mtp, its:ite, kts:kte))
        allocate(out_chem    (mtp, its:ite, kts:kte))
        allocate(out_chem_ref    (mtp, its:ite, kts:kte))

        ! Size of CNV_Tracers seem to match to mpt.  Note that this may change later
        allocate(CNV_Tracers(mtp))

        open(newunit=fileID, file=trim(dirName) // "/bc_meth_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) BC_METH
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/cum_ave_layer_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) CUM_AVE_LAYER
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/HEI_DOWN_LAND_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) HEI_DOWN_LAND
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/HEI_DOWN_OCEAN_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) HEI_DOWN_OCEAN
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/HEI_UPDF_LAND_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) HEI_UPDF_LAND
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/HEI_UPDF_OCEAN_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) HEI_UPDF_OCEAN
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MAX_EDT_LAND_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) MAX_EDT_LAND
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MAX_EDT_OCEAN_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) MAX_EDT_OCEAN
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/MAX_EDT_LAND_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) MAX_EDT_LAND
        ! print*,'nmp = ', nmp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // "/FADJ_MASSFLX_" // trim(rank_str) // ".in", &
        status='old', form="unformatted", action="read")
        read(fileID) FADJ_MASSFLX
        ! print*,'nmp = ', nmp
        close(fileID)

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

        open(newunit=fileID, file=trim(dirName) // "/nmp_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) nmp
        ! print*,'nmp = ', nmp
        close(fileID)

        allocate(mpqi(nmp, its:ite, kts:kte))
        allocate(mpqi_ref(nmp, its:ite, kts:kte))
        allocate(mpql(nmp, its:ite, kts:kte))
        allocate(mpql_ref(nmp, its:ite, kts:kte))
        allocate(mpcf(nmp, its:ite, kts:kte))
        allocate(mpcf_ref(nmp, its:ite, kts:kte))
        allocate(outmpqi(nmp, its:ite, kts:kte))
        allocate(outmpqi_ref(nmp, its:ite, kts:kte))
        allocate(outmpql(nmp, its:ite, kts:kte))
        allocate(outmpql_ref(nmp, its:ite, kts:kte))
        allocate(outmpcf(nmp, its:ite, kts:kte))
        allocate(outmpcf_ref(nmp, its:ite, kts:kte))

        open(newunit=fileID, file=trim(dirName) // "/plume_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) plume
        ! print*,'plume = ', plume
        close(fileID)

        ! print*,'cumulus_type(plume) = ', cumulus_type(plume)
        ! print*,'closure_choice(plume) = ', closure_choice(plume)
        ! print*,'cum_entr_rate(plume) = ', cum_entr_rate(plume)

        open(newunit=fileID, file=trim(dirName) // '/cum_use_excess_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cum_use_excess
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In cum_use_excess = ', cum_use_excess

        open(newunit=fileID, file=trim(dirName) // '/dx2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dx2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In dx2d = ', sum(dx2d)

        open(newunit=fileID, file=trim(dirName) // '/stochastic_sig_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) stochastic_sig
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In stochastic_sig = ', sum(stochastic_sig)

        open(newunit=fileID, file=trim(dirName) // '/col_sat_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) col_sat
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In col_sat = ', sum(col_sat)

        open(newunit=fileID, file=trim(dirName) // '/dt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dt
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In dt = ', dt

        open(newunit=fileID, file=trim(dirName) // '/kpbli_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) kpbli
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In kpbli = ', sum(kpbli)

        open(newunit=fileID, file=trim(dirName) // '/cum_ztexec_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cum_ztexec
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In cum_ztexec = ', sum(cum_ztexec)

        open(newunit=fileID, file=trim(dirName) // '/cum_zqexec_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cum_zqexec
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In cum_zqexec = ', sum(cum_zqexec)

        open(newunit=fileID, file=trim(dirName) // '/ccn_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ccn
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In ccn = ', sum(ccn)

        open(newunit=fileID, file=trim(dirName) // '/rhoi_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rhoi
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In rhoi = ', sum(rhoi)

        open(newunit=fileID, file=trim(dirName) // '/omeg_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) omeg
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In omeg = ', sum(omeg)

        open(newunit=fileID, file=trim(dirName) // '/temp_old_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) temp_old
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In temp_old = ', sum(temp_old)

        open(newunit=fileID, file=trim(dirName) // '/qv_old_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qv_old
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In qv_old = ', sum(qv_old)

        open(newunit=fileID, file=trim(dirName) // '/ter11_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ter11
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In ter11 = ', sum(ter11)

        open(newunit=fileID, file=trim(dirName) // '/h_sfc_flux_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) h_sfc_flux
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In h_sfc_flux = ', sum(h_sfc_flux)

        open(newunit=fileID, file=trim(dirName) // '/le_sfc_flux_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) le_sfc_flux
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In le_sfc_flux = ', sum(le_sfc_flux)

        open(newunit=fileID, file=trim(dirName) // '/xlons_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) xlons
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In xlons = ', sum(xlons)

        open(newunit=fileID, file=trim(dirName) // '/xlats_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) xlats
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In xlats = ', sum(xlats)

        open(newunit=fileID, file=trim(dirName) // '/xlandi_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) xlandi
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In xlandi = ', sum(xlandi)

        open(newunit=fileID, file=trim(dirName) // '/cnvfrc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cnvfrc
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In cnvfrc = ', sum(cnvfrc)

        open(newunit=fileID, file=trim(dirName) // '/srftype_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) srftype
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In srftype = ', sum(srftype)

        open(newunit=fileID, file=trim(dirName) // '/temp_new_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) temp_new
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In temp_new = ', sum(temp_new)

        open(newunit=fileID, file=trim(dirName) // '/qv_new_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qv_new
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In qv_new = ', sum(qv_new)

        open(newunit=fileID, file=trim(dirName) // '/temp_new_BL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) temp_new_BL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In temp_new_BL = ', sum(temp_new_BL)

        open(newunit=fileID, file=trim(dirName) // '/qv_new_BL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qv_new_BL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In qv_new_BL = ', sum(qv_new_BL)

        open(newunit=fileID, file=trim(dirName) // '/temp_new_ADV_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) temp_new_ADV
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In temp_new_ADV = ', sum(temp_new_ADV)

        open(newunit=fileID, file=trim(dirName) // '/qv_new_ADV_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qv_new_ADV
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In qv_new_ADV = ', sum(qv_new_ADV)

        open(newunit=fileID, file=trim(dirName) // '/zo_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) zo
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In zo = ', sum(zo)

        open(newunit=fileID, file=trim(dirName) // '/po_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) po
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In po = ', sum(po)

        open(newunit=fileID, file=trim(dirName) // '/tsur_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tsur
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In tsur = ', sum(tsur)

        open(newunit=fileID, file=trim(dirName) // '/psur_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) psur
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In psur = ', sum(psur)

        open(newunit=fileID, file=trim(dirName) // '/us_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) us
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In us = ', sum(us)

        open(newunit=fileID, file=trim(dirName) // '/vs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vs
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In vs = ', sum(vs)

        open(newunit=fileID, file=trim(dirName) // '/dm2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dm2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In dm2d = ', sum(dm2d)

        open(newunit=fileID, file=trim(dirName) // '/se_chem_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) se_chem
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In se_chem = ', sum(se_chem)

        open(newunit=fileID, file=trim(dirName) // '/zws_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) zws
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In zws = ', sum(zws)

        open(newunit=fileID, file=trim(dirName) // '/dhdt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dhdt
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In dhdt = ', sum(dhdt)

        open(newunit=fileID, file=trim(dirName) // '/buoy_exc2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) buoy_exc2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In buoy_exc2d = ', sum(buoy_exc2d)

        open(newunit=fileID, file=trim(dirName) // '/mpqi_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mpqi
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In mpqi = ', sum(mpqi)

        open(newunit=fileID, file=trim(dirName) // '/mpql_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mpql
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In mpql = ', sum(mpql)

        open(newunit=fileID, file=trim(dirName) // '/mpcf_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mpcf
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In mpcf = ', sum(mpcf)

        open(newunit=fileID, file=trim(dirName) // '/last_ierr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) last_ierr
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In last_ierr = ', sum(last_ierr)

        open(newunit=fileID, file=trim(dirName) // '/outt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outt
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outt = ', sum(outt)

        open(newunit=fileID, file=trim(dirName) // '/outq_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outq
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outq = ', sum(outq)

        open(newunit=fileID, file=trim(dirName) // '/outqc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outqc
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outqc = ', sum(outqc)

        open(newunit=fileID, file=trim(dirName) // '/outu_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outu
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outu = ', sum(outu)

        open(newunit=fileID, file=trim(dirName) // '/outv_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outv
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outv = ', sum(outv)

        open(newunit=fileID, file=trim(dirName) // '/outnliq_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outnliq
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outnliq = ', sum(outnliq)

        open(newunit=fileID, file=trim(dirName) // '/outnice_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outnice
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outnice = ', sum(outnice)

        open(newunit=fileID, file=trim(dirName) // '/outbuoy_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outbuoy
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outbuoy = ', sum(outbuoy)

        open(newunit=fileID, file=trim(dirName) // '/outmpqi_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outmpqi
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outmpqi = ', sum(outmpqi)

        open(newunit=fileID, file=trim(dirName) // '/outmpql_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outmpql
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outmpql = ', sum(outmpql)

        open(newunit=fileID, file=trim(dirName) // '/outmpcf_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) outmpcf
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In outmpcf = ', sum(outmpcf)

        open(newunit=fileID, file=trim(dirName) // '/out_chem_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) out_chem
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In out_chem = ', sum(out_chem)

        open(newunit=fileID, file=trim(dirName) // '/ierr4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ierr4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In ierr4d = ', sum(ierr4d)

        open(newunit=fileID, file=trim(dirName) // '/jmin4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) jmin4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In jmin4d = ', sum(jmin4d)

        open(newunit=fileID, file=trim(dirName) // '/klcl4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) klcl4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In klcl4d = ', sum(klcl4d)

        open(newunit=fileID, file=trim(dirName) // '/k224d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) k224d 
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In k224d = ', sum(k224d )

        open(newunit=fileID, file=trim(dirName) // '/kbcon4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) kbcon4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In kbcon4d = ', sum(kbcon4d)

        open(newunit=fileID, file=trim(dirName) // '/ktop4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ktop4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In ktop4d = ', sum(ktop4d)

        open(newunit=fileID, file=trim(dirName) // '/kstabi4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) kstabi4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In kstabi4d = ', sum(kstabi4d)

        open(newunit=fileID, file=trim(dirName) // '/kstabm4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) kstabm4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In kstabm4d = ', sum(kstabm4d)

        open(newunit=fileID, file=trim(dirName) // '/xmb4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) xmb4d 
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In xmb4d = ', sum(xmb4d )

        open(newunit=fileID, file=trim(dirName) // '/edt4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) edt4d 
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In edt4d = ', sum(edt4d )

        open(newunit=fileID, file=trim(dirName) // '/pwav4d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) pwav4d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In pwav4d = ', sum(pwav4d)

        open(newunit=fileID, file=trim(dirName) // '/pcup5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) pcup5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In pcup5d = ', sum(pcup5d)

        open(newunit=fileID, file=trim(dirName) // '/up_massentr5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) up_massentr5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In up_massentr5d = ', sum(up_massentr5d)

        open(newunit=fileID, file=trim(dirName) // '/up_massdetr5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) up_massdetr5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In up_massdetr5d = ', sum(up_massdetr5d)

        open(newunit=fileID, file=trim(dirName) // '/dd_massentr5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dd_massentr5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In dd_massentr5d = ', sum(dd_massentr5d)

        open(newunit=fileID, file=trim(dirName) // '/dd_massdetr5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dd_massdetr5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In dd_massdetr5d = ', sum(dd_massdetr5d)

        open(newunit=fileID, file=trim(dirName) // '/zup5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) zup5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In zup5d = ', sum(zup5d)

        open(newunit=fileID, file=trim(dirName) // '/zdn5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) zdn5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In zdn5d = ', sum(zdn5d)

        open(newunit=fileID, file=trim(dirName) // '/prup5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) prup5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In prup5d = ', sum(prup5d)

        open(newunit=fileID, file=trim(dirName) // '/prdn5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) prdn5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In prdn5d = ', sum(prdn5d)

        open(newunit=fileID, file=trim(dirName) // '/clwup5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) clwup5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In clwup5d = ', sum(clwup5d)

        open(newunit=fileID, file=trim(dirName) // '/tup5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tup5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In tup5d = ', sum(tup5d)

        open(newunit=fileID, file=trim(dirName) // '/conv_cld_fr5d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) conv_cld_fr5d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In conv_cld_fr5d = ', sum(conv_cld_fr5d)

        open(newunit=fileID, file=trim(dirName) // '/AA0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AA0
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AA0 = ', sum(AA0)

        open(newunit=fileID, file=trim(dirName) // '/AA1_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AA1
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AA1 = ', sum(AA1)

        open(newunit=fileID, file=trim(dirName) // '/AA2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AA2
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AA2 = ', sum(AA2)

        open(newunit=fileID, file=trim(dirName) // '/AA3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AA3
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AA3 = ', sum(AA3)

        open(newunit=fileID, file=trim(dirName) // '/AA1_BL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AA1_BL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AA1_BL = ', sum(AA1_BL)

        open(newunit=fileID, file=trim(dirName) // '/AA1_CIN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AA1_CIN
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AA1_CIN = ', sum(AA1_CIN)

        open(newunit=fileID, file=trim(dirName) // '/TAU_BL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) TAU_BL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In TAU_BL = ', sum(TAU_BL)

        open(newunit=fileID, file=trim(dirName) // '/TAU_EC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) TAU_EC
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In TAU_EC = ', sum(TAU_EC)

        open(newunit=fileID, file=trim(dirName) // '/TAU_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) TAU_DEEP
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In TAU_DEEP = ', sum(TAU_DEEP)

        open(newunit=fileID, file=trim(dirName) // '/revsu_gf_2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) revsu_gf_2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In revsu_gf_2d = ', sum(revsu_gf_2d)

        open(newunit=fileID, file=trim(dirName) // '/prfil_gf_2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) prfil_gf_2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In prfil_gf_2d = ', sum(prfil_gf_2d)

        open(newunit=fileID, file=trim(dirName) // '/var3d_agf_2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) var3d_agf_2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In var3d_agf_2d = ', sum(var3d_agf_2d)

        open(newunit=fileID, file=trim(dirName) // '/var3d_bgf_2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) var3d_bgf_2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In var3d_bgf_2d = ', sum(var3d_bgf_2d)

        open(newunit=fileID, file=trim(dirName) // '/Tpert_2d_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) Tpert_2d
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In Tpert_2d = ', sum(Tpert_2d)

        call start_timing()

        CALL CUP_gf_GF2020(its,ite,kts,kte, itf,ktf, mtp, nmp &
                  ,cumulus_type  (plume)            &
                  ,closure_choice(plume)            &
                  ,cum_entr_rate (plume)            &
                !   ,cum_use_excess(plume)            &
                  ,cum_use_excess            &
                  !- input data
                !   ,dx2d          (:,j)              &
                !   ,stochastic_sig(:,j)              &
                !   ,col_sat       (:,j)              &
                  ,dx2d          (:)              &
                  ,stochastic_sig(:)              &
                  ,col_sat       (:)              &
                  ,dt                               &
                  ,kpbli                            &
                  ,cum_ztexec                       &
                  ,cum_zqexec                       &
                  ,ccn                              &
                  ,rhoi                             &
                  ,omeg                             &
                  ,temp_old                         &
                  ,qv_old                           &
                  ,ter11                            &
                  , h_sfc_flux                      &
                  ,le_sfc_flux                      &
                  ,xlons                            &
                  ,xlats                            &
                  ,xlandi                           &
                !   ,cnvfrc(:,j)                      &
                !   ,srftype(:,j)                     &
                  ,cnvfrc(:)                      &
                  ,srftype(:)                     &
                  ,temp_new                         &
                  ,qv_new                           &
                  ,temp_new_BL                      &
                  ,qv_new_BL                        &
                  ,temp_new_ADV                     &
                  ,qv_new_ADV                       &
                  ,zo                               &
                  ,po                               &
                  ,tsur                             &
                  ,psur                             &
                  ,us                               &
                  ,vs                               &
                  ,dm2d                             &
                  ,se_chem                          &
                  ,zws                              &
                  ,dhdt                             &
                  ,buoy_exc2d                       &
                  ,mpqi                             &
                  ,mpql                             &
                  ,mpcf                             &
                  ,last_ierr            (:)         &
                  !output data
                !   ,outt                 (:,:,plume) &
                !   ,outq                 (:,:,plume) &
                !   ,outqc                (:,:,plume) &
                !   ,outu                 (:,:,plume) &
                !   ,outv                 (:,:,plume) &
                !   ,outnliq              (:,:,plume) &
                !   ,outnice              (:,:,plume) &
                !   ,outbuoy              (:,:,plume) &
                  ,outt                 (:,:) &
                  ,outq                 (:,:) &
                  ,outqc                (:,:) &
                  ,outu                 (:,:) &
                  ,outv                 (:,:) &
                  ,outnliq              (:,:) &
                  ,outnice              (:,:) &
                  ,outbuoy              (:,:) &
                  ,outmpqi            (:,:,:) &
                  ,outmpql            (:,:,:) &
                  ,outmpcf            (:,:,:) &
                  ,out_chem           (:,:,:) &
                  !- for convective transport
                !   ,ierr4d               (:,j,plume) &
                !   ,jmin4d               (:,j,plume) &
                !   ,klcl4d               (:,j,plume) &
                !   ,k224d                (:,j,plume) &
                !   ,kbcon4d              (:,j,plume) &
                !   ,ktop4d               (:,j,plume) &
                !   ,kstabi4d             (:,j,plume) &
                !   ,kstabm4d             (:,j,plume) &
                !   ,cprr4d               (:,j,plume) &
                !   ,xmb4d                (:,j,plume) &
                !   ,edt4d                (:,j,plume) &
                !   ,pwav4d               (:,j,plume) &
                !   ,sigma4d              (:,j,plume) &
                  ,ierr4d               (:) &
                  ,jmin4d               (:) &
                  ,klcl4d               (:) &
                  ,k224d                (:) &
                  ,kbcon4d              (:) &
                  ,ktop4d               (:) &
                  ,kstabi4d             (:) &
                  ,kstabm4d             (:) &
                  ,cprr4d               (:) &
                  ,xmb4d                (:) &
                  ,edt4d                (:) &
                  ,pwav4d               (:) &
                  ,sigma4d              (:) &
                !   ,pcup5d             (:,j,:,plume) &
                !   ,up_massentr5d      (:,j,:,plume) &
                !   ,up_massdetr5d      (:,j,:,plume) &
                !   ,dd_massentr5d      (:,j,:,plume) &
                !   ,dd_massdetr5d      (:,j,:,plume) &
                !   ,zup5d              (:,j,:,plume) &
                !   ,zdn5d              (:,j,:,plume) &
                !   ,prup5d             (:,j,:,plume) &
                !   ,prdn5d             (:,j,:,plume) &
                !   ,clwup5d            (:,j,:,plume) &
                !   ,tup5d              (:,j,:,plume) &
                !   ,conv_cld_fr5d      (:,j,:,plume) &
                  ,pcup5d             (:,:) &
                  ,up_massentr5d      (:,:) &
                  ,up_massdetr5d      (:,:) &
                  ,dd_massentr5d      (:,:) &
                  ,dd_massdetr5d      (:,:) &
                  ,zup5d              (:,:) &
                  ,zdn5d              (:,:) &
                  ,prup5d             (:,:) &
                  ,prdn5d             (:,:) &
                  ,clwup5d            (:,:) &
                  ,tup5d              (:,:) &
                  ,conv_cld_fr5d      (:,:) &
                  !-- for debug/diag
                !   ,AA0(:,j),AA1(:,j),AA2(:,j),AA3(:,j),AA1_BL(:,j),AA1_CIN(:,j),TAU_BL(:,j),TAU_EC(:,j) &
                  ,AA0(:),AA1(:),AA2(:),AA3(:),AA1_BL(:),AA1_CIN(:),TAU_BL(:),TAU_EC(:) &
                  !-- for diag
                !   ,lightn_dens  (:,j)               &
                !   ,var2d        (:,j)               &
                  ,lightn_dens  (:)               &
                  ,var2d        (:)               &
                  ,revsu_gf_2d                      &
                  ,prfil_gf_2d                      &
                  ,var3d_agf_2d                     &
                  ,var3d_bgf_2d                     &
                  ,Tpert_2d                         &
        )

        call end_timing()

        call print_timing()

        open(newunit=fileID, file=trim(dirName) // '/kpbli_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) kpbli_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out kpbli = ', sum(kpbli_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/cum_ztexec_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) cum_ztexec_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out cum_ztexec = ', sum(cum_ztexec_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/cum_zqexec_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) cum_zqexec_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out cum_zqexec = ', sum(cum_zqexec_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/ccn_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) ccn_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out ccn = ', sum(ccn_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/rhoi_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) rhoi_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out rhoi = ', sum(rhoi_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/omeg_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) omeg_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out omeg = ', sum(omeg_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/temp_old_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) temp_old_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out temp_old = ', sum(temp_old_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/qv_old_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) qv_old_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out qv_old = ', sum(qv_old_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/ter11_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) ter11_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out ter11 = ', sum(ter11_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/h_sfc_flux_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) h_sfc_flux_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out h_sfc_flux = ', sum(h_sfc_flux_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/le_sfc_flux_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) le_sfc_flux_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out le_sfc_flux = ', sum(le_sfc_flux_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/xlons_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) xlons_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out xlons = ', sum(xlons_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/xlats_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) xlats_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out xlats = ', sum(xlats_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/xlandi_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) xlandi_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out xlandi = ', sum(xlandi_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/temp_new_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) temp_new_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out temp_new = ', sum(temp_new_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/qv_new_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) qv_new_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out qv_new = ', sum(qv_new_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/temp_new_BL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) temp_new_BL_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out temp_new_BL = ', sum(temp_new_BL_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/qv_new_BL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) qv_new_BL_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out qv_new_BL = ', sum(qv_new_BL_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/temp_new_ADV_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) temp_new_ADV_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out temp_new_ADV = ', sum(temp_new_ADV_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/qv_new_ADV_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) qv_new_ADV_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out qv_new_ADV = ', sum(qv_new_ADV_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/zo_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) zo_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out zo = ', sum(zo_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/po_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) po_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out po = ', sum(po_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/tsur_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) tsur_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out tsur = ', sum(tsur_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/psur_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) psur_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out psur = ', sum(psur_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/us_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) us_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out us = ', sum(us_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/vs_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) vs_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out vs = ', sum(vs_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/dm2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) dm2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out dm2d = ', sum(dm2d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/se_chem_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) se_chem_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out se_chem = ', sum(se_chem_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/zws_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) zws_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out zws = ', sum(zws_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/dhdt_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) dhdt_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out dhdt = ', sum(dhdt_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/buoy_exc2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) buoy_exc2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out buoy_exc2d = ', sum(buoy_exc2d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/mpqi_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) mpqi_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out mpqi = ', sum(mpqi_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/mpql_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) mpql_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out mpql = ', sum(mpql_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/mpcf_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) mpcf_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out mpcf = ', sum(mpcf_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/last_ierr_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) last_ierr_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out last_ierr = ', sum(last_ierr_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outt_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outt_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outt = ', sum(outt_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outq_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outq_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outq = ', sum(outq_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outqc_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outqc_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outqc = ', sum(outqc_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outu_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outu_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outu = ', sum(outu_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outv_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outv_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outv = ', sum(outv_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outnliq_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outnliq_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outnliq = ', sum(outnliq_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outnice_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outnice_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outnice = ', sum(outnice_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outbuoy_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outbuoy_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outbuoy = ', sum(outbuoy_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outmpqi_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outmpqi_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outmpqi = ', sum(outmpqi_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outmpql_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outmpql_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outmpql = ', sum(outmpql_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/outmpcf_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) outmpcf_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out outmpcf = ', sum(outmpcf_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/out_chem_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) out_chem_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out out_chem = ', sum(out_chem_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/ierr4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) ierr4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out ierr4d = ', sum(ierr4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/jmin4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) jmin4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out jmin4d = ', sum(jmin4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/klcl4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) klcl4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out klcl4d = ', sum(klcl4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/k224d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) k224d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out k224d  = ', sum(k224d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/kbcon4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) kbcon4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out kbcon4d = ', sum(kbcon4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/ktop4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) ktop4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out ktop4d = ', sum(ktop4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/kstabi4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) kstabi4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out kstabi4d = ', sum(kstabi4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/kstabm4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) kstabm4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out kstabm4d = ', sum(kstabm4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/xmb4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) xmb4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out xmb4d  = ', sum(xmb4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/edt4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) edt4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out edt4d  = ', sum(edt4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/pwav4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) pwav4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out pwav4d = ', sum(pwav4d)
        
        open(newunit=fileID, file=trim(dirName) // '/pcup5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) pcup5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out pcup5d = ', sum(pcup5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/up_massentr5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) up_massentr5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out up_massentr5d = ', sum(up_massentr5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/up_massdetr5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) up_massdetr5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out up_massdetr5d = ', sum(up_massdetr5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/dd_massentr5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) dd_massentr5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out dd_massentr5d = ', sum(dd_massentr5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/dd_massdetr5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) dd_massdetr5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out dd_massdetr5d = ', sum(dd_massdetr5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/zup5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) zup5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out zup5d = ', sum(zup5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/zdn5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) zdn5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out zdn5d = ', sum(zdn5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/prup5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) prup5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out prup5d = ', sum(prup5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/prdn5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) prdn5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out prdn5d = ', sum(prdn5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/clwup5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) clwup5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out clwup5d = ', sum(clwup5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/tup5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) tup5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out tup5d = ', sum(tup5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/conv_cld_fr5d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) conv_cld_fr5d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out conv_cld_fr5d = ', sum(conv_cld_fr5d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/AA0_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) AA0_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out AA0 = ', sum(AA0_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/AA1_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) AA1_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out AA1 = ', sum(AA1_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/AA2_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) AA2_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out AA2 = ', sum(AA2_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/AA3_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) AA3_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out AA3 = ', sum(AA3_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/AA1_BL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) AA1_BL_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out AA1_BL = ', sum(AA1_BL_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/AA1_CIN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) AA1_CIN_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out AA1_CIN = ', sum(AA1_CIN_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/TAU_BL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) TAU_BL_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out TAU_BL = ', sum(TAU_BL_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/TAU_EC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) TAU_EC_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out TAU_EC = ', sum(TAU_EC_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/revsu_gf_2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) revsu_gf_2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out revsu_gf_2d = ', sum(revsu_gf_2d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/prfil_gf_2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) prfil_gf_2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out prfil_gf_2d = ', sum(prfil_gf_2d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/var3d_agf_2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) var3d_agf_2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out var3d_agf_2d = ', sum(var3d_agf_2d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/var3d_bgf_2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) var3d_bgf_2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out var3d_bgf_2d = ', sum(var3d_bgf_2d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/Tpert_2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) Tpert_2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out Tpert_2d = ', sum(Tpert_2d_ref)

        open(newunit=fileID, file=trim(dirName) // '/cprr4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) cprr4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out cprr4d = ', sum(cprr4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/sigma4d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) sigma4d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out sigma4d = ', sum(sigma4d_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/lightn_dens_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) lightn_dens_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out lightn_dens = ', sum(lightn_dens_ref)
        
        open(newunit=fileID, file=trim(dirName) // '/var2d_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) var2d_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out var2d = ', sum(var2d_ref)

        print*,'***'
        print*,'Compare sum(diff(kpbli)) = ',sum(kpbli_ref-kpbli)
        print*,'sum(diff(kpbli_ref)) = ',sum(kpbli_ref)
        print*,'sum(diff(kpbli)) = ',sum(kpbli)
        print*,'***'
        print*,'Compare sum(diff(cum_ztexec)) = ',sum(cum_ztexec_ref-cum_ztexec)
        print*,'sum(diff(cum_ztexec_ref)) = ',sum(cum_ztexec_ref)
        print*,'sum(diff(cum_ztexec)) = ',sum(cum_ztexec)
        print*,'***'
        print*,'Compare sum(diff(cum_zqexec)) = ',sum(cum_zqexec_ref-cum_zqexec)
        print*,'sum(diff(cum_zqexec_ref)) = ',sum(cum_zqexec_ref)
        print*,'sum(diff(cum_zqexec)) = ',sum(cum_zqexec)
        print*,'***'
        print*,'Compare sum(diff(ccn)) = ',sum(ccn_ref-ccn)
        print*,'sum(diff(ccn_ref)) = ',sum(ccn_ref)
        print*,'sum(diff(ccn)) = ',sum(ccn)
        print*,'***'
        print*,'Compare sum(diff(rhoi)) = ',sum(rhoi_ref-rhoi)
        print*,'sum(diff(rhoi_ref)) = ',sum(rhoi_ref)
        print*,'sum(diff(rhoi)) = ',sum(rhoi)
        print*,'***'
        print*,'Compare sum(diff(omeg)) = ',sum(omeg_ref-omeg)
        print*,'sum(diff(omeg_ref)) = ',sum(omeg_ref)
        print*,'sum(diff(omeg)) = ',sum(omeg)
        print*,'***'
        print*,'Compare sum(diff(temp_old)) = ',sum(temp_old_ref-temp_old)
        print*,'sum(diff(temp_old_ref)) = ',sum(temp_old_ref)
        print*,'sum(diff(temp_old)) = ',sum(temp_old)
        print*,'***'
        print*,'Compare sum(diff(qv_old)) = ',sum(qv_old_ref-qv_old)
        print*,'sum(diff(qv_old_ref)) = ',sum(qv_old_ref)
        print*,'sum(diff(qv_old)) = ',sum(qv_old)
        print*,'***'
        print*,'Compare sum(diff(ter11)) = ',sum(ter11_ref-ter11)
        print*,'sum(diff(ter11_ref)) = ',sum(ter11_ref)
        print*,'sum(diff(ter11)) = ',sum(ter11)
        print*,'***'
        print*,'Compare sum(diff(h_sfc_flux)) = ',sum(h_sfc_flux_ref-h_sfc_flux)
        print*,'sum(diff(h_sfc_flux_ref)) = ',sum(h_sfc_flux_ref)
        print*,'sum(diff(h_sfc_flux)) = ',sum(h_sfc_flux)
        print*,'***'
        print*,'Compare sum(diff(le_sfc_flux)) = ',sum(le_sfc_flux_ref-le_sfc_flux)
        print*,'sum(diff(le_sfc_flux_ref)) = ',sum(le_sfc_flux_ref)
        print*,'sum(diff(le_sfc_flux)) = ',sum(le_sfc_flux)
        print*,'***'
        print*,'Compare sum(diff(xlons)) = ',sum(xlons_ref-xlons)
        print*,'sum(diff(xlons_ref)) = ',sum(xlons_ref)
        print*,'sum(diff(xlons)) = ',sum(xlons)
        print*,'***'
        print*,'Compare sum(diff(xlats)) = ',sum(xlats_ref-xlats)
        print*,'sum(diff(xlats_ref)) = ',sum(xlats_ref)
        print*,'sum(diff(xlats)) = ',sum(xlats)
        print*,'***'
        print*,'Compare sum(diff(xlandi)) = ',sum(xlandi_ref-xlandi)
        print*,'sum(diff(xlandi_ref)) = ',sum(xlandi_ref)
        print*,'sum(diff(xlandi)) = ',sum(xlandi)
        print*,'***'
        print*,'Compare sum(diff(temp_new)) = ',sum(temp_new_ref-temp_new)
        print*,'sum(diff(temp_new_ref)) = ',sum(temp_new_ref)
        print*,'sum(diff(temp_new)) = ',sum(temp_new)
        print*,'***'
        print*,'Compare sum(diff(qv_new)) = ',sum(qv_new_ref-qv_new)
        print*,'sum(diff(qv_new_ref)) = ',sum(qv_new_ref)
        print*,'sum(diff(qv_new)) = ',sum(qv_new)
        print*,'***'
        print*,'Compare sum(diff(temp_new_BL)) = ',sum(temp_new_BL_ref-temp_new_BL)
        print*,'sum(diff(temp_new_BL_ref)) = ',sum(temp_new_BL_ref)
        print*,'sum(diff(temp_new_BL)) = ',sum(temp_new_BL)
        print*,'***'
        print*,'Compare sum(diff(qv_new_BL)) = ',sum(qv_new_BL_ref-qv_new_BL)
        print*,'sum(diff(qv_new_BL_ref)) = ',sum(qv_new_BL_ref)
        print*,'sum(diff(qv_new_BL)) = ',sum(qv_new_BL)
        print*,'***'
        print*,'Compare sum(diff(temp_new_ADV)) = ',sum(temp_new_ADV_ref-temp_new_ADV)
        print*,'sum(diff(temp_new_ADV_ref)) = ',sum(temp_new_ADV_ref)
        print*,'sum(diff(temp_new_ADV)) = ',sum(temp_new_ADV)
        print*,'***'
        print*,'Compare sum(diff(qv_new_ADV)) = ',sum(qv_new_ADV_ref-qv_new_ADV)
        print*,'sum(diff(qv_new_ADV_ref)) = ',sum(qv_new_ADV_ref)
        print*,'sum(diff(qv_new_ADV)) = ',sum(qv_new_ADV)
        print*,'***'
        print*,'Compare sum(diff(zo)) = ',sum(zo_ref-zo)
        print*,'sum(diff(zo_ref)) = ',sum(zo_ref)
        print*,'sum(diff(zo)) = ',sum(zo)
        print*,'***'
        print*,'Compare sum(diff(po)) = ',sum(po_ref-po)
        print*,'sum(diff(po_ref)) = ',sum(po_ref)
        print*,'sum(diff(po)) = ',sum(po)
        print*,'***'
        print*,'Compare sum(diff(tsur)) = ',sum(tsur_ref-tsur)
        print*,'sum(diff(tsur_ref)) = ',sum(tsur_ref)
        print*,'sum(diff(tsur)) = ',sum(tsur)
        print*,'***'
        print*,'Compare sum(diff(psur)) = ',sum(psur_ref-psur)
        print*,'sum(diff(psur_ref)) = ',sum(psur_ref)
        print*,'sum(diff(psur)) = ',sum(psur)
        print*,'***'
        print*,'Compare sum(diff(us)) = ',sum(us_ref-us)
        print*,'sum(diff(us_ref)) = ',sum(us_ref)
        print*,'sum(diff(us)) = ',sum(us)
        print*,'***'
        print*,'Compare sum(diff(vs)) = ',sum(vs_ref-vs)
        print*,'sum(diff(vs_ref)) = ',sum(vs_ref)
        print*,'sum(diff(vs)) = ',sum(vs)
        print*,'***'
        print*,'Compare sum(diff(dm2d)) = ',sum(dm2d_ref-dm2d)
        print*,'sum(diff(dm2d_ref)) = ',sum(dm2d_ref)
        print*,'sum(diff(dm2d)) = ',sum(dm2d)
        print*,'***'
        print*,'Compare sum(diff(se_chem)) = ',sum(se_chem_ref-se_chem)
        print*,'sum(diff(se_chem_ref)) = ',sum(se_chem_ref)
        print*,'sum(diff(se_chem)) = ',sum(se_chem)
        print*,'***'
        print*,'Compare sum(diff(zws)) = ',sum(zws_ref-zws)
        print*,'sum(diff(zws_ref)) = ',sum(zws_ref)
        print*,'sum(diff(zws)) = ',sum(zws)
        print*,'***'
        print*,'Compare sum(diff(dhdt)) = ',sum(dhdt_ref-dhdt)
        print*,'sum(diff(dhdt_ref)) = ',sum(dhdt_ref)
        print*,'sum(diff(dhdt)) = ',sum(dhdt)
        print*,'***'
        print*,'Compare sum(diff(buoy_exc2d)) = ',sum(buoy_exc2d_ref-buoy_exc2d)
        print*,'sum(diff(buoy_exc2d_ref)) = ',sum(buoy_exc2d_ref)
        print*,'sum(diff(buoy_exc2d)) = ',sum(buoy_exc2d)
        print*,'***'
        print*,'Compare sum(diff(mpqi)) = ',sum(mpqi_ref-mpqi)
        print*,'sum(diff(mpqi_ref)) = ',sum(mpqi_ref)
        print*,'sum(diff(mpqi)) = ',sum(mpqi)
        print*,'***'
        print*,'Compare sum(diff(mpql)) = ',sum(mpql_ref-mpql)
        print*,'sum(diff(mpql_ref)) = ',sum(mpql_ref)
        print*,'sum(diff(mpql)) = ',sum(mpql)
        print*,'***'
        print*,'Compare sum(diff(mpcf)) = ',sum(mpcf_ref-mpcf)
        print*,'sum(diff(mpcf_ref)) = ',sum(mpcf_ref)
        print*,'sum(diff(mpcf)) = ',sum(mpcf)
        print*,'***'
        print*,'Compare sum(diff(last_ierr)) = ',sum(last_ierr_ref-last_ierr)
        print*,'sum(diff(last_ierr_ref)) = ',sum(last_ierr_ref)
        print*,'sum(diff(last_ierr)) = ',sum(last_ierr)
        print*,'***'
        print*,'Compare sum(diff(outt)) = ',sum(outt_ref-outt)
        print*,'sum(diff(outt_ref)) = ',sum(outt_ref)
        print*,'sum(diff(outt)) = ',sum(outt)
        print*,'***'
        print*,'Compare sum(diff(outq)) = ',sum(outq_ref-outq)
        print*,'sum(diff(outq_ref)) = ',sum(outq_ref)
        print*,'sum(diff(outq)) = ',sum(outq)
        print*,'***'
        print*,'Compare sum(diff(outqc)) = ',sum(outqc_ref-outqc)
        print*,'sum(diff(outqc_ref)) = ',sum(outqc_ref)
        print*,'sum(diff(outqc)) = ',sum(outqc)
        print*,'***'
        print*,'Compare sum(diff(outu)) = ',sum(outu_ref-outu)
        print*,'sum(diff(outu_ref)) = ',sum(outu_ref)
        print*,'sum(diff(outu)) = ',sum(outu)
        print*,'***'
        print*,'Compare sum(diff(outv)) = ',sum(outv_ref-outv)
        print*,'sum(diff(outv_ref)) = ',sum(outv_ref)
        print*,'sum(diff(outv)) = ',sum(outv)
        print*,'***'
        print*,'Compare sum(diff(outnliq)) = ',sum(outnliq_ref-outnliq)
        print*,'sum(diff(outnliq_ref)) = ',sum(outnliq_ref)
        print*,'sum(diff(outnliq)) = ',sum(outnliq)
        print*,'***'
        print*,'Compare sum(diff(outnice)) = ',sum(outnice_ref-outnice)
        print*,'sum(diff(outnice_ref)) = ',sum(outnice_ref)
        print*,'sum(diff(outnice)) = ',sum(outnice)
        print*,'***'
        print*,'Compare sum(diff(outbuoy)) = ',sum(outbuoy_ref-outbuoy)
        print*,'sum(diff(outbuoy_ref)) = ',sum(outbuoy_ref)
        print*,'sum(diff(outbuoy)) = ',sum(outbuoy)
        print*,'***'
        print*,'Compare sum(diff(outmpqi)) = ',sum(outmpqi_ref-outmpqi)
        print*,'sum(diff(outmpqi_ref)) = ',sum(outmpqi_ref)
        print*,'sum(diff(outmpqi)) = ',sum(outmpqi)
        print*,'***'
        print*,'Compare sum(diff(outmpql)) = ',sum(outmpql_ref-outmpql)
        print*,'sum(diff(outmpql_ref)) = ',sum(outmpql_ref)
        print*,'sum(diff(outmpql)) = ',sum(outmpql)
        print*,'***'
        print*,'Compare sum(diff(outmpcf)) = ',sum(outmpcf_ref-outmpcf)
        print*,'sum(diff(outmpcf_ref)) = ',sum(outmpcf_ref)
        print*,'sum(diff(outmpcf)) = ',sum(outmpcf)
        print*,'***'
        print*,'Compare sum(diff(out_chem)) = ',sum(out_chem_ref-out_chem)
        print*,'sum(diff(out_chem_ref)) = ',sum(out_chem_ref)
        print*,'sum(diff(out_chem)) = ',sum(out_chem)
        print*,'***'
        print*,'Compare sum(diff(ierr4d)) = ',sum(ierr4d_ref-ierr4d)
        print*,'sum(diff(ierr4d_ref)) = ',sum(ierr4d_ref)
        print*,'sum(diff(ierr4d)) = ',sum(ierr4d)
        print*,'***'
        print*,'Compare sum(diff(jmin4d)) = ',sum(jmin4d_ref-jmin4d)
        print*,'sum(diff(jmin4d_ref)) = ',sum(jmin4d_ref)
        print*,'sum(diff(jmin4d)) = ',sum(jmin4d)
        print*,'***'
        print*,'Compare sum(diff(klcl4d)) = ',sum(klcl4d_ref-klcl4d)
        print*,'sum(diff(klcl4d_ref)) = ',sum(klcl4d_ref)
        print*,'sum(diff(klcl4d)) = ',sum(klcl4d)
        print*,'***'
        print*,'Compare sum(diff(k224d)) = ',sum(k224d_ref-k224d)
        print*,'sum(diff(k224d_ref)) = ',sum(k224d_ref)
        print*,'sum(diff(k224d)) = ',sum(k224d)
        print*,'***'
        print*,'Compare sum(diff(kbcon4d)) = ',sum(kbcon4d_ref-kbcon4d)
        print*,'sum(diff(kbcon4d_ref)) = ',sum(kbcon4d_ref)
        print*,'sum(diff(kbcon4d)) = ',sum(kbcon4d)
        print*,'***'
        print*,'Compare sum(diff(ktop4d)) = ',sum(ktop4d_ref-ktop4d)
        print*,'sum(diff(ktop4d_ref)) = ',sum(ktop4d_ref)
        print*,'sum(diff(ktop4d)) = ',sum(ktop4d)
        print*,'***'
        print*,'Compare sum(diff(kstabi4d)) = ',sum(kstabi4d_ref-kstabi4d)
        print*,'sum(diff(kstabi4d_ref)) = ',sum(kstabi4d_ref)
        print*,'sum(diff(kstabi4d)) = ',sum(kstabi4d)
        print*,'***'
        print*,'Compare sum(diff(kstabm4d)) = ',sum(kstabm4d_ref-kstabm4d)
        print*,'sum(diff(kstabm4d_ref)) = ',sum(kstabm4d_ref)
        print*,'sum(diff(kstabm4d)) = ',sum(kstabm4d)
        print*,'***'
        print*,'Compare sum(diff(xmb4d)) = ',sum(xmb4d_ref-xmb4d)
        print*,'sum(diff(xmb4d_ref)) = ',sum(xmb4d_ref)
        print*,'sum(diff(xmb4d)) = ',sum(xmb4d)
        print*,'***'
        print*,'Compare sum(diff(edt4d)) = ',sum(edt4d_ref-edt4d)
        print*,'sum(diff(edt4d_ref)) = ',sum(edt4d_ref)
        print*,'sum(diff(edt4d)) = ',sum(edt4d)
        print*,'***'
        print*,'Compare sum(diff(pwav4d)) = ',sum(pwav4d_ref-pwav4d)
        print*,'sum(diff(pwav4d_ref)) = ',sum(pwav4d_ref)
        print*,'sum(diff(pwav4d)) = ',sum(pwav4d)
        print*,'***'
        print*,'Compare sum(diff(pcup5d)) = ',sum(pcup5d_ref-pcup5d)
        print*,'sum(diff(pcup5d_ref)) = ',sum(pcup5d_ref)
        print*,'sum(diff(pcup5d)) = ',sum(pcup5d)
        print*,'***'
        print*,'Compare sum(diff(up_massentr5d)) = ',sum(up_massentr5d_ref-up_massentr5d)
        print*,'sum(diff(up_massentr5d_ref)) = ',sum(up_massentr5d_ref)
        print*,'sum(diff(up_massentr5d)) = ',sum(up_massentr5d)
        print*,'***'
        print*,'Compare sum(diff(up_massdetr5d)) = ',sum(up_massdetr5d_ref-up_massdetr5d)
        print*,'sum(diff(up_massdetr5d_ref)) = ',sum(up_massdetr5d_ref)
        print*,'sum(diff(up_massdetr5d)) = ',sum(up_massdetr5d)
        print*,'***'
        print*,'Compare sum(diff(dd_massentr5d)) = ',sum(dd_massentr5d_ref-dd_massentr5d)
        print*,'sum(diff(dd_massentr5d_ref)) = ',sum(dd_massentr5d_ref)
        print*,'sum(diff(dd_massentr5d)) = ',sum(dd_massentr5d)
        print*,'***'
        print*,'Compare sum(diff(dd_massdetr5d)) = ',sum(dd_massdetr5d_ref-dd_massdetr5d)
        print*,'sum(diff(dd_massdetr5d_ref)) = ',sum(dd_massdetr5d_ref)
        print*,'sum(diff(dd_massdetr5d)) = ',sum(dd_massdetr5d)
        print*,'***'
        print*,'Compare sum(diff(zup5d)) = ',sum(zup5d_ref-zup5d)
        print*,'sum(diff(zup5d_ref)) = ',sum(zup5d_ref)
        print*,'sum(diff(zup5d)) = ',sum(zup5d)
        print*,'***'
        print*,'Compare sum(diff(zdn5d)) = ',sum(zdn5d_ref-zdn5d)
        print*,'sum(diff(zdn5d_ref)) = ',sum(zdn5d_ref)
        print*,'sum(diff(zdn5d)) = ',sum(zdn5d)
        print*,'***'
        print*,'Compare sum(diff(prup5d)) = ',sum(prup5d_ref-prup5d)
        print*,'sum(diff(prup5d_ref)) = ',sum(prup5d_ref)
        print*,'sum(diff(prup5d)) = ',sum(prup5d)
        print*,'***'
        print*,'Compare sum(diff(prdn5d)) = ',sum(prdn5d_ref-prdn5d)
        print*,'sum(diff(prdn5d_ref)) = ',sum(prdn5d_ref)
        print*,'sum(diff(prdn5d)) = ',sum(prdn5d)
        print*,'***'
        print*,'Compare sum(diff(clwup5d)) = ',sum(clwup5d_ref-clwup5d)
        print*,'sum(diff(clwup5d_ref)) = ',sum(clwup5d_ref)
        print*,'sum(diff(clwup5d)) = ',sum(clwup5d)
        print*,'***'
        print*,'Compare sum(diff(tup5d)) = ',sum(tup5d_ref-tup5d)
        print*,'sum(diff(tup5d_ref)) = ',sum(tup5d_ref)
        print*,'sum(diff(tup5d)) = ',sum(tup5d)
        print*,'***'
        print*,'Compare sum(diff(conv_cld_fr5d)) = ',sum(conv_cld_fr5d_ref-conv_cld_fr5d)
        print*,'sum(diff(conv_cld_fr5d_ref)) = ',sum(conv_cld_fr5d_ref)
        print*,'sum(diff(conv_cld_fr5d)) = ',sum(conv_cld_fr5d)
        print*,'***'
        print*,'Compare sum(diff(AA0)) = ',sum(AA0_ref-AA0)
        print*,'sum(diff(AA0_ref)) = ',sum(AA0_ref)
        print*,'sum(diff(AA0)) = ',sum(AA0)
        print*,'***'
        print*,'Compare sum(diff(AA1)) = ',sum(AA1_ref-AA1)
        print*,'sum(diff(AA1_ref)) = ',sum(AA1_ref)
        print*,'sum(diff(AA1)) = ',sum(AA1)
        print*,'***'
        print*,'Compare sum(diff(AA2)) = ',sum(AA2_ref-AA2)
        print*,'sum(diff(AA2_ref)) = ',sum(AA2_ref)
        print*,'sum(diff(AA2)) = ',sum(AA2)
        print*,'***'
        print*,'Compare sum(diff(AA3)) = ',sum(AA3_ref-AA3)
        print*,'sum(diff(AA3_ref)) = ',sum(AA3_ref)
        print*,'sum(diff(AA3)) = ',sum(AA3)
        print*,'***'
        print*,'Compare sum(diff(AA1_BL)) = ',sum(AA1_BL_ref-AA1_BL)
        print*,'sum(diff(AA1_BL_ref)) = ',sum(AA1_BL_ref)
        print*,'sum(diff(AA1_BL)) = ',sum(AA1_BL)
        print*,'***'
        print*,'Compare sum(diff(AA1_CIN)) = ',sum(AA1_CIN_ref-AA1_CIN)
        print*,'sum(diff(AA1_CIN_ref)) = ',sum(AA1_CIN_ref)
        print*,'sum(diff(AA1_CIN)) = ',sum(AA1_CIN)
        print*,'***'
        print*,'Compare sum(diff(TAU_BL)) = ',sum(TAU_BL_ref-TAU_BL)
        print*,'sum(diff(TAU_BL_ref)) = ',sum(TAU_BL_ref)
        print*,'sum(diff(TAU_BL)) = ',sum(TAU_BL)
        print*,'***'
        print*,'Compare sum(diff(TAU_EC)) = ',sum(TAU_EC_ref-TAU_EC)
        print*,'sum(diff(TAU_EC_ref)) = ',sum(TAU_EC_ref)
        print*,'sum(diff(TAU_EC)) = ',sum(TAU_EC)
        print*,'***'
        print*,'Compare sum(diff(revsu_gf_2d)) = ',sum(revsu_gf_2d_ref-revsu_gf_2d)
        print*,'sum(diff(revsu_gf_2d_ref)) = ',sum(revsu_gf_2d_ref)
        print*,'sum(diff(revsu_gf_2d)) = ',sum(revsu_gf_2d)
        print*,'***'
        print*,'Compare sum(diff(prfil_gf_2d)) = ',sum(prfil_gf_2d_ref-prfil_gf_2d)
        print*,'sum(diff(prfil_gf_2d_ref)) = ',sum(prfil_gf_2d_ref)
        print*,'sum(diff(prfil_gf_2d)) = ',sum(prfil_gf_2d)
        print*,'***'
        print*,'Compare sum(diff(var3d_agf_2d)) = ',sum(var3d_agf_2d_ref-var3d_agf_2d)
        print*,'sum(diff(var3d_agf_2d_ref)) = ',sum(var3d_agf_2d_ref)
        print*,'sum(diff(var3d_agf_2d)) = ',sum(var3d_agf_2d)
        print*,'***'
        print*,'Compare sum(diff(var3d_bgf_2d)) = ',sum(var3d_bgf_2d_ref-var3d_bgf_2d)
        print*,'sum(diff(var3d_bgf_2d_ref)) = ',sum(var3d_bgf_2d_ref)
        print*,'sum(diff(var3d_bgf_2d)) = ',sum(var3d_bgf_2d)
        print*,'***'
        print*,'Compare sum(diff(Tpert_2d)) = ',sum(Tpert_2d_ref-Tpert_2d)
        print*,'sum(diff(Tpert_2d_ref)) = ',sum(Tpert_2d_ref)
        print*,'sum(diff(Tpert_2d)) = ',sum(Tpert_2d)
        print*,'***'
        print*,'Compare sum(diff(cprr4d)) = ',sum(cprr4d_ref-cprr4d)
        print*,'sum(diff(cprr4d_ref)) = ',sum(cprr4d_ref)
        print*,'sum(diff(cprr4d)) = ',sum(cprr4d)
        print*,'***'
        print*,'Compare sum(diff(sigma4d)) = ',sum(sigma4d_ref-sigma4d)
        print*,'sum(diff(sigma4d_ref)) = ',sum(sigma4d_ref)
        print*,'sum(diff(sigma4d)) = ',sum(sigma4d)
        print*,'***'
        print*,'Compare sum(diff(lightn_dens)) = ',sum(lightn_dens_ref-lightn_dens)
        print*,'sum(diff(lightn_dens_ref)) = ',sum(lightn_dens_ref)
        print*,'sum(diff(lightn_dens)) = ',sum(lightn_dens)
        print*,'***'
        print*,'Compare sum(diff(var2d)) = ',sum(var2d_ref-var2d)
        print*,'sum(diff(var2d_ref)) = ',sum(var2d_ref)
        print*,'sum(diff(var2d)) = ',sum(var2d)
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