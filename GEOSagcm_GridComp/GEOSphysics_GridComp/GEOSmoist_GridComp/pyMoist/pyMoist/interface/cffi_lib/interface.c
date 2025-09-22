#include <stdio.h>
#include <time.h>
#include "moist.h"

extern int pymoist_interface_c_init(moist_flags_t *flags)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, pymoist interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = pymoist_interface_py_init(flags);

    if (return_code < 0)
    {
        exit(return_code);
    }
}

extern int gfdl_1m_interface_c_init(gfdl_1m_flags_t *flags)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, gfdl_1m interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = gfdl_1m_interface_py_init(flags);

    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_run_GFDL1M(
    float dw_land, float dw_ocean, int PDFSHAPE, float TURNRHCRIT_PARAM,
    float DT_MOIST, float CCW_EVAP_EFF, float CCI_EVAP_EFF,
    int LMELTFRZ,
    float *AREA, float *CNV_FRC, float *SRF_TYPE, // 2D...
    int *KLCL,
    float *EIS, float *PLmb, float *PLEmb, float *NACTL, float *NACTI, float *QST, // 3D...
    float *T, float *Q, float *QLCN, float *QICN, float *QLLS, float *QILS, float *CLLS, float *CLCN,
    float *SUBLC, float *EVAPC, float *RHX)
{
    printf("> C LMELTFRZ = %d\n", LMELTFRZ);
    int return_code = pymoist_interface_py_run_GFDL1M(
        dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM,
        DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF,
        LMELTFRZ != 0,
        AREA, CNV_FRC, SRF_TYPE,
        KLCL,
        EIS, PLmb, PLEmb, NACTL, NACTI, QST,
        T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN,
        SUBLC, EVAPC, RHX);
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_run_GFDL_1M_driver(
    float *RAD_QV, float *RAD_QL, float *RAD_QR, float *RAD_QI, float *RAD_QS, float *RAD_QG, float *RAD_CF, float *NACTAll,
    float *DQVDTmic, float *DQLDTmic, float *DQRDTmic, float *DQIDTmic,
    float *DQSDTmic, float *DQGDTmic, float *DQADTmic, float *DTDTmic,
    float *T, float *W, float *U, float *V, float *DUDTmic, float *DVDTmic, float *DZ, float *DP,
    float *AREA, float *FRLAND, float *CNV_FRC, float *SRF_TYPE, float *EIS, float *RHCRIT3D,
    float DT_MOIST, float ANV_ICEFALL, float LS_ICEFALL,
    float *REV_LS, float *RSU_LS,
    float *PRCP_RAIN, float *PRCP_SNOW, float *PRCP_ICE, float *PRCP_GRAUPEL,
    float *PFL_LS, float *PFI_LS,
    int LHYDROSTATIC, int LPHYS_HYDROSTATIC)
{
    int return_code = pymoist_interface_py_run_GFDL_1M_driver(
        RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTAll,
        DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, DQSDTmic, DQGDTmic, DQADTmic, DTDTmic,
        T, W, U, V, DUDTmic, DVDTmic, DZ, DP, AREA, FRLAND, CNV_FRC, SRF_TYPE, EIS, RHCRIT3D,
        DT_MOIST, ANV_ICEFALL, LS_ICEFALL,
        REV_LS, RSU_LS, PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, PFL_LS, PFI_LS,
        LHYDROSTATIC, LPHYS_HYDROSTATIC);
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_run_AerActivation(
    // input
    float *aero_dgn, float *aero_num, float *aero_hygroscopicity, float *aero_sigma,
    float *frland, float nn_ocean, float nn_land,
    float *t, float *plo,
    float *qicn, float *qils, float *qlcn, float *qlls,
    float *vvel, float *tke,
    // output
    float *nacti, float *nwfa, float *nactl)
{
    int return_code = pymoist_interface_py_run_AerActivation(
        aero_dgn,
        aero_num,
        aero_hygroscopicity,
        aero_sigma,
        frland,
        nn_ocean,
        nn_land,
        t,
        plo,
        qicn,
        qils,
        qlcn,
        qlls,
        vvel,
        tke,
        nacti,
        nwfa,
        nactl);
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_finalize()
{
    int return_code = pymoist_interface_py_finalize();
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void compute_uwshcu_c_init(
    // inputs
    int ncnst, int k0, int windsrcavg
    // inputs-outputs

    // outputs

)
{
    int return_code = compute_uwshcu_py_init(ncnst, k0, windsrcavg);
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void compute_uwshcu_c_run_compute_uwshcu(
    // inputs
    int dotransport, int k0, int windsrcavg, float qtsrchgt, float qtsrc_fac,
    float thlsrc_fac, float frc_rasn, float rbuoy, float epsvarw,
    int use_CINcin, float mumin1, float rmaxfrac, float PGFc, float dt,
    int niter_xc, float criqc, float rle, int cridist_opt, float mixscale,
    float rdrag, float rkm, int use_self_detrain, float detrhgt,
    int use_cumpenent, float rpen, int use_momenflx, float rdrop, int iter_cin,
    const float *pifc0_inv, int *pifc0_inv_dim_sizes, int pifc0_inv_rank,
    const float *zifc0_inv, int *zifc0_inv_dim_sizes, int zifc0_inv_rank,
    const float *pmid0_inv, int *pmid0_inv_dim_sizes, int pmid0_inv_rank,
    const float *zmid0_inv, int *zmid0_inv_dim_sizes, int zmid0_inv_rank,
    const float *kpbl_inv, int *kpbl_inv_dim_sizes, int kpbl_inv_rank,
    const float *exnmid0_inv, int *exnmid0_inv_dim_sizes, int exnmid0_inv_rank,
    const float *exnifc0_inv, int *exnifc0_inv_dim_sizes, int exnifc0_inv_rank,
    const float *dp0_inv, int *dp0_inv_dim_sizes, int dp0_inv_rank,
    const float *u0_inv, int *u0_inv_dim_sizes, int u0_inv_rank,
    const float *v0_inv, int *v0_inv_dim_sizes, int v0_inv_rank,
    const float *qv0_inv, int *qv0_inv_dim_sizes, int qv0_inv_rank,
    const float *ql0_inv, int *ql0_inv_dim_sizes, int ql0_inv_rank,
    const float *qi0_inv, int *qi0_inv_dim_sizes, int qi0_inv_rank,
    const float *t0_inv, int *t0_inv_dim_sizes, int t0_inv_rank,
    const float *frland_in, int *frland_in_dim_sizes, int frland_in_rank,
    const float *tke_inv, int *tke_inv_dim_sizes, int tke_inv_rank,
    const float *rkfre, int *rkfre_dim_sizes, int rkfre_rank, const float *cush,
    int *cush_dim_sizes, int cush_rank, const float *shfx, int *shfx_dim_sizes,
    int shfx_rank, const float *evap, int *evap_dim_sizes, int evap_rank,
    const float *cnvtr, int *cnvtr_dim_sizes, int cnvtr_rank,
    const float *CNV_Tracers, int *CNV_Tracers_dim_sizes, int CNV_Tracers_rank,
    // inputs-outputs

    // outputs
    float *umf_inv, int *umf_inv_dim_sizes, int umf_inv_rank, float *dcm_inv,
    int *dcm_inv_dim_sizes, int dcm_inv_rank, float *qtflx_inv,
    int *qtflx_inv_dim_sizes, int qtflx_inv_rank, float *slflx_inv,
    int *slflx_inv_dim_sizes, int slflx_inv_rank, float *uflx_inv,
    int *uflx_inv_dim_sizes, int uflx_inv_rank, float *vflx_inv,
    int *vflx_inv_dim_sizes, int vflx_inv_rank, float *qvten_inv,
    int *qvten_inv_dim_sizes, int qvten_inv_rank, float *qlten_inv,
    int *qlten_inv_dim_sizes, int qlten_inv_rank, float *qiten_inv,
    int *qiten_inv_dim_sizes, int qiten_inv_rank, float *tten_inv,
    int *tten_inv_dim_sizes, int tten_inv_rank, float *uten_inv,
    int *uten_inv_dim_sizes, int uten_inv_rank, float *vten_inv,
    int *vten_inv_dim_sizes, int vten_inv_rank, float *qrten_inv,
    int *qrten_inv_dim_sizes, int qrten_inv_rank, float *qsten_inv,
    int *qsten_inv_dim_sizes, int qsten_inv_rank, float *cufrc_inv,
    int *cufrc_inv_dim_sizes, int cufrc_inv_rank, float *fer_inv,
    int *fer_inv_dim_sizes, int fer_inv_rank, float *fdr_inv,
    int *fdr_inv_dim_sizes, int fdr_inv_rank, float *ndrop_inv,
    int *ndrop_inv_dim_sizes, int ndrop_inv_rank, float *nice_inv,
    int *nice_inv_dim_sizes, int nice_inv_rank, float *qldet_inv,
    int *qldet_inv_dim_sizes, int qldet_inv_rank, float *qlsub_inv,
    int *qlsub_inv_dim_sizes, int qlsub_inv_rank, float *qidet_inv,
    int *qidet_inv_dim_sizes, int qidet_inv_rank, float *qisub_inv,
    int *qisub_inv_dim_sizes, int qisub_inv_rank, float *tpert_out,
    int *tpert_out_dim_sizes, int tpert_out_rank, float *qpert_out,
    int *qpert_out_dim_sizes, int qpert_out_rank)
{
    int return_code = compute_uwshcu_py_run(
        dotransport, k0, windsrcavg, qtsrchgt, qtsrc_fac, thlsrc_fac, frc_rasn,
        rbuoy, epsvarw, use_CINcin, mumin1, rmaxfrac, PGFc, dt, niter_xc, criqc,
        rle, cridist_opt, mixscale, rdrag, rkm, use_self_detrain, detrhgt,
        use_cumpenent, rpen, use_momenflx, rdrop, iter_cin, pifc0_inv,
        pifc0_inv_dim_sizes, pifc0_inv_rank, zifc0_inv, zifc0_inv_dim_sizes,
        zifc0_inv_rank, pmid0_inv, pmid0_inv_dim_sizes, pmid0_inv_rank, zmid0_inv,
        zmid0_inv_dim_sizes, zmid0_inv_rank, kpbl_inv, kpbl_inv_dim_sizes,
        kpbl_inv_rank, exnmid0_inv, exnmid0_inv_dim_sizes, exnmid0_inv_rank,
        exnifc0_inv, exnifc0_inv_dim_sizes, exnifc0_inv_rank, dp0_inv,
        dp0_inv_dim_sizes, dp0_inv_rank, u0_inv, u0_inv_dim_sizes, u0_inv_rank,
        v0_inv, v0_inv_dim_sizes, v0_inv_rank, qv0_inv, qv0_inv_dim_sizes,
        qv0_inv_rank, ql0_inv, ql0_inv_dim_sizes, ql0_inv_rank, qi0_inv,
        qi0_inv_dim_sizes, qi0_inv_rank, t0_inv, t0_inv_dim_sizes, t0_inv_rank,
        frland_in, frland_in_dim_sizes, frland_in_rank, tke_inv,
        tke_inv_dim_sizes, tke_inv_rank, rkfre, rkfre_dim_sizes, rkfre_rank, cush,
        cush_dim_sizes, cush_rank, shfx, shfx_dim_sizes, shfx_rank, evap,
        evap_dim_sizes, evap_rank, cnvtr, cnvtr_dim_sizes, cnvtr_rank,
        CNV_Tracers, CNV_Tracers_dim_sizes, CNV_Tracers_rank, umf_inv,
        umf_inv_dim_sizes, umf_inv_rank, dcm_inv, dcm_inv_dim_sizes, dcm_inv_rank,
        qtflx_inv, qtflx_inv_dim_sizes, qtflx_inv_rank, slflx_inv,
        slflx_inv_dim_sizes, slflx_inv_rank, uflx_inv, uflx_inv_dim_sizes,
        uflx_inv_rank, vflx_inv, vflx_inv_dim_sizes, vflx_inv_rank, qvten_inv,
        qvten_inv_dim_sizes, qvten_inv_rank, qlten_inv, qlten_inv_dim_sizes,
        qlten_inv_rank, qiten_inv, qiten_inv_dim_sizes, qiten_inv_rank, tten_inv,
        tten_inv_dim_sizes, tten_inv_rank, uten_inv, uten_inv_dim_sizes,
        uten_inv_rank, vten_inv, vten_inv_dim_sizes, vten_inv_rank, qrten_inv,
        qrten_inv_dim_sizes, qrten_inv_rank, qsten_inv, qsten_inv_dim_sizes,
        qsten_inv_rank, cufrc_inv, cufrc_inv_dim_sizes, cufrc_inv_rank, fer_inv,
        fer_inv_dim_sizes, fer_inv_rank, fdr_inv, fdr_inv_dim_sizes, fdr_inv_rank,
        ndrop_inv, ndrop_inv_dim_sizes, ndrop_inv_rank, nice_inv,
        nice_inv_dim_sizes, nice_inv_rank, qldet_inv, qldet_inv_dim_sizes,
        qldet_inv_rank, qlsub_inv, qlsub_inv_dim_sizes, qlsub_inv_rank, qidet_inv,
        qidet_inv_dim_sizes, qidet_inv_rank, qisub_inv, qisub_inv_dim_sizes,
        qisub_inv_rank, tpert_out, tpert_out_dim_sizes, tpert_out_rank, qpert_out,
        qpert_out_dim_sizes, qpert_out_rank

    );

    if (return_code < 0)
    {
        exit(return_code);
    }
}
