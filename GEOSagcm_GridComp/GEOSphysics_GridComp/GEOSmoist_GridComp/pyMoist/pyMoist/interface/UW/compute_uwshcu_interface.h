#pragma once

#include <stdbool.h>
#include <stdlib.h>

typedef union {
  int comm_int;
  void *comm_ptr;
} MPI_Comm_t;

extern void compute_uwshcu_py_init(
    // inputs
    int ncnst, int k0, int windsrcavg
    // inputs-outputs

    // outputs

);

extern void compute_uwshcu_py_run_compute_uwshcu(
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
    const int *kpbl_inv, int *kpbl_inv_dim_sizes, int kpbl_inv_rank,
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
    int *qpert_out_dim_sizes, int qpert_out_rank);
