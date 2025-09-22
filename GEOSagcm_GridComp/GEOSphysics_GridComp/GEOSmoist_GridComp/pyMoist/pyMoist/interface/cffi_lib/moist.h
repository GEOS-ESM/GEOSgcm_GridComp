#pragma once

/***
 * Moist physic parametrization configuration from GEOS (namelist and extra flags)
 ***/

#include <stdbool.h>
#include <stdlib.h>

// Fortran FlagStruct
typedef struct
{
	// Grid information
	int npx;
	int npy;
	int npz;
	int layout_x;
	int layout_y;
	int n_tiles;
	// Aer Activation
	int n_modes;
	// Magic number needs to be last item
	int mn_123456789;
} moist_flags_t;

typedef struct
{
	// GFDL_1M driver
	bool phys_hydrostatic;
	bool hydrostatic;
	bool do_qa;
	bool fix_negative;
	bool fast_sat_adj;
	bool const_vi;
	bool const_vs;
	bool const_vg;
	bool const_vr;
	bool use_ccn;
	bool do_bigg;
	bool do_evap;
	bool do_subl;
	bool z_slope_liq;
	bool z_slope_ice;
	bool prog_ccn;
	bool preciprad;
	bool use_ppm;
	bool mono_prof;
	bool do_sedi_heat;
	bool sedi_transport;
	bool do_sedi_w;
	bool de_ice;
	bool mp_print;
	float dt_moist;
	float mp_time;
	float t_min;
	float t_sub;
	float tau_r2g;
	float tau_smlt;
	float tau_g2r;
	float dw_land;
	float dw_ocean;
	float vi_fac;
	float vr_fac;
	float vs_fac;
	float vg_fac;
	float ql_mlt;
	float vi_max;
	float vs_max;
	float vg_max;
	float vr_max;
	float qs_mlt;
	float qs0_crt;
	float qi_gen;
	float ql0_max;
	float qi0_max;
	float qi0_crt;
	float qr0_crt;
	float rh_inc;
	float rh_ins;
	float rh_inr;
	float rthreshu;
	float rthreshs;
	float ccn_l;
	float ccn_o;
	float qc_crt;
	float tau_g2v;
	float tau_v2g;
	float tau_s2v;
	float tau_v2s;
	float tau_revp;
	float tau_frz;
	float sat_adj0;
	float c_piacr;
	float tau_imlt;
	float tau_v2l;
	float tau_l2v;
	float tau_i2v;
	float tau_i2s;
	float tau_l2r;
	float qi_lim;
	float ql_gen;
	float c_paut;
	float c_psaci;
	float c_pgacs;
	float c_pgaci;
	float c_cracw;
	float alin;
	float clin;
	float cld_min;
	int icloud_f;
	int irain_f;
	// Magic number needs to be last item
	int mn_123456789;
} gfdl_1m_flags_t;

// Carry wrapper for MPI
typedef union
{
	int comm_int;
	void *comm_ptr;
} MPI_Comm_t;

extern int pymoist_interface_py_init(moist_flags_t *flags);

extern int gfdl_1m_interface_py_init(gfdl_1m_flags_t *flags);

extern int pymoist_interface_py_run_AerActivation(
	// input
	float *aero_dgn, float *aero_num, float *aero_hygroscopicity, float *aero_sigma,
	float *frland, float nn_ocean, float nn_land,
	float *t, float *plo,
	float *qicn, float *qils, float *qlcn, float *qlls,
	float *vvel, float *tke,
	// output
	float *nacti, float *nwfa, float *nactl);

extern int pymoist_interface_py_run_GFDL1M(
	float dw_land, float dw_ocean, int PDFSHAPE, float TURNRHCRIT_PARAM,
	float DT_MOIST, float CCW_EVAP_EFF, float CCI_EVAP_EFF,
	int LMELTFRZ,
	float *AREA, float *CNV_FRC, float *SRF_TYPE, // 2D..
	int *KLCL,
	float *EIS, float *PLmb, float *PLEmb, float *NACTL, float *NACTI, float *QST, // 3D..
	float *T, float *Q, float *QLCN, float *QICN, float *QLLS, float *QILS, float *CLLS, float *CLCN,
	float *SUBLC, float *EVAPC, float *RHX);

extern int pymoist_interface_py_run_GFDL_1M_driver(
	float *RAD_QV, float *RAD_QL, float *RAD_QR, float *RAD_QI, float *RAD_QS, float *RAD_QG, float *RAD_CF, float *NACTAll,
	float *DQVDTmic, float *DQLDTmic, float *DQRDTmic, float *DQIDTmic,
	float *DQSDTmic, float *DQGDTmic, float *DQADTmic, float *DTDTmic,
	float *T, float *W, float *U, float *V, float *DUDTmic, float *DVDTmic, float *DZ, float *DP,
	float *AREA, float *FRLAND, float *CNV_FRC, float *SRF_TYPE, float *EIS, float *RHCRIT3D,
	float DT_MOIST, float ANV_ICEFALL, float LS_ICEFALL,
	float *REV_LS, float *RSU_LS,
	float *PRCP_RAIN, float *PRCP_SNOW, float *PRCP_ICE, float *PRCP_GRAUPEL,
	float *PFL_LS, float *PFI_LS,
	bool LHYDROSTATIC, bool LPHYS_HYDROSTATIC);

extern int pymoist_interface_py_finalize();

extern int compute_uwshcu_py_init(
	// inputs
	int ncnst, int k0, int windsrcavg
	// inputs-outputs

	// outputs

);

extern int compute_uwshcu_py_run(
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
	int *qpert_out_dim_sizes, int qpert_out_rank);
