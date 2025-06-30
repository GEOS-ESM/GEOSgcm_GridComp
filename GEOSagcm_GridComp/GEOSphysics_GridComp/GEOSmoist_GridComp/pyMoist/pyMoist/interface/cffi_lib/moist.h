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
	// GFDL_1M
	float DT_MOIST;
	float MP_TIME;
	float T_MIN;
	float T_SUB;
	float TAU_R2G;
	float TAU_SMLT;
	float TAU_G2R;
	float DW_LAND;
	float DW_OCEAN;
	float VI_FAC;
	float VR_FAC;
	float VS_FAC;
	float VG_FAC;
	float QL_MLT;
	bool DO_QA;
	bool FIX_NEGATIVE;
	float VI_MAX;
	float VS_MAX;
	float VG_MAX;
	float VR_MAX;
	float QS_MLT;
	float QS0_CRT;
	float QI_GEN;
	float QL0_MAX;
	float QI0_MAX;
	float QI0_CRT;
	float QR0_CRT;
	bool FAST_SAT_ADJ;
	float RH_INC;
	float RH_INS;
	float RH_INR;
	bool CONST_VI;
	bool CONST_VS;
	bool CONST_VG;
	bool CONST_VR;
	bool USE_CCN;
	float RTHRESHU;
	float RTHRESHS;
	float CCN_L;
	float CCN_O;
	float QC_CRT;
	float TAU_G2V;
	float TAU_V2G;
	float TAU_S2V;
	float TAU_V2S;
	float TAU_REVP;
	float TAU_FRZ;
	bool DO_BIGG;
	bool DO_EVAP;
	bool DO_SUBL;
	float SAT_ADJ0;
	float C_PIACR;
	float TAU_IMLT;
	float TAU_V2L;
	float TAU_L2V;
	float TAU_I2V;
	float TAU_I2S;
	float TAU_L2R;
	float QI_LIM;
	float QL_GEN;
	float C_PAUT;
	float C_PSACI;
	float C_PGACS;
	float C_PGACI;
	bool Z_SLOPE_LIQ;
	bool Z_SLOPE_ICE;
	bool PROG_CCN;
	float C_CRACW;
	float ALIN;
	float CLIN;
	bool PRECIPRAD;
	float CLD_MIN;
	bool USE_PPM;
	bool MONO_PROF;
	bool DO_SEDI_HEAT;
	bool SEDI_TRANSPORT;
	bool DO_SEDI_W;
	bool DE_ICE;
	int ICLOUD_F;
	int IRAIN_F;
	bool MP_PRINT;
	bool USE_BERGERON;
	// Magic number needs to be last item
	int mn_123456789;
} gfdl_1m_flags_t;

// Carry wrapper for MPI
typedef union
{
	int comm_int;
	void *comm_ptr;
} MPI_Comm_t;

extern int pymoist_interface_py_init(void *import, void *export, void *internal, void *mapl_comp, moist_flags_t *flags);

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

extern int pymoist_interface_py_run_GFDL_1M();

extern int pymoist_interface_py_finalize();
