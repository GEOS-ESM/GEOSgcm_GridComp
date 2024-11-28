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

// Carry wrapper for MPI
typedef union
{
	int comm_int;
	void *comm_ptr;
} MPI_Comm_t;

extern int pymoist_interface_py_init(moist_flags_t *flags);

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

extern int pymoist_interface_py_finalize();
