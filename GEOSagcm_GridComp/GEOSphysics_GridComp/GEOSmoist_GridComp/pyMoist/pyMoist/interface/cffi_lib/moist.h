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

extern int pymoist_interface_py_finalize();
