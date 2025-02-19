#pragma once

/***
 * C Header for the interface to python.
 * Define here any POD-strict structures and external functions
 * that will get exported by cffi from python (see interface.py)
 ***/

#include <stdbool.h>
#include <stdlib.h>

// POD-strict structure to pack options and flags efficiently
// Struct CANNOT hold pointers. The iso_c_binding does not allow for foolproof
// pointer memory packing.
// We use the low-embedded trick of the magic number to attempt to catch
// any type mismatch betweeen Fortran and C. This is not a foolproof method
// but it bring a modicum of check at the cost of a single integer.
typedef struct {
	int npx;
	int npy;
	int npz;
	// Magic number needs to be last item
	int mn_123456789;
} a_pod_struct_t;

// For complex type that can be exported with different
// types (like the MPI communication object), you can rely on C `union`
typedef union {
	int comm_int;
	void *comm_ptr;
} MPI_Comm_t;

// Python hook functions: defined as external so that the .so can link out ot them
// Though we define `in_buffer` as a `const float*` it is _not_ enforced
// by the interface. Treat as a developer hint only.

extern int pyMLINC_interface_run_py(a_pod_struct_t *options, const float *in_buffer, float *out_buffer);
extern int pyMLINC_interface_setservices_py();
