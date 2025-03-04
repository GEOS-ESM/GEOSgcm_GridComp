#pragma once

/***
 * C Header for the interface to python.
 * Define here any POD-strict structures and external functions
 * that will get exported by cffi from python (see interface.py)
 ***/

#include <stdbool.h>
#include <stdlib.h>

// For complex type that can be exported with different
// types (like the MPI communication object), one can rely on C `union`
typedef union {
	int comm_int;
	void *comm_ptr;
} MPI_Comm_t;

// Python hook functions: defined as external so that the .so can link out ot them
// Though we define `in_buffer` as a `const float*` it is _not_ enforced
// by the interface. Treat as a developer hint only.

extern int pyMLINC_interface_init_py(int magic_number);
extern int pyMLINC_interface_run_py(
    int xdim,
    int ydim,
    int zdim,
    const float *u,
    const float *v,
    const float *t,
    const float *qv,
    const float *ql,
    const float *qi,
    const float *qr,
    const float *qs,
    const float *qg,
    const float *ps,
    float *dtdt,
    int magic_number);
