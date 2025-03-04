#include <stdio.h>
#include <time.h>
#include "interface.h"

extern int pyMLINC_interface_init_c(int magic_number)
{
    // Check magic number
    if (magic_number != 123456789) {
        printf("[pyMLINC_interface_init_c] Magic number failed\n");
        exit(-1);
    }
    int rc = pyMLINC_interface_init_py(magic_number);
    if (rc != 0)
        exit(rc);
    return 0;
}

extern int pyMLINC_interface_run_c(
    // input
    int xdim, int ydim, int zdim,
    const float *u,
    const float *v,
    const float *t,
    const float *qv,
    const float *ql,
    const float *qi,
    const float *qr,
    const float *qs,
    const float *qg,
    // output
    float *dtdt,
    // LAST ARGUMENT
    int magic_number)
{
    // Check magic number
    if (magic_number != 123456789) {
        printf("[pyMLINC_interface_run_c] Magic number failed\n");
        exit(-1);
    }
    int rc = pyMLINC_interface_run_py(
        xdim, ydim, zdim,
        u, v, t,
        qv, ql, qi, qr, qs, qg,
        dtdt,
        magic_number);
    if (rc != 0)
        exit(rc);
    return 0;
}
