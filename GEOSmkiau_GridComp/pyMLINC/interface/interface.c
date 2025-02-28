#include <stdio.h>
#include <time.h>
#include "interface.h"

extern int pyMLINC_interface_init_c(int magic_number)
{
    // Check magic number
    if (magic_number != 123456789) {
        printf("Magic number failed on the C side\n");
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
    const float *qv,
    // output
    float *dtdt)
{
    int rc = pyMLINC_interface_run_py(xdim, ydim, zdim, qv, dtdt);
    if (rc != 0)
        exit(rc);
    return 0;
}
