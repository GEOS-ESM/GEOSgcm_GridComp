#include <stdio.h>
#include <time.h>
#include "interface.h"

extern int pyMKIAU_interface_c_setservice()
{
    // Check magic number
    int return_code = pyMKIAU_interface_py_setservices();

    if (return_code < 0)
    {
        exit(return_code);
    }
}

extern int pyMKIAU_interface_c_run(a_pod_struct_t *options, const float *in_buffer, float *out_buffer)
{
    // Check magic number
    if (options->mn_123456789 != 123456789)
    {
        printf("Magic number failed, pyMKIAU interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = pyMKIAU_interface_py_run(options, in_buffer, out_buffer);

    if (return_code < 0)
    {
        exit(return_code);
    }
}
