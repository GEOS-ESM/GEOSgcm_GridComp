#include <stdio.h>
#include <time.h>
#include "interface.h"

extern int pyMLINC_interface_init_c() {
    int rc = pyMLINC_interface_init_py();
    if (rc != 0) {
        exit(rc);
    }
    return 0;
}

extern int pyMLINC_interface_run_c(a_pod_struct_t *options, const float *in_buffer, float *out_buffer) {
    // Check magic number
    if (options->mn_123456789 != 123456789) {
        printf("Magic number failed, pyMLINC interface is broken on the C side\n");
        exit(-1);
    }
    int rc = pyMLINC_interface_run_py(options, in_buffer, out_buffer);
    if (rc != 0) {
        exit(rc);
    }
    return 0;
}
