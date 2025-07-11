#include <stdio.h>
#include <time.h>
#include "moist.h"

void pymoist_interface_c_init(void *import, void *export, void *mapl, moist_flags_t *flags)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, pymoist interface is broken on the C side\n");
        exit(-1);
    }
    int return_code = pymoist_interface_py_init(import, export, mapl, flags);

    if (return_code < 0)
    {
        exit(return_code);
    }
}

void gfdl_1m_interface_c_init(gfdl_1m_flags_t *flags, void *internal)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, gfdl_1m interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = gfdl_1m_interface_py_init(flags, internal);

    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_run_AerActivation(
    // input
    float *aero_dgn, float *aero_num, float *aero_hygroscopicity, float *aero_sigma,
    float *frland, float nn_ocean, float nn_land,
    float *t, float *plo,
    float *qicn, float *qils, float *qlcn, float *qlls,
    float *vvel, float *tke,
    // output
    float *nacti, float *nwfa, float *nactl)
{
    int return_code = pymoist_interface_py_run_AerActivation(
        aero_dgn,
        aero_num,
        aero_hygroscopicity,
        aero_sigma,
        frland,
        nn_ocean,
        nn_land,
        t,
        plo,
        qicn,
        qils,
        qlcn,
        qlls,
        vvel,
        tke,
        nacti,
        nwfa,
        nactl);
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_run_GFDL_1M()
{
    int return_code = pymoist_interface_py_run_GFDL_1M();
    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_finalize()
{
    int return_code = pymoist_interface_py_finalize();
    if (return_code < 0)
    {
        exit(return_code);
    }
}
