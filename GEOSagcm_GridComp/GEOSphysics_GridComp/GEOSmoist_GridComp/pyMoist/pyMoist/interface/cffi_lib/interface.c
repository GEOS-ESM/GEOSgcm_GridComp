#include <stdio.h>
#include <time.h>
#include "moist.h"

extern int pymoist_interface_c_init(moist_flags_t *flags)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, pymoist interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = pymoist_interface_py_init(flags);

    if (return_code < 0)
    {
        exit(return_code);
    }
}

void pymoist_interface_c_run_GFDL1M(
    float dw_land, float dw_ocean, int PDFSHAPE, float TURNRHCRIT_PARAM,
    float DT_MOIST, float CCW_EVAP_EFF, float CCI_EVAP_EFF,
    int LMELTFRZ,
    float *AREA, float *CNV_FRC, float *SRF_TYPE, // 2D...
    int *KLCL,
    float *EIS, float *PLmb, float *PLEmb, float *NACTL, float *NACTI, float *QST, // 3D...
    float *T, float *Q, float *QLCN, float *QICN, float *QLLS, float *QILS, float *CLLS, float *CLCN,
    float *SUBLC, float *EVAPC, float *RHX)
{
    printf("> C LMELTFRZ = %d\n", LMELTFRZ);
    int return_code = pymoist_interface_py_run_GFDL1M(
        dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM,
        DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF,
        LMELTFRZ != 0,
        AREA, CNV_FRC, SRF_TYPE,
        KLCL,
        EIS, PLmb, PLEmb, NACTL, NACTI, QST,
        T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN,
        SUBLC, EVAPC, RHX);
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

void pymoist_interface_c_finalize()
{
    int return_code = pymoist_interface_py_finalize();
    if (return_code < 0)
    {
        exit(return_code);
    }
}
