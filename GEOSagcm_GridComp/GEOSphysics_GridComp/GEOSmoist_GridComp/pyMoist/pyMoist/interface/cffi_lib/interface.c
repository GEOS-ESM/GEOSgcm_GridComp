#include <stdio.h>
#include <time.h>
#include "moist.h"

extern int pymoist_interface_c_init(void *import, void *export, void *internal, void *mapl_comp, moist_flags_t *flags)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, pymoist interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = pymoist_interface_py_init(import, export, internal, mapl_comp, flags);

    if (return_code < 0)
    {
        exit(return_code);
    }
}

extern int gfdl_1m_interface_c_init(gfdl_1m_flags_t *flags)
{
    // Check magic number
    if (flags->mn_123456789 != 123456789)
    {
        printf("Magic number failed, gfdl_1m interface is broken on the C side\n");
        exit(-1);
    }

    int return_code = gfdl_1m_interface_py_init(flags);

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

void pymoist_interface_c_run_GFDL_1M_driver(
    float *RAD_QV, float *RAD_QL, float *RAD_QR, float *RAD_QI, float *RAD_QS, float *RAD_QG, float *RAD_CF, float *NACTAll,
    float *DQVDTmic, float *DQLDTmic, float *DQRDTmic, float *DQIDTmic,
    float *DQSDTmic, float *DQGDTmic, float *DQADTmic, float *DTDTmic,
    float *T, float *W, float *U, float *V, float *DUDTmic, float *DVDTmic, float *DZ, float *DP,
    float *AREA, float *FRLAND, float *CNV_FRC, float *SRF_TYPE, float *EIS, float *RHCRIT3D,
    float DT_MOIST, float ANV_ICEFALL, float LS_ICEFALL,
    float *REV_LS, float *RSU_LS,
    float *PRCP_RAIN, float *PRCP_SNOW, float *PRCP_ICE, float *PRCP_GRAUPEL,
    float *PFL_LS, float *PFI_LS,
    int LHYDROSTATIC, int LPHYS_HYDROSTATIC)
{
    int return_code = pymoist_interface_py_run_GFDL_1M_driver(
        RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTAll,
        DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, DQSDTmic, DQGDTmic, DQADTmic, DTDTmic,
        T, W, U, V, DUDTmic, DVDTmic, DZ, DP, AREA, FRLAND, CNV_FRC, SRF_TYPE, EIS, RHCRIT3D,
        DT_MOIST, ANV_ICEFALL, LS_ICEFALL,
        REV_LS, RSU_LS, PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, PFL_LS, PFI_LS,
        LHYDROSTATIC, LPHYS_HYDROSTATIC);
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
