#include <stdbool.h>
#include <stdio.h>

extern void c_held_suarez_oacc(
    float *, float *, float *, float *, float *, // CPHI2, DISS, DTDT, DUDT, DVDT
    float *, float *, float *, float *, float *, float *, float *, // HFCN, P_I, PLE, SPHI2, TAUX, TAUY, T,
    float *, float *, float *, float *, // THEQ, T_EQ, U, V,
    float, float, float, float, int, int, // DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp,
    int, float, float, int, int, int, float, float, float, //FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX
    float, float, float, float, float, float, int, int); // SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank

void c_call_hs_oacc(
    float * CPHI2, float * DISS, float * DTDT, float * DUDT, float * DVDT,
    float * HFCN,  float * P_I, float * PLE, float * SPHI2, float * TAUX, float * TAUY, float * T,
    float * THEQ, float * T_EQ, float * U, float * V,
    float DAYLEN, float DELH, float DELV1, float DT, int FRICQ, int FriendlyTemp,
    int FriendlyWind, float GAM_D, float GAM_I, int IM, int JM, int LM, float P_1, float P_D, float QMAX,
    float SIG1, float TAUA, float TAUF, float TAUS, float TSTRT, float T0, int compType, int rank) {

    //printf("calling c_held_suarez_oacc\n");
    c_held_suarez_oacc(
        CPHI2, DISS, DTDT, DUDT, DVDT,
        HFCN, P_I, PLE, SPHI2, TAUX, TAUY, T,
        THEQ, T_EQ, U, V,
        DAYLEN, DELH, DELV1, DT, FRICQ, FriendlyTemp, 
        FriendlyWind, GAM_D, GAM_I, IM, JM, LM, P_1, P_D, QMAX, 
        SIG1, TAUA, TAUF, TAUS, TSTRT, T0, compType, rank);
}
