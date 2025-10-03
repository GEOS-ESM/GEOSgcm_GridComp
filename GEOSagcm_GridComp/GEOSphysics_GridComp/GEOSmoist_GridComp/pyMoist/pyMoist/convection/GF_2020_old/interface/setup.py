from ndsl.dsl.gt4py import PARALLEL, interval, computation, FORWARD, sqrt, max, K, BACKWARD
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, Int
import pyMoist.constants as constants
from pyMoist.convective_parameterization.GF_2020.field_types import FloatField_nmp


def reset_to_zero(
    WQT_DC: FloatField,
    CNV_MFC: FloatField,
    PRFIL: FloatField,
    CNV_MF0: FloatField,
    CNV_PRC3: FloatField,
    CNV_MFD: FloatField,
    CNV_DQCDT: FloatField,
    CNV_UPDF: FloatField,
    CNV_CVW: FloatField,
    CNV_QC: FloatField,
    ENTLAM: FloatField,
    REVSU: FloatField,
    DQDT_GF: FloatField,
    DTDT_GF: FloatField,
    DUDT_GF: FloatField,
    DVDT_GF: FloatField,
    MUPDP: FloatField,
    MDNDP: FloatField,
    MUPSH: FloatField,
    MUPMD: FloatField,
    CNPCPRATE: FloatFieldIJ,
    LIGHTN_DENS: FloatFieldIJ,
    SIGMA_DEEP: FloatFieldIJ,
    SIGMA_MID: FloatFieldIJ,
    MFDP: FloatFieldIJ,
    MFSH: FloatFieldIJ,
    MFMD: FloatFieldIJ,
    ERRDP: FloatFieldIJ,
    ERRSH: FloatFieldIJ,
    ERRMD: FloatFieldIJ,
    AA0: FloatFieldIJ,
    AA1: FloatFieldIJ,
    AA2: FloatFieldIJ,
    AA3: FloatFieldIJ,
    AA1_BL: FloatFieldIJ,
    AA1_CIN: FloatFieldIJ,
    TAU_BL: FloatFieldIJ,
    TAU_EC: FloatFieldIJ,
):
    """
    Reset a bunch of fields to zero. They will be filled (in part or in full) later,
    either deeper in the driver or at the conclusion of the driver interface.

    This stencil MUST be built using Z_INTERFACE_DIM to function properly.
    """

    # 3D Z_INTERFACE_DIM fields
    with computation(PARALLEL), interval(...):
        WQT_DC = 0.0
        CNV_MFC = 0.0
        PRFIL = 0.0

    # 3D Z_DIM fields
    with computation(PARALLEL), interval(0, -1):
        CNV_MF0 = 0.0
        CNV_PRC3 = 0.0
        CNV_MFD = 0.0
        CNV_DQCDT = 0.0
        CNV_UPDF = 0.0
        CNV_CVW = 0.0
        CNV_QC = 0.0
        ENTLAM = 0.0
        REVSU = 0.0
        DQDT_GF = 0.0
        DTDT_GF = 0.0
        DUDT_GF = 0.0
        DVDT_GF = 0.0
        MUPDP = 0.0
        MDNDP = 0.0
        MUPSH = 0.0
        MUPMD = 0.0

    # 2D fields
    with computation(FORWARD), interval(0, 1):
        CNPCPRATE = 0.0
        LIGHTN_DENS = 0.0
        SIGMA_DEEP = 0.0
        SIGMA_MID = 0.0
        MFDP = 0.0
        MFSH = 0.0
        MFMD = 0.0
        ERRDP = 0.0
        ERRSH = 0.0
        ERRMD = 0.0
        AA0 = 0.0
        AA1 = 0.0
        AA2 = 0.0
        AA3 = 0.0
        AA1_BL = 0.0
        AA1_CIN = 0.0
        TAU_BL = 0.0
        TAU_EC = 0.0


def driver_interface(
    # inputs
    PLE: FloatField,
    PLO: FloatField,
    ZLE: FloatField,
    ZLO: FloatField,
    PK: FloatField,
    MASS: FloatField,
    KH: FloatField,
    T1: FloatField,
    Q1: FloatField,
    U1: FloatField,
    V1: FloatField,
    W1: FloatField,
    BYNCY: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    QLLS: FloatField,
    QILS: FloatField,
    CLCN: FloatField,
    CLLS: FloatField,
    QV_DYN_IN: FloatField,
    PLE_DYN_IN: FloatField,
    U_DYN_IN: FloatField,
    V_DYN_IN: FloatField,
    T_DYN_IN: FloatField,
    RADSW: FloatField,
    RADLW: FloatField,
    DQDT_BL: FloatField,
    DTDT_BL: FloatField,
    FRLAND: FloatFieldIJ,
    AREA: FloatFieldIJ,
    T2M: FloatFieldIJ,
    SH: FloatFieldIJ,
    EVAP: FloatFieldIJ,
    PHIS: FloatFieldIJ,
    KPBLIN: FloatFieldIJ,
    DTDTDYN: FloatField,
    DQVDTDYN: FloatField,
    TPWI: FloatFieldIJ,
    TPWI_star: FloatFieldIJ,
    CNV_TR: FloatField,
    # fields computed here and passed to the driver
    aot500: FloatFieldIJ,
    temp2m: FloatFieldIJ,
    sflux_r: FloatFieldIJ,
    sflux_t: FloatFieldIJ,
    topt: FloatFieldIJ,
    xland: FloatFieldIJ,
    dx2d: FloatFieldIJ,
    kpbl: FloatFieldIJ,
    temp: FloatField,
    press: FloatField,
    rvap: FloatField,
    up: FloatField,
    vp: FloatField,
    wp: FloatField,
    zt3d: FloatField,
    zm3d: FloatField,
    dm3d: FloatField,
    khloc: FloatField,
    curr_rvap: FloatField,
    mp_ice: FloatField_nmp,
    mp_liq: FloatField_nmp,
    mp_cf: FloatField_nmp,
    buoy_exc: FloatField,
    # fields computed here and used immediately after driver conclusion
    DZ: FloatField,
    AIR_DEN: FloatField,
    # fields computed here and used elsewhere in the model
    entr3d: FloatField,
    # workarounds
    maximum_t2m: Float,  # because max of entire field cannot be determined within a stencil
    i_end: Int,  # i_end cannot be used like k_end, gtscript error
    j_end: Int,  # j_end cannot be used like k_end, gtscript error
    ZLE_N: FloatField,
    ZLE_N_surface: FloatFieldIJ,
):
    from __externals__ import k_end, GF_ENV_SETTING, ENTRVERSION, CONVECTION_TRACER

    with computation(PARALLEL), interval(...):
        if i_end == 1 and j_end == 1 and maximum_t2m < 1e-6:
            single_column_stop = True
        else:
            single_column_stop = False

    with computation(PARALLEL), interval(...):
        if single_column_stop == False:
            # define the vector "flip" to invert the z-axis orientation
            flip = k_end - K

    # 2d input data
    with computation(FORWARD), interval(0, 1):
        if single_column_stop == False:
            aot500 = 0.1

            if maximum_t2m < 1.0e-6:
                temp2m = T1.at(K=k_end)  # Kelvin
            else:
                temp2m = T2M  # Kelvin

            # moisture flux from sfc
            sflux_r = EVAP  # kg m-2 s-1
            # sensible heat flux (sh) comes in W m-2, below it is converted to K m s-1
            sflux_t = SH / (
                1004.0 * PLE.at(K=k_end) / (287.04 * T1.at(K=k_end) * (1.0 + 0.608 * Q1.at(K=k_end)))
            )  # K m s-1
            # topography height  (m)
            topt = PHIS / constants.MAPL_GRAV
            # land/ocean fraction: land if < 1 ,ocean if = 1
            xland = 1.0 - FRLAND
            # grid length for the scale awareness
            if i_end == 1 and j_end == 1:
                # set area for single column runs
                dx2d = 100000.0
            else:
                dx2d = sqrt(AREA)  # meters

            # TODO not confident this is correc
            if KPBLIN != 0.0:
                kpbl = flip.at(K=int(round(KPBLIN)))
            else:
                kpbl = 1

    # 3-d input data
    # any var with index "1" (and w and pk) are already updated with dynamics
    # tendencies and everything else (from physics) that was called before moist
    with computation(PARALLEL), interval(...):
        if single_column_stop == False:
            if GF_ENV_SETTING == 0:
                # 1st setting: enviromental state is the one already modified by dyn + physics
                DZ = -(ZLE[0, 0, 1] - ZLE)
                AIR_DEN = PLO / (287.04 * T1 * (1.0 + 0.608 * Q1))
                temp = T1.at(K=flip)
                press = PLO.at(K=flip)  # Pa
                rvap = Q1.at(K=flip)
                up = U1.at(K=flip)  # already @ A-grid (m/s)
                vp = V1.at(K=flip)  # already @ A-grid (m/s)
                wp = W1.at(K=flip)  # m/s
                zt3d = ZLO.at(K=flip)  # mid -layer level
                zm3d = ZLE.at(K=flip)  # edge-layer level
                dm3d = MASS.at(K=flip)
                khloc = KH.at(K=flip)
                curr_rvap = Q1.at(K=flip)  # current rvap

                if ENTRVERSION == 0:
                    # eq 6 of https://doi.org/10.1029/2021JD034881
                    ec3d = 0.71 * max(0.5, W1.at(K=flip)) ** (-1.17) * max(0.1, BYNCY.at(K=flip)) ** (-0.36)
                else:
                    ec3d = 1.0

    with computation(PARALLEL), interval(...):
        if single_column_stop == False:
            if GF_ENV_SETTING == 0:
                # NOTE NEED A WAY TO WRITE TO A FLIPPED INDEXs
                entr3d = ec3d.at(K=flip)

                mp_ice[0, 0, 0][0] = QILS.at(K=flip)
                mp_liq[0, 0, 0][0] = QLLS.at(K=flip)
                mp_cf[0, 0, 0][0] = CLLS.at(K=flip)
                mp_ice[0, 0, 0][1] = QICN.at(K=flip)
                mp_liq[0, 0, 0][1] = QLCN.at(K=flip)
                mp_cf[0, 0, 0][1] = CLCN.at(K=flip)

                # sfc pressure (Pa)
                sfc_press = PLE.at(K=k_end)
                # Grid and sub-grid scale forcings for convection
                gsf_q = 0.0
                gsf_t = 0.0
                sgsf_t = 0.0
                sgsf_q = 0.0
                advf_t = 0.0
            elif GF_ENV_SETTING == 1:
                # 2nd setting: environmental state is that one before any tendency
                # is applied (i.e, at begin of each time step).
                # Get back the model state, heights and others variables at time N
                # (or at the beggining of current time step)
                # In physics, the state vars (T,U,V,PLE) are untouched and represent the
                # model state after dynamics phase 1. But, "Q" is modified by physics, so
                # depending on what was called before this subroutine, "Q" may be already
                # changed from what it was just after dynamics phase 1. To solve this issue,
                # "Q" just after dynamics is saved in the var named "QV_DYN_IN" in "GEOS_AgcmGridComp.F90".
                MASS_N = (PLE_DYN_IN[0, 0, 1] - PLE_DYN_IN) * (1.0 / constants.MAPL_GRAV)
                PLO_N = 0.5 * (PLE_DYN_IN + PLE_DYN_IN[0, 0, 1])
                PKE_N = (PLE_DYN_IN / constants.MAPL_P00) ** (constants.MAPL_RGAS / constants.MAPL_CP)
                if K == k_end:
                    PKE_N_surface = (PLE_DYN_IN[0, 0, 1] / constants.MAPL_P00) ** (
                        constants.MAPL_RGAS / constants.MAPL_CP
                    )
                PK_N = (PLO_N / constants.MAPL_P00) ** (constants.MAPL_RGAS / constants.MAPL_CP)
                ZLE_N = (T_DYN_IN / PK_N) * (1.0 + constants.MAPL_VIREPS * QV_DYN_IN)

    with computation(FORWARD), interval(0, 1):
        if single_column_stop == False:
            if GF_ENV_SETTING == 1:
                ZLE_N_surface = 0.0

    with computation(BACKWARD), interval(...):
        if single_column_stop == False:
            if GF_ENV_SETTING == 1:
                if K == k_end:
                    ZLO_N = (
                        ZLE_N_surface
                        + (constants.MAPL_CP / constants.MAPL_GRAV) * (PKE_N_surface - PK_N) * ZLE_N
                    )
                else:
                    ZLO_N = (
                        ZLE_N[0, 0, 1]
                        + (constants.MAPL_CP / constants.MAPL_GRAV) * (PKE_N[0, 0, 1] - PK_N) * ZLE_N
                    )
                ZLE_N = ZLO_N + (constants.MAPL_CP / constants.MAPL_GRAV) * (PK_N - PKE_N) * ZLE_N

    with computation(PARALLEL), interval(...):
        if single_column_stop == False:
            if GF_ENV_SETTING == 1:
                if K == k_end:
                    DZ = -(ZLE_N_surface - ZLE_N)
                else:
                    DZ = -(ZLE_N[0, 0, 1] - ZLE_N)
                AIR_DEN = PLO_N / (287.04 * T_DYN_IN * (1.0 + 0.608 * QV_DYN_IN))
                temp = T_DYN_IN.at(K=flip)  # (K)
                press = PLO_N.at(K=flip)  # (Pa) @ mid-layer level
                rvap = QV_DYN_IN.at(K=flip)  # water vapor mix ratio
                up = U_DYN_IN.at(K=flip)  # already @ A-grid (m/s)
                vp = V_DYN_IN.at(K=flip)  # already @ A-grid (m/s)
                wp = W1.at(K=flip)  # (m/s)
                zt3d = ZLO_N.at(K=flip)  # mid -layer level (m)
                zm3d = ZLE_N.at(K=flip + 1)  # edge-layer level (m)
                dm3d = MASS_N.at(K=flip)
                khloc = KH.at(K=flip)
                curr_rvap = Q1.at(K=flip)  # current rvap (dyn+phys)

                if ENTRVERSION == 0:
                    # eq 6 of https://doi.org/10.1029/2021JD034881
                    ec3d = 0.71 * max(0.5, W1.at(K=flip)) ** (-1.17) * max(0.1, BYNCY.at(K=flip)) ** (-0.36)
                else:
                    ec3d = 1.0

                mp_ice[0, 0, 0][0] = QILS.at(K=flip)
                mp_liq[0, 0, 0][0] = QLLS.at(K=flip)
                mp_cf[0, 0, 0][0] = CLLS.at(K=flip)
                mp_ice[0, 0, 0][1] = QICN.at(K=flip)
                mp_liq[0, 0, 0][1] = QLCN.at(K=flip)
                mp_cf[0, 0, 0][1] = CLCN.at(K=flip)

                # sfc pressure (Pa)
                sfc_press = PLE_DYN_IN.at(K=k_end)
                # Grid and sub-grid scale forcings for convection
                gsf_t = DTDTDYN.at(K=flip) + RADSW.at(K=flip) + RADLW.at(K=flip)
                gsf_q = DQVDTDYN.at(K=flip)
                sgsf_t = DTDT_BL.at(K=flip)
                sgsf_q = DQDT_BL.at(K=flip)
                advf_t = DTDTDYN.at(K=flip)

    with computation(PARALLEL), interval(...):
        if single_column_stop == False:
            if GF_ENV_SETTING == 1:
                entr3d = ec3d.at(K=flip)
            if CONVECTION_TRACER == 1:
                buoy_exc = CNV_TR.at(K=flip)
            else:
                buoy_exc = 0.0

    # NOTE need to deal with tracer stuff later
    # mtp = size(CNV_Tracers)
    # IF(.not.allocated(SRC_CHEM)) THEN
    #    allocate(SRC_CHEM(mtp, mzp, mxp, myp),stat=alloc_stat) !- tendency from convection
    #    IF(alloc_stat==0) SRC_CHEM=0.0
    # ENDIF
