from dataclasses import dataclass

import numpy as np
from ndsl.dsl.gt4py import int32

from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.constants import constants


@dataclass
class ConfigConstants:
    C_AIR: Float
    C_VAP: Float
    P_NONHYDRO: bool
    D0_VAP: Float
    LV00: Float
    DO_SEDI_W: bool
    LATV: Float
    LATI: Float
    LATS: Float
    LAT2: Float
    LCP: Float
    ICP: Float
    TCP: Float
    MPDT: Float
    RDT: Float
    NTIMES: Float
    DTS: Float
    CPAUT: Float
    RDTS: Float
    FAC_IMLT: Float
    FAC_I2S: Float
    FAC_V2L: Float
    FAC_L2V: Float
    FAC_I2V: Float
    FAC_S2V: Float
    FAC_V2S: Float
    FAC_G2V: Float
    FAC_V2G: Float
    FAC_FRZ: Float
    CGACS: Float
    CSACW: Float
    CRACI: Float
    CSACI: Float
    CGACW: Float
    CGACI: Float
    CRACW: Float
    CGFR_0: Float
    CGFR_1: Float
    CSSUB_0: Float
    CSSUB_1: Float
    CSSUB_2: Float
    CSSUB_3: Float
    CSSUB_4: Float
    CGSUB_0: Float
    CGSUB_1: Float
    CGSUB_2: Float
    CGSUB_3: Float
    CGSUB_4: Float
    CREVP_0: Float
    CREVP_1: Float
    CREVP_2: Float
    CREVP_3: Float
    CREVP_4: Float
    CSMLT_0: Float
    CSMLT_1: Float
    CSMLT_2: Float
    CSMLT_3: Float
    CSMLT_4: Float
    CGMLT_0: Float
    CGMLT_1: Float
    CGMLT_2: Float
    CGMLT_3: Float
    CGMLT_4: Float

    @classmethod
    def make(cls, GFDL_1M_config: GFDL1MConfig):
        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vap based on hydrostatical property
        # -----------------------------------------------------------------------

        if GFDL_1M_config.PHYS_HYDROSTATIC or GFDL_1M_config.HYDROSTATIC:
            C_AIR = constants.CP_AIR
            C_VAP = constants.CP_VAP
            P_NONHYDRO = False
        else:
            C_AIR = constants.CV_AIR
            C_VAP = constants.CV_VAP
            P_NONHYDRO = True
        D0_VAP = C_VAP - constants.C_LIQ
        LV00 = constants.HLV0 - D0_VAP * constants.T_ICE

        if GFDL_1M_config.HYDROSTATIC:
            DO_SEDI_W = False
        else:
            DO_SEDI_W = True

        # -----------------------------------------------------------------------
        # define latent heat coefficient used in wet bulb and bigg mechanism
        # -----------------------------------------------------------------------

        LATV = constants.HLV
        LATI = constants.HLF
        LATS = LATV + LATI
        LAT2 = LATS * LATS

        LCP = LATV / constants.CP_AIR
        ICP = LATI / constants.CP_AIR
        TCP = (LATV + LATI) / constants.CP_AIR

        # -----------------------------------------------------------------------
        # define cloud microphysics sub time step
        # -----------------------------------------------------------------------

        MPDT = min(GFDL_1M_config.DT_MOIST, GFDL_1M_config.MP_TIME)
        RDT = Float(1.0) / GFDL_1M_config.DT_MOIST
        NTIMES = int32(np.round(GFDL_1M_config.DT_MOIST / MPDT))

        # small time step:
        DTS = GFDL_1M_config.DT_MOIST / Float(NTIMES)

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        CPAUT = GFDL_1M_config.C_PAUT * Float(0.104) * constants.GRAV / Float(1.717e-5)

        # -----------------------------------------------------------------------
        # define conversion scalar / factor for icloud
        # -----------------------------------------------------------------------
        RDTS = Float(1.0) / DTS
        FAC_IMLT = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_IMLT, dtype=Float)
        FAC_I2S = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_I2S, dtype=Float)
        FAC_V2L = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_V2L, dtype=Float)
        FAC_L2V = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_L2V, dtype=Float)
        FAC_I2V = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_I2V, dtype=Float)
        FAC_S2V = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_S2V, dtype=Float)
        FAC_V2S = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_V2S, dtype=Float)
        FAC_G2V = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_G2V, dtype=Float)
        FAC_V2G = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_V2G, dtype=Float)
        FAC_FRZ = Float(1.0) - np.exp(-DTS / GFDL_1M_config.TAU_FRZ, dtype=Float)

        # -----------------------------------------------------------------------
        # constatns from setupm
        # -----------------------------------------------------------------------

        CGACS = constants.PISQ * constants.RNZG * constants.RNZS * constants.RHOS
        CGACS = CGACS * GFDL_1M_config.C_PGACS

        CSACW = (
            constants.PIE
            * constants.RNZS
            * GFDL_1M_config.CLIN
            * constants.GAM325
            / (Float(4.0) * np.power(constants.ACT[0], 0.8125, dtype=Float))
        )
        # decreasing csacw to reduce cloud water --- > snow

        CRACI = (
            constants.PIE
            * constants.RNZR
            * GFDL_1M_config.ALIN
            * constants.GAM380
            / (Float(4.0) * np.power(constants.ACT[1], 0.95, dtype=Float))
        )
        CSACI = CSACW * GFDL_1M_config.C_PSACI

        CGACW = (
            constants.PIE
            * constants.RNZG
            * constants.GAM350
            * constants.GCON
            / (Float(4.0) * np.power(constants.ACT[5], 0.875, dtype=Float))
        )

        CGACI = CGACW * GFDL_1M_config.C_PGACI

        CRACW = CRACI  # cracw = 3.27206196043822
        CRACW = GFDL_1M_config.C_CRACW * CRACW

        CSSUB = np.zeros(5)
        CGSUB = np.zeros(5)
        CREVP = np.zeros(5)

        CSSUB[0] = (
            Float(2.0) * constants.PIE * constants.VDIFU * constants.TCOND * constants.RVGAS * constants.RNZS
        )
        CGSUB[0] = (
            Float(2.0) * constants.PIE * constants.VDIFU * constants.TCOND * constants.RVGAS * constants.RNZG
        )
        CREVP[0] = (
            Float(2.0) * constants.PIE * constants.VDIFU * constants.TCOND * constants.RVGAS * constants.RNZR
        )
        CSSUB[1] = Float(0.78) / np.sqrt(constants.ACT[0], dtype=Float)
        CGSUB[1] = Float(0.78) / np.sqrt(constants.ACT[5], dtype=Float)
        CREVP[1] = Float(0.78) / np.sqrt(constants.ACT[1], dtype=Float)
        CSSUB[2] = (
            Float(0.31)
            * constants.SCM3
            * constants.GAM263
            * np.sqrt(GFDL_1M_config.CLIN / constants.VISK, dtype=Float)
            / np.power(constants.ACT[0], Float(0.65625), dtype=Float)
        )
        CGSUB[2] = (
            Float(0.31)
            * constants.SCM3
            * constants.GAM275
            * np.sqrt(constants.GCON / constants.VISK, dtype=Float)
            / np.power(constants.ACT[5], Float(0.6875), dtype=Float)
        )
        CREVP[2] = (
            Float(0.31)
            * constants.SCM3
            * constants.GAM209
            * np.sqrt(GFDL_1M_config.ALIN / constants.VISK, dtype=Float)
            / np.power(constants.ACT[1], Float(0.725), dtype=Float)
        )
        CSSUB[3] = constants.TCOND * constants.RVGAS
        CSSUB[4] = np.power(constants.HLTS, 2, dtype=Float) * constants.VDIFU
        CGSUB[3] = CSSUB[3]
        CREVP[3] = CSSUB[3]
        CGSUB[4] = CSSUB[4]
        CREVP[4] = np.power(constants.HLTC, 2, dtype=Float) * constants.VDIFU

        CGFR_0 = (
            Float(20.0e2)
            * constants.PISQ
            * constants.RNZR
            * constants.RHOR
            / np.power(constants.ACT[1], Float(1.75), dtype=Float)
        )
        CGFR_1 = Float(0.66)

        CSSUB_0 = CSSUB[0]
        CSSUB_1 = CSSUB[1]
        CSSUB_2 = CSSUB[2]
        CSSUB_3 = CSSUB[3]
        CSSUB_4 = CSSUB[4]

        CGSUB_0 = CGSUB[0]
        CGSUB_1 = CGSUB[1]
        CGSUB_2 = CGSUB[2]
        CGSUB_3 = CGSUB[3]
        CGSUB_4 = CGSUB[4]

        CREVP_0 = CREVP[0]
        CREVP_1 = CREVP[1]
        CREVP_2 = CREVP[2]
        CREVP_3 = CREVP[3]
        CREVP_4 = CREVP[4]

        # smlt: five constants (lin et al. 1983)

        CSMLT = np.zeros(5)
        CSMLT[0] = Float(2.0) * constants.PIE * constants.TCOND * constants.RNZS / constants.HLTF
        CSMLT[1] = (
            Float(2.0) * constants.PIE * constants.VDIFU * constants.RNZS * constants.HLTC / constants.HLTF
        )
        CSMLT[2] = CSSUB[1]
        CSMLT[3] = CSSUB[2]
        CSMLT[4] = constants.CH2O / constants.HLTF

        CSMLT_0 = CSMLT[0]
        CSMLT_1 = CSMLT[1]
        CSMLT_2 = CSMLT[2]
        CSMLT_3 = CSMLT[3]
        CSMLT_4 = CSMLT[4]

        # gmlt: five constants

        CGMLT = np.zeros(5)
        CGMLT[0] = Float(2.0) * constants.PIE * constants.TCOND * constants.RNZG / constants.HLTF
        CGMLT[1] = (
            Float(2.0) * constants.PIE * constants.VDIFU * constants.RNZG * constants.HLTC / constants.HLTF
        )
        CGMLT[2] = CGSUB[1]
        CGMLT[3] = CGSUB[2]
        CGMLT[4] = constants.CH2O / constants.HLTF

        CGMLT_0 = CGMLT[0]
        CGMLT_1 = CGMLT[1]
        CGMLT_2 = CGMLT[2]
        CGMLT_3 = CGMLT[3]
        CGMLT_4 = CGMLT[4]

        return cls(
            C_AIR,
            C_VAP,
            P_NONHYDRO,
            D0_VAP,
            LV00,
            DO_SEDI_W,
            LATV,
            LATI,
            LATS,
            LAT2,
            LCP,
            ICP,
            TCP,
            MPDT,
            RDT,
            NTIMES,
            DTS,
            CPAUT,
            RDTS,
            FAC_IMLT,
            FAC_I2S,
            FAC_V2L,
            FAC_L2V,
            FAC_I2V,
            FAC_S2V,
            FAC_V2S,
            FAC_G2V,
            FAC_V2G,
            FAC_FRZ,
            CGACS,
            CSACW,
            CRACI,
            CSACI,
            CGACW,
            CGACI,
            CRACW,
            CGFR_0,
            CGFR_1,
            CSSUB_0,
            CSSUB_1,
            CSSUB_2,
            CSSUB_3,
            CSSUB_4,
            CGSUB_0,
            CGSUB_1,
            CGSUB_2,
            CGSUB_3,
            CGSUB_4,
            CREVP_0,
            CREVP_1,
            CREVP_2,
            CREVP_3,
            CREVP_4,
            CSMLT_0,
            CSMLT_1,
            CSMLT_2,
            CSMLT_3,
            CSMLT_4,
            CGMLT_0,
            CGMLT_1,
            CGMLT_2,
            CGMLT_3,
            CGMLT_4,
        )
