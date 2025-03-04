import numpy as np
from gt4py.cartesian.gtscript import i32

from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.driver.config import MicrophysicsConfiguration
from pyMoist.GFDL_1M.driver.constants import constants


class ConfigConstants:
    def __init__(self, GFDL_1M_config: MicrophysicsConfiguration):
        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vap based on hydrostatical property
        # -----------------------------------------------------------------------

        if GFDL_1M_config.PHYS_HYDROSTATIC or GFDL_1M_config.HYDROSTATIC:
            self.C_AIR = constants.CP_AIR
            self.C_VAP = constants.CP_VAP
            self.P_NONHYDRO = False
        else:
            self.C_AIR = constants.CV_AIR
            self.C_VAP = constants.CV_VAP
            self.P_NONHYDRO = True
        self.D0_VAP = self.C_VAP - constants.C_LIQ
        self.LV00 = constants.HLV0 - self.D0_VAP * constants.T_ICE

        if GFDL_1M_config.HYDROSTATIC:
            self.DO_SEDI_W = False

        # -----------------------------------------------------------------------
        # define latent heat coefficient used in wet bulb and bigg mechanism
        # -----------------------------------------------------------------------

        self.LATV = constants.HLV
        self.LATI = constants.HLF
        self.LATS = self.LATV + self.LATI
        self.LAT2 = self.LATS * self.LATS

        self.LCP = self.LATV / constants.CP_AIR
        self.ICP = self.LATI / constants.CP_AIR
        self.TCP = (self.LATV + self.LATI) / constants.CP_AIR

        # -----------------------------------------------------------------------
        # define cloud microphysics sub time step
        # -----------------------------------------------------------------------

        self.MPDT = min(GFDL_1M_config.DT_MOIST, GFDL_1M_config.MP_TIME)
        self.RDT = Float(1.0) / GFDL_1M_config.DT_MOIST
        self.NTIMES = i32(np.round(GFDL_1M_config.DT_MOIST / self.MPDT))

        # small time step:
        self.DTS = GFDL_1M_config.DT_MOIST / Float(self.NTIMES)

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        self.CPAUT = (
            GFDL_1M_config.C_PAUT * Float(0.104) * constants.GRAV / Float(1.717e-5)
        )

        # -----------------------------------------------------------------------
        # define conversion scalar / factor for icloud
        # -----------------------------------------------------------------------
        self.RDTS = Float(1.0) / self.DTS
        self.FAC_IMLT = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_IMLT, dtype=Float
        )
        self.FAC_I2S = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_I2S, dtype=Float
        )
        self.FAC_V2L = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_V2L, dtype=Float
        )
        self.FAC_L2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_L2V, dtype=Float
        )
        self.FAC_I2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_I2V, dtype=Float
        )
        self.FAC_S2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_S2V, dtype=Float
        )
        self.FAC_V2S = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_V2S, dtype=Float
        )
        self.FAC_G2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_G2V, dtype=Float
        )
        self.FAC_V2G = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_V2G, dtype=Float
        )
        self.FAC_FRZ = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.TAU_FRZ, dtype=Float
        )

        # -----------------------------------------------------------------------
        # constatns from setupm
        # -----------------------------------------------------------------------

        self.CGACS = constants.PISQ * constants.RNZG * constants.RNZS * constants.RHOS
        self.CGACS = self.CGACS * GFDL_1M_config.C_PGACS

        self.CSACW = (
            constants.PIE
            * constants.RNZS
            * GFDL_1M_config.CLIN
            * constants.GAM325
            / (Float(4.0) * np.power(constants.ACT[0], 0.8125, dtype=Float))
        )
        # decreasing csacw to reduce cloud water --- > snow

        self.CRACI = (
            constants.PIE
            * constants.RNZR
            * GFDL_1M_config.ALIN
            * constants.GAM380
            / (Float(4.0) * np.power(constants.ACT[1], 0.95, dtype=Float))
        )
        self.CSACI = self.CSACW * GFDL_1M_config.C_PSACI

        self.CGACW = (
            constants.PIE
            * constants.RNZG
            * constants.GAM350
            * constants.GCON
            / (Float(4.0) * np.power(constants.ACT[5], 0.875, dtype=Float))
        )

        self.CGACI = self.CGACW * GFDL_1M_config.C_PGACI

        self.CRACW = self.CRACI  # cracw = 3.27206196043822
        self.CRACW = GFDL_1M_config.C_CRACW * self.CRACW

        self.CSSUB = np.zeros(5)
        self.CGSUB = np.zeros(5)
        self.CREVP = np.zeros(5)

        self.CSSUB[0] = (
            Float(2.0)
            * constants.PIE
            * constants.VDIFU
            * constants.TCOND
            * constants.RVGAS
            * constants.RNZS
        )
        self.CGSUB[0] = (
            Float(2.0)
            * constants.PIE
            * constants.VDIFU
            * constants.TCOND
            * constants.RVGAS
            * constants.RNZG
        )
        self.CREVP[0] = (
            Float(2.0)
            * constants.PIE
            * constants.VDIFU
            * constants.TCOND
            * constants.RVGAS
            * constants.RNZR
        )
        self.CSSUB[1] = Float(0.78) / np.sqrt(constants.ACT[0], dtype=Float)
        self.CGSUB[1] = Float(0.78) / np.sqrt(constants.ACT[5], dtype=Float)
        self.CREVP[1] = Float(0.78) / np.sqrt(constants.ACT[1], dtype=Float)
        self.CSSUB[2] = (
            Float(0.31)
            * constants.SCM3
            * constants.GAM263
            * np.sqrt(GFDL_1M_config.CLIN / constants.VISK, dtype=Float)
            / np.power(constants.ACT[0], Float(0.65625), dtype=Float)
        )
        self.CGSUB[2] = (
            Float(0.31)
            * constants.SCM3
            * constants.GAM275
            * np.sqrt(constants.GCON / constants.VISK, dtype=Float)
            / np.power(constants.ACT[5], Float(0.6875), dtype=Float)
        )
        self.CREVP[2] = (
            Float(0.31)
            * constants.SCM3
            * constants.GAM209
            * np.sqrt(GFDL_1M_config.ALIN / constants.VISK, dtype=Float)
            / np.power(constants.ACT[1], Float(0.725), dtype=Float)
        )
        self.CSSUB[3] = constants.TCOND * constants.RVGAS
        self.CSSUB[4] = np.power(constants.HLTS, 2, dtype=Float) * constants.VDIFU
        self.CGSUB[3] = self.CSSUB[3]
        self.CREVP[3] = self.CSSUB[3]
        self.CGSUB[4] = self.CSSUB[4]
        self.CREVP[4] = np.power(constants.HLTC, 2, dtype=Float) * constants.VDIFU

        self.CGFR_0 = (
            Float(20.0e2)
            * constants.PISQ
            * constants.RNZR
            * constants.RHOR
            / np.power(constants.ACT[1], Float(1.75), dtype=Float)
        )
        self.CGFR_1 = Float(0.66)

        self.CSSUB_0 = self.CSSUB[0]
        self.CSSUB_1 = self.CSSUB[1]
        self.CSSUB_2 = self.CSSUB[2]
        self.CSSUB_3 = self.CSSUB[3]
        self.CSSUB_4 = self.CSSUB[4]

        self.CGSUB_0 = self.CGSUB[0]
        self.CGSUB_1 = self.CGSUB[1]
        self.CGSUB_2 = self.CGSUB[2]
        self.CGSUB_3 = self.CGSUB[3]
        self.CGSUB_4 = self.CGSUB[4]

        self.CREVP_0 = self.CREVP[0]
        self.CREVP_1 = self.CREVP[1]
        self.CREVP_2 = self.CREVP[2]
        self.CREVP_3 = self.CREVP[3]
        self.CREVP_4 = self.CREVP[4]

        # smlt: five constants (lin et al. 1983)

        self.CSMLT = np.zeros(5)
        self.CSMLT[0] = (
            Float(2.0)
            * constants.PIE
            * constants.TCOND
            * constants.RNZS
            / constants.HLTF
        )
        self.CSMLT[1] = (
            Float(2.0)
            * constants.PIE
            * constants.VDIFU
            * constants.RNZS
            * constants.HLTC
            / constants.HLTF
        )
        self.CSMLT[2] = self.CSSUB[1]
        self.CSMLT[3] = self.CSSUB[2]
        self.CSMLT[4] = constants.CH2O / constants.HLTF

        self.CSMLT_0 = self.CSMLT[0]
        self.CSMLT_1 = self.CSMLT[1]
        self.CSMLT_2 = self.CSMLT[2]
        self.CSMLT_3 = self.CSMLT[3]
        self.CSMLT_4 = self.CSMLT[4]

        # gmlt: five constants

        self.CGMLT = np.zeros(5)
        self.CGMLT[0] = (
            Float(2.0)
            * constants.PIE
            * constants.TCOND
            * constants.RNZG
            / constants.HLTF
        )
        self.CGMLT[1] = (
            Float(2.0)
            * constants.PIE
            * constants.VDIFU
            * constants.RNZG
            * constants.HLTC
            / constants.HLTF
        )
        self.CGMLT[2] = self.CGSUB[1]
        self.CGMLT[3] = self.CGSUB[2]
        self.CGMLT[4] = constants.CH2O / constants.HLTF

        self.CGMLT_0 = self.CGMLT[0]
        self.CGMLT_1 = self.CGMLT[1]
        self.CGMLT_2 = self.CGMLT[2]
        self.CGMLT_3 = self.CGMLT[3]
        self.CGMLT_4 = self.CGMLT[4]
