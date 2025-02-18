import numpy as np
from gt4py.cartesian.gtscript import i32
from ndsl.dsl.typing import Float
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
import pyMoist.GFDL_1M.GFDL_1M_driver.constants as constants
from pyMoist.GFDL_1M.GFDL_1M_driver.config import config
from ndsl.dsl.typing import Int


def check_flags(
    GFDL_1M_config: config,
    dts: Float,
):
    failed_keywords = []
    if not GFDL_1M_config.phys_hydrostatic:
        failed_keywords.append("phys_hydrostatic")
    if GFDL_1M_config.hydrostatic:
        failed_keywords.append("hydrostatic")
    if GFDL_1M_config.const_vi:
        failed_keywords.append("const_vi")
    if GFDL_1M_config.const_vs:
        failed_keywords.append("const_vs")
    if GFDL_1M_config.const_vg:
        failed_keywords.append("const_vg")
    if GFDL_1M_config.const_vr:
        failed_keywords.append("const_vr")
    if GFDL_1M_config.use_ppm:
        failed_keywords.append("use_ppm")
    if not GFDL_1M_config.use_ccn:
        failed_keywords.append("use_ccn")
    if GFDL_1M_config.do_qa:
        failed_keywords.append("do_qa")
    if not GFDL_1M_config.fix_negative:
        failed_keywords.append("fix_negative")
    if GFDL_1M_config.fast_sat_adj:
        failed_keywords.append("fast_sat_adj")
    if GFDL_1M_config.do_bigg:
        failed_keywords.append("do_bigg")
    if GFDL_1M_config.do_evap:
        failed_keywords.append("do_evap")
    if GFDL_1M_config.do_subl:
        failed_keywords.append("do_subl")
    if not GFDL_1M_config.z_slope_liq:
        failed_keywords.append("z_slope_liq")
    if not GFDL_1M_config.z_slope_ice:
        failed_keywords.append("z_slope_ice")
    if not GFDL_1M_config.prog_ccn:
        failed_keywords.append("prog_ccn")
    if not GFDL_1M_config.preciprad:
        failed_keywords.append("preciprad")
    if not GFDL_1M_config.mono_prof:
        failed_keywords.append("mono_prof")
    if GFDL_1M_config.do_sedi_heat:
        failed_keywords.append("do_sedi_heat")
    if not GFDL_1M_config.sedi_transport:
        failed_keywords.append("sedi_transport")
    if GFDL_1M_config.do_sedi_w:
        failed_keywords.append("do_sedi_w")
    if GFDL_1M_config.de_ice:
        failed_keywords.append("de_ice")
    if GFDL_1M_config.mp_print:
        failed_keywords.append("mp_print")
    if dts >= 300:
        failed_keywords.append("dts")

    if len(failed_keywords) > 0:
        raise ValueError(
            "One or more namelist parameters do not meet \
                expected values. Failing parameters: ",
            failed_keywords,
        )


class config_constants:
    def __init__(self, GFDL_1M_config: config):
        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vap based on hydrostatical property
        # -----------------------------------------------------------------------

        if GFDL_1M_config.phys_hydrostatic or GFDL_1M_config.hydrostatic:
            self.C_AIR = constants.CP_AIR
            self.C_VAP = constants.CP_VAP
            self.P_NONHYDRO = False
        else:
            self.C_AIR = constants.CV_AIR
            self.C_VAP = constants.CV_VAP
            self.P_NONHYDRO = True
        self.D0_VAP = self.C_VAP - constants.C_LIQ
        self.LV00 = constants.HLV0 - self.D0_VAP * constants.T_ICE

        if GFDL_1M_config.hydrostatic:
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

        self.MPDT = min(GFDL_1M_config.dt_moist, GFDL_1M_config.mp_time)
        self.RDT = Float(1.0) / GFDL_1M_config.dt_moist
        self.NTIMES = i32(np.round(GFDL_1M_config.dt_moist / self.MPDT))

        # small time step:
        self.DTS = GFDL_1M_config.dt_moist / Float(self.NTIMES)

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        self.CPAUT = (
            GFDL_1M_config.c_paut * Float(0.104) * constants.GRAV / Float(1.717e-5)
        )

        # -----------------------------------------------------------------------
        # define conversion scalar / factor for icloud
        # -----------------------------------------------------------------------
        self.RDTS = Float(1.0) / self.DTS
        self.FAC_IMLT = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_imlt, dtype=Float
        )
        self.FAC_I2S = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_i2s, dtype=Float
        )
        self.FAC_V2L = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_v2l, dtype=Float
        )
        self.FAC_L2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_l2v, dtype=Float
        )
        self.FAC_I2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_i2v, dtype=Float
        )
        self.FAC_S2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_s2v, dtype=Float
        )
        self.FAC_V2S = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_v2s, dtype=Float
        )
        self.FAC_G2V = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_g2v, dtype=Float
        )
        self.FAC_V2G = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_v2g, dtype=Float
        )
        self.FAC_FRZ = Float(1.0) - np.exp(
            -self.DTS / GFDL_1M_config.tau_frz, dtype=Float
        )

        # -----------------------------------------------------------------------
        # constatns from setupm
        # -----------------------------------------------------------------------

        self.CGACS = constants.PISQ * constants.RNZG * constants.RNZS * constants.RHOS
        self.CGACS = self.CGACS * GFDL_1M_config.c_pgacs

        self.CSACW = (
            constants.PIE
            * constants.RNZS
            * GFDL_1M_config.clin
            * constants.GAM325
            / (Float(4.0) * np.power(constants.ACT[0], 0.8125, dtype=Float))
        )
        # decreasing csacw to reduce cloud water --- > snow

        self.CRACI = (
            constants.PIE
            * constants.RNZR
            * GFDL_1M_config.alin
            * constants.GAM380
            / (Float(4.0) * np.power(constants.ACT[1], 0.95, dtype=Float))
        )
        self.CSACI = self.CSACW * GFDL_1M_config.c_psaci

        self.CGACW = (
            constants.PIE
            * constants.RNZG
            * constants.GAM350
            * constants.GCON
            / (Float(4.0) * np.power(constants.ACT[5], 0.875, dtype=Float))
        )

        self.CGACI = self.CGACW * GFDL_1M_config.c_pgaci

        self.CRACW = self.CRACI  # cracw = 3.27206196043822
        self.CRACW = GFDL_1M_config.c_cracw * self.CRACW

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
            * np.sqrt(GFDL_1M_config.clin / constants.VISK, dtype=Float)
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
            * np.sqrt(GFDL_1M_config.alin / constants.VISK, dtype=Float)
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


class outputs:
    def __init__(self, quantity_factory):
        # -----------------------------------------------------------------------
        # initialize precipitation outputs
        # -----------------------------------------------------------------------

        self.rain = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.ice = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.graupel = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.m2_rain = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m2_sol = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.revap = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.isubl = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")


class temporaries:
    def __init__(self, quantity_factory):
        self.t1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.dp1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.omq = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qv0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ql0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qr0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qi0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qs0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qg0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qa0 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qv1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ql1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qr1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qi1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qs1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qg1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qa1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.dz1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.den = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.den1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.denfac = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.p_dry = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.u1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.v1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.w1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.onemsig = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.ccn = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.c_praut = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.rh_limited = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ze = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.zt = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.lhi = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.icpk = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.hold_data = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vti = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vts = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vtg = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vtr = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m1_sol = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m1_rain = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.rain1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.graupel1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.snow1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.ice1 = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.evap1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.subl1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")


class masks:
    def __init__(self, quantity_factory):
        self.is_frozen = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.precip_fall = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        # TODO: temporary code requires mask used within a double k loop
        # will be removed once a proper feature for double k loop is introduces
        self.melting_mask_1 = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.melting_mask_2 = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.current_k_level = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=Int
        )
        for k in range(self.current_k_level.view[:].shape[2]):
            self.current_k_level.view[:, :, k] = k
