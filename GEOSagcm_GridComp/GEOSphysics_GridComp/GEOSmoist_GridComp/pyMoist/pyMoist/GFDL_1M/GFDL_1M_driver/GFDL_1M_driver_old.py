"""GFDL_1M driver"""

from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, Int
import numpy as np
import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_core import (
    gfdl_1m_driver_preloop,
    gfdl_1m_driver_loop,
    gfdl_1m_driver_postloop,
    create_temporaries,
)
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_tables import get_tables


class GFDL_1M_driver:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        phys_hydrostatic: bool,
        hydrostatic: bool,
        dt_moist: Float,
        mp_time: Float,
        t_min: Float,
        t_sub: Float,
        tau_r2g: Float,
        tau_smlt: Float,
        tau_g2r: Float,
        dw_land: Float,
        dw_ocean: Float,
        vi_fac: Float,
        vr_fac: Float,
        vs_fac: Float,
        vg_fac: Float,
        ql_mlt: Float,
        do_qa: bool,
        fix_negative: bool,
        vi_max: Float,
        vs_max: Float,
        vg_max: Float,
        vr_max: Float,
        qs_mlt: Float,
        qs0_crt: Float,
        qi_gen: Float,
        ql0_max: Float,
        qi0_max: Float,
        qi0_crt: Float,
        qr0_crt: Float,
        fast_sat_adj: bool,
        rh_inc: Float,
        rh_ins: Float,
        rh_inr: Float,
        const_vi: bool,
        const_vs: bool,
        const_vg: bool,
        const_vr: bool,
        use_ccn: bool,
        rthreshu: Float,
        rthreshs: Float,
        ccn_l: Float,
        ccn_o: Float,
        qc_crt: Float,
        tau_g2v: Float,
        tau_v2g: Float,
        tau_s2v: Float,
        tau_v2s: Float,
        tau_revp: Float,
        tau_frz: Float,
        do_bigg: bool,
        do_evap: bool,
        do_subl: bool,
        sat_adj0: Float,
        c_piacr: Float,
        tau_imlt: Float,
        tau_v2l: Float,
        tau_l2v: Float,
        tau_i2v: Float,
        tau_i2s: Float,
        tau_l2r: Float,
        qi_lim: Float,
        ql_gen: Float,
        c_paut: Float,
        c_psaci: Float,
        c_pgacs: Float,
        c_pgaci: Float,
        z_slope_liq: bool,
        z_slope_ice: bool,
        prog_ccn: bool,
        c_cracw: Float,
        alin: Float,
        clin: Float,
        preciprad: bool,
        cld_min: Float,
        use_ppm: bool,
        mono_prof: bool,
        do_sedi_heat: bool,
        sedi_transport: bool,
        do_sedi_w: bool,
        de_ice: bool,
        icloud_f: Float,
        irain_f: Float,
        mp_print: bool,
    ):
        """Wrapper for the GFDL_1M driver"""

        self.check_flags(
            phys_hydrostatic,
            hydrostatic,
            const_vi,
            const_vs,
            const_vg,
            const_vr,
            use_ppm,
            use_ccn,
            do_qa,
            fast_sat_adj,
            do_bigg,
            do_evap,
            do_subl,
            z_slope_liq,
            z_slope_ice,
            prog_ccn,
            preciprad,
            mono_prof,
            do_sedi_heat,
            sedi_transport,
            do_sedi_w,
            de_ice,
            mp_print,
        )
        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vapor based on hydrostatical property
        # -----------------------------------------------------------------------

        if phys_hydrostatic or hydrostatic:
            c_air = driver_constants.cp_air
            c_vap = driver_constants.cp_vap
            p_nonhydro = False
        else:
            c_air = driver_constants.cv_air
            c_vap = driver_constants.cv_vap
            p_nonhydro = True
        d0_vap = c_vap - driver_constants.c_liq
        lv00 = driver_constants.hlv0 - d0_vap * driver_constants.t_ice

        if hydrostatic:
            do_sedi_w = False

        # -----------------------------------------------------------------------
        # define latent heat coefficient used in wet bulb and bigg mechanism
        # -----------------------------------------------------------------------

        latv = driver_constants.hlv
        lati = driver_constants.hlf
        lats = latv + lati
        lat2 = lats * lats

        lcp = latv / driver_constants.cp_air
        icp = lati / driver_constants.cp_air
        tcp = (latv + lati) / driver_constants.cp_air

        # tendency zero out for am moist processes should be done outside the driver

        # -----------------------------------------------------------------------
        # define cloud microphysics sub time step
        # -----------------------------------------------------------------------

        mpdt = min(dt_moist, mp_time)
        rdt = 1.0 / dt_moist
        ntimes = int(
            dt_moist / mpdt
        )  # need to implement a round function here. for now int does it, but this is temporary\
        self.ntimes = ntimes  # save so that this can be referenced from within __call__

        # small time step:
        dts = (
            dt_moist / ntimes
        )  # make sure dts is a Float, if not, need to cast ntimes to float before calcualtion

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        cpaut = c_paut * 0.104 * driver_constants.grav / 1.717e-5

        # -----------------------------------------------------------------------
        # define conversion scalar / factor for icloud
        # -----------------------------------------------------------------------
        rdts = 1.0 / dts
        fac_imlt = 1.0 - np.exp(-dts / tau_imlt)
        fac_i2s = 1.0 - np.exp(-dts / tau_i2s)
        fac_v2l = 1.0 - np.exp(-dts / tau_v2l)
        fac_l2v = 1.0 - np.exp(-dts / tau_l2v)
        fac_i2v = 1.0 - np.exp(-dts / tau_i2v)
        fac_s2v = 1.0 - np.exp(-dts / tau_s2v)
        fac_v2s = 1.0 - np.exp(-dts / tau_v2s)
        fac_g2v = 1.0 - np.exp(-dts / tau_g2v)
        fac_v2g = 1.0 - np.exp(-dts / tau_v2g)
        fac_frz = 1.0 - np.exp(-dts / tau_frz)

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

        # -----------------------------------------------------------------------
        # constatns from setupm
        # -----------------------------------------------------------------------

        cgacs = (
            driver_constants.pisq
            * driver_constants.rnzg
            * driver_constants.rnzs
            * driver_constants.rhos
        )
        cgacs = cgacs * c_pgacs

        csacw = (
            driver_constants.pie
            * driver_constants.rnzs
            * clin
            * driver_constants.gam325
            / (4.0 * driver_constants.act[0] ** 0.8125)
        )
        # decreasing csacw to reduce cloud water --- > snow

        craci = (
            driver_constants.pie
            * driver_constants.rnzr
            * alin
            * driver_constants.gam380
            / (4.0 * driver_constants.act[1] ** 0.95)
        )
        csaci = csacw * c_psaci

        cgacw = (
            driver_constants.pie
            * driver_constants.rnzg
            * driver_constants.gam350
            * driver_constants.gcon
            / (4.0 * driver_constants.act[5] ** 0.875)
        )

        cgaci = cgacw * c_pgaci

        cracw = craci  # cracw = 3.27206196043822
        cracw = c_cracw * cracw

        cssub = np.zeros(5)
        cgsub = np.zeros(5)
        crevp = np.zeros(5)

        cssub[0] = (
            2.0
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.tcond
            * driver_constants.rvgas
            * driver_constants.rnzs
        )
        cgsub[0] = (
            2.0
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.tcond
            * driver_constants.rvgas
            * driver_constants.rnzg
        )
        crevp[0] = (
            2.0
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.tcond
            * driver_constants.rvgas
            * driver_constants.rnzr
        )
        cssub[1] = 0.78 / np.sqrt(driver_constants.act[0])
        cgsub[1] = 0.78 / np.sqrt(driver_constants.act[5])
        crevp[1] = 0.78 / np.sqrt(driver_constants.act[1])
        cssub[2] = (
            0.31
            * driver_constants.scm3
            * driver_constants.gam263
            * np.sqrt(clin / driver_constants.visk)
            / driver_constants.act[0] ** 0.65625
        )
        cgsub[2] = (
            0.31
            * driver_constants.scm3
            * driver_constants.gam275
            * np.sqrt(driver_constants.gcon / driver_constants.visk)
            / driver_constants.act[5] ** 0.6875
        )
        crevp[2] = (
            0.31
            * driver_constants.scm3
            * driver_constants.gam290
            * np.sqrt(alin / driver_constants.visk)
            / driver_constants.act[1] ** 0.725
        )
        cssub[3] = driver_constants.tcond * driver_constants.rvgas
        cssub[4] = driver_constants.hlts**2 * driver_constants.vdifu
        cgsub[3] = cssub[3]
        crevp[3] = cssub[3]
        cgsub[4] = cssub[4]
        crevp[4] = driver_constants.hltc**2 * driver_constants.vdifu

        cgfr_0 = (
            20.0e2
            * driver_constants.pisq
            * driver_constants.rnzr
            * driver_constants.rhor
            / driver_constants.act[1] ** 1.75
        )
        cgfr_1 = 0.66

        cssub_0 = cssub[0]
        cssub_1 = cssub[1]
        cssub_2 = cssub[2]
        cssub_3 = cssub[3]
        cssub_4 = cssub[4]

        cgsub_0 = cgsub[0]
        cgsub_1 = cgsub[1]
        cgsub_2 = cgsub[2]
        cgsub_3 = cgsub[3]
        cgsub_4 = cgsub[4]

        crevp_0 = crevp[0]
        crevp_1 = crevp[1]
        crevp_2 = crevp[2]
        crevp_3 = crevp[3]
        crevp_4 = crevp[4]

        # smlt: five constants (lin et al. 1983)

        csmlt = np.zeros(5)
        csmlt[0] = (
            2.0
            * driver_constants.pie
            * driver_constants.tcond
            * driver_constants.rnzs
            / driver_constants.hltf
        )
        csmlt[1] = (
            2.0
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.rnzs
            * driver_constants.hltc
            / driver_constants.hltf
        )
        csmlt[2] = cssub[1]
        csmlt[3] = cssub[2]
        csmlt[4] = driver_constants.ch2o / driver_constants.hltf

        csmlt_0 = csmlt[0]
        csmlt_1 = csmlt[1]
        csmlt_2 = csmlt[2]
        csmlt_3 = csmlt[3]
        csmlt_4 = csmlt[4]

        # gmlt: five constants

        cgmlt = np.zeros(5)
        cgmlt[0] = (
            2.0
            * driver_constants.pie
            * driver_constants.tcond
            * driver_constants.rnzg
            / driver_constants.hltf
        )
        cgmlt[1] = (
            2.0
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.rnzg
            * driver_constants.hltc
            / driver_constants.hltf
        )
        cgmlt[2] = cgsub[1]
        cgmlt[3] = cgsub[2]
        cgmlt[4] = driver_constants.ch2o / driver_constants.hltf

        cgmlt_0 = cgmlt[0]
        cgmlt_1 = cgmlt[1]
        cgmlt_2 = cgmlt[2]
        cgmlt_3 = cgmlt[3]
        cgmlt_4 = cgmlt[4]

        es0 = 6.107799961e2  # ~6.1 mb
        ces0 = driver_constants.eps * es0

        # -----------------------------------------------------------------------
        # initialize temporaries
        # -----------------------------------------------------------------------

        self.t1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.dp1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.omq = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
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
        self.p_dry = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.u1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.v1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.w1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ccn = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.c_praut = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.rh_limited = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ze = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.zt = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.lhi = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.icpk = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.is_frozen = quantity_factory.ones([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.precip_fall = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.hold_data = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # -----------------------------------------------------------------------
        # generate saturation specific humidity tables
        # -----------------------------------------------------------------------

        self.sat_tables = get_tables()

        # -----------------------------------------------------------------------
        # initialize stencils
        # -----------------------------------------------------------------------

        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._create_temporaries = stencil_factory.from_dims_halo(
            func=create_temporaries,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "cpaut": cpaut,
            },
        )
        self._gfdl_1m_driver_preloop = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_preloop,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": c_air,
                "c_vap": c_vap,
                "p_nonhydro": p_nonhydro,
                "d0_vap": d0_vap,
                "lv00": lv00,
                "latv": latv,
                "lati": lati,
                "lats": lats,
                "lat2": lat2,
                "lcp": lcp,
                "icp": icp,
                "tcp": tcp,
                "mpdt": mpdt,
                "rdt": rdt,
                "ntimes": ntimes,
                "dts": dts,
                "do_sedi_w": do_sedi_w,
                "cpaut": cpaut,
                "hydrostatic": hydrostatic,
                "phys_hydrostatic": phys_hydrostatic,
                "fix_negative": fix_negative,
                "sedi_transport": sedi_transport,
                "const_vi": const_vi,
                "const_vs": const_vs,
                "const_vg": const_vg,
            },
        )

        self._gfdl_1m_driver_loop = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_loop,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": c_air,
                "c_vap": c_vap,
                "p_nonhydro": p_nonhydro,
                "d0_vap": d0_vap,
                "lv00": lv00,
                "latv": latv,
                "lati": lati,
                "lats": lats,
                "lat2": lat2,
                "lcp": lcp,
                "icp": icp,
                "tcp": tcp,
                "mpdt": mpdt,
                "rdt": rdt,
                "ntimes": ntimes,
                "dts": dts,
                "rdts": rdts,
                "do_sedi_w": do_sedi_w,
                "cpaut": cpaut,
                "hydrostatic": hydrostatic,
                "phys_hydrostatic": phys_hydrostatic,
                "fix_negative": fix_negative,
                "sedi_transport": sedi_transport,
                "const_vi": const_vi,
                "const_vs": const_vs,
                "const_vg": const_vg,
                "use_ppm": use_ppm,
                "fac_imlt": fac_imlt,
                "fac_i2s": fac_i2s,
                "fac_v2l": fac_v2l,
                "fac_l2v": fac_l2v,
                "fac_i2v": fac_i2v,
                "fac_s2v": fac_s2v,
                "fac_v2s": fac_v2s,
                "fac_g2v": fac_g2v,
                "fac_v2g": fac_v2g,
                "fac_frz": fac_frz,
                "ql_mlt": ql_mlt,
                "qi0_crt": qi0_crt,
                "vi_fac": vi_fac,
                "vi_max": vi_max,
                "vs_fac": vs_fac,
                "vs_max": vs_max,
                "vg_fac": vg_fac,
                "vg_max": vg_max,
                "cgacs": cgacs,
                "csacw": csacw,
                "craci": craci,
                "csaci": csaci,
                "cgacw": cgacw,
                "cgaci": cgaci,
                "cracw": cracw,
                "cgfr_0": cgfr_0,
                "cgfr_1": cgfr_1,
                "cssub_0": cssub_0,
                "cssub_1": cssub_1,
                "cssub_2": cssub_2,
                "cssub_3": cssub_3,
                "cssub_4": cssub_4,
                "cgsub_0": cgsub_0,
                "cgsub_1": cgsub_1,
                "cgsub_2": cgsub_2,
                "cgsub_3": cgsub_3,
                "cgsub_4": cgsub_4,
                "crevp_0": crevp_0,
                "crevp_1": crevp_1,
                "crevp_2": crevp_2,
                "crevp_3": crevp_3,
                "crevp_4": crevp_4,
                "csmlt_0": csmlt_0,
                "csmlt_1": csmlt_1,
                "csmlt_2": csmlt_2,
                "csmlt_3": csmlt_3,
                "csmlt_4": csmlt_4,
                "cgmlt_0": cgmlt_0,
                "cgmlt_1": cgmlt_1,
                "cgmlt_2": cgmlt_2,
                "cgmlt_3": cgmlt_3,
                "cgmlt_4": cgmlt_4,
                "tau_imlt": tau_imlt,
                "z_slope_ice": z_slope_ice,
                "do_qa": do_qa,
                "do_evap": do_evap,
                "do_bigg": do_bigg,
                "do_subl": do_subl,
                "ces0": ces0,
                "qc_crt": qc_crt,
                "qi0_crt": qi0_crt,
                "qs0_crt": qs0_crt,
                "qr0_crt": qr0_crt,
                "qi_gen": qi_gen,
                "ql_gen": ql_gen,
                "qi_lim": qi_lim,
                "qi0_max": qi0_max,
                "ql0_max": ql0_max,
                "ql_mlt": ql_mlt,
                "qs_mlt": qs_mlt,
                "rh_inc": rh_inc,
                "rh_ins": rh_ins,
                "rh_inr": rh_inr,
                "t_min": t_min,
                "t_sub": t_sub,
                "preciprad": preciprad,
                "icloud_f": icloud_f,
            },
        )

        self._gfdl_1m_driver_postloop = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_postloop,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": c_air,
                "c_vap": c_vap,
                "p_nonhydro": p_nonhydro,
                "d0_vap": d0_vap,
                "lv00": lv00,
                "latv": latv,
                "lati": lati,
                "lats": lats,
                "lat2": lat2,
                "lcp": lcp,
                "icp": icp,
                "tcp": tcp,
                "mpdt": mpdt,
                "rdt": rdt,
                "ntimes": ntimes,
                "dts": dts,
                "do_sedi_w": do_sedi_w,
                "cpaut": cpaut,
                "hydrostatic": hydrostatic,
                "phys_hydrostatic": phys_hydrostatic,
                "fix_negative": fix_negative,
                "sedi_transport": sedi_transport,
                "const_vi": const_vi,
                "const_vs": const_vs,
                "const_vg": const_vg,
            },
        )

        self.TESTVAR_1 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_2 = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

    def __call__(
        self,
        qv: FloatField,
        ql: FloatField,
        qr: FloatField,
        qi: FloatField,
        qs: FloatField,
        qg: FloatField,
        qa: FloatField,
        qn: FloatField,  # NACTL + NACTI
        qv_dt: FloatField,
        ql_dt: FloatField,
        qr_dt: FloatField,
        qi_dt: FloatField,
        qs_dt: FloatField,
        qg_dt: FloatField,
        qa_dt: FloatField,
        t_dt: FloatField,
        t: FloatField,
        w: FloatField,
        uin: FloatField,
        vin: FloatField,
        udt: FloatField,
        vdt: FloatField,
        dz: FloatField,
        dp: FloatField,
        area: FloatFieldIJ,
        dt_moist: Float,
        fr_land: FloatFieldIJ,
        cnv_frc: FloatFieldIJ,
        srf_type: FloatFieldIJ,
        eis: FloatFieldIJ,
        rhcrit3d: FloatField,
        anv_icefall: Float,
        ls_icefall: Float,
        hydrostatic: bool,
        phys_hydrostatic: bool,
        kmin: Int,
        kmax: Int,
        fix_negative: bool,
        sedi_transport: bool,
    ):

        # The driver modifies a number of variables (t, p, qX) but does not pass the changes back
        # to the rest of the model. To replicate this behavior, temporary copies of these variables
        # are used throughout the driver.
        self._create_temporaries(
            t,
            dp,
            rhcrit3d,
            qv,
            ql,
            qi,
            qr,
            qs,
            qg,
            qa,
            qn,
            dz,
            uin,
            vin,
            w,
            area,
            self.t1,
            self.dp1,
            self.omq,
            self.qv1,
            self.ql1,
            self.qr1,
            self.qi1,
            self.qs1,
            self.qg1,
            self.qa1,
            self.den,
            self.p_dry,
            self.m1,
            self.u1,
            self.v1,
            self.w1,
            self.ccn,
            self.c_praut,
            self.rh_limited,
        )

        self._gfdl_1m_driver_preloop(
            self.t1,
            self.qv1,
            self.ql1,
            self.qr1,
            self.qi1,
            self.qs1,
            self.qg1,
            self.dp1,
        )

        for n in range(1):  # range(self.ntimes):
            self._gfdl_1m_driver_loop(
                self.qv1,
                self.ql1,
                self.qr1,
                self.qi1,
                self.qs1,
                self.qg1,
                self.qa1,
                qn,
                qv_dt,
                ql_dt,
                qr_dt,
                qi_dt,
                qs_dt,
                qg_dt,
                qa_dt,
                t_dt,
                t,
                self.t1,
                w,
                self.w1,
                uin,
                self.u1,
                vin,
                self.v1,
                udt,
                vdt,
                dz,
                self.dz1,
                dp,
                self.dp1,
                self.den,
                self.den1,
                self.p_dry,
                area,
                dt_moist,
                fr_land,
                cnv_frc,
                srf_type,
                eis,
                self.rh_limited,
                anv_icefall,
                ls_icefall,
                self.revap,
                self.isubl,
                self.rain,
                self.snow,
                self.ice,
                self.graupel,
                self.m2_rain,
                self.m2_sol,
                self.ze,
                self.zt,
                self.lhi,
                self.icpk,
                self.is_frozen,
                self.precip_fall,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
                self.hold_data,
                kmin,
                kmax,
                self.TESTVAR_1,
                self.TESTVAR_2,
            )

        self._gfdl_1m_driver_postloop(
            self.qv1,
            self.ql1,
            self.qr1,
            self.qi1,
            self.qs1,
            self.qg1,
            self.qa1,
            self.ccn,
            qv_dt,
            ql_dt,
            qr_dt,
            qi_dt,
            qs_dt,
            qg_dt,
            qa_dt,
            t_dt,
            t,
            self.t1,
            self.w1,
            self.u1,
            self.v1,
            udt,
            vdt,
            dz,
            self.dp1,
            self.den,
            self.p_dry,
            area,
            dt_moist,
            fr_land,
            cnv_frc,
            srf_type,
            eis,
            self.rh_limited,
            anv_icefall,
            ls_icefall,
            self.revap,
            self.isubl,
            self.rain,
            self.snow,
            self.ice,
            self.graupel,
            self.m2_rain,
            self.m2_sol,
            self.TESTVAR_1,
            self.TESTVAR_2,
        )

        # print("TESTVAR_1: ", self.TESTVAR_1.view[:])
        # print("TESTVAR_2: ", self.TESTVAR_2.view[:])
        # print(self.snow.view[:])

    def check_flags(
        self,
        phys_hydrostatic: bool,
        hydrostatic: bool,
        const_vi: bool,
        const_vs: bool,
        const_vg: bool,
        const_vr: bool,
        use_ppm: bool,
        use_ccn: bool,
        do_qa: bool,
        fast_set_adj: bool,
        do_bigg: bool,
        do_evap: bool,
        do_subl: bool,
        z_slope_liq: bool,
        z_slope_ice: bool,
        prog_ccn: bool,
        preciprad: bool,
        mono_prof: bool,
        do_sedi_heat: bool,
        sedi_transport: bool,
        do_sedi_w: bool,
        de_ice: bool,
        mp_print: bool,
    ):
        if not phys_hydrostatic:
            raise ValueError(
                "Untested phys_hydrostatic option, disable this error manually to proceed."
            )
        if hydrostatic:
            raise ValueError(
                "Untested hydrostatic option, disable this error manually to proceed."
            )
        if const_vi:
            raise ValueError(
                "Untested const_vi option, disable this error manually to proceed."
            )
        if const_vs:
            raise ValueError(
                "Untested const_vs option, disable this error manually to proceed."
            )
        if const_vg:
            raise ValueError(
                "Untested const_vg option, disable this error manually to proceed."
            )
        if const_vr:
            raise ValueError(
                "const_vr option not implemented, contact dev team for further assistance."
            )
        if use_ppm:
            raise NotImplementedError(
                "use_ppm option not implemented, contact dev team for further assistance."
            )
        if not use_ccn:
            raise NotImplementedError(
                "use_ccn option not implemented, contact dev team for further assistance."
            )
        if do_qa:
            raise NotImplementedError(
                "do_qa option not implemented, contact dev team for further assistance."
            )
        if fast_set_adj:
            raise NotImplementedError(
                "fast_set_adj option not implemented, contact dev team for further assistance."
            )
        if do_bigg:
            raise NotImplementedError(
                "do_bigg option not implemented, contact dev team for further assistance."
            )
        if do_evap:
            raise NotImplementedError(
                "do_evap option not implemented, contact dev team for further assistance."
            )
        if do_subl:
            raise NotImplementedError(
                "do_subl option not implemented, contact dev team for further assistance."
            )
        if not z_slope_liq:
            raise NotImplementedError(
                "z_slope_liq option not implemented, contact dev team for further assistance."
            )
        if not z_slope_ice:
            raise NotImplementedError(
                "z_slope_ice option not implemented, contact dev team for further assistance."
            )
        if not prog_ccn:
            raise NotImplementedError(
                "prog_ccn option not implemented, contact dev team for further assistance."
            )
        if not preciprad:
            raise NotImplementedError(
                "preciprad option not implemented, contact dev team for further assistance."
            )
        if not mono_prof:
            raise NotImplementedError(
                "mono_prof option not implemented, contact dev team for further assistance."
            )
        if do_sedi_heat:
            raise NotImplementedError(
                "do_sedi_heat option not implemented, contact dev team for further assistance."
            )
        if not sedi_transport:
            raise NotImplementedError(
                "sedi_transport option not implemented, contact dev team for further assistance."
            )
        if do_sedi_w:
            raise NotImplementedError(
                "do_sedi_w option not implemented, contact dev team for further assistance."
            )
        if de_ice:
            raise NotImplementedError(
                "de_ice option not implemented, contact dev team for further assistance."
            )
        if mp_print:
            raise NotImplementedError(
                "mp_print option not implemented, contact dev team for further assistance."
            )
