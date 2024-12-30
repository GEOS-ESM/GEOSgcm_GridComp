"""GFDL_1M driver"""

import numpy as np
from gt4py.cartesian.gtscript import i32

import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_core import (
    create_temporaries,
    gfdl_1m_driver_loop_1,
    gfdl_1m_driver_loop_2,
    gfdl_1m_driver_loop_3,
    gfdl_1m_driver_loop_4,
    gfdl_1m_driver_postloop,
    gfdl_1m_driver_preloop,
)
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_tables import get_tables
from pyMoist.GFDL_1M.GFDL_1M_driver.icloud import icloud
from pyMoist.GFDL_1M.GFDL_1M_driver.terminal_fall import terminal_fall
from pyMoist.GFDL_1M.GFDL_1M_driver.warm_rain import warm_rain


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

        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vap based on hydrostatical property
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
        rdt = Float(1.0) / dt_moist
        ntimes = i32(np.round(dt_moist / mpdt))
        self.ntimes = ntimes  # save so that this can be referenced from within __call__

        # small time step:
        dts = dt_moist / Float(ntimes)

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        cpaut = c_paut * Float(0.104) * driver_constants.grav / Float(1.717e-5)

        # -----------------------------------------------------------------------
        # define conversion scalar / factor for icloud
        # -----------------------------------------------------------------------
        rdts = Float(1.0) / dts
        fac_imlt = Float(1.0) - np.exp(-dts / tau_imlt, dtype=Float)
        fac_i2s = Float(1.0) - np.exp(-dts / tau_i2s, dtype=Float)
        fac_v2l = Float(1.0) - np.exp(-dts / tau_v2l, dtype=Float)
        fac_l2v = Float(1.0) - np.exp(-dts / tau_l2v, dtype=Float)
        fac_i2v = Float(1.0) - np.exp(-dts / tau_i2v, dtype=Float)
        fac_s2v = Float(1.0) - np.exp(-dts / tau_s2v, dtype=Float)
        fac_v2s = Float(1.0) - np.exp(-dts / tau_v2s, dtype=Float)
        fac_g2v = Float(1.0) - np.exp(-dts / tau_g2v, dtype=Float)
        fac_v2g = Float(1.0) - np.exp(-dts / tau_v2g, dtype=Float)
        fac_frz = Float(1.0) - np.exp(-dts / tau_frz, dtype=Float)

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
            / (Float(4.0) * np.power(driver_constants.act[0], 0.8125, dtype=Float))
        )
        # decreasing csacw to reduce cloud water --- > snow

        craci = (
            driver_constants.pie
            * driver_constants.rnzr
            * alin
            * driver_constants.gam380
            / (Float(4.0) * np.power(driver_constants.act[1], 0.95, dtype=Float))
        )
        csaci = csacw * c_psaci

        cgacw = (
            driver_constants.pie
            * driver_constants.rnzg
            * driver_constants.gam350
            * driver_constants.gcon
            / (Float(4.0) * np.power(driver_constants.act[5], 0.875, dtype=Float))
        )

        cgaci = cgacw * c_pgaci

        cracw = craci  # cracw = 3.27206196043822
        cracw = c_cracw * cracw

        cssub = np.zeros(5)
        cgsub = np.zeros(5)
        crevp = np.zeros(5)

        cssub[0] = (
            Float(2.0)
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.tcond
            * driver_constants.rvgas
            * driver_constants.rnzs
        )
        cgsub[0] = (
            Float(2.0)
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.tcond
            * driver_constants.rvgas
            * driver_constants.rnzg
        )
        crevp[0] = (
            Float(2.0)
            * driver_constants.pie
            * driver_constants.vdifu
            * driver_constants.tcond
            * driver_constants.rvgas
            * driver_constants.rnzr
        )
        cssub[1] = Float(0.78) / np.sqrt(driver_constants.act[0], dtype=Float)
        cgsub[1] = Float(0.78) / np.sqrt(driver_constants.act[5], dtype=Float)
        crevp[1] = Float(0.78) / np.sqrt(driver_constants.act[1], dtype=Float)
        cssub[2] = (
            Float(0.31)
            * driver_constants.scm3
            * driver_constants.gam263
            * np.sqrt(clin / driver_constants.visk)
            / driver_constants.act[0] ** Float(0.65625)
        )
        cgsub[2] = (
            Float(0.31)
            * driver_constants.scm3
            * driver_constants.gam275
            * np.sqrt(driver_constants.gcon / driver_constants.visk)
            / driver_constants.act[5] ** Float(0.6875)
        )
        crevp[2] = (
            Float(0.31)
            * driver_constants.scm3
            * driver_constants.gam290
            * np.sqrt(alin / driver_constants.visk)
            / driver_constants.act[1] ** Float(0.725)
        )
        cssub[3] = driver_constants.tcond * driver_constants.rvgas
        cssub[4] = driver_constants.hlts ** 2 * driver_constants.vdifu
        cgsub[3] = cssub[3]
        crevp[3] = cssub[3]
        cgsub[4] = cssub[4]
        crevp[4] = driver_constants.hltc ** 2 * driver_constants.vdifu

        cgfr_0 = (
            Float(20.0e2)
            * driver_constants.pisq
            * driver_constants.rnzr
            * driver_constants.rhor
            / driver_constants.act[1] ** Float(1.75)
        )
        cgfr_1 = Float(0.66)

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
            Float(2.0)
            * driver_constants.pie
            * driver_constants.tcond
            * driver_constants.rnzs
            / driver_constants.hltf
        )
        csmlt[1] = (
            Float(2.0)
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
            Float(2.0)
            * driver_constants.pie
            * driver_constants.tcond
            * driver_constants.rnzg
            / driver_constants.hltf
        )
        cgmlt[1] = (
            Float(2.0)
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

        es0 = Float(6.107799961e2)  # ~6.1 mb
        ces0 = driver_constants.eps * es0

        # Check values for untested code paths
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
            dts,
        )

        # -----------------------------------------------------------------------
        # initialize temporaries
        # -----------------------------------------------------------------------

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

        # -----------------------------------------------------------------------
        # initialize masks
        # -----------------------------------------------------------------------
        self.is_frozen = quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.precip_fall = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
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

        # -----------------------------------------------------------------------
        # generate saturation specific humidity tables
        # -----------------------------------------------------------------------

        self.sat_tables = get_tables(stencil_factory.backend)

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

        self._gfdl_1m_driver_loop_1 = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_loop_1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "p_nonhydro": p_nonhydro,
                "const_vi": const_vi,
                "const_vs": const_vs,
                "const_vg": const_vg,
                "vi_fac": vi_fac,
                "vi_max": vi_max,
                "vs_fac": vs_fac,
                "vs_max": vs_max,
                "vg_fac": vg_fac,
                "vg_max": vg_max,
            },
        )

        self._terminal_fall = stencil_factory.from_dims_halo(
            func=terminal_fall,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": dts,
                "tau_imlt": tau_imlt,
                "ql_mlt": ql_mlt,
                "vi_fac": vi_fac,
                "do_sedi_w": do_sedi_w,
                "use_ppm": use_ppm,
                "tau_smlt": tau_smlt,
                "tau_g2r": tau_g2r,
                "c_air": c_air,
                "c_vap": c_vap,
                "d0_vap": d0_vap,
                "lv00": lv00,
            },
        )

        self._gfdl_1m_driver_loop_2 = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_loop_2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._warm_rain = stencil_factory.from_dims_halo(
            func=warm_rain,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": dts,
                "do_qa": do_qa,
                "rthreshs": rthreshs,
                "rthreshu": rthreshu,
                "irain_f": irain_f,
                "ql0_max": ql0_max,
                "z_slope_liq": z_slope_liq,
                "vr_fac": vr_fac,
                "const_vr": const_vr,
                "vr_max": vr_max,
                "tau_revp": tau_revp,
                "lv00": lv00,
                "d0_vap": d0_vap,
                "c_air": c_air,
                "c_vap": c_vap,
                "crevp_0": crevp_0,
                "crevp_1": crevp_1,
                "crevp_2": crevp_2,
                "crevp_3": crevp_3,
                "crevp_4": crevp_4,
                "cracw": cracw,
                "do_sedi_w": do_sedi_w,
                "use_ppm": use_ppm,
            },
        )

        self._gfdl_1m_driver_loop_3 = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_loop_3,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._icloud = stencil_factory.from_dims_halo(
            func=icloud,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": c_air,
                "c_vap": c_vap,
                "dts": dts,
                "rdts": rdts,
                "const_vi": const_vi,
                "fac_g2v": fac_g2v,
                "fac_i2s": fac_i2s,
                "fac_imlt": fac_imlt,
                "fac_frz": fac_frz,
                "fac_l2v": fac_l2v,
                "fac_s2v": fac_s2v,
                "fac_v2s": fac_v2s,
                "fac_v2g": fac_v2g,
                "cgacs": cgacs,
                "csacw": csacw,
                "csaci": csaci,
                "cgacw": cgacw,
                "cgaci": cgaci,
                "cgfr_0": cgfr_0,
                "cgfr_1": cgfr_1,
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
                "cssub_0": cssub_0,
                "cssub_1": cssub_1,
                "cssub_2": cssub_2,
                "cssub_3": cssub_3,
                "cssub_4": cssub_4,
                "qi0_crt": qi0_crt,
                "qs0_crt": qs0_crt,
                "qs_mlt": qs_mlt,
                "ql_mlt": ql_mlt,
                "z_slope_ice": z_slope_ice,
                "lv00": lv00,
                "d0_vap": d0_vap,
                "lat2": lat2,
                "do_qa": do_qa,
                "do_evap": do_evap,
                "do_bigg": do_bigg,
                "qc_crt": qc_crt,
                "qi_lim": qi_lim,
                "rh_inc": rh_inc,
                "rh_inr": rh_inr,
                "t_min": t_min,
                "t_sub": t_sub,
                "preciprad": preciprad,
                "icloud_f": icloud_f,
            },
        )

        self._gfdl_1m_driver_loop_4 = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_loop_4,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._gfdl_1m_driver_postloop = stencil_factory.from_dims_halo(
            func=gfdl_1m_driver_postloop,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": c_air,
                "c_vap": c_vap,
                "rdt": rdt,
                "do_sedi_w": do_sedi_w,
                "sedi_transport": sedi_transport,
                "do_qa": do_qa,
            },
        )

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

        # The driver modifies a number of variables (t, p, qX) but does not pass
        # the changes back to the rest of the model. To replicate this behavior,
        # temporary copies of these variables are used throughout the driver.
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
            self.qv0,
            self.ql0,
            self.qr0,
            self.qi0,
            self.qs0,
            self.qg0,
            self.qa0,
            self.qv1,
            self.ql1,
            self.qr1,
            self.qi1,
            self.qs1,
            self.qg1,
            self.qa1,
            dz,
            uin,
            vin,
            w,
            area,
            self.t1,
            self.dp1,
            self.omq,
            self.den,
            self.p_dry,
            self.m1,
            self.u1,
            self.v1,
            self.w1,
            self.onemsig,
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

        for n in range(self.ntimes):
            self._gfdl_1m_driver_loop_1(
                self.ql1,
                self.qi1,
                self.qs1,
                self.qg1,
                t,
                self.t1,
                dz,
                self.dz1,
                self.den,
                self.den1,
                self.denfac,
                self.p_dry,
                self.vti,
                self.vts,
                self.vtg,
                cnv_frc,
                anv_icefall,
                ls_icefall,
            )

            self._terminal_fall(
                self.t1,
                self.qv1,
                self.ql1,
                self.qr1,
                self.qg1,
                self.qs1,
                self.qi1,
                self.dz1,
                self.dp1,
                self.den1,
                self.vtg,
                self.vts,
                self.vti,
                self.rain1,
                self.graupel1,
                self.snow1,
                self.ice1,
                self.m1_sol,
                self.w1,
                self.ze,
                self.zt,
                self.is_frozen,
                self.precip_fall,
                self.melting_mask_1,
                self.melting_mask_2,
                self.current_k_level,
            )

            self._gfdl_1m_driver_loop_2(
                self.rain,
                self.graupel,
                self.snow,
                self.ice,
                self.rain1,
                self.graupel1,
                self.snow1,
                self.ice1,
            )

            self._warm_rain(
                self.dp1,
                self.dz1,
                self.t1,
                self.qv1,
                self.ql1,
                self.qr1,
                self.qi1,
                self.qs1,
                self.qg1,
                self.qa1,
                self.ccn,
                self.den,
                self.denfac,
                self.c_praut,
                self.vtr,
                self.evap1,
                self.m1_rain,
                self.w1,
                self.rh_limited,
                eis,
                self.onemsig,
                self.rain1,
                self.ze,
                self.zt,
                self.precip_fall,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
            )

            self._gfdl_1m_driver_loop_3(
                self.rain,
                self.rain1,
                self.evap1,
                self.revap,
                self.m1_rain,
                self.m2_rain,
                self.m1_sol,
                self.m2_sol,
                self.m1,
            )

            self._icloud(
                self.t1,
                self.p_dry,
                self.dp1,
                self.qv1,
                self.ql1,
                self.qr1,
                self.qi1,
                self.qs1,
                self.qg1,
                self.qa1,
                self.den1,
                self.denfac,
                self.vts,
                self.vtg,
                self.vtr,
                self.subl1,
                self.rh_limited,
                self.ccn,
                cnv_frc,
                srf_type,
                self.sat_tables.table1,
                self.sat_tables.table2,
                self.sat_tables.table3,
                self.sat_tables.table4,
                self.sat_tables.des1,
                self.sat_tables.des2,
                self.sat_tables.des3,
                self.sat_tables.des4,
            )

            self._gfdl_1m_driver_loop_4(
                self.isubl,
                self.subl1,
            )

        self._gfdl_1m_driver_postloop(
            self.qv0,
            self.ql0,
            self.qr0,
            self.qi0,
            self.qs0,
            self.qg0,
            self.qa0,
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
            t,
            self.t1,
            t_dt,
            w,
            self.w1,
            uin,
            self.u1,
            udt,
            vin,
            self.v1,
            vdt,
            dz,
            dp,
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
            self.m1,
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
        )

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
        dts: Float,
    ):
        failed_keywords = []
        if not phys_hydrostatic:
            failed_keywords.append("phys_hydrostatic")
        if hydrostatic:
            failed_keywords.append("hydrostatic")
        if const_vi:
            failed_keywords.append("const_vi")
        if const_vs:
            failed_keywords.append("const_vs")
        if const_vg:
            failed_keywords.append("const_vg")
        if const_vr:
            failed_keywords.append("const_vr")
        if use_ppm:
            failed_keywords.append("use_ppm")
        if not use_ccn:
            failed_keywords.append("use_ccn")
        if do_qa:
            failed_keywords.append("do_qa")
        if fast_set_adj:
            failed_keywords.append("fast_set_adj")
        if do_bigg:
            failed_keywords.append("do_bigg")
        if do_evap:
            failed_keywords.append("do_evap")
        if do_subl:
            failed_keywords.append("do_subl")
        if not z_slope_liq:
            failed_keywords.append("z_slope_liq")
        if not z_slope_ice:
            failed_keywords.append("z_slope_ice")
        if not prog_ccn:
            failed_keywords.append("prog_ccn")
        if not preciprad:
            failed_keywords.append("preciprad")
        if not mono_prof:
            failed_keywords.append("mono_prof")
        if do_sedi_heat:
            failed_keywords.append("do_sedi_heat")
        if not sedi_transport:
            failed_keywords.append("sedi_transport")
        if do_sedi_w:
            failed_keywords.append("do_sedi_w")
        if de_ice:
            failed_keywords.append("de_ice")
        if mp_print:
            failed_keywords.append("mp_print")
        if dts >= 300:
            failed_keywords.append("dts")

        if len(failed_keywords) > 0:
            raise ValueError(
                "One or more namelist parameters do not meet \
                    expected values. Failing parameters: ",
                failed_keywords,
            )
