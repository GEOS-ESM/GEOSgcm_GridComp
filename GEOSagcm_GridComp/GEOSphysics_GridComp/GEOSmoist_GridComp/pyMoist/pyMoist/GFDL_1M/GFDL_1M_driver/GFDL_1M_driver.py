"""GFDL_1M driver"""

import numpy as np
from gt4py.cartesian.gtscript import i32

import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from ndsl import QuantityFactory, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_core import (
    fall_speed,
    fix_negative_values,
    icloud_update,
    init_temporaries,
    terminal_fall_update,
    update_tendencies,
    warm_rain_update,
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
        self.dt_moist = dt_moist
        self.fix_negative = fix_negative
        self.sedi_transport = sedi_transport

        self.namelist_constants = namelist_setup(
            phys_hydrostatic,
            hydrostatic,
            dt_moist,
            mp_time,
            t_min,
            t_sub,
            tau_r2g,
            tau_smlt,
            tau_g2r,
            dw_land,
            dw_ocean,
            vi_fac,
            vr_fac,
            vs_fac,
            vg_fac,
            ql_mlt,
            do_qa,
            fix_negative,
            vi_max,
            vs_max,
            vg_max,
            vr_max,
            qs_mlt,
            qs0_crt,
            qi_gen,
            ql0_max,
            qi0_max,
            qi0_crt,
            qr0_crt,
            fast_sat_adj,
            rh_inc,
            rh_ins,
            rh_inr,
            const_vi,
            const_vs,
            const_vg,
            const_vr,
            use_ccn,
            rthreshu,
            rthreshs,
            ccn_l,
            ccn_o,
            qc_crt,
            tau_g2v,
            tau_v2g,
            tau_s2v,
            tau_v2s,
            tau_revp,
            tau_frz,
            do_bigg,
            do_evap,
            do_subl,
            sat_adj0,
            c_piacr,
            tau_imlt,
            tau_v2l,
            tau_l2v,
            tau_i2v,
            tau_i2s,
            tau_l2r,
            qi_lim,
            ql_gen,
            c_paut,
            c_psaci,
            c_pgacs,
            c_pgaci,
            z_slope_liq,
            z_slope_ice,
            prog_ccn,
            c_cracw,
            alin,
            clin,
            preciprad,
            cld_min,
            use_ppm,
            mono_prof,
            do_sedi_heat,
            sedi_transport,
            do_sedi_w,
            de_ice,
            icloud_f,
            irain_f,
            mp_print,
        )

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
            fix_negative,
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
            self.namelist_constants.DTS,
        )

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
        self._init_temporaries = stencil_factory.from_dims_halo(
            func=init_temporaries,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "cpaut": self.namelist_constants.CPAUT,
            },
        )

        self._gfdl_1m_driver_preloop = stencil_factory.from_dims_halo(
            func=fix_negative_values,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.namelist_constants.C_AIR,
                "c_vap": self.namelist_constants.C_VAP,
                "p_nonhydro": self.namelist_constants.P_NONHYDRO,
                "d0_vap": self.namelist_constants.D0_VAP,
                "lv00": self.namelist_constants.LV00,
                "latv": self.namelist_constants.LATV,
                "lati": self.namelist_constants.LATI,
                "lats": self.namelist_constants.LATS,
                "lat2": self.namelist_constants.LAT2,
                "lcp": self.namelist_constants.LCP,
                "icp": self.namelist_constants.ICP,
                "tcp": self.namelist_constants.TCP,
                "mpdt": self.namelist_constants.MPDT,
                "rdt": self.namelist_constants.RDT,
                "ntimes": self.namelist_constants.NTIMES,
                "dts": self.namelist_constants.DTS,
                "do_sedi_w": do_sedi_w,
                "cpaut": self.namelist_constants.CPAUT,
                "hydrostatic": hydrostatic,
                "phys_hydrostatic": phys_hydrostatic,
                "fix_negative": fix_negative,
                "sedi_transport": sedi_transport,
                "const_vi": const_vi,
                "const_vs": const_vs,
                "const_vg": const_vg,
            },
        )

        self._fall_speed = stencil_factory.from_dims_halo(
            func=fall_speed,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "p_nonhydro": self.namelist_constants.P_NONHYDRO,
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
                "dts": self.namelist_constants.DTS,
                "tau_imlt": tau_imlt,
                "ql_mlt": ql_mlt,
                "vi_fac": vi_fac,
                "do_sedi_w": do_sedi_w,
                "use_ppm": use_ppm,
                "tau_smlt": tau_smlt,
                "tau_g2r": tau_g2r,
                "c_air": self.namelist_constants.C_AIR,
                "c_vap": self.namelist_constants.C_VAP,
                "d0_vap": self.namelist_constants.D0_VAP,
                "lv00": self.namelist_constants.LV00,
            },
        )

        self._terminal_fall_update = stencil_factory.from_dims_halo(
            func=terminal_fall_update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._warm_rain = stencil_factory.from_dims_halo(
            func=warm_rain,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": self.namelist_constants.DTS,
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
                "lv00": self.namelist_constants.LV00,
                "d0_vap": self.namelist_constants.D0_VAP,
                "c_air": self.namelist_constants.C_AIR,
                "c_vap": self.namelist_constants.C_VAP,
                "crevp_0": self.namelist_constants.CREVP_0,
                "crevp_1": self.namelist_constants.CREVP_1,
                "crevp_2": self.namelist_constants.CREVP_2,
                "crevp_3": self.namelist_constants.CREVP_3,
                "crevp_4": self.namelist_constants.CREVP_4,
                "cracw": self.namelist_constants.CRACW,
                "do_sedi_w": do_sedi_w,
                "use_ppm": use_ppm,
            },
        )

        self._warm_rain_update = stencil_factory.from_dims_halo(
            func=warm_rain_update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._icloud = stencil_factory.from_dims_halo(
            func=icloud,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.namelist_constants.C_AIR,
                "c_vap": self.namelist_constants.C_VAP,
                "dts": self.namelist_constants.DTS,
                "rdts": self.namelist_constants.RDTS,
                "const_vi": const_vi,
                "fac_g2v": self.namelist_constants.FAC_G2V,
                "fac_i2s": self.namelist_constants.FAC_I2S,
                "fac_imlt": self.namelist_constants.FAC_IMLT,
                "fac_frz": self.namelist_constants.FAC_FRZ,
                "fac_l2v": self.namelist_constants.FAC_L2V,
                "fac_s2v": self.namelist_constants.FAC_S2V,
                "fac_v2s": self.namelist_constants.FAC_V2S,
                "fac_v2g": self.namelist_constants.FAC_V2G,
                "cgacs": self.namelist_constants.CGACS,
                "csacw": self.namelist_constants.CSACW,
                "csaci": self.namelist_constants.CSACI,
                "cgacw": self.namelist_constants.CGACW,
                "cgaci": self.namelist_constants.CGACI,
                "cgfr_0": self.namelist_constants.CGFR_0,
                "cgfr_1": self.namelist_constants.CGFR_1,
                "csmlt_0": self.namelist_constants.CSMLT_0,
                "csmlt_1": self.namelist_constants.CSMLT_1,
                "csmlt_2": self.namelist_constants.CSMLT_2,
                "csmlt_3": self.namelist_constants.CSMLT_3,
                "csmlt_4": self.namelist_constants.CSMLT_4,
                "cgmlt_0": self.namelist_constants.CGMLT_0,
                "cgmlt_1": self.namelist_constants.CGMLT_1,
                "cgmlt_2": self.namelist_constants.CGMLT_2,
                "cgmlt_3": self.namelist_constants.CGMLT_3,
                "cgmlt_4": self.namelist_constants.CGMLT_4,
                "cssub_0": self.namelist_constants.CSSUB_0,
                "cssub_1": self.namelist_constants.CSSUB_1,
                "cssub_2": self.namelist_constants.CSSUB_2,
                "cssub_3": self.namelist_constants.CSSUB_3,
                "cssub_4": self.namelist_constants.CSSUB_4,
                "qi0_crt": qi0_crt,
                "qs0_crt": qs0_crt,
                "qs_mlt": qs_mlt,
                "ql_mlt": ql_mlt,
                "z_slope_ice": z_slope_ice,
                "lv00": self.namelist_constants.LV00,
                "d0_vap": self.namelist_constants.D0_VAP,
                "lat2": self.namelist_constants.LAT2,
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

        self._icloud_update = stencil_factory.from_dims_halo(
            func=icloud_update,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": self.namelist_constants.C_AIR,
                "c_vap": self.namelist_constants.C_VAP,
                "rdt": self.namelist_constants.RDT,
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
        fr_land: FloatFieldIJ,
        cnv_frc: FloatFieldIJ,
        srf_type: FloatFieldIJ,
        eis: FloatFieldIJ,
        rhcrit3d: FloatField,
        anv_icefall: Float,
        ls_icefall: Float,
    ):
        rain_x = 7
        rain_y = 11
        print("Should be zero (at call) --> ", sum(sum(self.rain.view[:])))

        # The driver modifies a number of variables (t, p, qX) but does not pass
        # the changes back to the rest of the model. To replicate this behavior,
        # temporary copies of these variables are used throughout the driver.
        self._init_temporaries(
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
            self.rain,
            self.snow,
            self.graupel,
            self.ice,
            self.m2_rain,
            self.m2_sol,
            self.revap,
            self.isubl,
        )

        print(
            "Should be zero (post init_temporaries) --> ", sum(sum(self.rain.view[:]))
        )

        print("rain 11 7 post init: ", self.rain.view[rain_x, rain_y])

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

        for n in range(self.namelist_constants.NTIMES):
            self._fall_speed(
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

            print("rain 11 7 post terminal_fall: ", self.rain.view[rain_x, rain_y])
            print("rain1 11 7 post terminal_fall: ", self.rain1.view[rain_x, rain_y])

            self._terminal_fall_update(
                self.rain,
                self.graupel,
                self.snow,
                self.ice,
                self.rain1,
                self.graupel1,
                self.snow1,
                self.ice1,
            )

            print(
                "rain 11 7 post terminal_fall_update: ", self.rain.view[rain_x, rain_y]
            )
            print(
                "rain1 11 7 post terminal_fall_update: ",
                self.rain1.view[rain_x, rain_y],
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

            print("rain 11 7 post warm_rain: ", self.rain.view[rain_x, rain_y])
            print("rain1 11 7 post warm_rain: ", self.rain1.view[rain_x, rain_y])

            self._warm_rain_update(
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

            print("rain 11 7 post warm_rain_update: ", self.rain.view[rain_x, rain_y])
            print("rain1 11 7 post warm_rain_update: ", self.rain1.view[rain_x, rain_y])

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

            self._icloud_update(
                self.isubl,
                self.subl1,
            )

        self._update_tendencies(
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
            self.dt_moist,
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
        fix_negative: bool,
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
        if not fix_negative:
            failed_keywords.append("fix_negative")
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


class namelist_setup:
    def __init__(
        self,
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
        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vap based on hydrostatical property
        # -----------------------------------------------------------------------

        if phys_hydrostatic or hydrostatic:
            self.C_AIR = driver_constants.CP_AIR
            self.C_VAP = driver_constants.CP_VAP
            self.P_NONHYDRO = False
        else:
            self.C_AIR = driver_constants.CV_AIR
            self.C_VAP = driver_constants.CV_VAP
            self.P_NONHYDRO = True
        self.D0_VAP = self.C_VAP - driver_constants.C_LIQ
        self.LV00 = driver_constants.HLV0 - self.D0_VAP * driver_constants.T_ICE

        if hydrostatic:
            self.DO_SEDI_W = False

        # -----------------------------------------------------------------------
        # define latent heat coefficient used in wet bulb and bigg mechanism
        # -----------------------------------------------------------------------

        self.LATV = driver_constants.HLV
        self.LATI = driver_constants.HLF
        self.LATS = self.LATV + self.LATI
        self.LAT2 = self.LATS * self.LATS

        self.LCP = self.LATV / driver_constants.CP_AIR
        self.ICP = self.LATI / driver_constants.CP_AIR
        self.TCP = (self.LATV + self.LATI) / driver_constants.CP_AIR

        # -----------------------------------------------------------------------
        # define cloud microphysics sub time step
        # -----------------------------------------------------------------------

        self.MPDT = min(dt_moist, mp_time)
        self.RDT = Float(1.0) / dt_moist
        self.NTIMES = i32(np.round(dt_moist / self.MPDT))

        # small time step:
        self.DTS = dt_moist / Float(self.NTIMES)

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        self.CPAUT = c_paut * Float(0.104) * driver_constants.GRAV / Float(1.717e-5)

        # -----------------------------------------------------------------------
        # define conversion scalar / factor for icloud
        # -----------------------------------------------------------------------
        self.RDTS = Float(1.0) / self.DTS
        self.FAC_IMLT = Float(1.0) - np.exp(-self.DTS / tau_imlt, dtype=Float)
        self.FAC_I2S = Float(1.0) - np.exp(-self.DTS / tau_i2s, dtype=Float)
        self.FAC_V2L = Float(1.0) - np.exp(-self.DTS / tau_v2l, dtype=Float)
        self.FAC_L2V = Float(1.0) - np.exp(-self.DTS / tau_l2v, dtype=Float)
        self.FAC_I2V = Float(1.0) - np.exp(-self.DTS / tau_i2v, dtype=Float)
        self.FAC_S2V = Float(1.0) - np.exp(-self.DTS / tau_s2v, dtype=Float)
        self.FAC_V2S = Float(1.0) - np.exp(-self.DTS / tau_v2s, dtype=Float)
        self.FAC_G2V = Float(1.0) - np.exp(-self.DTS / tau_g2v, dtype=Float)
        self.FAC_V2G = Float(1.0) - np.exp(-self.DTS / tau_v2g, dtype=Float)
        self.FAC_FRZ = Float(1.0) - np.exp(-self.DTS / tau_frz, dtype=Float)

        # -----------------------------------------------------------------------
        # constatns from setupm
        # -----------------------------------------------------------------------

        self.CGACS = (
            driver_constants.PISQ
            * driver_constants.RNZG
            * driver_constants.RNZS
            * driver_constants.RHOS
        )
        self.CGACS = self.CGACS * c_pgacs

        self.CSACW = (
            driver_constants.PIE
            * driver_constants.RNZS
            * clin
            * driver_constants.GAM325
            / (Float(4.0) * np.power(driver_constants.ACT[0], 0.8125, dtype=Float))
        )
        # decreasing csacw to reduce cloud water --- > snow

        self.CRACI = (
            driver_constants.PIE
            * driver_constants.RNZR
            * alin
            * driver_constants.GAM380
            / (Float(4.0) * np.power(driver_constants.ACT[1], 0.95, dtype=Float))
        )
        self.CSACI = self.CSACW * c_psaci

        self.CGACW = (
            driver_constants.PIE
            * driver_constants.RNZG
            * driver_constants.GAM350
            * driver_constants.GCON
            / (Float(4.0) * np.power(driver_constants.ACT[5], 0.875, dtype=Float))
        )

        self.CGACI = self.CGACW * c_pgaci

        self.CRACW = self.CRACI  # cracw = 3.27206196043822
        self.CRACW = c_cracw * self.CRACW

        self.CSSUB = np.zeros(5)
        self.CGSUB = np.zeros(5)
        self.CREVP = np.zeros(5)

        self.CSSUB[0] = (
            Float(2.0)
            * driver_constants.PIE
            * driver_constants.VDIFU
            * driver_constants.TCOND
            * driver_constants.RVGAS
            * driver_constants.RNZS
        )
        self.CGSUB[0] = (
            Float(2.0)
            * driver_constants.PIE
            * driver_constants.VDIFU
            * driver_constants.TCOND
            * driver_constants.RVGAS
            * driver_constants.RNZG
        )
        self.CREVP[0] = (
            Float(2.0)
            * driver_constants.PIE
            * driver_constants.VDIFU
            * driver_constants.TCOND
            * driver_constants.RVGAS
            * driver_constants.RNZR
        )
        self.CSSUB[1] = Float(0.78) / np.sqrt(driver_constants.ACT[0], dtype=Float)
        self.CGSUB[1] = Float(0.78) / np.sqrt(driver_constants.ACT[5], dtype=Float)
        self.CREVP[1] = Float(0.78) / np.sqrt(driver_constants.ACT[1], dtype=Float)
        self.CSSUB[2] = (
            Float(0.31)
            * driver_constants.SCM3
            * driver_constants.GAM263
            * np.sqrt(clin / driver_constants.VISK, dtype=Float)
            / np.power(driver_constants.ACT[0], Float(0.65625), dtype=Float)
        )
        self.CGSUB[2] = (
            Float(0.31)
            * driver_constants.SCM3
            * driver_constants.GAM275
            * np.sqrt(driver_constants.GCON / driver_constants.VISK, dtype=Float)
            / np.power(driver_constants.ACT[5], Float(0.6875), dtype=Float)
        )
        self.CREVP[2] = (
            Float(0.31)
            * driver_constants.SCM3
            * driver_constants.GAM209
            * np.sqrt(alin / driver_constants.VISK, dtype=Float)
            / np.power(driver_constants.ACT[1], Float(0.725), dtype=Float)
        )
        self.CSSUB[3] = driver_constants.TCOND * driver_constants.RVGAS
        self.CSSUB[4] = (
            np.power(driver_constants.HLTS, 2, dtype=Float) * driver_constants.VDIFU
        )
        self.CGSUB[3] = self.CSSUB[3]
        self.CREVP[3] = self.CSSUB[3]
        self.CGSUB[4] = self.CSSUB[4]
        self.CREVP[4] = (
            np.power(driver_constants.HLTC, 2, dtype=Float) * driver_constants.VDIFU
        )

        self.CGFR_0 = (
            Float(20.0e2)
            * driver_constants.PISQ
            * driver_constants.RNZR
            * driver_constants.RHOR
            / np.power(driver_constants.ACT[1], Float(1.75), dtype=Float)
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
            * driver_constants.PIE
            * driver_constants.TCOND
            * driver_constants.RNZS
            / driver_constants.HLTF
        )
        self.CSMLT[1] = (
            Float(2.0)
            * driver_constants.PIE
            * driver_constants.VDIFU
            * driver_constants.RNZS
            * driver_constants.HLTC
            / driver_constants.HLTF
        )
        self.CSMLT[2] = self.CSSUB[1]
        self.CSMLT[3] = self.CSSUB[2]
        self.CSMLT[4] = driver_constants.CH2O / driver_constants.HLTF

        self.CSMLT_0 = self.CSMLT[0]
        self.CSMLT_1 = self.CSMLT[1]
        self.CSMLT_2 = self.CSMLT[2]
        self.CSMLT_3 = self.CSMLT[3]
        self.CSMLT_4 = self.CSMLT[4]

        # gmlt: five constants

        self.CGMLT = np.zeros(5)
        self.CGMLT[0] = (
            Float(2.0)
            * driver_constants.PIE
            * driver_constants.TCOND
            * driver_constants.RNZG
            / driver_constants.HLTF
        )
        self.CGMLT[1] = (
            Float(2.0)
            * driver_constants.PIE
            * driver_constants.VDIFU
            * driver_constants.RNZG
            * driver_constants.HLTC
            / driver_constants.HLTF
        )
        self.CGMLT[2] = self.CGSUB[1]
        self.CGMLT[3] = self.CGSUB[2]
        self.CGMLT[4] = driver_constants.CH2O / driver_constants.HLTF

        self.CGMLT_0 = self.CGMLT[0]
        self.CGMLT_1 = self.CGMLT[1]
        self.CGMLT_2 = self.CGMLT[2]
        self.CGMLT_3 = self.CGMLT[3]
        self.CGMLT_4 = self.CGMLT[4]
