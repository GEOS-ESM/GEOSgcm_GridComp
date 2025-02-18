from ndsl.dsl.typing import Float


class config:
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
        self.phys_hydrostatic = phys_hydrostatic
        self.hydrostatic = hydrostatic
        self.dt_moist = dt_moist
        self.mp_time = mp_time
        self.t_min = t_min
        self.t_sub = t_sub
        self.tau_r2g = tau_r2g
        self.tau_smlt = tau_smlt
        self.tau_g2r = tau_g2r
        self.dw_land = dw_land
        self.dw_ocean = dw_ocean
        self.vi_fac = vi_fac
        self.vr_fac = vr_fac
        self.vs_fac = vs_fac
        self.vg_fac = vg_fac
        self.ql_mlt = ql_mlt
        self.do_qa = do_qa
        self.fix_negative = fix_negative
        self.vi_max = vi_max
        self.vs_max = vs_max
        self.vg_max = vg_max
        self.vr_max = vr_max
        self.qs_mlt = qs_mlt
        self.qs0_crt = qs0_crt
        self.qi_gen = qi_gen
        self.ql0_max = ql0_max
        self.qi0_max = qi0_max
        self.qi0_crt = qi0_crt
        self.qr0_crt = qr0_crt
        self.fast_sat_adj = fast_sat_adj
        self.rh_inc = rh_inc
        self.rh_ins = rh_ins
        self.rh_inr = rh_inr
        self.const_vi = const_vi
        self.const_vs = const_vs
        self.const_vg = const_vg
        self.const_vr = const_vr
        self.use_ccn = use_ccn
        self.rthreshu = rthreshu
        self.rthreshs = rthreshs
        self.ccn_l = ccn_l
        self.ccn_o = ccn_o
        self.qc_crt = qc_crt
        self.tau_g2v = tau_g2v
        self.tau_v2g = tau_v2g
        self.tau_s2v = tau_s2v
        self.tau_v2s = tau_v2s
        self.tau_revp = tau_revp
        self.tau_frz = tau_frz
        self.do_bigg = do_bigg
        self.do_evap = do_evap
        self.do_subl = do_subl
        self.sat_adj0 = sat_adj0
        self.c_piacr = c_piacr
        self.tau_imlt = tau_imlt
        self.tau_v2l = tau_v2l
        self.tau_l2v = tau_l2v
        self.tau_i2v = tau_i2v
        self.tau_i2s = tau_i2s
        self.tau_l2r = tau_l2r
        self.qi_lim = qi_lim
        self.ql_gen = ql_gen
        self.c_paut = c_paut
        self.c_psaci = c_psaci
        self.c_pgacs = c_pgacs
        self.c_pgaci = c_pgaci
        self.z_slope_liq = z_slope_liq
        self.z_slope_ice = z_slope_ice
        self.prog_ccn = prog_ccn
        self.c_cracw = c_cracw
        self.alin = alin
        self.clin = clin
        self.preciprad = preciprad
        self.cld_min = cld_min
        self.use_ppm = use_ppm
        self.mono_prof = mono_prof
        self.do_sedi_heat = do_sedi_heat
        self.sedi_transport = sedi_transport
        self.do_sedi_w = do_sedi_w
        self.de_ice = de_ice
        self.icloud_f = icloud_f
        self.irain_f = irain_f
        self.mp_print = mp_print
