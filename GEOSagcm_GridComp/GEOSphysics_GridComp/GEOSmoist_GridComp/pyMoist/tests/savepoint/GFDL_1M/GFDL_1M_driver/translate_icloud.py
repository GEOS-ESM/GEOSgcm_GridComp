import numpy as np

import pyMoist.GFDL_1M.GFDL_1M_driver.constants as constants
from ndsl import Namelist, Quantity, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M_driver.sat_tables import get_tables
from pyMoist.GFDL_1M.GFDL_1M_driver.icloud import icloud


class Translateicloud(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "t1_icloud": {},
            "p1_icloud": {},
            "dp1_icloud": {},
            "qv_icloud": {},
            "ql_icloud": {},
            "qr_icloud": {},
            "qg_icloud": {},
            "qs_icloud": {},
            "qi_icloud": {},
            "qa_icloud": {},
            "den1_icloud": {},
            "denfac_icloud": {},
            "vtg_icloud": {},
            "vts_icloud": {},
            "vtr_icloud": {},
            "ccn_icloud": {},
            "rh_limited_icloud": {},
            "cnv_frc_icloud": {},
            "srf_type_icloud": {},
            "mp_time_icloud": {},
            "t_min_icloud": {},
            "t_sub_icloud": {},
            "tau_r2g_icloud": {},
            "tau_smlt_icloud": {},
            "tau_g2r_icloud": {},
            "dw_land_icloud": {},
            "dw_ocean_icloud": {},
            "vi_fac_icloud": {},
            "vr_fac_icloud": {},
            "vs_fac_icloud": {},
            "vg_fac_icloud": {},
            "ql_mlt_icloud": {},
            "do_qa_icloud": {},
            "fix_negative_icloud": {},
            "vi_max_icloud": {},
            "vs_max_icloud": {},
            "vg_max_icloud": {},
            "vr_max_icloud": {},
            "qs_mlt_icloud": {},
            "qs0_crt_icloud": {},
            "qi_gen_icloud": {},
            "ql0_max_icloud": {},
            "qi0_max_icloud": {},
            "qi0_crt_icloud": {},
            "qr0_crt_icloud": {},
            "fast_sat_adj_icloud": {},
            "rh_inc_icloud": {},
            "rh_ins_icloud": {},
            "rh_inr_icloud": {},
            "const_vi_icloud": {},
            "const_vs_icloud": {},
            "const_vg_icloud": {},
            "const_vr_icloud": {},
            "use_ccn_icloud": {},
            "rthreshu_icloud": {},
            "rthreshs_icloud": {},
            "ccn_l_icloud": {},
            "ccn_o_icloud": {},
            "qc_crt_icloud": {},
            "tau_g2v_icloud": {},
            "tau_v2g_icloud": {},
            "tau_s2v_icloud": {},
            "tau_v2s_icloud": {},
            "tau_revp_icloud": {},
            "tau_frz_icloud": {},
            "do_bigg_icloud": {},
            "do_evap_icloud": {},
            "do_subl_icloud": {},
            "sat_adj0_icloud": {},
            "c_piacr_icloud": {},
            "tau_imlt_icloud": {},
            "tau_v2l_icloud": {},
            "tau_l2v_icloud": {},
            "tau_i2v_icloud": {},
            "tau_i2s_icloud": {},
            "tau_l2r_icloud": {},
            "qi_lim_icloud": {},
            "ql_gen_icloud": {},
            "c_paut_icloud": {},
            "c_psaci_icloud": {},
            "c_pgacs_icloud": {},
            "c_pgaci_icloud": {},
            "z_slope_liq_icloud": {},
            "z_slope_ice_icloud": {},
            "prog_ccn_icloud": {},
            "c_cracw_icloud": {},
            "alin_icloud": {},
            "clin_icloud": {},
            "preciprad_icloud": {},
            "cld_min_icloud": {},
            "use_ppm_icloud": {},
            "mono_prof_icloud": {},
            "do_sedi_heat_icloud": {},
            "sedi_transport_icloud": {},
            "do_sedi_w_icloud": {},
            "de_ice_icloud": {},
            "icloud_f_icloud": {},
            "irain_f_icloud": {},
            "mp_print_icloud": {},
            "hydrostatic_icloud": {},
            "dts_icloud": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "t1_icloud": self.grid.compute_dict(),
            "qv_icloud": self.grid.compute_dict(),
            "ql_icloud": self.grid.compute_dict(),
            "qr_icloud": self.grid.compute_dict(),
            "qg_icloud": self.grid.compute_dict(),
            "qs_icloud": self.grid.compute_dict(),
            "qi_icloud": self.grid.compute_dict(),
            "qa_icloud": self.grid.compute_dict(),
            "subl1_icloud": self.grid.compute_dict(),
        }

    def make_ij_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM],
            "n/a",
        )
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty

    def make_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM],
            "n/a",
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty

    def compute(self, inputs):
        # FloatField Variables
        t1_icloud = self.make_ijk_field(inputs["t1_icloud"])
        p1_icloud = self.make_ijk_field(inputs["p1_icloud"])
        dp1_icloud = self.make_ijk_field(inputs["dp1_icloud"])
        qv_icloud = self.make_ijk_field(inputs["qv_icloud"])
        ql_icloud = self.make_ijk_field(inputs["ql_icloud"])
        qr_icloud = self.make_ijk_field(inputs["qr_icloud"])
        qg_icloud = self.make_ijk_field(inputs["qg_icloud"])
        qs_icloud = self.make_ijk_field(inputs["qs_icloud"])
        qi_icloud = self.make_ijk_field(inputs["qi_icloud"])
        qa_icloud = self.make_ijk_field(inputs["qa_icloud"])
        den1_icloud = self.make_ijk_field(inputs["den1_icloud"])
        denfac_icloud = self.make_ijk_field(inputs["denfac_icloud"])
        vtg_icloud = self.make_ijk_field(inputs["vtg_icloud"])
        vts_icloud = self.make_ijk_field(inputs["vts_icloud"])
        vtr_icloud = self.make_ijk_field(inputs["vtr_icloud"])
        ccn_icloud = self.make_ijk_field(inputs["ccn_icloud"])
        rh_limited_icloud = self.make_ijk_field(inputs["rh_limited_icloud"])
        cnv_frc_icloud = self.make_ij_field(inputs["cnv_frc_icloud"])
        srf_type_icloud = self.make_ij_field(inputs["srf_type_icloud"])

        # Float Variables
        # Namelist options
        mp_time = Float(inputs["mp_time_icloud"][0])
        t_min = Float(inputs["t_min_icloud"][0])
        t_sub = Float(inputs["t_sub_icloud"][0])
        tau_r2g = Float(inputs["tau_r2g_icloud"][0])
        tau_smlt = Float(inputs["tau_smlt_icloud"][0])
        tau_g2r = Float(inputs["tau_g2r_icloud"][0])
        dw_land = Float(inputs["dw_land_icloud"][0])
        dw_ocean = Float(inputs["dw_ocean_icloud"][0])
        vi_fac = Float(inputs["vi_fac_icloud"][0])
        vr_fac = Float(inputs["vr_fac_icloud"][0])
        vs_fac = Float(inputs["vs_fac_icloud"][0])
        vg_fac = Float(inputs["vg_fac_icloud"][0])
        ql_mlt = Float(inputs["ql_mlt_icloud"][0])
        do_qa = Float(inputs["do_qa_icloud"][0])
        fix_negative = Float(inputs["fix_negative_icloud"][0])
        vi_max = Float(inputs["vi_max_icloud"][0])
        vs_max = Float(inputs["vs_max_icloud"][0])
        vg_max = Float(inputs["vg_max_icloud"][0])
        vr_max = Float(inputs["vr_max_icloud"][0])
        qs_mlt = Float(inputs["qs_mlt_icloud"][0])
        qs0_crt = Float(inputs["qs0_crt_icloud"][0])
        qi_gen = Float(inputs["qi_gen_icloud"][0])
        ql0_max = Float(inputs["ql0_max_icloud"][0])
        qi0_max = Float(inputs["qi0_max_icloud"][0])
        qi0_crt = Float(inputs["qi0_crt_icloud"][0])
        qr0_crt = Float(inputs["qr0_crt_icloud"][0])
        fast_sat_adj = Float(inputs["fast_sat_adj_icloud"][0])
        rh_inc = Float(inputs["rh_inc_icloud"][0])
        rh_ins = Float(inputs["rh_ins_icloud"][0])
        rh_inr = Float(inputs["rh_inr_icloud"][0])
        const_vi = Float(inputs["const_vi_icloud"][0])
        const_vs = Float(inputs["const_vs_icloud"][0])
        const_vg = Float(inputs["const_vg_icloud"][0])
        const_vr = Float(inputs["const_vr_icloud"][0])
        use_ccn = Float(inputs["use_ccn_icloud"][0])
        rthreshu = Float(inputs["rthreshu_icloud"][0])
        rthreshs = Float(inputs["rthreshs_icloud"][0])
        ccn_l = Float(inputs["ccn_l_icloud"][0])
        ccn_o = Float(inputs["ccn_o_icloud"][0])
        qc_crt = Float(inputs["qc_crt_icloud"][0])
        tau_g2v = Float(inputs["tau_g2v_icloud"][0])
        tau_v2g = Float(inputs["tau_v2g_icloud"][0])
        tau_s2v = Float(inputs["tau_s2v_icloud"][0])
        tau_v2s = Float(inputs["tau_v2s_icloud"][0])
        tau_revp = Float(inputs["tau_revp_icloud"][0])
        tau_frz = Float(inputs["tau_frz_icloud"][0])
        do_bigg = Float(inputs["do_bigg_icloud"][0])
        do_evap = Float(inputs["do_evap_icloud"][0])
        do_subl = Float(inputs["do_subl_icloud"][0])
        sat_adj0 = Float(inputs["sat_adj0_icloud"][0])
        c_piacr = Float(inputs["c_piacr_icloud"][0])
        tau_imlt = Float(inputs["tau_imlt_icloud"][0])
        tau_v2l = Float(inputs["tau_v2l_icloud"][0])
        tau_l2v = Float(inputs["tau_l2v_icloud"][0])
        tau_i2v = Float(inputs["tau_i2v_icloud"][0])
        tau_i2s = Float(inputs["tau_i2s_icloud"][0])
        tau_l2r = Float(inputs["tau_l2r_icloud"][0])
        qi_lim = Float(inputs["qi_lim_icloud"][0])
        ql_gen = Float(inputs["ql_gen_icloud"][0])
        c_paut = Float(inputs["c_paut_icloud"][0])
        c_psaci = Float(inputs["c_psaci_icloud"][0])
        c_pgacs = Float(inputs["c_pgacs_icloud"][0])
        c_pgaci = Float(inputs["c_pgaci_icloud"][0])
        z_slope_liq = Float(inputs["z_slope_liq_icloud"][0])
        z_slope_ice = Float(inputs["z_slope_ice_icloud"][0])
        prog_ccn = Float(inputs["prog_ccn_icloud"][0])
        c_cracw = Float(inputs["c_cracw_icloud"][0])
        alin = Float(inputs["alin_icloud"][0])
        clin = Float(inputs["clin_icloud"][0])
        preciprad = Float(inputs["preciprad_icloud"][0])
        cld_min = Float(inputs["cld_min_icloud"][0])
        use_ppm = Float(inputs["use_ppm_icloud"][0])
        mono_prof = Float(inputs["mono_prof_icloud"][0])
        do_sedi_heat = Float(inputs["do_sedi_heat_icloud"][0])
        sedi_transport = Float(inputs["sedi_transport_icloud"][0])
        do_sedi_w = Float(inputs["do_sedi_w_icloud"][0])
        de_ice = Float(inputs["de_ice_icloud"][0])
        icloud_f = Float(inputs["icloud_f_icloud"][0])
        irain_f = Float(inputs["irain_f_icloud"][0])
        mp_print = Float(inputs["mp_print_icloud"][0])
        hydrostatic = Float(inputs["hydrostatic_icloud"][0])
        dts = Float(inputs["dts_icloud"][0])

        # Revert to boolean
        if hydrostatic == 1:
            hydrostatic = True
        else:
            hydrostatic = False
        if fix_negative == 1:
            fix_negative = True
        else:
            fix_negative = False
        if sedi_transport == 1:
            sedi_transport = True
        else:
            sedi_transport = False
        if const_vi == 1:
            const_vi = True
        else:
            const_vi = False
        if const_vs == 1:
            const_vs = True
        else:
            const_vs = False
        if const_vg == 1:
            const_vg = True
        else:
            const_vg = False
        if use_ppm == 1:
            use_ppm = True
        else:
            use_ppm = False
        if do_qa == 1:
            do_qa = True
        else:
            do_qa = False
        if fast_sat_adj == 1:
            fast_sat_adj = True
        else:
            fast_sat_adj = False

        # Calculate additional constants
        # -----------------------------------------------------------------------
        # define heat capacity of dry air and water vap based on hydrostatical property
        # -----------------------------------------------------------------------

        phys_hydrostatic = (
            True  # not about to try to serialize this, it should always be true
        )
        if phys_hydrostatic or hydrostatic:
            c_air = constants.CP_AIR
            c_vap = constants.CP_VAP
            p_nonhydro = False
        else:
            c_air = constants.CV_AIR
            c_vap = constants.CV_VAP
            p_nonhydro = True
        d0_vap = c_vap - constants.C_LIQ
        lv00 = constants.HLV0 - d0_vap * constants.T_ICE

        if hydrostatic:
            do_sedi_w = False

        # -----------------------------------------------------------------------
        # define latent heat coefficient used in wet bulb and bigg mechanism
        # -----------------------------------------------------------------------

        latv = constants.HLV
        lati = constants.HLF
        lats = latv + lati
        lat2 = lats * lats

        lcp = latv / constants.CP_AIR
        icp = lati / constants.CP_AIR
        tcp = (latv + lati) / constants.CP_AIR

        # -----------------------------------------------------------------------
        # calculate cloud condensation nuclei (ccn)
        # the following is based on klein eq. 15
        # -----------------------------------------------------------------------

        cpaut = c_paut * 0.104 * constants.GRAV / 1.717e-5

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
        # constatns from setupm
        # -----------------------------------------------------------------------

        cgacs = (
            constants.PISQ
            * constants.RNZG
            * constants.RNZS
            * constants.RHOS
        )
        cgacs = cgacs * c_pgacs

        csacw = (
            constants.PIE
            * constants.RNZS
            * clin
            * constants.GAM325
            / (4.0 * constants.ACT[0] ** 0.8125)
        )
        # decreasing csacw to reduce cloud water --- > snow

        craci = (
            constants.PIE
            * constants.RNZR
            * alin
            * constants.GAM380
            / (4.0 * constants.ACT[1] ** 0.95)
        )
        csaci = csacw * c_psaci

        cgacw = (
            constants.PIE
            * constants.RNZG
            * constants.GAM350
            * constants.GCON
            / (4.0 * constants.ACT[5] ** 0.875)
        )

        cgaci = cgacw * c_pgaci

        cracw = craci  # cracw = 3.27206196043822
        cracw = c_cracw * cracw

        cssub = np.zeros(5)
        cgsub = np.zeros(5)
        crevp = np.zeros(5)

        cssub[0] = (
            2.0
            * constants.PIE
            * constants.VDIFU
            * constants.TCOND
            * constants.RVGAS
            * constants.RNZS
        )
        cgsub[0] = (
            2.0
            * constants.PIE
            * constants.VDIFU
            * constants.TCOND
            * constants.RVGAS
            * constants.RNZG
        )
        crevp[0] = (
            2.0
            * constants.PIE
            * constants.VDIFU
            * constants.TCOND
            * constants.RVGAS
            * constants.RNZR
        )
        cssub[1] = 0.78 / np.sqrt(constants.ACT[0])
        cgsub[1] = 0.78 / np.sqrt(constants.ACT[5])
        crevp[1] = 0.78 / np.sqrt(constants.ACT[1])
        cssub[2] = (
            0.31
            * constants.SCM3
            * constants.GAM263
            * np.sqrt(clin / constants.VISK)
            / constants.ACT[0] ** 0.65625
        )
        cgsub[2] = (
            0.31
            * constants.SCM3
            * constants.GAM275
            * np.sqrt(constants.GCON / constants.VISK)
            / constants.ACT[5] ** 0.6875
        )
        crevp[2] = (
            0.31
            * constants.SCM3
            * constants.GAM209
            * np.sqrt(alin / constants.VISK)
            / constants.ACT[1] ** 0.725
        )
        cssub[3] = constants.TCOND * constants.RVGAS
        cssub[4] = constants.HLTS**2 * constants.VDIFU
        cgsub[3] = cssub[3]
        crevp[3] = cssub[3]
        cgsub[4] = cssub[4]
        crevp[4] = constants.HLTC**2 * constants.VDIFU

        cgfr_0 = (
            20.0e2
            * constants.PISQ
            * constants.RNZR
            * constants.RHOR
            / constants.ACT[1] ** 1.75
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
            * constants.PIE
            * constants.TCOND
            * constants.RNZS
            / constants.HLTF
        )
        csmlt[1] = (
            2.0
            * constants.PIE
            * constants.VDIFU
            * constants.RNZS
            * constants.HLTC
            / constants.HLTF
        )
        csmlt[2] = cssub[1]
        csmlt[3] = cssub[2]
        csmlt[4] = constants.CH2O / constants.HLTF

        csmlt_0 = csmlt[0]
        csmlt_1 = csmlt[1]
        csmlt_2 = csmlt[2]
        csmlt_3 = csmlt[3]
        csmlt_4 = csmlt[4]

        # gmlt: five constants

        cgmlt = np.zeros(5)
        cgmlt[0] = (
            2.0
            * constants.PIE
            * constants.TCOND
            * constants.RNZG
            / constants.HLTF
        )
        cgmlt[1] = (
            2.0
            * constants.PIE
            * constants.VDIFU
            * constants.RNZG
            * constants.HLTC
            / constants.HLTF
        )
        cgmlt[2] = cgsub[1]
        cgmlt[3] = cgsub[2]
        cgmlt[4] = constants.CH2O / constants.HLTF

        cgmlt_0 = cgmlt[0]
        cgmlt_1 = cgmlt[1]
        cgmlt_2 = cgmlt[2]
        cgmlt_3 = cgmlt[3]
        cgmlt_4 = cgmlt[4]

        es0 = 6.107799961e2  # ~6.1 mb
        ces0 = constants.EPS * es0

        # make temporaries
        self.TESTVAR_1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_2 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.csaci = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.csaci.view[:] = csaci
        self.csacw = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.csacw.view[:] = csacw
        self.c_psaci = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.c_psaci.view[:] = c_psaci

        # make outputs
        self.subl1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        orchestrate(obj=self, config=self.stencil_factory.config.dace_config)
        self._stencil = self.stencil_factory.from_dims_halo(
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

        self.sat_tables = get_tables(self.stencil_factory.backend)

        self._stencil(
            t1_icloud,
            p1_icloud,
            dp1_icloud,
            qv_icloud,
            ql_icloud,
            qr_icloud,
            qi_icloud,
            qs_icloud,
            qg_icloud,
            qa_icloud,
            den1_icloud,
            denfac_icloud,
            vts_icloud,
            vtg_icloud,
            vtr_icloud,
            self.subl1,
            rh_limited_icloud,
            ccn_icloud,
            cnv_frc_icloud,
            srf_type_icloud,
            self.sat_tables.table1,
            self.sat_tables.table2,
            self.sat_tables.table3,
            self.sat_tables.table4,
            self.sat_tables.des1,
            self.sat_tables.des2,
            self.sat_tables.des3,
            self.sat_tables.des4,
        )

        return {
            "t1_icloud": t1_icloud.view[:],
            "qv_icloud": qv_icloud.view[:],
            "ql_icloud": ql_icloud.view[:],
            "qr_icloud": qr_icloud.view[:],
            "qi_icloud": qi_icloud.view[:],
            "qs_icloud": qs_icloud.view[:],
            "qg_icloud": qg_icloud.view[:],
            "qa_icloud": qa_icloud.view[:],
            "subl1_icloud": self.subl1.view[:],
        }
