from ndsl import Namelist, Quantity, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_tables import get_tables
from pyMoist.GFDL_1M.GFDL_1M_driver.warm_rain import warm_rain_stencil


class Translatewarm_rain(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "dp1_warm_rain": {},
            "dz1_warm_rain": {},
            "t1_warm_rain": {},
            "qv_warm_rain": {},
            "ql_warm_rain": {},
            "qr_warm_rain": {},
            "qi_warm_rain": {},
            "qs_warm_rain": {},
            "qg_warm_rain": {},
            "qa_warm_rain": {},
            "ccn_warm_rain": {},
            "den_warm_rain": {},
            "denfac_warm_rain": {},
            "c_praut_warm_rain": {},
            "vtr_warm_rain": {},
            "evap1_warm_rain": {},
            "m1_rain_warm_rain": {},
            "w1_warm_rain": {},
            "rh_limited_warm_rain": {},
            "eis_warm_rain": {},
            "onemsig_warm_rain": {},
            "dts_warm_rain": {},
            "do_qa_warm_rain": {},
            "rthreshs_warm_rain": {},
            "rthreshu_warm_rain": {},
            "irain_f_warm_rain": {},
            "ql0_max_warm_rain": {},
            "z_slope_liq_warm_rain": {},
            "vr_fac_warm_rain": {},
            "const_vr_warm_rain": {},
            "vr_max_warm_rain": {},
            "vr_min_warm_rain": {},
            "tau_revp_warm_rain": {},
            "lv00_warm_rain": {},
            "d0_vap_warm_rain": {},
            "c_air_warm_rain": {},
            "c_vap_warm_rain": {},
            "crevp_0_warm_rain": {},
            "crevp_1_warm_rain": {},
            "crevp_2_warm_rain": {},
            "crevp_3_warm_rain": {},
            "crevp_4_warm_rain": {},
            "cracw_warm_rain": {},
            "do_sedi_w_warm_rain": {},
            "use_ppm_warm_rain": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "t1_warm_rain": self.grid.compute_dict(),
            "qv_warm_rain": self.grid.compute_dict(),
            "ql_warm_rain": self.grid.compute_dict(),
            "qr_warm_rain": self.grid.compute_dict(),
            "qi_warm_rain": self.grid.compute_dict(),
            "qs_warm_rain": self.grid.compute_dict(),
            "qg_warm_rain": self.grid.compute_dict(),
            "qa_warm_rain": self.grid.compute_dict(),
            "vtr_warm_rain": self.grid.compute_dict(),
            "evap1_warm_rain": self.grid.compute_dict(),
            "m1_rain_warm_rain": self.grid.compute_dict(),
            "r1_warm_rain": self.grid.compute_dict(),
            # "qsat_warm_rain": self.grid.compute_dict(),
            # "dqsdt_warm_rain": self.grid.compute_dict(),
        }

        self.sat_tables = get_tables()

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
        dp1_warm_rain = self.make_ijk_field(inputs["dp1_warm_rain"])
        dp1_original = dp1_warm_rain
        dz1_warm_rain = self.make_ijk_field(inputs["dz1_warm_rain"])
        dz1_original = dz1_warm_rain
        t1_warm_rain = self.make_ijk_field(inputs["t1_warm_rain"])
        t1_original = t1_warm_rain
        qv_warm_rain = self.make_ijk_field(inputs["qv_warm_rain"])
        qv_original = qv_warm_rain
        ql_warm_rain = self.make_ijk_field(inputs["ql_warm_rain"])
        ql_original = ql_warm_rain
        qr_warm_rain = self.make_ijk_field(inputs["qr_warm_rain"])
        qr_original = qr_warm_rain
        qi_warm_rain = self.make_ijk_field(inputs["qi_warm_rain"])
        qi_original = qi_warm_rain
        qs_warm_rain = self.make_ijk_field(inputs["qs_warm_rain"])
        qs_original = qs_warm_rain
        qg_warm_rain = self.make_ijk_field(inputs["qg_warm_rain"])
        qg_original = qg_warm_rain
        qa_warm_rain = self.make_ijk_field(inputs["qa_warm_rain"])
        qa_original = qa_warm_rain
        ccn_warm_rain = self.make_ijk_field(inputs["ccn_warm_rain"])
        ccn_original = ccn_warm_rain
        den_warm_rain = self.make_ijk_field(inputs["den_warm_rain"])
        den_original = den_warm_rain
        denfac_warm_rain = self.make_ijk_field(inputs["denfac_warm_rain"])
        denfac_original = denfac_warm_rain
        c_praut_warm_rain = self.make_ijk_field(inputs["c_praut_warm_rain"])
        c_praut_original = c_praut_warm_rain
        vtr_warm_rain = self.make_ijk_field(inputs["vtr_warm_rain"])
        vtr_original = vtr_warm_rain
        evap1_warm_rain = self.make_ijk_field(inputs["evap1_warm_rain"])
        evap1_original = evap1_warm_rain
        m1_rain_warm_rain = self.make_ijk_field(inputs["m1_rain_warm_rain"])
        m1_rain_original = m1_rain_warm_rain
        w1_warm_rain = self.make_ijk_field(inputs["w1_warm_rain"])
        w1_original = w1_warm_rain
        rh_limited_warm_rain = self.make_ijk_field(inputs["rh_limited_warm_rain"])
        rh_limited_original = rh_limited_warm_rain
        eis_warm_rain = self.make_ij_field(inputs["eis_warm_rain"])
        eis_original = eis_warm_rain
        onemsig_warm_rain = self.make_ij_field(inputs["onemsig_warm_rain"])
        onemsig_original = onemsig_warm_rain

        # Float Variables
        dts = Float(inputs["dts_warm_rain"][0])
        do_qa = Float(inputs["do_qa_warm_rain"][0])
        rthreshs = Float(inputs["rthreshs_warm_rain"][0])
        rthreshu = Float(inputs["rthreshu_warm_rain"][0])
        irain_f = Float(inputs["irain_f_warm_rain"][0])
        ql0_max = Float(inputs["ql0_max_warm_rain"][0])
        z_slope_liq = Float(inputs["z_slope_liq_warm_rain"][0])
        vr_fac = Float(inputs["vr_fac_warm_rain"][0])
        const_vr = Float(inputs["const_vr_warm_rain"][0])
        vr_max = Float(inputs["vr_max_warm_rain"][0])
        vr_min = Float(inputs["vr_min_warm_rain"][0])
        tau_revp = Float(inputs["tau_revp_warm_rain"][0])
        lv00 = Float(inputs["lv00_warm_rain"][0])
        d0_vap = Float(inputs["d0_vap_warm_rain"][0])
        c_air = Float(inputs["c_air_warm_rain"][0])
        c_vap = Float(inputs["c_vap_warm_rain"][0])
        crevp_0 = Float(inputs["crevp_0_warm_rain"][0])
        crevp_1 = Float(inputs["crevp_1_warm_rain"][0])
        crevp_2 = Float(inputs["crevp_2_warm_rain"][0])
        crevp_3 = Float(inputs["crevp_3_warm_rain"][0])
        crevp_4 = Float(inputs["crevp_4_warm_rain"][0])
        cracw = Float(inputs["cracw_warm_rain"][0])
        do_sedi_w = Float(inputs["do_sedi_w_warm_rain"][0])
        use_ppm = Float(inputs["use_ppm_warm_rain"][0])

        # Revert to boolean
        if do_qa == 1:
            do_qa = True
        else:
            do_qa = False
        if z_slope_liq == 1:
            z_slope_liq = True
        else:
            z_slope_liq = False
        if const_vr == 1:
            const_vr = True
        else:
            const_vr = False
        if do_sedi_w == 1:
            do_sedi_w = True
        else:
            do_sedi_w = False

        # make masks
        self.is_frozen = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.precip_fall = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.melting_mask_1 = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.melting_mask_2 = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )
        self.current_k_level = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=Int
        )
        for k in range(self.current_k_level.view[:].shape[2]):
            self.current_k_level.view[:, :, k] = k

        # make temporaries
        self.TESTVAR_1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.TESTVAR_2 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ze = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.zt = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")

        # make outputs
        self.rain = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.evap1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        orchestrate(obj=self, config=self.stencil_factory.config.dace_config)
        self._stencil = self.stencil_factory.from_dims_halo(
            func=warm_rain_stencil,
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
                "vr_min": vr_min,
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

        self._stencil(
            dp1_warm_rain,
            dz1_warm_rain,
            t1_warm_rain,
            qv_warm_rain,
            ql_warm_rain,
            qr_warm_rain,
            qi_warm_rain,
            qs_warm_rain,
            qg_warm_rain,
            qa_warm_rain,
            ccn_warm_rain,
            den_warm_rain,
            denfac_warm_rain,
            c_praut_warm_rain,
            vtr_warm_rain,
            self.evap1,
            m1_rain_warm_rain,
            w1_warm_rain,
            rh_limited_warm_rain,
            eis_warm_rain,
            onemsig_warm_rain,
            self.rain,
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
            self.TESTVAR_1,
            self.TESTVAR_2,
        )

        return {
            "t1_warm_rain": t1_warm_rain.view[:],
            "qv_warm_rain": qv_warm_rain.view[:],
            "ql_warm_rain": ql_warm_rain.view[:],
            "qr_warm_rain": qr_warm_rain.view[:],
            "qi_warm_rain": qi_warm_rain.view[:],
            "qs_warm_rain": qs_warm_rain.view[:],
            "qg_warm_rain": qg_warm_rain.view[:],
            "qa_warm_rain": qa_warm_rain.view[:],
            "vtr_warm_rain": vtr_warm_rain.view[:],
            "evap1_warm_rain": self.evap1.view[:],
            "m1_rain_warm_rain": m1_rain_warm_rain.view[:],
            "r1_warm_rain": self.rain.view[:],
            # "qsat_warm_rain": self.TESTVAR_1.view[:],
            # "dqsdt_warm_rain": self.TESTVAR_2.view[:],
        }
