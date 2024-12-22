from ndsl import Namelist, Quantity, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver import GFDL_1M_driver
import xarray as xr
from pyMoist.GFDL_1M.GFDL_1M_driver.terminal_fall import terminal_fall
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)


class Translateterminal_fall(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "t_terminal_fall": {},
            "qv_terminal_fall": {},
            "ql_terminal_fall": {},
            "qr_terminal_fall": {},
            "qg_terminal_fall": {},
            "qs_terminal_fall": {},
            "qi_terminal_fall": {},
            "dz1_terminal_fall": {},
            "dp1_terminal_fall": {},
            "den1_terminal_fall": {},
            "vtg_terminal_fall": {},
            "vts_terminal_fall": {},
            "vti_terminal_fall": {},
            "m1_sol_terminal_fall": {},
            "w1_terminal_fall": {},
            "dts_terminal_fall": {},
            "tau_imlt_terminal_fall": {},
            "ql_mlt_terminal_fall": {},
            "c_air_terminal_fall": {},
            "c_vap_terminal_fall": {},
            "d0_vap_terminal_fall": {},
            "lv00_terminal_fall": {},
            "vi_fac_terminal_fall": {},
            "do_sedi_w_terminal_fall": {},
            "use_ppm_terminal_fall": {},
            "tau_smlt_terminal_fall": {},
            "tau_g2r_terminal_fall": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "t_terminal_fall": self.grid.compute_dict(),
            "qv_terminal_fall": self.grid.compute_dict(),
            "ql_terminal_fall": self.grid.compute_dict(),
            "qr_terminal_fall": self.grid.compute_dict(),
            "qg_terminal_fall": self.grid.compute_dict(),
            "qs_terminal_fall": self.grid.compute_dict(),
            "qi_terminal_fall": self.grid.compute_dict(),
            "dz1_terminal_fall": self.grid.compute_dict(),
            "dp1_terminal_fall": self.grid.compute_dict(),
            "den1_terminal_fall": self.grid.compute_dict(),
            "m1_sol_terminal_fall": self.grid.compute_dict(),
            "r1_terminal_fall": self.grid.compute_dict(),
            "g1_terminal_fall": self.grid.compute_dict(),
            "s1_terminal_fall": self.grid.compute_dict(),
            "i1_terminal_fall": self.grid.compute_dict(),
            # "precip_fall_snow": self.grid.compute_dict(),
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
        t_terminal_fall = self.make_ijk_field(inputs["t_terminal_fall"])
        t_original = t_terminal_fall
        qv_terminal_fall = self.make_ijk_field(inputs["qv_terminal_fall"])
        qv_original = qv_terminal_fall
        ql_terminal_fall = self.make_ijk_field(inputs["ql_terminal_fall"])
        ql_original = ql_terminal_fall
        qr_terminal_fall = self.make_ijk_field(inputs["qr_terminal_fall"])
        qr_original = qr_terminal_fall
        qg_terminal_fall = self.make_ijk_field(inputs["qg_terminal_fall"])
        qg_original = qg_terminal_fall
        qs_terminal_fall = self.make_ijk_field(inputs["qs_terminal_fall"])
        qs_original = qs_terminal_fall
        qi_terminal_fall = self.make_ijk_field(inputs["qi_terminal_fall"])
        qi_original = qi_terminal_fall
        dz1_terminal_fall = self.make_ijk_field(inputs["dz1_terminal_fall"])
        dz1_original = dz1_terminal_fall
        dp1_terminal_fall = self.make_ijk_field(inputs["dp1_terminal_fall"])
        dp1_original = dp1_terminal_fall
        den1_terminal_fall = self.make_ijk_field(inputs["den1_terminal_fall"])
        den1_original = den1_terminal_fall
        vtg_terminal_fall = self.make_ijk_field(inputs["vtg_terminal_fall"])
        vtg_original = vtg_terminal_fall
        vts_terminal_fall = self.make_ijk_field(inputs["vts_terminal_fall"])
        vts_original = vts_terminal_fall
        vti_terminal_fall = self.make_ijk_field(inputs["vti_terminal_fall"])
        vti_original = vti_terminal_fall
        m1_sol_terminal_fall = self.make_ijk_field(inputs["m1_sol_terminal_fall"])
        m1_sol_original = m1_sol_terminal_fall
        w1_terminal_fall = self.make_ijk_field(inputs["w1_terminal_fall"])
        w1_original = w1_terminal_fall

        # Float Variables
        dts = Float(inputs["dts_terminal_fall"][0])
        tau_imlt = Float(inputs["tau_imlt_terminal_fall"][0])
        ql_mlt = Float(inputs["ql_mlt_terminal_fall"][0])
        c_air = Float(inputs["c_air_terminal_fall"][0])
        c_vap = Float(inputs["c_vap_terminal_fall"][0])
        d0_vap = Float(inputs["d0_vap_terminal_fall"][0])
        lv00 = Float(inputs["lv00_terminal_fall"][0])
        vi_fac = Float(inputs["vi_fac_terminal_fall"][0])
        do_sedi_w = Float(inputs["do_sedi_w_terminal_fall"][0])
        use_ppm = Float(inputs["use_ppm_terminal_fall"][0])
        tau_smlt = Float(inputs["tau_smlt_terminal_fall"][0])
        tau_g2r = Float(inputs["tau_g2r_terminal_fall"][0])

        # Revert to boolean
        if do_sedi_w == 1:
            do_sedi_w = True
        else:
            do_sedi_w = False
        if use_ppm == 1:
            use_ppm = True
        else:
            use_ppm = False

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
        self.m1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ze = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.zt = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")

        # make outputs
        self.rain = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.snow = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        self.ice = self.quantity_factory.ones([X_DIM, Y_DIM], "n/a")
        self.graupel = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        orchestrate(obj=self, config=self.stencil_factory.config.dace_config)
        self._stencil = self.stencil_factory.from_dims_halo(
            func=terminal_fall,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "dts": dts,
                "tau_imlt": tau_imlt,
                "ql_mlt": ql_mlt,
                "c_air": c_air,
                "c_vap": c_vap,
                "d0_vap": d0_vap,
                "lv00": lv00,
                "vi_fac": vi_fac,
                "do_sedi_w": do_sedi_w,
                "use_ppm": use_ppm,
                "tau_smlt": tau_smlt,
                "tau_g2r": tau_g2r,
            },
        )

        self._stencil(
            t_terminal_fall,
            qv_terminal_fall,
            ql_terminal_fall,
            qr_terminal_fall,
            qg_terminal_fall,
            qs_terminal_fall,
            qi_terminal_fall,
            dz1_terminal_fall,
            dp1_terminal_fall,
            den1_terminal_fall,
            vtg_terminal_fall,
            vts_terminal_fall,
            vti_terminal_fall,
            self.rain,
            self.graupel,
            self.snow,
            self.ice,
            m1_sol_terminal_fall,
            w1_terminal_fall,
            self.ze,
            self.zt,
            self.is_frozen,
            self.precip_fall,
            self.melting_mask_1,
            self.melting_mask_2,
            self.current_k_level,
        )

        return {
            "t_terminal_fall": t_terminal_fall.view[:],
            "qv_terminal_fall": qv_terminal_fall.view[:],
            "ql_terminal_fall": ql_terminal_fall.view[:],
            "qr_terminal_fall": qr_terminal_fall.view[:],
            "qg_terminal_fall": qg_terminal_fall.view[:],
            "qs_terminal_fall": qs_terminal_fall.view[:],
            "qi_terminal_fall": qi_terminal_fall.view[:],
            "dz1_terminal_fall": dz1_terminal_fall.view[:],
            "dp1_terminal_fall": dp1_terminal_fall.view[:],
            "den1_terminal_fall": den1_terminal_fall.view[:],
            "m1_sol_terminal_fall": m1_sol_terminal_fall.view[:],
            "r1_terminal_fall": self.rain.view[:],
            "g1_terminal_fall": self.graupel.view[:],
            "s1_terminal_fall": self.snow.view[:],
            "i1_terminal_fall": self.ice.view[:],
            # "precip_fall_snow": self.precip_fall.view[:],
        }
