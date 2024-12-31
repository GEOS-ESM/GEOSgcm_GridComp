import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from ndsl import Namelist, Quantity, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_core import (
    init_temporaries,
    fix_negative_values,
)


class TranslateGFDL_1M_driver_preloop(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "u_GFDL_1M_driver_preloop": {},
            "v_GFDL_1M_driver_preloop": {},
            "w_GFDL_1M_driver_preloop": {},
            "area_GFDL_1M_driver_preloop": {},
            "qs_GFDL_1M_driver_preloop": {},
            "qi_GFDL_1M_driver_preloop": {},
            "qg_GFDL_1M_driver_preloop": {},
            "ql_GFDL_1M_driver_preloop": {},
            "qr_GFDL_1M_driver_preloop": {},
            "qa_GFDL_1M_driver_preloop": {},
            "qn_GFDL_1M_driver_preloop": {},
            "qv_GFDL_1M_driver_preloop": {},
            "t_GFDL_1M_driver_preloop": {},
            "dp_GFDL_1M_driver_preloop": {},
            "dz_GFDL_1M_driver_preloop": {},
            "rhcrit_GFDL_1M_driver_preloop": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "t_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "p_dry_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "ql_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "qs_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "qg_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "qi_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "qr_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "qv_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "dp_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "den_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "c_praut_GFDL_1M_driver_preloop": self.grid.compute_dict(),
            "omq_GFDL_1M_driver_preloop": self.grid.compute_dict(),
        }

        # initialize temporaries
        self.t1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.dp1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.omq = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qv1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ql1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qr1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qi1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qs1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qg1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.qa1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.den = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.p_dry = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.m1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.u1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.v1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.w1 = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ccn = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.c_praut = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.rh_limited = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.ze = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.zt = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        self.cvm = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.lhi = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.icpk = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.is_frozen = self.quantity_factory.ones([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.precip_fall = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.hold_data = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # set externals manually to the need for avoid another workaround
        c_paut = 1
        cpaut = c_paut * 0.104 * driver_constants.grav / 1.717e-5
        c_air = driver_constants.cp_air
        c_vap = driver_constants.cp_vap
        d0_vap = c_vap - driver_constants.c_liq
        lv00 = driver_constants.hlv0 - d0_vap * driver_constants.t_ice

        # initalize stencils
        orchestrate(obj=self, config=stencil_factory.config.dace_config)
        self._create_temporaries = stencil_factory.from_dims_halo(
            func=init_temporaries,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "cpaut": cpaut,
            },
        )

        self._gfdl_1m_driver_preloop = stencil_factory.from_dims_halo(
            func=fix_negative_values,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "c_air": c_air,
                "c_vap": c_vap,
                "d0_vap": d0_vap,
                "lv00": lv00,
            },
        )

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
        u = self.make_ijk_field(inputs["u_GFDL_1M_driver_preloop"])
        v = self.make_ijk_field(inputs["v_GFDL_1M_driver_preloop"])
        w = self.make_ijk_field(inputs["w_GFDL_1M_driver_preloop"])
        area = self.make_ij_field(inputs["area_GFDL_1M_driver_preloop"])
        qs = self.make_ijk_field(inputs["qs_GFDL_1M_driver_preloop"])
        qi = self.make_ijk_field(inputs["qi_GFDL_1M_driver_preloop"])
        qg = self.make_ijk_field(inputs["qg_GFDL_1M_driver_preloop"])
        ql = self.make_ijk_field(inputs["ql_GFDL_1M_driver_preloop"])
        qr = self.make_ijk_field(inputs["qr_GFDL_1M_driver_preloop"])
        qa = self.make_ijk_field(inputs["qa_GFDL_1M_driver_preloop"])
        qn = self.make_ijk_field(inputs["qn_GFDL_1M_driver_preloop"])
        qv = self.make_ijk_field(inputs["qv_GFDL_1M_driver_preloop"])
        t = self.make_ijk_field(inputs["t_GFDL_1M_driver_preloop"])
        dp = self.make_ijk_field(inputs["dp_GFDL_1M_driver_preloop"])
        dz = self.make_ijk_field(inputs["dz_GFDL_1M_driver_preloop"])
        rhcrit = self.make_ijk_field(inputs["rhcrit_GFDL_1M_driver_preloop"])

        # make outputs
        self.vti = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vts = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vtg = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.onemsig = self.quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        self._create_temporaries(
            t,
            dp,
            rhcrit,
            qv,
            ql,
            qi,
            qr,
            qs,
            qg,
            qa,
            qn,
            dz,
            u,
            v,
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

        return {
            "t_GFDL_1M_driver_preloop": self.t1.view[:],
            "p_dry_GFDL_1M_driver_preloop": self.p_dry.view[:],
            "ql_GFDL_1M_driver_preloop": self.ql1.view[:],
            "qs_GFDL_1M_driver_preloop": self.qs1.view[:],
            "qg_GFDL_1M_driver_preloop": self.qg1.view[:],
            "qi_GFDL_1M_driver_preloop": self.qi1.view[:],
            "qr_GFDL_1M_driver_preloop": self.qr1.view[:],
            "qv_GFDL_1M_driver_preloop": self.qv1.view[:],
            "dp_GFDL_1M_driver_preloop": self.dp1.view[:],
            "den_GFDL_1M_driver_preloop": self.den.view[:],
            "c_praut_GFDL_1M_driver_preloop": self.c_praut.view[:],
            "omq_GFDL_1M_driver_preloop": self.omq.view[:],
        }
