from ndsl import Namelist, Quantity, StencilFactory, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, Int
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_core import fall_speed
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)


class Translatefall_speed(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "p1_fall_speed": {},
            "cnv_frc_fall_speed": {},
            "anv_icefall_fall_speed": {},
            "lsc_icefall_fall_speed": {},
            "den_fall_speed": {},
            "qs_fall_speed": {},
            "qi_fall_speed": {},
            "qg_fall_speed": {},
            "ql_fall_speed": {},
            "t_fall_speed": {},
            "const_vi_fall_speed": {},
            "const_vs_fall_speed": {},
            "const_vg_fall_speed": {},
            "vi_fac_fall_speed": {},
            "vi_max_fall_speed": {},
            "vs_fac_fall_speed": {},
            "vs_max_fall_speed": {},
            "vg_fac_fall_speed": {},
            "vg_max_fall_speed": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "vti_fall_speed": self.grid.compute_dict(),
            "vts_fall_speed": self.grid.compute_dict(),
            "vtg_fall_speed": self.grid.compute_dict(),
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
        p1 = self.make_ijk_field(inputs["p1_fall_speed"])
        qs = self.make_ijk_field(inputs["qs_fall_speed"])
        qi = self.make_ijk_field(inputs["qi_fall_speed"])
        qg = self.make_ijk_field(inputs["qg_fall_speed"])
        ql = self.make_ijk_field(inputs["ql_fall_speed"])
        den = self.make_ijk_field(inputs["den_fall_speed"])
        t = self.make_ijk_field(inputs["t_fall_speed"])
        cnv_frc = self.make_ij_field(inputs["cnv_frc_fall_speed"])

        # Float Variables
        anv_icefall = Float(inputs["anv_icefall_fall_speed"][0])
        lsc_icefall = Float(inputs["lsc_icefall_fall_speed"][0])
        const_vi = Float(inputs["const_vi_fall_speed"][0])
        const_vs = Float(inputs["const_vs_fall_speed"][0])
        const_vg = Float(inputs["const_vg_fall_speed"][0])
        vi_fac = Float(inputs["vi_fac_fall_speed"][0])
        vi_max = Float(inputs["vi_max_fall_speed"][0])
        vs_fac = Float(inputs["vs_fac_fall_speed"][0])
        vs_max = Float(inputs["vs_max_fall_speed"][0])
        vg_fac = Float(inputs["vg_fac_fall_speed"][0])
        vg_max = Float(inputs["vg_max_fall_speed"][0])

        # Revert to boolean
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

        # make outputs
        self.vti = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vts = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self.vtg = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        orchestrate(obj=self, config=self.stencil_factory.config.dace_config)
        self._stencil = self.stencil_factory.from_dims_halo(
            func=stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
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

        self._stencil(
            p1,
            cnv_frc,
            anv_icefall,
            lsc_icefall,
            den,
            qs,
            qi,
            qg,
            ql,
            t,
            self.vti,
            self.vts,
            self.vtg,
        )

        return {
            "vti_fall_speed": self.vti.view[:],
            "vts_fall_speed": self.vts.view[:],
            "vtg_fall_speed": self.vtg.view[:],
        }


def stencil(
    p1: FloatField,
    cnv_frc: FloatFieldIJ,
    anv_icefall: Float,
    lsc_icefall: Float,
    den: FloatField,
    qs: FloatField,
    qi: FloatField,
    qg: FloatField,
    ql: FloatField,
    t: FloatField,
    vti: FloatField,
    vts: FloatField,
    vtg: FloatField,
):
    with computation(PARALLEL), interval(...):
        vti, vts, vtg = fall_speed(
            p1, cnv_frc, anv_icefall, lsc_icefall, den, qs, qi, qg, ql, t
        )
