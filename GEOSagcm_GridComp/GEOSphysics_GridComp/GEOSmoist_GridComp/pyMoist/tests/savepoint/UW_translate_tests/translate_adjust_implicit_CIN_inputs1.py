from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import adjust_implicit_CIN_inputs1
from pyMoist.UW.config import UWConfiguration


class TranslateAdjustImplicitCINInputs1(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
        # UW_config: UWConfiguration,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        # self.UW_config = UW_config

        self._adjust_implicit_CIN_inputs1 = self.stencil_factory.from_dims_halo(
            func=adjust_implicit_CIN_inputs1,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "qi0": {},
            "qiten": {},
            "ql0": {},
            "qlten": {},
            "qv0": {},
            "qvten": {},
            "s0": {},
            "sten": {},
            "t0": {},
            "tr0_AdjustCIN": {},
            "trten": {},
            "u0": {},
            "uten": {},
            "v0": {},
            "vten": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = [
            "dotransport",
            "ncnst",
            "k0",
            "tr0",
            "windsrcavg",
            "qtsrchgt",
            "qtsrc_fac",
            "thlsrc_fac",
            "frc_rasn",
            "rbuoy",
            "epsvarw",
            "use_CINcin",
            "mumin1",
            "rmaxfrac",
            "PGFc",
            "niter_xc",
            "criqc",
            "rle",
            "cridist_opt",
            "mixscale",
            "rkm",
            "dt",
            "detrhgt",
            "rdrag",
            "use_self_detrain",
            "detrhgt",
            "use_cumpenent",
            "rpen",
            "use_momenflx",
            "rdrop",
            "iter_cin",
        ]

        # FloatField Outputs
        self.out_vars = {
            "qi0_s": self.grid.compute_dict(),
            "qiten_s": self.grid.compute_dict(),
            "ql0_s": self.grid.compute_dict(),
            "qlten_s": self.grid.compute_dict(),
            "qv0_s": self.grid.compute_dict(),
            "qvten_s": self.grid.compute_dict(),
            "s0_s": self.grid.compute_dict(),
            "sten_s": self.grid.compute_dict(),
            "t0_s": self.grid.compute_dict(),
            "tr0_s": self.grid.compute_dict(),
            "u0_s": self.grid.compute_dict(),
            "uten_s": self.grid.compute_dict(),
            "v0_s": self.grid.compute_dict(),
            "vten_s": self.grid.compute_dict(),
        }

    def compute(self, inputs):
        self.UW_config = UWConfiguration(Int(inputs["ncnst"]), Int(inputs["k0"]), Int(inputs["windsrcavg"]))

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": constants.NCNST,
            }
        )

        # Float/Int Inputs
        dotransport = Int(inputs["dotransport"])
        k0 = Int(inputs["k0"])
        windsrcavg = Int(inputs["windsrcavg"])
        qtsrchgt = Float(inputs["qtsrchgt"])
        qtsrc_fac = Float(inputs["qtsrc_fac"])
        thlsrc_fac = Float(inputs["thlsrc_fac"])
        frc_rasn = Float(inputs["frc_rasn"])
        rbuoy = Float(inputs["rbuoy"])
        epsvarw = Float(inputs["epsvarw"])
        use_CINcin = Int(inputs["use_CINcin"])
        mumin1 = Float(inputs["mumin1"])
        rmaxfrac = Float(inputs["rmaxfrac"])
        PGFc = Float(inputs["PGFc"])
        dt = Float(inputs["dt"])
        niter_xc = Int(inputs["niter_xc"])
        criqc = Float(inputs["criqc"])
        rle = Float(inputs["rle"])
        cridist_opt = Int(inputs["cridist_opt"])
        mixscale = Float(inputs["mixscale"])
        rdrag = Float(inputs["rdrag"])
        rkm = Float(inputs["rkm"])
        use_self_detrain = Int(inputs["use_self_detrain"])
        detrhgt = Float(inputs["detrhgt"])
        use_cumpenent = Int(inputs["use_cumpenent"])
        rpen = Float(inputs["rpen"])
        use_momenflx = Int(inputs["use_momenflx"])
        rdrop = Float(inputs["rdrop"])
        iter_cin = Int(inputs["iter_cin"])

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        qi0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qi0.view[:], inputs["qi0"])
        qiten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qiten.view[:], inputs["qiten"])
        ql0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ql0.view[:], inputs["ql0"])
        qlten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qlten.view[:], inputs["qlten"])
        qv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qv0.view[:], inputs["qv0"])
        qvten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qvten.view[:], inputs["qvten"])
        s0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(s0.view[:], inputs["s0"])
        sten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(sten.view[:], inputs["sten"])
        t0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(t0.view[:], inputs["t0"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_AdjustCIN"])
        trten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(trten.view[:], inputs["trten"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        uten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(uten.view[:], inputs["uten"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])
        vten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(vten.view[:], inputs["vten"])

        # Outputs
        qi0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qiten_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ql0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qlten_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qv0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qvten_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        s0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        sten_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        t0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        tr0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        u0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        uten_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        v0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        vten_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._adjust_implicit_CIN_inputs1(
            condensation=condensation,
            qv0=qv0,
            qvten=qvten,
            dt=dt,
            ql0=ql0,
            qlten=qlten,
            qi0=qi0,
            qiten=qiten,
            s0=s0,
            sten=sten,
            u0=u0,
            uten=uten,
            v0=v0,
            vten=vten,
            t0=t0,
            dotransport=dotransport,
            tr0_s=tr0_s,
            tr0=tr0,
            trten=trten,
            qv0_s=qv0_s,
            ql0_s=ql0_s,
            qi0_s=qi0_s,
            s0_s=s0_s,
            t0_s=t0_s,
            u0_s=u0_s,
            v0_s=v0_s,
            qvten_s=qvten_s,
            qlten_s=qlten_s,
            qiten_s=qiten_s,
            sten_s=sten_s,
            uten_s=uten_s,
            vten_s=vten_s,
        )

        return {
            "qi0_s": qi0_s.view[:],
            "qiten_s": qiten_s.view[:],
            "ql0_s": ql0_s.view[:],
            "qlten_s": qlten_s.view[:],
            "qv0_s": qv0_s.view[:],
            "qvten_s": qvten_s.view[:],
            "s0_s": s0_s.view[:],
            "sten_s": sten_s.view[:],
            "t0_s": t0_s.view[:],
            "tr0_s": tr0_s.view[:],
            "u0_s": u0_s.view[:],
            "uten_s": uten_s.view[:],
            "v0_s": v0_s.view[:],
            "vten_s": vten_s.view[:],
        }
