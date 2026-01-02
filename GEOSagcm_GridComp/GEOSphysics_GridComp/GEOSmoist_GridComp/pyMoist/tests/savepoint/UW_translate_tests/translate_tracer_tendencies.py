from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import calc_tracer_tendencies
from pyMoist.UW.config import UWConfiguration


class TranslateTracerTendencies(TranslateFortranData2Py):
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

        self._calc_tracer_tendencies = self.stencil_factory.from_dims_halo(
            func=calc_tracer_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "dp0": {},
            "tr0_TracerTendencies": {},
            "trflx": {},
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
            "trten": self.grid.compute_dict(),
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

        dp0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_TracerTendencies"])
        trflx = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        safe_assign_array(trflx.view[:], inputs["trflx"])

        # Outputs
        trten = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        trflx_d = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        trflx_u = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM, "ntracers"], units="n/a")
        trmin = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, "ntracers"], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_tracer_tendencies(
            condensation=condensation,
            dotransport=dotransport,
            k0=k0,
            dt=dt,
            dp0=dp0,
            trflx_d=trflx_d,
            trflx_u=trflx_u,
            trmin=trmin,
            tr0=tr0,
            trflx=trflx,
            trten=trten,
            iteration=iter_test,
        )

        return {
            "trten": trten.view[:],
        }
