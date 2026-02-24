from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import calc_ppen
from pyMoist.UW.config import UWConfiguration


# Dev NOTE: The data for this translate test comes from combining two files in
#           a single nc file using the following NCO tool:
#           `ncks -A ComputeUwshcu-In.nc CalcPpen-In.nc`


class TranslateCalcPpen(TranslateFortranData2Py):
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

        self._calc_ppen = self.stencil_factory.from_dims_halo(
            func=calc_ppen,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "bogbot": {},
            "bogtop": {},
            "condensation": {},
            "dp0": {},
            "drage": {},
            "kpen": {},
            "pifc0": {},
            "rhomid0j": {},
            "wtwb": {},
            "wu": {},
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
            "ppen": self.grid.compute_dict(),
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

        # Field inputs
        condensation = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        kpen_IJ = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen_IJ.view[:], inputs["kpen"] - 1)
        bogbot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(bogbot.view[:], inputs["bogbot"])
        bogtop = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(bogtop.view[:], inputs["bogtop"])
        dp0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dp0.view[:], inputs["dp0"])
        drage = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(drage.view[:], inputs["drage"])
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        rhomid0j = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(rhomid0j.view[:], inputs["rhomid0j"])
        wtwb = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(wtwb.view[:], inputs["wtwb"])
        wu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(wu.view[:], inputs["wu"])

        # Outputs
        ppen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        kpen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_ppen(
            condensation=condensation,
            drage=drage,
            bogbot=bogbot,
            bogtop=bogtop,
            pifc0=pifc0,
            kpen_IJ=kpen_IJ,
            kpen=kpen,
            wu=wu,
            rhomid0j=rhomid0j,
            dp0=dp0,
            wtwb=wtwb,
            ppen=ppen,
            iteration=iter_test,
        )

        return {
            "ppen": ppen.view[:],
        }
