from f90nml import Namelist
from gt4py.cartesian.gtscript import int32

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table
from pyMoist.UW.compute_uwshcu import recalc_environmental_variables
from pyMoist.UW.config import UWConfiguration


class TranslateRecalcEnvVariables(TranslateFortranData2Py):
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

        self._recalc_env_vars = self.stencil_factory.from_dims_halo(
            func=recalc_environmental_variables,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "pifc0": {},
            "pmid0": {},
            "exnmid0": {},
            "qi0_s": {},
            "ql0_s": {},
            "qv0_s": {},
            "s0_s": {},
            "t0_s": {},
            "tr0_RecalcEnv": {},
            "u0": {},
            "v0": {},
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
            "qi0": self.grid.compute_dict(),
            "ql0": self.grid.compute_dict(),
            "qt0": self.grid.compute_dict(),
            "qv0": self.grid.compute_dict(),
            "s0": self.grid.compute_dict(),
            "ssqt0": self.grid.compute_dict(),
            "ssthl0": self.grid.compute_dict(),
            "sstr0": self.grid.compute_dict(),
            "ssu0": self.grid.compute_dict(),
            "ssv0": self.grid.compute_dict(),
            "t0": self.grid.compute_dict(),
            "thl0": self.grid.compute_dict(),
            "thv0bot": self.grid.compute_dict(),
            "thv0top": self.grid.compute_dict(),
            "thvl0top": self.grid.compute_dict(),
            "thvl0bot": self.grid.compute_dict(),
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

        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        pmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(pmid0.view[:], inputs["pmid0"])
        exnmid0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(exnmid0.view[:], inputs["exnmid0"])
        qi0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qi0_s.view[:], inputs["qi0_s"])
        ql0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ql0_s.view[:], inputs["ql0_s"])
        qv0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qv0_s.view[:], inputs["qv0_s"])
        s0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(s0_s.view[:], inputs["s0_s"])
        t0_s = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(t0_s.view[:], inputs["t0_s"])
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0.view[:], inputs["tr0_RecalcEnv"])
        u0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(u0.view[:], inputs["u0"])
        v0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(v0.view[:], inputs["v0"])

        # Outputs
        qi0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ql0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        s0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ssqt0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ssthl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        sstr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        ssu0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        ssv0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        t0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvl0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thv0bot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thv0top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvl0top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        thvl0bot = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        tr0_temp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._recalc_env_vars(
            condensation=condensation,
            qv0_s=qv0_s,
            ql0_s=ql0_s,
            qi0_s=qi0_s,
            s0_s=s0_s,
            t0_s=t0_s,
            exnmid0=exnmid0,
            pmid0=pmid0,
            dotransport=dotransport,
            sstr0=sstr0,
            tr0=tr0,
            u0=u0,
            v0=v0,
            pifc0=pifc0,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            thvl0bot=thvl0bot,
            thv0bot=thv0bot,
            thvl0top=thvl0top,
            thv0top=thv0top,
            thl0=thl0,
            qt0=qt0,
            thvl0=thvl0,
            ssthl0=ssthl0,
            ssu0=ssu0,
            ssv0=ssv0,
            ssqt0=ssqt0,
            qv0=qv0,
            ql0=ql0,
            qi0=qi0,
            s0=s0,
            t0=t0,
            tr0_temp=tr0_temp,
            cush_inout=cush_inout,
            iteration=iter_test,
        )

        return {
            "qi0": qi0.view[:],
            "ql0": ql0.view[:],
            "qt0": qt0.view[:],
            "qv0": qv0.view[:],
            "s0": s0.view[:],
            "ssqt0": ssqt0.view[:],
            "ssthl0": ssthl0.view[:],
            "sstr0": sstr0.view[:],
            "ssu0": ssu0.view[:],
            "ssv0": ssv0.view[:],
            "t0": t0.view[:],
            "thl0": thl0.view[:],
            "thv0bot": thv0bot.view[:],
            "thv0top": thv0top.view[:],
            "thvl0top": thvl0top.view[:],
            "thvl0bot": thvl0bot.view[:],
        }
