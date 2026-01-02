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
from pyMoist.UW.compute_uwshcu import calc_cumulus_condensate_at_interface
from pyMoist.UW.config import UWConfiguration


class TranslateCalcCumulusCondensate(TranslateFortranData2Py):
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

        self._calc_cumulus_condensate = self.stencil_factory.from_dims_halo(
            func=calc_cumulus_condensate_at_interface,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "kpen": {},
            "krel": {},
            "pifc0": {},
            "ppen": {},
            "qcubelow": {},
            "qiubelow": {},
            "qlubelow": {},
            "qtu": {},
            "prel": {},
            "qtu_top": {},
            "rcwp": {},
            "riwp": {},
            "rlwp": {},
            "thlu": {},
            "thlu_top": {},
            "ufrc": {},
            "ufrclcl": {},
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
            "qcubelow": self.grid.compute_dict(),
            "qiubelow": self.grid.compute_dict(),
            "qlubelow": self.grid.compute_dict(),
            "rcwp": self.grid.compute_dict(),
            "riwp": self.grid.compute_dict(),
            "rlwp": self.grid.compute_dict(),
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

        prel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(prel.view[:], inputs["prel"])
        qtu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtu.view[:], inputs["qtu"])
        thlu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(thlu.view[:], inputs["thlu"])
        krel = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(krel.view[:], inputs["krel"] - 1)
        kpen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        pifc0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(pifc0.view[:], inputs["pifc0"])
        ppen = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(ppen.view[:], inputs["ppen"])
        qtu_top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(qtu_top.view[:], inputs["qtu_top"])
        thlu_top = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(thlu_top.view[:], inputs["thlu_top"])
        ufrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(ufrc.view[:], inputs["ufrc"])
        ufrclcl = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ufrclcl.view[:], inputs["ufrclcl"])

        # Outputs
        qcubelow = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(qcubelow.view[:], inputs["qcubelow"])
        qiubelow = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(qiubelow.view[:], inputs["qiubelow"])
        qlubelow = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(qlubelow.view[:], inputs["qlubelow"])
        rcwp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(rcwp.view[:], inputs["rcwp"])
        riwp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(riwp.view[:], inputs["riwp"])
        rlwp = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(rlwp.view[:], inputs["rlwp"])

        saturation_vapor_pressure_table = get_saturation_vapor_pressure_table(self.stencil_factory.backend)
        self.ese = saturation_vapor_pressure_table.ese
        self.esx = saturation_vapor_pressure_table.esx

        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        qcu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qlu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qiu = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cufrc = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._calc_cumulus_condensate(
            condensation=condensation,
            krel=krel,
            kpen=kpen,
            pifc0=pifc0,
            ppen=ppen,
            thlu_top=thlu_top,
            qtu_top=qtu_top,
            thlu=thlu,
            qtu=qtu,
            ese=self.ese,
            esx=self.esx,
            umf_out=umf_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            ufrc=ufrc,
            ufrclcl=ufrclcl,
            prel=prel,
            criqc=criqc,
            qcu=qcu,
            qlu=qlu,
            qiu=qiu,
            qcubelow=qcubelow,
            qlubelow=qlubelow,
            qiubelow=qiubelow,
            rcwp=rcwp,
            rlwp=rlwp,
            riwp=riwp,
            cufrc=cufrc,
            cush_inout=cush_inout,
            iteration=iter_test,
        )

        return {
            "qcubelow": qcubelow.view[:],
            "qiubelow": qiubelow.view[:],
            "qlubelow": qlubelow.view[:],
            "rcwp": rcwp.view[:],
            "riwp": riwp.view[:],
            "rlwp": rlwp.view[:],
        }
