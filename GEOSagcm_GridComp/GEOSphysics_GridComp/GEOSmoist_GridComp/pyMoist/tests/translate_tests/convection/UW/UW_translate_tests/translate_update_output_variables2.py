from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

import pyMoist.constants as constants
from pyMoist.convection.UW.compute_uwshcu import update_output_variables2
from pyMoist.convection.UW.config import UWConfiguration


class TranslateUpdateOutputVars2(TranslateFortranData2Py):
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

        self._update_output_vars2 = self.stencil_factory.from_dims_halo(
            func=update_output_variables2,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "fdr": {},
            "fer": {},
            "qiten_det": {},
            "qiten_sink": {},
            "qlten_det": {},
            "qlten_sink": {},
            "qtflx": {},
            "slflx": {},
            "tr0_inout": {},
            "trten": {},
            "uflx": {},
            "vflx": {},
            "kpen": {},
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
            "fdr_out": self.grid.compute_dict(),
            "fer_out": self.grid.compute_dict(),
            "qidet_out": self.grid.compute_dict(),
            "qisub_out": self.grid.compute_dict(),
            "qldet_out": self.grid.compute_dict(),
            "qlsub_out": self.grid.compute_dict(),
            "qtflx_out": self.grid.compute_dict(),
            "slflx_out": self.grid.compute_dict(),
            "tr0_inout": self.grid.compute_dict(),
            "uflx_out": self.grid.compute_dict(),
            "vflx_out": self.grid.compute_dict(),
            "ndrop_out": self.grid.compute_dict(),
            "nice_out": self.grid.compute_dict(),
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
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        fdr = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(fdr.view[:], inputs["fdr"])
        fer = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(fer.view[:], inputs["fer"])
        kpen = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kpen.view[:], inputs["kpen"] - 1)
        qiten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qiten_det.view[:], inputs["qiten_det"])
        qiten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qiten_sink.view[:], inputs["qiten_sink"])
        qlten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qlten_det.view[:], inputs["qlten_det"])
        qlten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qlten_sink.view[:], inputs["qlten_sink"])
        qtflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtflx.view[:], inputs["qtflx"])
        slflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(slflx.view[:], inputs["slflx"])
        tr0_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0_inout.view[:], inputs["tr0_inout"])
        trten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        safe_assign_array(trten.view[:], inputs["trten"])
        uflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(uflx.view[:], inputs["uflx"])
        vflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(vflx.view[:], inputs["vflx"])

        # Outputs
        fdr_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fer_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qidet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qisub_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qldet_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlsub_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qtflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        tr0_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM, "ntracers"], units="n/a")
        uflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        ndrop_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        nice_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._update_output_vars2(
            condensation=condensation,
            dotransport=dotransport,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qlsub_out=qlsub_out,
            qisub_out=qisub_out,
            ndrop_out=ndrop_out,
            nice_out=nice_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            fer=fer,
            fdr=fdr,
            kpen=kpen,
            fer_out=fer_out,
            fdr_out=fdr_out,
            dt=dt,
            rdrop=rdrop,
            qlten_det=qlten_det,
            qiten_det=qiten_det,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            qtflx=qtflx,
            slflx=slflx,
            uflx=uflx,
            vflx=vflx,
            tr0_inout=tr0_inout,
            trten=trten,
        )

        return {
            "fdr_out": fdr_out.view[:],
            "fer_out": fer_out.view[:],
            "qidet_out": qidet_out.view[:],
            "qisub_out": qisub_out.view[:],
            "qldet_out": qldet_out.view[:],
            "qlsub_out": qlsub_out.view[:],
            "qtflx_out": qtflx_out.view[:],
            "slflx_out": slflx_out.view[:],
            "tr0_inout": tr0_inout.view[:],
            "uflx_out": uflx_out.view[:],
            "vflx_out": vflx_out.view[:],
            "ndrop_out": ndrop_out.view[:],
            "nice_out": nice_out.view[:],
        }
