from f90nml import Namelist

import pyMoist.constants as constants
from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.dsl.typing import Float, Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array
from pyMoist.UW.compute_uwshcu import compute_uwshcu_invert_after
from pyMoist.UW.config import UWConfiguration


class TranslateComputeUwshcuInvertAfter(TranslateFortranData2Py):
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

        self._invert_after = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu_invert_after,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"ncnst": 23},
        )

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "cnvtr": {},
            "cufrc_out": {},
            "cush_inout": {},
            "dcm_out": {},
            "evap": {},
            "fdr_out": {},
            "fer_out": {},
            "ndrop_out": {},
            "nice_out": {},
            "qidet_out": {},
            "qisub_out": {},
            "qiten_out": {},
            "qlten_out": {},
            "qldet_out": {},
            "qlsub_out": {},
            "qtflx_out": {},
            "slflx_out": {},
            "tr0_inout": {},
            "uflx_out": {},
            "vflx_out": {},
            "qpert_out": {},
            "qrten_out": {},
            "qsten_out": {},
            "qvten_out": {},
            "shfx": {},
            "sten_out": {},
            "tpert_out": {},
            "uten_out": {},
            "vten_out": {},
            "umf_out": {},
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
            "cufrc_inv": self.grid.compute_dict(),
            "cush": self.grid.compute_dict(),
            "dcm_inv": self.grid.compute_dict(),
            "dotransport": self.grid.compute_dict(),
            "fdr_inv": self.grid.compute_dict(),
            "fer_inv": self.grid.compute_dict(),
            "qidet_inv": self.grid.compute_dict(),
            "qisub_inv": self.grid.compute_dict(),
            "qiten_inv": self.grid.compute_dict(),
            "qldet_inv": self.grid.compute_dict(),
            "qlsub_inv": self.grid.compute_dict(),
            "qlten_inv": self.grid.compute_dict(),
            "qrten_inv": self.grid.compute_dict(),
            "qsten_inv": self.grid.compute_dict(),
            "qtflx_inv": self.grid.compute_dict(),
            "slflx_inv": self.grid.compute_dict(),
            "tten_inv": self.grid.compute_dict(),
            "qvten_inv": self.grid.compute_dict(),
            "uflx_inv": self.grid.compute_dict(),
            "vflx_inv": self.grid.compute_dict(),
            "ndrop_inv": self.grid.compute_dict(),
            "nice_inv": self.grid.compute_dict(),
            "umf_inv": self.grid.compute_dict(),
            "uten_inv": self.grid.compute_dict(),
            "vten_inv": self.grid.compute_dict(),
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
        fdr_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(fdr_out.view[:], inputs["fdr_out"])
        fer_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(fer_out.view[:], inputs["fer_out"])
        qiten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qiten_out.view[:], inputs["qiten_out"])
        qlten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qlten_out.view[:], inputs["qlten_out"])
        qtflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtflx_out.view[:], inputs["qtflx_out"])
        slflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(slflx_out.view[:], inputs["slflx_out"])
        tr0_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        safe_assign_array(tr0_inout.view[:], inputs["tr0_inout"])
        uflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(uflx_out.view[:], inputs["uflx_out"])
        vflx_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(vflx_out.view[:], inputs["vflx_out"])
        cush_inout = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        safe_assign_array(cush_inout.view[:], inputs["cush_inout"])
        ndrop_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(ndrop_out.view[:], inputs["ndrop_out"])
        nice_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(nice_out.view[:], inputs["nice_out"])
        umf_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf_out.view[:], inputs["umf_out"])
        dcm_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(dcm_out.view[:], inputs["dcm_out"])
        qvten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qvten_out.view[:], inputs["qvten_out"])
        sten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(sten_out.view[:], inputs["sten_out"])
        uten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(uten_out.view[:], inputs["uten_out"])
        vten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(vten_out.view[:], inputs["vten_out"])
        qrten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qrten_out.view[:], inputs["qrten_out"])
        qsten_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qsten_out.view[:], inputs["qsten_out"])
        cufrc_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(cufrc_out.view[:], inputs["cufrc_out"])
        qldet_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qldet_out.view[:], inputs["qldet_out"])
        qidet_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qidet_out.view[:], inputs["qidet_out"])
        qlsub_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qlsub_out.view[:], inputs["qlsub_out"])
        qisub_out = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        safe_assign_array(qisub_out.view[:], inputs["qisub_out"])

        # Outputs
        fdr_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        fer_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qidet_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qisub_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qldet_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qlsub_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qtflx_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        slflx_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        tr0 = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        uflx_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        vflx_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        ndrop_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        nice_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        umf_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM], units="n/a")
        dcm_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qvten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qlten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qiten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        tten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        uten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        vten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qrten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        qsten_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cufrc_inv = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="n/a")
        cnvtr = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM, "ntracers"], units="n/a")
        cush = self.quantity_factory.zeros(dims=[X_DIM, Y_DIM], units="n/a")
        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._invert_after(
            k0=k0,
            umf_outvar=umf_out,
            qtflx_outvar=qtflx_out,
            slflx_outvar=slflx_out,
            uflx_outvar=uflx_out,
            vflx_outvar=vflx_out,
            dcm_outvar=dcm_out,
            qvten_outvar=qvten_out,
            qlten_outvar=qlten_out,
            qiten_outvar=qiten_out,
            sten_outvar=sten_out,
            uten_outvar=uten_out,
            vten_outvar=vten_out,
            qrten_outvar=qrten_out,
            qsten_outvar=qsten_out,
            cufrc_outvar=cufrc_out,
            qldet_outvar=qldet_out,
            qidet_outvar=qidet_out,
            qlsub_outvar=qlsub_out,
            qisub_outvar=qisub_out,
            fer_outvar=fer_out,
            fdr_outvar=fdr_out,
            ndrop_out=ndrop_out,
            nice_out=nice_out,
            tr0=tr0,
            tr0_inoutvar=tr0_inout,
            cush_inoutvar=cush_inout,
            # Outputs
            umf_inv=umf_inv,
            dcm_inv=dcm_inv,
            qtflx_inv=qtflx_inv,
            slflx_inv=slflx_inv,
            uflx_inv=uflx_inv,
            vflx_inv=vflx_inv,
            qvten_inv=qvten_inv,
            qlten_inv=qlten_inv,
            qiten_inv=qiten_inv,
            tten_inv=tten_inv,
            uten_inv=uten_inv,
            vten_inv=vten_inv,
            qrten_inv=qrten_inv,
            qsten_inv=qsten_inv,
            cufrc_inv=cufrc_inv,
            fer_inv=fer_inv,
            fdr_inv=fdr_inv,
            ndrop_inv=ndrop_inv,
            nice_inv=nice_inv,
            qldet_inv=qldet_inv,
            qlsub_inv=qlsub_inv,
            qidet_inv=qidet_inv,
            qisub_inv=qisub_inv,
            CNV_Tracers=cnvtr,
            dotransport=dotransport,
            cush=cush,
        )

        return {
            "cufrc_inv": cufrc_inv.view[:],
            "cush": cush.view[:],
            "dcm_inv": dcm_inv.view[:],
            "dotransport": dotransport,
            "fdr_inv": fdr_inv.view[:],
            "fer_inv": fer_inv.view[:],
            "qidet_inv": qidet_inv.view[:],
            "qisub_inv": qisub_inv.view[:],
            "qiten_inv": qiten_inv.view[:],
            "qldet_inv": qldet_inv.view[:],
            "qlsub_inv": qlsub_inv.view[:],
            "qlten_inv": qlten_inv.view[:],
            "qrten_inv": qrten_inv.view[:],
            "qsten_inv": qsten_inv.view[:],
            "qtflx_inv": qtflx_inv.view[:],
            "slflx_inv": slflx_inv.view[:],
            "tten_inv": tten_inv.view[:],
            "qvten_inv": qvten_inv.view[:],
            "uflx_inv": uflx_inv.view[:],
            "vflx_inv": vflx_inv.view[:],
            "ndrop_inv": ndrop_inv.view[:],
            "nice_inv": nice_inv.view[:],
            "umf_inv": umf_inv.view[:],
            "uten_inv": uten_inv.view[:],
            "vten_inv": vten_inv.view[:],
        }
