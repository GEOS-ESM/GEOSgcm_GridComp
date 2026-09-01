from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Int
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import update_output_variables1
from pyMoist.convection.UW.config import UWConfiguration


class TranslateUpdateOutputVars1(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "condensation": {},
            "cufrc": {},
            "cush": {},
            "dcm": {},
            "kinv": {},
            "qiten": {},
            "qlten": {},
            "qrten": {},
            "qsten": {},
            "qvten": {},
            "sten": {},
            "umf": {},
            "uten": {},
            "vten": {},
            "zifc0": {},
        }

        # Float/Int Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "cufrc_out": self.grid.compute_dict(),
            "cush_inout": self.grid.compute_dict(),
            "dcm_out": self.grid.compute_dict(),
            "qiten_out": self.grid.compute_dict(),
            "qlten_out": self.grid.compute_dict(),
            "qrten_out": self.grid.compute_dict(),
            "qsten_out": self.grid.compute_dict(),
            "qvten_out": self.grid.compute_dict(),
            "sten_out": self.grid.compute_dict(),
            "umf_out": self.grid.compute_dict(),
            "uten_out": self.grid.compute_dict(),
            "vten_out": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")
        self.constants["JASON"] = True

    def compute(self, inputs):
        self.UW_config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": self.UW_config.NCNST,
            }
        )

        self._update_output_vars1 = self.stencil_factory.from_dims_halo(
            func=update_output_variables1,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        # Inputs
        condensation = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a", dtype=bool)

        for i in range(0, 24):
            for j in range(0, 24):
                if inputs["condensation"][i, j] == 1:
                    condensation.view[i, j] = False
                else:
                    condensation.view[i, j] = True

        cufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(cufrc.view[:], inputs["cufrc"])
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        safe_assign_array(cush.view[:], inputs["cush"])
        dcm = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(dcm.view[:], inputs["dcm"])
        kinv = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a", dtype=Int)
        safe_assign_array(kinv.view[:], inputs["kinv"])
        qiten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qiten.view[:], inputs["qiten"])
        qlten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qlten.view[:], inputs["qlten"])
        qrten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qrten.view[:], inputs["qrten"])
        qsten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qsten.view[:], inputs["qsten"])
        qvten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(qvten.view[:], inputs["qvten"])
        sten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(sten.view[:], inputs["sten"])
        umf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        uten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(uten.view[:], inputs["uten"])
        vten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        safe_assign_array(vten.view[:], inputs["vten"])
        zifc0 = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(zifc0.view[:], inputs["zifc0"])

        # Outputs
        cufrc_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush_inout = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        dcm_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qvten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        sten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        umf_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        vten_out = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # Call stencils
        self._update_output_vars1(
            condensation=condensation,
            umf_zint=umf,
            kinv=kinv,
            zifc0=zifc0,
            dcm=dcm,
            qvten=qvten,
            qlten=qlten,
            qiten=qiten,
            sten=sten,
            uten=uten,
            vten=vten,
            qrten=qrten,
            qsten=qsten,
            cufrc=cufrc,
            cush=cush,
            umf_out=umf_out,
            dcm_out=dcm_out,
            qvten_out=qvten_out,
            qlten_out=qlten_out,
            qiten_out=qiten_out,
            sten_out=sten_out,
            uten_out=uten_out,
            vten_out=vten_out,
            qrten_out=qrten_out,
            qsten_out=qsten_out,
            cufrc_out=cufrc_out,
            cush_inout=cush_inout,
        )

        # For some reason, `cush_inout` is a 3d field in the translate test
        # data. We thus just copy the lowest level into all other levels.
        cush_inout_3d = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        for k in range(self.grid.npz):
            cush_inout_3d[:, :, k] = cush_inout[:, :]

        return {
            "cufrc_out": cufrc_out.view[:],
            "cush_inout": cush_inout_3d.view[:],
            "dcm_out": dcm_out.view[:],
            "qiten_out": qiten_out.view[:],
            "qlten_out": qlten_out.view[:],
            "qrten_out": qrten_out.view[:],
            "qsten_out": qsten_out.view[:],
            "qvten_out": qvten_out.view[:],
            "sten_out": sten_out.view[:],
            "umf_out": umf_out.view[:],
            "uten_out": uten_out.view[:],
            "vten_out": vten_out.view[:],
        }
