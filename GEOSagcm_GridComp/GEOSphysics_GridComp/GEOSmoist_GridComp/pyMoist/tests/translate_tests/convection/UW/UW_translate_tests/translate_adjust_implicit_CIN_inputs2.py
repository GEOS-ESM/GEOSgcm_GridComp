from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import int32
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW.compute_uwshcu import adjust_implicit_CIN_inputs2
from pyMoist.convection.UW.config import UWConfiguration


class TranslateAdjustImplicitCINInputs2(TranslateFortranData2Py):
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
            "qtflx": {},
            "slflx": {},
            "uflx": {},
            "ufrc": {},
            "umf": {},
            "vflx": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "qtflx_s": self.grid.compute_dict(),
            "slflx_s": self.grid.compute_dict(),
            "uflx_s": self.grid.compute_dict(),
            "ufrc_s": self.grid.compute_dict(),
            "umf_s": self.grid.compute_dict(),
            "vflx_s": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")
        self.constants["JASON"] = True

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)

        self.quantity_factory.add_data_dimensions(
            {
                "ntracers": config.NCNST,
            }
        )

        self._adjust_implicit_CIN_inputs2 = self.stencil_factory.from_dims_halo(
            func=adjust_implicit_CIN_inputs2,
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

        qtflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(qtflx.view[:], inputs["qtflx"])
        slflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(slflx.view[:], inputs["slflx"])
        uflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(uflx.view[:], inputs["uflx"])
        ufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(ufrc.view[:], inputs["ufrc"])
        umf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(umf.view[:], inputs["umf"])
        vflx = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        safe_assign_array(vflx.view[:], inputs["vflx"])

        dcm = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cufrc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qcu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiu = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fer = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        xco = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cin = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cinlcl = self.quantity_factory.zeros(dims=[I_DIM, J_DIM], units="n/a")
        cbmf = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qc = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_det = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qiten_sink = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # Outputs
        qtflx_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        slflx_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        uflx_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        ufrc_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        umf_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")
        vflx_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_INTERFACE_DIM], units="n/a")

        dcm_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qrten_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qsten_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qldet_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qidet_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qlsub_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        qisub_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cush_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        cufrc_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fer_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")
        fdr_s = self.quantity_factory.zeros(dims=[I_DIM, J_DIM, K_DIM], units="n/a")

        # The iteration you want to test
        iter_test = int32(0)

        # # Call stencils
        self._adjust_implicit_CIN_inputs2(
            condensation=condensation,
            umf_s=umf_s,
            umf_zint=umf,
            dcm=dcm,
            qrten=qrten,
            qsten=qsten,
            cush=cush,
            cufrc=cufrc,
            slflx_s=slflx_s,
            slflx=slflx,
            qtflx_s=qtflx_s,
            qtflx=qtflx,
            uflx_s=uflx_s,
            uflx=uflx,
            vflx_s=vflx_s,
            vflx=vflx,
            qcu=qcu,
            qlu=qlu,
            qiu=qiu,
            fer=fer,
            fdr=fdr,
            xco=xco,
            cin_IJ=cin,
            cinlcl_IJ=cinlcl,
            cbmf=cbmf,
            qc=qc,
            qlten_det=qlten_det,
            qiten_det=qiten_det,
            qlten_sink=qlten_sink,
            qiten_sink=qiten_sink,
            ufrc_s=ufrc_s,
            ufrc=ufrc,
            dcm_s=dcm_s,
            qrten_s=qrten_s,
            qsten_s=qsten_s,
            qldet_s=qldet_s,
            qidet_s=qidet_s,
            qlsub_s=qlsub_s,
            qisub_s=qisub_s,
            cush_s=cush_s,
            cufrc_s=cufrc_s,
            fer_s=fer_s,
            fdr_s=fdr_s,
            iteration=iter_test,
        )

        return {
            "qtflx_s": qtflx_s.view[:],
            "slflx_s": slflx_s.view[:],
            "uflx_s": uflx_s.view[:],
            "ufrc_s": ufrc_s.view[:],
            "umf_s": umf_s.view[:],
            "vflx_s": vflx_s.view[:],
        }
