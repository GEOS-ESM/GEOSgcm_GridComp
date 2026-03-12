from f90nml import Namelist

from ndsl import Quantity, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.UW.compute_uwshcu import ComputeUwshcuInv
from pyMoist.UW.config import UWConfiguration
from pyMoist.UW.state import UWState


# Dev NOTE: The data for this translate test comes from combining several ncfiles
# Run ComputeUwshcuInv_postprocess_netcdfs.py before running translate test.


class TranslateComputeUwshcuInv(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "PLE": {},
            "ZLE": {},
            "QLLS": {},
            "QILS": {},
            "QLCN": {},
            "QICN": {},
            "CLCN": {},
            "kpbl_inv": {},
            "u0_inv": {},
            "v0_inv": {},
            "qv0_inv": {},
            "t0_inv": {},
            "frland": {},
            "tke_inv": {},
            "cush": {},
            "shfx": {},
            "evap": {},
            "cnvtr": {},
            "CNV_Tracers": {},
        }

        # FloatField Outputs
        self.out_vars = {
            "CNV_Tracers": self.grid.compute_dict(),
            "CNPCRATE": self.grid.compute_dict(),
            "RKFRE": self.grid.compute_dict(),
            "DQADT_SC": self.grid.compute_dict(),
            "MFD_SC": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
            "QIDET_SC": self.grid.compute_dict(),
            "QIENT_SC": self.grid.compute_dict(),
            "QLDET_SC": self.grid.compute_dict(),
            "QLENT_SC": self.grid.compute_dict(),
            "CNV_MFC": self.grid.compute_dict(),
            "CNV_MFD": self.grid.compute_dict(),
            "SHLW_PRC3": self.grid.compute_dict(),
            "SHLW_SNO3": self.grid.compute_dict(),
            "QITOT": self.grid.compute_dict(),
            "QLTOT": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "U": self.grid.compute_dict(),
            "V": self.grid.compute_dict(),
            "vten_inv": self.grid.compute_dict(),
            "cufrc_inv": self.grid.compute_dict(),
            "qlten_inv": self.grid.compute_dict(),
            "fer_inv": self.grid.compute_dict(),
            "qrten_inv": self.grid.compute_dict(),
            "qisub_inv": self.grid.compute_dict(),
            "tten_inv": self.grid.compute_dict(),
            "qiten_inv": self.grid.compute_dict(),
            "uten_inv": self.grid.compute_dict(),
            "dcm_inv": self.grid.compute_dict(),
            "vflx_inv": self.grid.compute_dict(),
            "qldet_inv": self.grid.compute_dict(),
            "umf_inv": self.grid.compute_dict(),
            "nice_inv": self.grid.compute_dict(),
            "fdr_inv": self.grid.compute_dict(),
            "qidet_inv": self.grid.compute_dict(),
            "qlsub_inv": self.grid.compute_dict(),
            "qtflx_inv": self.grid.compute_dict(),
            "uflx_inv": self.grid.compute_dict(),
            "slflx_inv": self.grid.compute_dict(),
            "dotransport": self.grid.compute_dict(),
            "ndrop_inv": self.grid.compute_dict(),
            "qvten_inv": self.grid.compute_dict(),
            "qsten_inv": self.grid.compute_dict(),
            "tpert_out": self.grid.compute_dict(),
            "qpert_out": self.grid.compute_dict(),
            "cush": self.grid.compute_dict(),
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("ComputeUwshcuInv-constants")

    def compute(self, inputs):
        config = UWConfiguration(**self.constants)
        state = UWState.zeros(
            self.quantity_factory,
            data_dimensions={
                "ntracers": config.NCNST,
            },
        )

        compute_uwshcu = ComputeUwshcuInv(
            self.stencil_factory,
            self.quantity_factory,
            config,
        )

        # Inputs
        state.input.PLE.field[:] = inputs["PLE"]
        state.input.ZLE.field[:] = inputs["ZLE"]
        state.input.QLLS.field[:] = inputs["QLLS"]
        state.input.QILS.field[:] = inputs["QILS"]
        state.input.QLCN.field[:] = inputs["QLCN"]
        state.input.QICN.field[:] = inputs["QICN"]
        state.input.kpbl_inv.field[:] = inputs["kpbl_inv"]
        state.input.frland.field[:] = inputs["frland"]
        state.input.tke_inv.field[:] = inputs["tke_inv"]
        state.input.shfx.field[:] = inputs["shfx"]
        state.input.evap.field[:] = inputs["evap"]

        # In/outs
        state.input_output.u0_inv.field[:] = inputs["u0_inv"]
        state.input_output.v0_inv.field[:] = inputs["v0_inv"]
        state.input_output.qv0_inv.field[:] = inputs["qv0_inv"]
        state.input_output.t0_inv.field[:] = inputs["t0_inv"]
        state.input_output.cush.field[:] = inputs["cush"]
        state.input_output.CNV_Tracers.field[:] = inputs["CNV_Tracers"]
        state.input_output.cnvtr.field[:] = inputs["cnvtr"]
        state.input_output.CLCN.field[:] = inputs["CLCN"]

        compute_uwshcu(
            state,
        )

        return {
            "CNV_Tracers": state.input_output.CNV_Tracers.field[:],
            "CNPCRATE": state.input_output.cnvtr.field[:],
            "RKFRE": state.output.RKFRE.field[:],
            "CLCN": state.input_output.CLCN.field[:],
            "Q": state.input_output.qv0_inv.view[:],
            "T": state.input_output.t0_inv.view[:],
            "U": state.input_output.u0_inv.view[:],
            "V": state.input_output.v0_inv.view[:],
            "QLTOT": state.output.ql0_inv.view[:],
            "QITOT": state.output.qi0_inv.view[:],
            "QIDET_SC": state.output.qidet_inv.view[:],
            "QLDET_SC": state.output.qldet_inv.view[:],
            "MFD_SC": state.output.MFD_SC.view[:],
            "DQADT_SC": state.output.DQADT_SC.view[:],
            "QLENT_SC": state.output.QLENT_SC.view[:],
            "QIENT_SC": state.output.QIENT_SC.view[:],
            "CNV_MFC": state.output.CNV_MFC.view[:],
            "CNV_MFD": state.output.CNV_MFD.view[:],
            "SHLW_PRC3": state.output.SHLW_PRC3.view[:],
            "SHLW_SNO3": state.output.SHLW_SNO3.view[:],
            "umf_inv": state.output.umf_inv.view[:],
            "dcm_inv": state.output.dcm_inv.view[:],
            "qtflx_inv": state.output.qtflx_inv.view[:],
            "slflx_inv": state.output.slflx_inv.view[:],
            "uflx_inv": state.output.uflx_inv.view[:],
            "vflx_inv": state.output.vflx_inv.view[:],
            "qvten_inv": state.output.qvten_inv.view[:],
            "qlten_inv": state.output.qlten_inv.view[:],
            "qiten_inv": state.output.qiten_inv.view[:],
            "tten_inv": state.output.tten_inv.view[:],
            "uten_inv": state.output.uten_inv.view[:],
            "vten_inv": state.output.vten_inv.view[:],
            "qrten_inv": state.output.qrten_inv.view[:],
            "qsten_inv": state.output.qsten_inv.view[:],
            "cufrc_inv": state.output.cufrc_inv.view[:],
            "fer_inv": state.output.fer_inv.view[:],
            "fdr_inv": state.output.fdr_inv.view[:],
            "ndrop_inv": state.output.ndrop_inv.view[:],
            "nice_inv": state.output.nice_inv.view[:],
            "qldet_inv": state.output.qldet_inv.view[:],
            "qlsub_inv": state.output.qlsub_inv.view[:],
            "qidet_inv": state.output.qidet_inv.view[:],
            "qisub_inv": state.output.qisub_inv.view[:],
            "tpert_out": state.output.tpert_out.view[:],
            "qpert_out": state.output.qpert_out.view[:],
            "dotransport": config.dotransport,
            "cush": state.input_output.cush.view[:],
        }
