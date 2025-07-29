from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.PhaseChange.melt_freeze import melt_freeze


class TranslateGFDL_1M_melt_freeze(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "t": grid.compute_dict() | {"serialname": "T"},
            "convection_fraction": grid.compute_dict() | {"serialname": "CNV_FRC"},
            "surface_type": grid.compute_dict() | {"serialname": "SRF_TYPE"},
            "liquid": grid.compute_dict() | {"serialname": "QLCN"},
            "ice": grid.compute_dict() | {"serialname": "QICN"},
        }

        self.in_vars["parameters"] = [
            "DT_MOIST",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["convection_fraction"],
            self.out_vars["surface_type"],
        )

    def compute_from_storage(self, inputs):
        _meltfrz = self.stencil_factory.from_dims_halo(
            func=melt_freeze,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": inputs.pop("DT_MOIST"),
            },
        )

        _meltfrz(**inputs)

        return inputs
