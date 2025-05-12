from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from ndsl.stencils.testing.grid import Grid
from pyMoist.GFDL_1M.PhaseChange.sublimate import sublimate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


class Translatesublimate(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "p_mb": grid.compute_dict() | {"serialname": "PLmb"},
            "t": grid.compute_dict() | {"serialname": "T"},
            "q": grid.compute_dict() | {"serialname": "Q"},
            "qlcn": grid.compute_dict() | {"serialname": "QLCN"},
            "qicn": grid.compute_dict() | {"serialname": "QICN"},
            "clcn": grid.compute_dict() | {"serialname": "CLCN"},
            "nactl": grid.compute_dict() | {"serialname": "NACTL"},
            "nacti": grid.compute_dict() | {"serialname": "NACTI"},
            "qst": grid.compute_dict() | {"serialname": "QST3"},
            "sublc": grid.compute_dict() | {"serialname": "SUBLC"},
        }

        self.in_vars["parameters"] = [
            "CCI_EVAP_EFF_s",
            "DT_MOIST_s",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["p_mb"],
            self.out_vars["nactl"],
            self.out_vars["nacti"],
            self.out_vars["qst"],
        )

    def compute_from_storage(self, inputs):
        _sublimate = self.stencil_factory.from_dims_halo(
            func=sublimate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "CCI_EVAP_EFF": inputs.pop("CCI_EVAP_EFF_s"),
                "DT_MOIST": inputs.pop("DT_MOIST_s"),
            },
        )

        _sublimate(**inputs)

        return inputs
