from ndsl import Namelist, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.PhaseChange.evaporate import evaporate


class TranslateGFDL_1M_evaporate(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "p_mb": grid.compute_dict() | {"serialname": "PLmb"},
            "t": grid.compute_dict() | {"serialname": "T"},
            "vapor": grid.compute_dict() | {"serialname": "Q"},
            "convective_liquid": grid.compute_dict() | {"serialname": "QLCN"},
            "convective_ice": grid.compute_dict() | {"serialname": "QICN"},
            "convective_cloud_fraction": grid.compute_dict() | {"serialname": "CLCN"},
            "nactl": grid.compute_dict() | {"serialname": "NACTL"},
            "nacti": grid.compute_dict() | {"serialname": "NACTI"},
            "qsat": grid.compute_dict() | {"serialname": "QST3"},
            "evapc": grid.compute_dict() | {"serialname": "EVAPC"},
        }

        self.in_vars["parameters"] = [
            "CCW_EVAP_EFF_s",
            "DT_MOIST_s",
        ]

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["p_mb"],
            self.out_vars["nactl"],
            self.out_vars["nacti"],
            self.out_vars["qsat"],
        )

    def compute_from_storage(self, inputs):
        _evaporate = self.stencil_factory.from_dims_halo(
            func=evaporate,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "CCW_EVAP_EFF": inputs.pop("CCW_EVAP_EFF_s"),
                "DT_MOIST": inputs.pop("DT_MOIST_s"),
            },
        )

        _evaporate(**inputs)

        return inputs
