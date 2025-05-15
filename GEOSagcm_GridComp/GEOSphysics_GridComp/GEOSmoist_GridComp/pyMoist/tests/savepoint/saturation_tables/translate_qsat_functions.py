from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.saturation_tables.qsat_functions import QSat_Float_Liquid
import numpy as np


class Translateqsat_functions(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        self.in_vars["data_vars"] = {
            "t": grid.compute_dict() | {"serialname": "T"},
            "p": grid.compute_dict() | {"serialname": "PLmb"},
        }

        # Set Up Outputs
        self.out_vars = {
            "SER_QSATLQU": {},
            "SER_QSATICE": {},
            "SER_DQSL": {},
            "SER_DQSI": {},
        }

    def compute(self, inputs):
        tables = SaturationVaporPressureTable()
        domain = self._grid.compute_dict()
        halo = self._grid.halo
        out_qsatlqu = np.zeros(
            (domain["iend"] - halo + 1, domain["jend"] - halo + 1, domain["kend"] + 1)
        )
        out_qsatice = np.zeros(
            (domain["iend"] - halo + 1, domain["jend"] - halo + 1, domain["kend"] + 1)
        )
        out_dqsl = np.zeros(
            (domain["iend"] - halo + 1, domain["jend"] - halo + 1, domain["kend"] + 1)
        )
        out_dqsi = np.zeros(
            (domain["iend"] - halo + 1, domain["jend"] - halo + 1, domain["kend"] + 1)
        )
        t = inputs.pop("T")
        p = inputs.pop("PLmb")
        print(f"README {type(tables.ese)}")
        # for i in range(0, domain["iend"] - halo):
        #     for j in range(0, domain["jend"] - halo):
        #         for k in range(0, domain["kend"]):
        #             out_qsatlqu, out_dqsl = QSat_Float_Liquid(
        #                 tables.esw, tables.lqu, t[i, j, k], p[i, j, k], True
        #             )

        inputs.update(
            {
                "SER_QSATLQU": None,
                "SER_QSATICE": None,
                "SER_DQSL": None,
                "SER_DQSI": None,
            }
        )

        return inputs
