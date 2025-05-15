from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class Translatesaturation_tables(TranslateFortranData2Py):
    def __init__(self, grid: Grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        # Set Up Outputs
        self.out_vars = {
            "ESTBLE": {},
            "ESTBLW": {},
            "ESTBLX": {},
            "ESTFRZ": {},
            "ESTLQU": {},
        }

    def compute(self, inputs):
        tables = SaturationVaporPressureTable()
        inputs.update(
            {
                "ESTBLE": tables.ese,
                "ESTBLX": tables.esx,
                "ESTBLW": tables.esw,
                "ESTFRZ": tables.frz,
                "ESTLQU": tables.lqu,
            }
        )

        return inputs
