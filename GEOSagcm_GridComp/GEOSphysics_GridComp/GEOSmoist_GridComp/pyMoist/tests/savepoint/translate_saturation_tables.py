from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation import QSat
from pyMoist.saturation.table import get_table


class Translatesaturation_tables(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        self.nmodes_quantity_factory = QSat.make_extra_dim_quantity_factory(
            self.quantity_factory
        )

        # Get Inputs
        self.in_vars["parameters"] = ["TABLESIZE"]

        # Set Up Outputs
        self.out_vars = {
            "ESTBLE": {},
            "ESTBLW": {},
            "ESTBLX": {},
            "ESTFRZ": {},
            "ESTLQU": {},
        }

    def compute(self, inputs):
        table = get_table()

        inputs.update(
            {
                "ESTBLE": table.ese,
                "ESTBLX": table.esx,
                "ESTBLW": table.esw,
                "ESTFRZ": table.frz,
                "ESTLQU": table.lqu,
            }
        )

        return inputs
