from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.driver.sat_tables import get_tables


class TranslateGFDL_driver_tables(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid

        # FloatField Inputs
        self.in_vars["data_vars"] = {}

        # FloatField Outputs
        self.out_vars = {
            "table1_driver": self.grid.compute_dict(),
            "table2_driver": self.grid.compute_dict(),
            "table3_driver": self.grid.compute_dict(),
            "table4_driver": self.grid.compute_dict(),
            "des1_driver": self.grid.compute_dict(),
            "des2_driver": self.grid.compute_dict(),
            "des3_driver": self.grid.compute_dict(),
            "des4_driver": self.grid.compute_dict(),
        }

    def compute(self, inputs):

        self.sat_tables = get_tables(self.stencil_factory.backend)

        return {
            "table1_driver": self.sat_tables.table1,
            "table2_driver": self.sat_tables.table2,
            "table3_driver": self.sat_tables.table3,
            "table4_driver": self.sat_tables.table4,
            "des1_driver": self.sat_tables.des1,
            "des2_driver": self.sat_tables.des2,
            "des3_driver": self.sat_tables.des3,
            "des4_driver": self.sat_tables.des4,
        }
