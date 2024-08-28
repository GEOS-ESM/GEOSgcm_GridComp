from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation.qsat import QSat


class TranslateQSat(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        self.nmodes_quantity_factory = QSat.make_extra_dim_quantity_factory(
            self.quantity_factory
        )

        #FloatField Inputs
        self.in_vars["data_vars"] = {
            "PL": {},
            "T": {},
        }

        #Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "QSAT": self.grid.compute_dict(),
        }

    def make_ij_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM],
            "n/a",
        )
        qty.view[:, :] = qty.np.asarray(data[:, :])
        return qty
    
    def make_ijk_field(self, data) -> Quantity:
        qty = self.quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM],
            "n/a",
        )
        qty.view[:, :, :] = qty.np.asarray(data[:, :, :])
        return qty
    
    def make_extra_dim_field(self, data) -> Quantity:
        qty = self.nmodes_quantity_factory.empty(
            [Z_DIM, "table_axis"],
            "n/a",
        )
        qty.view[:] = qty.np.asarray(data[:])
        return qty
    
    def compute(self, inputs):
        code = QSat(
            self.stencil_factory,
            self.quantity_factory,
        )

        # FloatField Variables
        T = self.make_ijk_field(inputs["T"])
        PL = self.make_ijk_field(inputs["PL"])


        code(T,
             PL,
        )

        return {
            "QSAT": code.QSat.view[:],
        }