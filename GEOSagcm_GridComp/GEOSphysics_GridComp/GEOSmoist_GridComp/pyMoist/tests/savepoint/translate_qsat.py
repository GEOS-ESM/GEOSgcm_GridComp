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
            "QSAT": {},
            "IT": {},
            "ESTBLE_TEST": {},
            "ESTBLW_TEST": {},
            "ESTBLX_TEST": {},
            "TEMP_IN_QSAT": {},
        }

        #Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "QSAT": self.grid.compute_dict(),
            "IT": self.grid.compute_dict(),
            "ESTBLE_TEST": self.grid.compute_dict(),
            "ESTBLW_TEST": self.grid.compute_dict(),
            "ESTBLX_TEST": self.grid.compute_dict(),
            "TEMP_IN_QSAT": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
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
        QSAT = self.make_ijk_field(inputs["QSAT"])
        IT = self.make_ijk_field(inputs["IT"])

        # Extra Dim Variables
        ESTBLE_TEST = self.make_extra_dim_field(inputs["ESTBLE_TEST"])
        ESTBLW_TEST = self.make_extra_dim_field(inputs["ESTBLW_TEST"])
        ESTBLX_TEST = self.make_extra_dim_field(inputs["ESTBLX_TEST"])
        TEMP_IN_QSAT = self.make_extra_dim_field(inputs["TEMP_IN_QSAT"])


        code(T,
             PL,
             QSAT,
             IT,
             ESTBLE_TEST,
             ESTBLW_TEST,
             ESTBLX_TEST,
             TEMP_IN_QSAT,
        )

        return {
            "QSAT": QSAT.view[:],
            "IT": IT.view[:],
            "ESTBLE_TEST": ESTBLE_TEST.view[0],
            "ESTBLW_TEST": ESTBLW_TEST.view[0],
            "ESTBLX_TEST": ESTBLX_TEST.view[0],
            "TEMP_IN_QSAT": TEMP_IN_QSAT.view[0],
            "T": T.view[:],
        }