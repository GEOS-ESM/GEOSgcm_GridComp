import xarray as xr

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

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "PL": {},
            "T": {},
            "ESTBLE_TEST": {},
            "ESTBLW_TEST": {},
            "ESTBLX_TEST": {},
            "IT": {},
            "IFELSE": {},
            "QSAT_HALFWAY": {},
            "TI": {},
        }

        # Float Inputs
        self.in_vars["parameters"] = []

        # FloatField Outputs
        self.out_vars = {
            "QSAT": self.grid.compute_dict(),
            "ESTBLE_TEST": self.grid.compute_dict(),
            "ESTBLW_TEST": self.grid.compute_dict(),
            "ESTBLX_TEST": self.grid.compute_dict(),
            # "IT": self.grid.compute_dict(),
            # "IFELSE": self.grid.compute_dict(),
            # "QSAT_HALFWAY": self.grid.compute_dict(),
            # "TI": self.grid.compute_dict(),
            # "LOC": self.grid.compute_dict(),
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

        with xr.open_dataset("/home/charleskrop/netcdfs/QSat-Out.nc") as ds:
            ese_array = ds.data_vars["ESTBLE_TEST"].values[0, 0, :]
            esw_array = ds.data_vars["ESTBLW_TEST"].values[0, 0, :]
            esx_array = ds.data_vars["ESTBLX_TEST"].values[0, 0, :]

        # FloatField_Extra_Dim Variables
        ese = self.make_extra_dim_field(ese_array)
        esw = self.make_extra_dim_field(esw_array)
        esx = self.make_extra_dim_field(esx_array)

        code(
            T,
            PL,
            # ese=ese,
            # esw=esw,
            # esx=esx,
        )

        return {
            "QSAT": code.QSat.view[:],
            "ESTBLE_TEST": code.ese.view[0],
            "ESTBLW_TEST": code.esw.view[0],
            "ESTBLX_TEST": code.esx.view[0],
            # "IT": code._IT.view[:],
            # "IFELSE": code._IFELSE.view[:],
            # "QSAT_HALFWAY": code._QSAT_HALFWAY.view[:],
            # "TI": code._TI.view[:],
            # "LOC": code.table._LOC,
        }
