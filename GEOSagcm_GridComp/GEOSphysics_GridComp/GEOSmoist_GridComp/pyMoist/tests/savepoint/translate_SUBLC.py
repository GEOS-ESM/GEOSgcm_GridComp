from ndsl import Namelist, Quantity, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.GFDL_1M.GFDL_1M import GFDL_1M


class TranslateSUBLC(TranslateFortranData2Py):
    def __init__(self, grid, namelist: Namelist, stencil_factory: StencilFactory):
        super().__init__(grid, stencil_factory)
        self.stencil_factory = stencil_factory
        self.quantity_factory = grid.quantity_factory
        self._grid = grid
        self.max_error = 1e-9

        # FloatField Inputs
        self.in_vars["data_vars"] = {
            "PLmb": {},
            "CLCN": {},
            "T": {},
            "QST3": {},
            "QICN": {},
            "QLCN": {},
            "NACTI": {},
            "NACTL": {},
            "Q": {},
        }

        # Float Inputs
        self.in_vars["parameters"] = [
            "DT_MOIST_s",
            "CCI_EVAP_EFF_s",
        ]

        # FloatField Outputs
        self.out_vars = {
            "T": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "CLCN": self.grid.compute_dict(),
            "SUBLC": self.grid.compute_dict(),
        }

        self._GFDL_1M = GFDL_1M(self.stencil_factory, self.quantity_factory)

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

    def compute(self, inputs):
        # FloatField Variables
        PLmb = self.make_ijk_field(inputs["PLmb"])
        T = self.make_ijk_field(inputs["T"])
        Q = self.make_ijk_field(inputs["Q"])
        QLCN = self.make_ijk_field(inputs["QLCN"])
        QICN = self.make_ijk_field(inputs["QICN"])
        CLCN = self.make_ijk_field(inputs["CLCN"])
        NACTL = self.make_ijk_field(inputs["NACTL"])
        NACTI = self.make_ijk_field(inputs["NACTI"])
        QST = self.make_ijk_field(inputs["QST3"])

        # Float Variables
        DT_MOIST = Float(inputs["DT_MOIST_s"])
        CCI_EVAP_EFF = Float(inputs["CCI_EVAP_EFF_s"])

        self._GFDL_1M._subl(
            DT_MOIST,
            CCI_EVAP_EFF,
            PLmb,
            T,
            Q,
            QLCN,
            QICN,
            CLCN,
            NACTL,
            NACTI,
            QST,
            self._GFDL_1M.sublc,
        )

        return {
            "T": T.view[:],
            "Q": Q.view[:],
            "QLCN": QLCN.view[:],
            "QICN": QICN.view[:],
            "CLCN": CLCN.view[:],
            "SUBLC": self._GFDL_1M.sublc.view[:],
        }
