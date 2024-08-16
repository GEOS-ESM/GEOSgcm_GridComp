from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation.qsat import QSat
from ndsl.constants import X_DIM, Y_DIM, Z_DIM

class TranslateQSat(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist, 
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = QSat(
            self.stencil_factory,
            self.grid.quantity_factory,
        )
        print(grid.__dict__)
        #FloatField Inputs
        self.in_vars["data_vars"] = {
            "PLmb": {},
            "T": {},
        }

        #Float Inputs
        self.in_vars["parameters"] = []

        #FloatField Outputs
        self.out_vars = {
            "QSAT": self.grid.compute_dict(),
        }

    #Calculated Outputs
    def compute_from_storage(self, inputs):
        self.compute_func(**inputs)
        return {"T": inputs["T"], "QSAT": inputs["QSAT"]}
