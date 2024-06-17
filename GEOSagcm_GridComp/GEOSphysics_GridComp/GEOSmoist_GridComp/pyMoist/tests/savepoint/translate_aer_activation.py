from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.aer_activation import AerActivation
from ndsl.constants import X_DIM, Y_DIM, Z_DIM

class TranslateAerActivation(TranslateFortranData2Py):
    def __init__(
        self,
        grid,
        namelist: Namelist, 
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.compute_func = AerActivation(  # type: ignore
            self.stencil_factory,
            self.grid.quantity_factory,
        )
        self._grid = grid
        self.max_error = 1e-9

        #FloatField Inputs
        self.in_vars["data_vars"] = {
            #"Q": self.grid.compute_dict(),
        }

        #Float Inputs
        self.in_vars["parameters"] = []

        #FloatField Outputs
        self.out_vars = {
            #"Q": self.grid.compute_dict(),
        }

    #Calculated Outputs
    def compute_from_storage(self, inputs):
        outputs = {
            #"RAD_QV": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
        }
        self.compute_func(**inputs, **outputs)
        #return {**outputs}