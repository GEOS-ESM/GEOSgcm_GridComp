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
            "NACTL": {},
            "ZL0": {},
            "CCN_OCN": {},
            "T": {},
            "PLmb": {},
            "EUAP": {},
            "FRLAND": {},
            "QLCN": {},
            "Q": {},
            "TMP3D": {},
            "QILS": {},
            "QLLS": {},
            "KPBL": {},
            "ZLEO": {},
            "PLE": {},
            "NACTI": {},
            "CCN_LND": {},
            "TKE": {},
            "QICN": {},
            "SH": {},
            "NWFA": {},
        }

        #Float/Int Inputs
        self.in_vars["parameters"] = ["n_modes"]

        #FloatField Outputs
        self.out_vars["data_vars"] = {
            "NACTL": {},
            "NACTI": {},
            "NWFA": {},
        }

     #Calculated Outputs
    def compute_from_storage(self, inputs):
        outputs = {
            "NACTL": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "NACTI": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
            "NWFA": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown")
        }
        self.compute_func(**inputs, **outputs)
        return {**outputs}