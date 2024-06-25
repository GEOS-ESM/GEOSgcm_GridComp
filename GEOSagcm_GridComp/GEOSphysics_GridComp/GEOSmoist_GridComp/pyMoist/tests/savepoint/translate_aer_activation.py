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
            "NACTL": self.grid.compute_dict(),
            "ZL0": self.grid.compute_dict(),
            "CCN_OCN": self.grid.compute_dict(),
            "T": self.grid.compute_dict(),
            "PLmb": self.grid.compute_dict(),
            "IM": self.grid.compute_dict(),
            "EUAP": self.grid.compute_dict(),
            "FRLAND": self.grid.compute_dict(),
            "QLCN": self.grid.compute_dict(),
            "Q": self.grid.compute_dict(),
            "JM": self.grid.compute_dict(),
            "TMP3D": self.grid.compute_dict(),
            "QILS": self.grid.compute_dict(),
            "QLLS": self.grid.compute_dict(),
            "LM": self.grid.compute_dict(),
            "KPBL": self.grid.compute_dict(),
            "ZLEO": self.grid.compute_dict(),
            "USE_AERO_BUFFER": self.grid.compute_dict(),
            "PLE": self.grid.compute_dict(),
            "NACTI": self.grid.compute_dict(),
            "CCN_LND": self.grid.compute_dict(),
            "TKE": self.grid.compute_dict(),
            "QICN": self.grid.compute_dict(),
            "SH": self.grid.compute_dict(),
            "NWFA": self.grid.compute_dict()
        }

        #Float/Int Inputs
        self.in_vars["parameters"] = ["n_modes"]

        #FloatField Outputs
        self.out_vars = {
            "NACTL": self.grid.compute_dict(),
            "NACTI": self.grid.compute_dict(),
            "NWFA": self.grid.compute_dict()
        }

    #Calculated Outputs
    def compute_from_storage(self, inputs):
        outputs = {
            #"RAD_QV": self._grid.quantity_factory.zeros(dims=[X_DIM, Y_DIM, Z_DIM], units="unknown"),
        }
        self.compute_func(**inputs, **outputs)
        #return {**outputs}